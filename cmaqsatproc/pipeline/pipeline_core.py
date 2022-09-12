__all__ = ['Pipeline']
from ..utils import weight_vars, getcmrlinks, rootremover, csp_formatwarnings


class Pipeline:
    def __init__(
        self, cmaqgrid, short_name, reader, link_filter=None,
        varkeys2d=None, renamer2d=None, outtmpl2d=None,
        varkeys3d=None, renamer3d=None, outtmpl3d=None,
        **cmr_kw
    ):
        """
        The basic process for creating CMAQ gridded satellite data is:
        1. Query the NASA Common Metadata Repository for valid granules.
          * Usually based on a date and spatial overlap with domain
          * requires a CMAQ grid and product short_name.
        2. Choose the links that will provide access
          * typically opendap, but could be a derivative to point to local
            copies
          * requires a function that will filter links
        3. For each link, calculate pixel contributions to grid cells.
          * This process is specific to each satellite reader.
          * requires a reader object with the weights method
        4. For each link, weight pixel values and produces grid cell values.
          * This process is specific to each satellite reader.
          * requires a reader object with the weighted method
        5. Combine weighted pixel values to produce overall value
        6. Rename output variables for IOAPI compliance.
          * requires a custom dictionary
        6. Save to a file on disk
          * requires a string template

        Arguments
        ---------
        cmaqgrid : cmaqsatproc.cmaq.CMAQGrid or str
            If a str, then a CMAQGrid is created using the built in GRIDDESC
        short_name : str
            Recognized by NASA Common Metadata Repository as a unique product.
        reader : cmaqsatproc.readers.satellite subclass
            Usually one of:
                cmaqsatproc.readers.modis.MOD04
                cmaqsatproc.readers.modis.OMNO2 .OMNO2d or .OMHCHO
        link_filter : func
            Takes all links returned by cmaqsatproc.utils.getcmrlinks and
            returns the subset to be used.
        varkeys2d : list
            List of keys for output in 2-dimensions (grid no layer)
        renamer2d : mappable
            Translates satellite variables to IOAPI compliant names
        outtmpl2d : str
            String to be used to make output file paths using
            outtmpl2d.format(date=date)
        varkesy3d : list
            List of keys for output in 3-dimensions (grid with layers)
        renamer3d : mappable
            Translates satellite variables to IOAPI compliant names
        outtmpl3d : str
            String to be used to make output file paths using
            outtmpl3d.format(date=date)
        """
        if isinstance(cmaqgrid, str):
            from .. import cmaq
            cmaqgrid = cmaq.CMAQGrid(gdpath=None, GDNAM=cmaqgrid)
        self.cmaqgrid = cmaqgrid
        self.approxpoly = cmaqgrid.exterior.to_crs(
            4326
        ).geometry.iloc[0].simplify(0.01)
        self.cmr_kw = cmr_kw.copy()
        self.cmr_kw['short_name'] = short_name
        self.link_filter = link_filter
        self.reader = reader
        if varkeys2d is None and renamer2d is not None:
            varkeys2d = list(renamer2d)
        self.varkeys2d = varkeys2d
        self.renamer2d = renamer2d
        if outtmpl2d is None:
            outtmpl2d = (
                '{date:%Y-%m}/{short_name}_{date:%F}_' + cmaqgrid.GDNAM + '.nc'
            )
        self.outtmpl2d = outtmpl2d
        if varkeys3d is None and renamer3d is not None:
            varkeys3d = list(renamer3d)

        self.varkeys3d = varkeys3d
        self.renamer3d = renamer3d
        if outtmpl3d is None:
            outtmpl3d = (
                '{date:%Y-%m}/{short_name}_{date:%F}_' + cmaqgrid.GDNAM + '.nc'
            )
        self.outtmpl3d = outtmpl3d

    def get_links(self, date):
        """
        Arguments
        ---------
        date : str
            Anything that can be interpreted by pandas.to_datetime as a date

        Returns
        -------
        links : list
            Links from utils.getcmrlinks after custom processing by link_filter
        """
        import pandas as pd

        date = pd.to_datetime(date)
        links = getcmrlinks(
            temporal=f'{date:%F}T00:00:00Z/{date:%F}T23:59:59Z',
            poly=self.approxpoly, **self.cmr_kw
        )
        finallinks = self.link_filter(links)
        return finallinks

    def process_dates(self, date_range, verbose=0, output=False, makedirs=True):
        """
        Arguments
        ---------
        date_range : iterable
            Iterable of dates
        """
        import warnings
        import pandas as pd

        if output:
            outputs = []
        for date in date_range:
            date = pd.to_datetime(date)
            try:
                out = self.process_date(
                    date, verbose=verbose, makedirs=makedirs
                )
                if output:
                    outputs.append(out)
                else:
                    for outf in out:
                        outf.close()
                    del out
            except ValueError:
                if output:
                    outputs.append({})
                warnings.warn(
                    f'Processing failed for {date:%F}. No valid pixels'
                )

        if output:
            return outputs

    def process_date(self, date, verbose=0, links=None, **io_kw):
        """
        Arguments
        ---------
        date : str
            Anything that can be interpreted by pandas.to_datetime as a date

        Returns
        -------
        out : list
            List of outputs either dataframes or ioapi-like files. Type depends
            on whether the outtmpl was provided. If yes, then ioapi.
        """
        import warnings
        import pandas as pd
        import os

        with warnings.catch_warnings(record=True) as warning_list:
            date = pd.to_datetime(date)
            cg = self.cmaqgrid

            if verbose > 0:
                print('Start query', flush=True)

            if links is None:
                links = self.get_links(date)
            if verbose > 0:
                print('Start read')
                print(links, flush=True)

            outputs = {
                'links': links,
                'links_with_pixels': [],
                '2d': [],
                '3d': []
            }
            for link in links:
                if verbose > 0:
                    print(f'Starting {link}', flush=True)
                self.sat = sat = self.reader(link)

                if verbose > 1:
                    print('Calculate weights', flush=True)

                self.wgts = wgts = sat.weights(
                    cg.geodf, clip=self.cmaqgrid.exterior
                )
                if wgts.shape[0] == 0:
                    slink = os.path.basename(link)
                    warnings.warn(
                        f'No valid pixels for {slink}'
                    )
                    continue

                outputs['links_with_pixels'].append(link)
                if self.varkeys2d is not None:
                    if verbose > 1:
                        print('Process 2d data', flush=True)

                    # Produce 2d output
                    wgtd2d = sat.weighted(
                        *self.varkeys2d, groupkeys=['ROW', 'COL'], wgtdf=wgts
                    )
                    outputs['2d'].append(wgtd2d)

                if self.varkeys3d is not None:
                    if verbose > 1:
                        print('Process 3d data', flush=True)

                    # Produce 2d output
                    wgtd3d = sat.weighted(
                        *self.varkeys3d, groupkeys=['ROW', 'COL'], wgtdf=wgts
                    )
                    outputs['3d'].append(wgtd3d)

        outputs['history'] = csp_formatwarnings(warning_list)
        outputio = self.to_ioapi(date, outputs, verbose=verbose, **io_kw)
        if len(outputio) == 0:
            # Must be not perist, so return raw output
            return outputs
        else:
            outputio['links'] = outputs['links']
            outputio['links_with_pixels'] = outputs['links_with_pixels']
            outputio['history'] = outputs['history']
            return outputio

    def to_ioapi(
        self, date, inputs, makedirs=True, laykey='nLevels', verbose=0
    ):
        """
        Arguments
        ---------
        date : str
            Anything that can be interpreted by pandas.to_datetime as a date
        inputs : mappable
            If contains 2d and/or 3d and outtmpls are not False, then files
            will be written to disk and replace standard output. Must contain
            links.
        """
        import os
        import pandas as pd

        date = pd.to_datetime(date)
        cg = self.cmaqgrid
        links = inputs.get('links', [])
        short_name = self.cmr_kw['short_name']
        outputs = {}
        if len(inputs.get('2d', [])) > 0 and self.outtmpl2d is not False:
            outpath2d = self.outtmpl2d.format(
                date=date, short_name=short_name
            )
            if verbose > 0:
                print(f'Make {outpath2d}', flush=True)
            wgtd2d = weight_vars(pd.concat(inputs['2d']))
            outf2d = cg.to_ioapi(wgtd2d, rename=self.renamer2d)
            outf2d.FILEDESC = (
                f'{short_name} files for {date:%F} oversampled\n'
                + '\n'.join(rootremover(links, insert=True)[1])
            )[:60 * 80]
            outf2d.HISTORY = inputs['history'][:60 * 80]
            outdir2d = os.path.dirname(outpath2d)
            if makedirs and outdir2d != '':
                os.makedirs(outdir2d, exist_ok=True)
            diskf = outf2d.save(outpath2d, complevel=1, verbose=verbose)
            outputs['2d'] = diskf

        if len(inputs.get('3d', [])) > 0 and self.outtmpl3d is not False:
            outpath3d = self.outtmpl3d.format(
                date=date, short_name=short_name
            )
            if verbose > 0:
                print(f'Make {outpath2d}', flush=True)
            wgtd3d = weight_vars(
                pd.concat(inputs['3d']), groupkeys=['ROW', 'COL', laykey]
            )
            outf3d = cg.to_ioapi(wgtd3d, rename=self.renamer3d)
            outf3d.FILEDESC = (
                f'{short_name} files for {date:%F} oversampled\n'
                + '\n'.join(rootremover(links, insert=True)[1])
            )[:60 * 80]
            outf3d.HISTORY = inputs['history'][:60 * 80]
            outdir3d = os.path.dirname(outpath3d)
            if makedirs and outdir3d != '':
                os.makedirs(outdir3d, exist_ok=True)
            diskf = outf3d.save(outpath3d, complevel=1, verbose=verbose)
            outputs['3d'] = diskf

        return outputs
