#!./venv/bin/python
import PseudoNetCDF as pnc
import numpy as np
import cftime
import os
import gc


class TropOMIv001(pnc.PseudoNetCDFFile):
    @classmethod
    def boundingbox(self, path, keys=['time']):
        tmpf = pnc.pncopen(path, format='netcdf')
        out = {}
        if 'time' in keys:
            rtf = pnc.PseudoNetCDFFile()
            rtf.createDimension('time', 1)
            rtf.copyVariable(tmpf['PRODUCT/time'], key='time')
            refdate = rtf.getTimes()[0]
            tunit = refdate.strftime('milliseconds since %F %H:%M:%S+0000')
            
            tf = pnc.PseudoNetCDFFile()
            tf.createDimension('time', 1)
            tf.copyDimension(tmpf['PRODUCT'].dimensions['scanline'])
            tf.copyVariable(tmpf['PRODUCT/delta_time'], key='time')
            tf.variables['time'].units = tunit
            tf = tf.removeSingleton()
            del tmpf
            times = tf.getTimes()
            out['time'] = times.min(), times.max()

        if 'longitude' in keys:
            longitude = tmpf['PRODUCT/longitude'][:]
            out['longitude'] = longitude.min(), longitude.max()

        if 'longitude' in keys:
            latitude = tmpf['PRODUCT/latitude'][:]
            out['latitude'] = latitude.min(), latitude.max()

        return  out

    def __init__(self, path):
        tmpf = pnc.pncopen(path, format='netcdf')
        geogrpk = 'PRODUCT/SUPPORT_DATA/GEOLOCATIONS/'
        outkeys = dict(
            time='PRODUCT/delta_time',
            qa_value='PRODUCT/qa_value',
            latitude='PRODUCT/latitude',
            longitude='PRODUCT/longitude',
            level='PRODUCT/layer',
            hyai='PRODUCT/tm5_constant_a',
            hybi='PRODUCT/tm5_constant_b',
            tropopause_level_index='PRODUCT/tm5_tropopause_layer_index',
            averaging_kernel='PRODUCT/averaging_kernel',
            nitrogendioxide_tropospheric_column='PRODUCT/nitrogendioxide_tropospheric_column',
            air_mass_factor_troposphere='PRODUCT/air_mass_factor_troposphere',
            air_mass_factor_total='PRODUCT/air_mass_factor_total',
            surface_pressure='PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure',
            longitude_bounds=geogrpk + 'longitude_bounds',
            latitude_bounds=geogrpk + 'latitude_bounds',
            viewing_zenith_angle=geogrpk + 'viewing_zenith_angle',
            solar_zenith_angle=geogrpk + 'solar_zenith_angle'
        )
        f = pnc.PseudoNetCDFFile()
        for ok, ik in outkeys.items():
            iv = tmpf[ik]
            for dk, dl in zip(iv.dimensions, iv.shape):
                if dk not in f.dimensions:
                    f.createDimension(dk, dl)
            f.copyVariable(iv, key=ok)
        tf = pnc.PseudoNetCDFFile()
        tf.createDimension('time', 1)
        tf.copyVariable(tmpf['PRODUCT/time'], key='time')
        refdate = tf.getTimes()[0]
        x = np.arange(len(f.dimensions['scanline']))
        y = np.arange(len(f.dimensions['ground_pixel']))
        X, Y = np.meshgrid(x, y)
        outf = f.removeSingleton().slice(
            scanline=X.ravel(),
            ground_pixel=Y.ravel(),
            newdims=('retrieval',)
        ).slice(
            scanline=X.ravel(),
            newdims=('retrieval',)
        ).slice(
            ground_pixel=Y.ravel(),
            newdims=('retrieval',)
        )
        outf.renameDimensions(scanline='retrieval', inplace=True)
        outf.renameDimensions(ground_pixel='retrieval', inplace=True)
        tunit = refdate.strftime('milliseconds since %F %H:%M:%S+0000')
        outf.variables['time'].units = tunit
        self.variables = outf.variables
        self.dimensions = outf.dimensions
        self.setncatts(outf.getncatts())
        del tmpf

    def clip(self, func):
        """
        Arguments
        ---------
        func : function
            function takes longitude and latitude and returns True
            if in domain and False if out of domain

        Returns
        -------
        outf : PseudoNetCDF-file 
            only scanlines in the domain
        """
        longitude=self.variables['longitude'][:]
        latitude=self.variables['latitude'][:]
        result = func(longitude, latitude)
        return self.slice(retrieval=np.where(result))

    def clip_to_pnc(self, pncf):
        def clipper(lon, lat):
            I, J = pncf.ll2ij(lon, lat, clean='mask')
            indomain = (I.mask | J.mask) == False
            return indomain

        return self.clip(clipper)

    def filter(self, minqa=0.7):
        return self.slice(
            retrieval=np.where(
                np.ma.filled(
                    self.variables['qa_value'][:] >= minqa, False
                )
            )
        )

    def approxSigma(self):
        sp = self.variables['surface_pressure'][:].mean()
        pvals = (self.variables['hyai'] +  self.variables['hybi'] * sp)
        pedges = np.append(pvals[:, 0], pvals[-1, -1])
        sigma = (pedges - pedges[-1]) / (pedges[0] - pedges[-1])
        return dict(vglvls=sigma, vgtop=pedges[-1])

    def averaging_kernel_troposphere(self):
        tak = eval(
            'air_mass_factor_total[:][:, None]' +
            '/ air_mass_factor_troposphere[:][:, None]' +
            '* averaging_kernel[:]', None, self.variables
        )
        Kabove = (
            (np.arange(self.variables['level'].size)[None, :] + 1) >
            self.variables['tropopause_level_index'][:][:, None]
        )
        tak[Kabove] = 0
        return tak

    def process_cmaq(self, pncf, key='NO2'):
        ffc = self.clip_to_pnc(pncf)
        if len(ffc.dimensions['retrieval']) == 0:
            raise ValueError('No valid pixels')

        times = ffc.getTimes()
        H = np.array([t.hour for t in times])
        I, J = pncf.ll2ij(ffc.variables['longitude'][:], ffc.variables['latitude'][:])
        pixelf = pncf.slice(
            TSTEP=H, ROW=J, COL=I, newdims=('retrieval',)
        )
        # pixelf.renameDimensions(TSTEP='retrieval', inplace=True)
        tmf = pixelf.interpSigma(
            **ffc.approxSigma(), interptype='conserve'
        ).mask(invalid=True)
        outf = tmf.mask(invalid=True)
        outf.createDimension('layer', len(ffc.dimensions['layer']))
        tak = ffc.averaging_kernel_troposphere()
        ak = ffc.variables['averaging_kernel'][:]
        akv = outf.createVariable(
            'averaging_kernel', 'f', ('retrieval', 'layer')
        )
        akv[:] = ak[:]
        takv = outf.createVariable(
            'averaging_kernel_troposphere', 'f', ('retrieval', 'layer')
        )
        takv[:] = tak[:]
        vcdv = outf.createVariable(
            'nitrogendioxide_total_column', 'f', ('retrieval',)
        )
        vcdv[:] = 0
        tvcdv = outf.createVariable(
            'nitrogendioxide_tropospheric_column', 'f', ('retrieval',)
        )
        tvcdv[:] = 0
        mcdv = outf.createVariable(
            'nitrogendioxide_model_column', 'f', ('retrieval',)
        )
        mcdv[:] = 0
        pedges = eval('surface_pressure[:][:, None, None] * hybi[None, :] + hyai[None, :]', None, ffc.variables)
        # surface_pressure in Pa, need dp in hPa
        dp = (pedges[..., 0] - pedges[..., 1]) / 100

        # https://aura.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level2/OMO3PR.003/doc/README.OMO3PR.pdf
        # http://www.temis.nl/data/conversions.pdf
        # assumes mixing ratio in PPM and dp in hPa
        hPa_to_du = (
            10 * 1.3807e-23 * 6.022e23 / 0.02894 * 273.15 / 9.80665 / 101325.
        )

        ppm = outf.variables[key][:]
        factor = hPa_to_du * 2.69e16 / 6.022e23 * 1e4

        for li in range(len(outf.dimensions['layer'])):
            lc = np.ma.filled(ppm[:, li] * dp[:, li] * factor, 0)
            mcdv[:] += lc
            vcdv[:] += ak[:, li] * lc
            tvcdv[:] += tak[:, li] * lc

        return ffc, outf
