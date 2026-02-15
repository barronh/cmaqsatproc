class s3access:
    def __init__(self, credurl, buffer='2min', verbose=0):
        """
        Simple AWS credentialing interface.
        Assumes credentials are in ~/.netrc

        Arguments
        ---------
        credurl : str
            Url to use to refresh credentials.
        buffer : str
            If the credentials will expire in less than buffer, refresh before
            next request.
        verbose : int
            Level of verbosity.
        """
        import pandas as pd
        self._fs = None
        self._credurl = credurl
        self._buffer = buffer
        self._expiration = pd.to_datetime('1970-01-01 00:00:00+00:00')
        self._verbose = verbose

    def get_cred(self):
        """
        Get AWS credentials
        """
        import requests
        # Get credentials using ~/.netrc
        if self._verbose > 0:
            print('INFO:: Requesting AWS credentials...', flush=True)
        credsr = requests.get(self._credurl)
        credsr.raise_for_status()
        try:
            creds = credsr.json()
        except Exception:
            print('ERROR::', credsr.text)
        if self._verbose > 0:
            print('INFO:: Requested AWS credentials', flush=True)
        return creds

    def ls(self, path):
        """
        List the contents of path

        Arguments
        ---------
        path : str
            Path on aws
        Returns
        -------
        paths : list
            List of paths on aws.
        """
        fs = self.get_fs()
        return fs.ls(path)

    def get_fs(self, refresh=False):
        """
        Arguments
        ---------
        refresh : bool
            Force a refresh of credentials.

        Returns
        -------
        fs : s3fs.FileSystem
            File system object
        """
        import s3fs
        import pandas as pd
        if self._verbose > 0:
            print('INFO:: Requesting AWS Filesystem...', flush=True)
        now = pd.to_datetime('now', utc=True)
        refreshtime = self._expiration - pd.to_timedelta(self._buffer)
        refreshnow = refresh or (now > refreshtime)
        if self._fs is None or refreshnow:
            # Build s3fs filesystem (fs) object with credentials for remote
            # file operations
            creds = self.get_cred()
            self._expiration = pd.to_datetime(creds['expiration'])
            fs = self._fs = s3fs.S3FileSystem(
                key=creds['accessKeyId'], secret=creds['secretAccessKey'],
                anon=False, token=creds['sessionToken']
            )
        else:
            fs = self._fs
        if self._verbose > 0:
            print('INFO:: Requested AWS Filesystem.', flush=True)
        return fs

    def download(self, path, dest=None, outdir=None):
        """
        Arguments
        ---------
        path : str
            Path on AWS to download.
        dest : str
            Path to download to on disk (default to path).
        outdir : str
            If dest is None, dest = os.path.join(outdir, path)

        Returns
        -------
        localpath : str
            Path on disk (same as dest)
        """
        import os
        import botocore.exceptions
        if dest is None:
            if outdir is None:
                outdir = '.'
            dest = os.path.join(outdir, path)

        # If not available, download it
        if os.path.exists(dest):
            if self._verbose > 0:
                print('INFO::', dest, 'exists; keeping cached', flush=True)
        else:
            if self._verbose > 0:
                print('INFO:: downloading', path, flush=True)
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            try:
                fs = self.get_fs()
                fs.download(path, dest)
            except botocore.exceptions.ClientError:
                # If fails to read, then assume credentials timed out and
                # force a refresh of the fs
                if self._verbose > 0:
                    print('INFO:: retrying download', path, flush=True)
                fs = self.get_fs(refresh=True)
                fs.download(path, dest)
            if self._verbose > 0:
                print('INFO:: download success!', flush=True)

        return dest
