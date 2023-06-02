Instructions For Colab
======================

These directions explain how to get cmaqsatproc working on Google Colab.
Google Colab is a commonly used web-based Jupyter Hub system. It has a
set of pre-installed packages that makes it easy to work with. The
first step is to open a new Notebook on https://colab.research.google.com.
After that, Colab requires two library installations and credential
configuration for OpenDAP. Then, you can process data on Colab.

Library Installations
---------------------
The first two are installations. You must install an older version of `netcdf`
that works with OpenDAP, and you must install `cmaqsatproc`. Both are done
using the commands below in their own Jupyter cell.


```
!pip install 'netcdf4<=1.5.3' h5netcdf
!pip install -qqq https://github.com/barronh/cmaqsatproc/archive/refs/heads/main.zip
```

Credential Configuration for OpenDAP
------------------------------------

Next you must configure the machine to work with URS Earthdata credentials. The
instructions here assume you already have an account. If not, first go to
NASA's https://earthdata.nasa.gov and register for an account. Once you have an
account, configuring Colab to use the account in OpenDAP requires two steps.
First, install a `.netrc` with your credentials using the code below. Copy the
code into a Notebook cell and run that cell. It will prompt you for your
username and password.

```
import getpass
with open('/root/.netrc', 'w') as nrc:
  nrc.write(f"""machine urs.earthdata.nasa.gov
  login {getpass.getpass('URS Username')}
  password {getpass.getpass('URS Password')}
""")
!chmod 600 /root/.netrc
```

Second, you must configure a `.dodsrc` to use the `.netrc` and a `cookie jar`.
Copy the three lines of code below into a Notebook cell and run that cell.

```
%%writefile /root/.dodsrc
HTTP.NETRC=/root/.netrc
HTTP.COOKIEJAR=/root/.urs_cookies
```

Process Data
------------

Now any of the [examples](README.md) will run on Colab. Just cut and past the
example code into a Notebook cell and run it. Note that Colab works in scratch
space, so your results will not be retained unless you download them.
automatically saved an