import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("cmaqsatproc/__init__.py", "r") as fh:
    for l in fh:
        if l.startswith('__version__'):
            exec(l)
            break
    else:
        __version__ = 'x.y.z'

setuptools.setup(
    name="cmaqsatproc",
    version=__version__,
    author="Barron H. Henderson",
    author_email="barronh@gmail.com",
    description="Processor to grid satellite data for comparison to CMAQ.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/barronh/cmaqsatproc",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Development Status :: 2 - Pre-Alpha",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "numpy", "matplotlib", "pandas", "geopandas", "xarray", "pyproj",
        "shapely", "pygeos", "pycno"
    ],
    extras_require={
        "gdal":  ["gdal"],
    }
)
