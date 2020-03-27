import setuptools

with open('/Users/gsharman/Box/detritalPy/pip_installation/README.md', "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="detritalpy",
    packages=['detritalpy'],    
    version="1.3.8",
    author="Glenn Sharman",
    author_email="gsharman@uark.edu",
    description="A Python-based toolset for visualizing and analyzing detrital geo-thermochronologic data",
    keywords = 'detrital, zircon, provenance, stratigraphy, geochronology, thermochronology',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/grsharman/detritalpy",
    download_url="https://github.com/grsharman/detritalPy/archive/v1.3.8.tar.gz",
    install_requires=['numpy','matplotlib','pandas','xlrd','folium','vincent','simplekml',
        'scipy','sklearn','statsmodels','peakutils'],
    classifiers=[
        "Programming Language :: Python :: 3",
        'Intended Audience :: Science/Research',
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
)
