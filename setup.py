import setuptools

with open('README.md', "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="detritalpy",
    packages=['detritalpy'],    
    version="1.4.3",
    author="Glenn Sharman",
    author_email="gsharman@uark.edu",
    description="A Python-based toolset for visualizing and analyzing detrital geo-thermochronologic data",
    keywords = 'detrital, zircon, provenance, stratigraphy, geochronology, thermochronology',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/grsharman/detritalpy",
    download_url="https://github.com/grsharman/detritalPy/archive/v1.4.3.tar.gz",
    install_requires=['numpy','matplotlib','pandas','xlrd','folium','vincent','simplekml',
        'scipy','scikit-learn','KDEpy','peakutils','openpyxl','xlsxwriter'],
    classifiers=[
        "Programming Language :: Python :: 3",
        'Intended Audience :: Science/Research',
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
)
