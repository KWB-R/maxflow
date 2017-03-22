# Maxflow

## Requirement: installation of FloPy (modified by KWB to work for Maxflow)

Python versions:

FloPy requires Python 2.7 or Python 3.3 (or higher)

Dependencies:

FloPy requires NumPy 1.9 (or higher) and matplotlib 1.4 (or higher). The mapping and cross-section capabilities in the flopy.plot submodule and shapefile export capabilities (to_shapefile()) require Pyshp 1.2 (or higher). The NetCDF export capabilities in the flopy.export submodule require python-dateutil 2.4 (or higher), netcdf4 1.1 (or higher), and pyproj 1.9 (or higher). Other NetCDF dependencies are detailed on the UniData website. The get_dataframes method in the ListBudget class in the flopy.utils submodule require pandas 0.15 (or higher).

For base Python distributions to install FloPy type:

For Phython 2.7:

pip install https://github.com/mrustl/flopy/archive/3.2.5_kwb.zip

For Python 3.3 (or higher):

pip3 install https://github.com/mrustl/flopy/archive/3.2.5_kwb.zip
