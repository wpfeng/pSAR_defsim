# pSAR_defsim
a python package to simulate surface deformation due to co- and post-seismic faulting processes in E/N/U and InSAR line of sight (LOS) directions.
the script is released along with the publication below:

Zhang, Z., Feng, W., Xu, X., & Samsonov, S. (2023). Performance of Common Scene Stacking Atmospheric Correction on Nonlinear InSAR Deformation Retrieval. Remote Sensing, 15(22), 5399.

To run the script for surface deformation, you need to add MY_PYTHON into your personal python path. Following steps can help configure your python enviroments. 
1) please add MY_PYTHON into your PYTHONPATH
2) add pSAR_defSIM.py into your system path.
3) run pSAR_defSIM.py to see help information...

There may still have several dependencies required by our own modules, e.g. pSAR and pokada. Please check the dependencies list,
# Name                    Version                   Build  Channel
gdal                      3.7.3           py311h815a124_5    conda-forge
geojson                   3.1.0                    pypi_0    pypi
geomet                    1.1.0                    pypi_0    pypi
geopy                     2.4.0              pyhd8ed1ab_0    conda-forge
geos                      3.12.0               h59595ed_0    conda-forge
geotiff                   1.7.1               hf074850_14    conda-forge
gettext                   0.21.1               h27087fc_0    conda-forge
matplotlib                3.8.1           py311h38be061_0    conda-forge
netcdf4                   1.6.5           nompi_py311he8ad708_100    conda-forge
numpy                     1.26.0          py311h64a7726_0    conda-forge
python                    3.11.6          hab00c5b_0_cpython    conda-forge
scipy                     1.11.3          py311h64a7726_1    conda-forge
shapely                   2.0.2                    pypi_0    pypi

If you use conda to manage external modules, you can modify defsim_requirements.txt to keep name only. The version information may also be important as some modules can change quicely. SO you can modify, for example, a line as 

gdal=3.7.3

Then run below script to update your python enviroment,

conda install -y -c conda-forge --file <modified_defsim_requirements.txt>

When we release the codes (latest updated on 2023/11/18), we keep the latest versions of required modules. SO no speical requiremnts were found. 

Good luck for your application!
