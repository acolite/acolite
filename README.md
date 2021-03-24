## About ACOLITE
ACOLITE combines the atmospheric correction algorithms for aquatic applications of Landsat and Sentinel-2 developed at RBINS. This repository hosts the (more) generic version of ACOLITE with the aim of bringing together the processing of all different sensors. This new generic version was started 4 February 2021.

It allows simple and fast processing of Landsat (5/7/8) and Sentinel-2 (A/B) images for coastal and inland water applications. Features include generation of RGB images before and after atmospheric correction, atmospheric correction of water bodies and extraction of rectangular regions of interest (defined by bounding coordinates). Level 2 outputs are surface reflectance (ρs=Rrs⋅π) and derived products that can be saved as PNG maps and geolocated datasets in a NCDF (NetCDF) file. The atmospheric correction is image based and needs no external inputs.

The Dark Spectrum Fitting (DSF) algorithm was presented in:
Vanhellemont and Ruddick 2018, [Atmospheric correction of metre-scale optical satellite data for inland and coastal water applications](https://www.sciencedirect.com/science/article/pii/S0034425718303481)
Vanhellemont 2019a, [Adaptation of the dark spectrum fitting atmospheric correction for aquatic applications of the Landsat and Sentinel-2 archives](https://doi.org/10.1016/j.rse.2019.03.010).
Vanhellemont 2019b, [Daily metre-scale mapping of water turbidity using CubeSat imagery.](https://doi.org/10.1364/OE.27.0A1372).

New settings were suggested in:
Vanhellemont 2020, [Sensitivity analysis of the dark spectrum fitting atmospheric correction for metre- and decametre-scale satellite imagery using autonomous hyperspectral radiometry](https://doi.org/10.1364/OE.397456)


ACOLITE development was funded under various projects, a.o. by the Belgian Science Policy Office STEREO program under contracts SR/37/135 (JELLYFOR project) and SR/00/325 (PONDER project), and by the European Community's Seventh Framework Programme (FP7/2007-2013) under grant agreement n° 606797 (HIGHROC project).

**ACOLITE is provided by RBINS as an experimental tool, without explicit or implied warranty. Use of the program is at your own discretion and risk.**

## Distribution
This repository contains all the source code to run ACOLITE in a Python environment. ACOLITE is also distributed as a binary package on the [REMSEM page](http://odnature.naturalsciences.be/remsem/software-and-data/acolite) and supported on the [ACOLITE forum](http://odnature.naturalsciences.be/remsem/acolite-forum/). Currently the previous Landsat/Sentinel-2 version (v20210114.0) is available.

## Dependencies
ACOLITE is coded in Python 3, and requires the following Python packages to run with all functionality:`matplotlib scipy pyproj gdal netcdf4 pyhdf requests statsmodels basemap pillow scikit-image pyresample`

A suitable Python environment can for example be set up using conda and the packages on conda-forge:

            conda create -n acolite -c conda-forge python=3
            conda activate acolite
            conda install -c conda-forge matplotlib scipy pyproj gdal netcdf4 pyhdf requests statsmodels basemap pillow scikit-image pyresample

## Installation
* cd into a suitable directory and clone the git repository: `git clone https://github.com/acolite/acolite`
* cd into the new acolite directory `cd acolite`
* run `python launch_acolite.py`
