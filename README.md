## About ACOLITE
ACOLITE combines the atmospheric correction algorithms for aquatic applications of Landsat and Sentinel-2 developed at RBINS. This repository hosts the (more) generic version of ACOLITE with the aim of bringing together the processing of all different sensors. This new generic version was started 4 February 2021.

**On 21 April 2021 this code was released as public beta. Please contact Quinten via email/ACOLITE Forum/GitHub if you find any issues. The settings files are largely compatible with the previous version, but it is recommended to create a new one.**

ACOLITE allows simple and fast processing of imagery from various satellites, including Landsat (5/7/8) and Sentinel-2/MSI (A/B), PlanetScope and RapidEye, Venµs, SPOT and Pléiades, WorldView-2 and -3, and Sentinel-3/OLCI (A/B) for coastal and inland water applications. The Dark Spectrum Fitting atmospheric correction algorithm works especially well for turbid and productive waters, but can also be applied over clear waters and land with reasonable success.

Features include generation of RGB images before and after atmospheric correction, atmospheric correction of water bodies and extraction of regions of interest (defined by bounding coordinates or a polygon file). Level 2 outputs are surface reflectance (ρs=Rrs⋅π) and derived products that are saved as geolocated datasets in a NetCDF file, and can be exported as PNG maps. The atmospheric correction is image based and needs no external inputs.

ACOLITE development was funded under various projects by the Belgian Science Policy Office STEREO program under contracts SR/37/135 (JELLYFOR project) and SR/00/325 (PONDER project), and by the European Community's Seventh Framework Programme (FP7/2007-2013) under grant agreement n° 606797 (HIGHROC project). TACT development was funded by the Belgian Science Policy Office BRAIN-be program under contract BR/165/A1/MICROBIAN.

**ACOLITE is provided by RBINS as an experimental tool, without explicit or implied warranty. Use of the program is at your own discretion and risk.**

## References
The Dark Spectrum Fitting (DSF) algorithm was presented in:

* Vanhellemont and Ruddick 2018, [Atmospheric correction of metre-scale optical satellite data for inland and coastal water applications](https://www.sciencedirect.com/science/article/pii/S0034425718303481)

* Vanhellemont 2019a, [Adaptation of the dark spectrum fitting atmospheric correction for aquatic applications of the Landsat and Sentinel-2 archives](https://doi.org/10.1016/j.rse.2019.03.010)

* Vanhellemont 2019b, [Daily metre-scale mapping of water turbidity using CubeSat imagery.](https://doi.org/10.1364/OE.27.0A1372)

New settings were suggested in:

* Vanhellemont 2020c, [Sensitivity analysis of the dark spectrum fitting atmospheric correction for metre- and decametre-scale satellite imagery using autonomous hyperspectral radiometry](https://doi.org/10.1364/OE.397456)

The adaptation to Sentinel-3/OLCI was presented in:

* Vanhellemont and Ruddick 2021, [Atmospheric correction of Sentinel-3/OLCI data for mapping of suspended particulate matter and chlorophyll-a concentration in Belgian turbid coastal waters](https://doi.org/10.1016/j.rse.2021.112284)

The Thermal Atmospheric Correction Tool (TACT) is now integrated in ACOLITE and allows for processing of Landsat thermal band data to surface temperatures. TACT was presented in:

* Vanhellemont 2020a, [Automated water surface temperature retrieval from Landsat 8/TIRS](https://doi.org/10.1016/j.rse.2019.111518)

* Vanhellemont 2020b, [Combined land surface emissivity and temperature estimation from Landsat 8 OLI and TIRS](https://doi.org/10.1016/j.isprsjprs.2020.06.007)

## Distribution
This repository contains all the source code to run ACOLITE in a Python environment. ACOLITE is also distributed as a binary package on the [REMSEM page](http://odnature.naturalsciences.be/remsem/software-and-data/acolite) and supported on the [ACOLITE forum](http://odnature.naturalsciences.be/remsem/acolite-forum/). Currently only the previous Landsat/Sentinel-2 version (v20210114.0) is available as binary on the REMSEM page.

## Dependencies
ACOLITE is coded in Python 3, and requires the following Python packages to run with all functionality:`numpy matplotlib scipy gdal pyproj astropy cartopy scikit-image pyhdf pyresample`

A suitable Python environment can for example be set up using conda and the packages on conda-forge:

            conda create -n acolite -c conda-forge python=3
            conda activate acolite
            conda install -c conda-forge numpy matplotlib scipy gdal pyproj astropy cartopy scikit-image pyhdf pyresample

## Installation
* cd into a suitable directory and clone the git repository: `git clone https://github.com/acolite/acolite`
* cd into the new acolite directory `cd acolite`
* run `python launch_acolite.py`

## Ancillary and SRTM DEM data download
ACOLITE can automatically retrieve ancillary data (ozone, water vapour, pressure and wind speed) from the servers of the Ocean Biology Processing Group (OBPG) at NASA (https://oceandata.sci.gsfc.nasa.gov) and SRTM DEM data from NASA's Land Processes Distributed Active Archive Center (LP DAAC) at the USGS Earth Resources Observation and
Science (EROS) Center (https://e4ftl01.cr.usgs.gov/MEASURES/) for both resources an EarthData account is requred: https://urs.earthdata.nasa.gov/users/new

ACOLITE tries and retrieve the EARTHDATA_u and EARTHDATA_p (for username and password) from the system environment variables. Alternatively they can be specified in the processing settings file or the general ACOLITE config file at config/config.txt.

## TACT installation
TACT is now integrated into ACOLITE, and allows for processing of the Landsat thermal data to surface temperature. TACT needs libRadtran to be installed to perform simulations of the atmospheric down and upwelling radiances and transmittance: http://libradtran.org/doku.php

TACT needs the user to have an account at the Research Data Archive (RDA) at the University Corporation for Atmospheric Research (UCAR) to retrieve atmospheric profile data: https://rda.ucar.edu/

## TACT configuration
* Edit the configuration file `nano config/config.txt` and add the full path to the libRadtran directory on your system to the libradtran_dir= setting
* Edit your .netrc file to add your RDA UCAR credentials: `nano $HOME/.netrc`, with $l and $p your login and password respectively for the RDA:

            machine rda.ucar.edu
            login $l
            password $p

* Edit your .dodsrc file to point to your .netrc file: `nano $HOME/.dodsrc`. Write the full path explicitly:
            HTTP.NETRC=/path/to/.netrc
