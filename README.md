## About ACOLITE
ACOLITE combines the atmospheric correction algorithms for aquatic applications of Landsat and Sentinel-2 developed at RBINS. This repository hosts the (more) generic version of ACOLITE with the aim of bringing together the processing of all different sensors. This new generic version was started 4 February 2021, and was released to the public on 21 April 2021. Binary releases from 20210802.0 onward are based on this codebase. Please contact Quinten via email/ACOLITE Forum/GitHub if you find any issues. The settings files are largely compatible with previous versions, but it is recommended to create a new settings file configuring only the settings you want to change.

ACOLITE allows simple and fast processing of imagery from various satellites, including Landsat (5/7/8/9) and Sentinel-2/MSI (A/B), PlanetScope (Dove and SuperDove) and RapidEye, Venµs, SPOT and Pléiades, QuickBird2, WorldView-2 and -3, Sentinel-3/OLCI (A/B) and ENVISAT/MERIS, as well as several hyperspectral sensors (CHRIS, HYPERION, HICO, PRISMA, DESIS) for coastal and inland water applications. The Dark Spectrum Fitting atmospheric correction algorithm works especially well for turbid and productive waters, but can also be applied over clear waters and land with reasonable success.

Features include generation of RGB images before and after atmospheric correction, atmospheric correction of water bodies and extraction of regions of interest (defined by bounding coordinates or a polygon file). Level 2 outputs are surface level reflectance (ρs=Rrs⋅π) and derived products that are saved as geolocated datasets in a NetCDF file, and can be exported as PNG maps. The atmospheric correction is image based and needs no external inputs.

To reduce initial download volume, no LUTs are provided in this repository. (LUTs are still in the history of the repository so it is advised to clone with a low depth, e.g. using --depth 1, to avoid pulling those in.) Required LUTs will be retrieved automatically from a separate repository (https://github.com/acolite/acolite_luts) the first time they are needed.

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

TACT performance for Antarctic mountain sites and near shore waters was evaluated in:

* Vanhellemont et al. 2021a, [Towards physical habitat characterisation in the Antarctic Sør Rondane Mountains using satellite remote sensing](https://doi.org/10.1016/j.rsase.2021.100529)

* Vanhellemont et al. 2022, [Validation of Landsat 8 high resolution Sea Surface Temperature using surfers](https://doi.org/10.1016/j.ecss.2021.107650)

## Distribution
This repository contains all the source code to run ACOLITE in a Python environment. ACOLITE is also distributed as a binary package on the [releases page](https://github.com/acolite/acolite/releases) and is supported on the [ACOLITE forum](http://odnature.naturalsciences.be/remsem/acolite-forum/).

## Dependencies
ACOLITE is coded in Python 3, and requires the following Python packages to run: `numpy matplotlib scipy gdal pyproj scikit-image pyhdf pyresample netcdf4 h5py requests pygrib astropy cartopy`

A suitable Python environment can for example be set up using conda and the packages on conda-forge:

            conda create -n acolite -c conda-forge python=3
            conda activate acolite
            conda install -c conda-forge numpy matplotlib scipy gdal pyproj scikit-image pyhdf pyresample netcdf4 h5py requests pygrib astropy cartopy

## Installation
* cd into a suitable directory and clone the git repository: `git clone --depth 1 https://github.com/acolite/acolite`
* cd into the new acolite directory `cd acolite`
* run `python launch_acolite.py`

## Ancillary and SRTM DEM data download
ACOLITE can automatically retrieve ancillary data (ozone, water vapour, pressure and wind speed) from the servers of the Ocean Biology Processing Group (OBPG) at NASA (https://oceandata.sci.gsfc.nasa.gov) and SRTM DEM data from NASA's Land Processes Distributed Active Archive Center (LP DAAC) at the USGS Earth Resources Observation and
Science (EROS) Center (https://e4ftl01.cr.usgs.gov/MEASURES/) for both resources an EarthData account is required: https://urs.earthdata.nasa.gov/users/new

Your EarthData account needs to have approval for OB.DAAC Data Access and LP DAAC Data Pool to access the ancillary and SRTM DEM datasets. Login to your EarthData account to check Authorised Apps, and click Approve More Applications if necessary: https://urs.earthdata.nasa.gov/profile

ACOLITE tries and retrieve the EARTHDATA_u and EARTHDATA_p (for username and password) from the system environment variables. Alternatively they can be specified in the processing settings file or the general ACOLITE config file at config/config.txt.

## TACT installation
TACT is now integrated into ACOLITE, and allows for processing of the Landsat thermal data to surface temperature. TACT needs libRadtran to be installed to perform simulations of the atmospheric down and upwelling radiances and transmittance: http://libradtran.org/doku.php

By default ACOLITE expects the libRadtran (v2.0.2) to be in external/libRadtran-2.0.2 directory, but this can be changed by editing the configuration file `nano config/config.txt` and providing the full path to the libRadtran directory at a different location to the libradtran_dir= setting. To install libRadtran in the external directory (when in the main acolite directory):

            mkdir external
            cd external
            wget http://www.libradtran.org/download/history/libRadtran-2.0.2.tar.gz
            gzip -d libRadtran-2.0.2.tar.gz
            tar -xvf libRadtran-2.0.2.tar
            cd libRadtran-2.0.2
            ./configure
            make
            make check

## TACT configuration
**2021-11-22 QV I believe this next section is no longer required**

TACT needs the user to have an account at the Research Data Archive (RDA) at the University Corporation for Atmospheric Research (UCAR) to retrieve atmospheric profile data: https://rda.ucar.edu/

* Edit your .netrc file to add your RDA UCAR credentials: `nano $HOME/.netrc`, with $l and $p your login and password respectively for the RDA:

            machine rda.ucar.edu
            login $l
            password $p

* Edit your .dodsrc file to point to your .netrc file: `nano $HOME/.dodsrc`. Write the full path explicitly:
            HTTP.NETRC=/path/to/.netrc
