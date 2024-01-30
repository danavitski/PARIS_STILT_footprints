# README file for STILT model with added sparse footprint saving multi-core execution functionality
### Id: STILT_sparse_README.md, 29-01-2024 D. Kivits $
---
DISCLAIMER: THIS IS THE README FILE FOR THE STILT MODEL WITH ADDED SPARSE FOOTPRINT SAVING MULTI-CORE EXECUTION FUNCTIONALITY. 
For *only* installing the "regular" version of STILT follow steps 1-6 of the instructions given in the STILT_README.md file. To install STILT *including* the added functionalities, follow this README file! 

This directory contains all R scripts needed to run STILT (the receptor oriented modelling package), as well as added functionality to save the footprint in a sparse format and run the model in parallel on multiple cores. 

- Instructions to install STILT are given under the chapter 'INSTALLING STILT'. 
- A few examples on how to run STILT with the different functionalities are given below under 'RUNNING STILT'.
- A listing of all functions used and what they do is given under 'DIRECTORY LISTING'.
---

## INSTALLING STILT
Fetch the STILT + sparse functionality from <[this repository on Github](https://github.com/DaanKivits/STILT_PARIS.git)> into a user-defined STILT parent directory (on Snellius currenty, this is </projects/0/ctdas/PARIS/transport_models/STILT_Model/>). This can be done by using either git clone, svn checkout, or by downloading the zip file from the repository.

This version of STILT was provided by Santiago Botia and is *not* the most recent version of STILT. A tutorial on how to install the most recent version of STILT can be found on the STILT model development website, which seems to currently be offline: <[www.stilt-model.org]>(www.stilt-model.org). You will need to contact the developers to apply for a user account to access the Trac system, which is the STILT model development website.

## USING STILT + added functionalities
To prepare for a STILT run with the added sparse and multi-core functionality, the followings files are needed (beside the default STILT R files and libraries):

1. a bash submit script, that can communicate with the SLURM scheduler on our HPC cluster (Snellius) and set most of the neccessary STILT parameters before submitting STILT to the cluster. An example of such a file can be found at <batch_scripts/stilt-multistation.sh>. 
    Using this submit script, you can set the following STILT-specific run parameters:
    - the neccessary directories: the parent directory (basedir), the run directory (rundir), the source directory that contains all the STILT-specific R functions (sourcepath), a directory that contains all template files that are later put into the run directory (bdyfiles_dir), and the output directory (path).
    - the station file name (FILENAME) that contains all the station-specific metadata (i.e. station name, station code, latitude, longitude, STILT-corrected elevation), in which each row represents a different station. An example of this file is located over at <stationfiles/stationfile_all.csv>.
    - **the --sparse and --dense options: these options are used to save the footprint in either a dense or sparse format, the former of which means that the influence field is saved for each grid cell (even if the influence is zero); and the latter of which means that only the nonzero values of the simualted influence field are saved. This is done to save disk space, as the footprint can be quite large.**
    - **the option to --filter-times: this option is used to control whether STILT is ran for every hour of the day or only during well-mixed conditions.** For lowland sites (i.e. station with a STILT-corrected elevation below 1000m) this means STILT is ran only between 11 and 16 UTC, but for mountaineous sites (i.e. station with a STILT-corrected elevation above 1000m) STILT is ran only between 23 and 04 UTC to avoid any unwanted stratification effects in the nocturnal boundary layer (that wouldn't represent the 'average' conditions).
    - the option to define --nhrs: this option is used to define the number of hours that STILT is ran backwards in time for. After some tests we found 240 hours to be sufficient, but the group at University of Bristol uses a hybrid scheme with 720 hours (30 days) backwards in time in their RHIME/InTEM inversions.
    - the option to define --npars: this option is used to define the number of particle sets that are released in each STILT simulation timestep. The default value is 250, as this proved to be sufficient after some tests (mostly reduced to save computational time). The group at University of Bristol uses 10000 particles per timestep in their RHIME/InTEM inversions.
    - the option to --overwrite-localization: the option to overwrite .RData localization files if they already exist. This is useful (for instance) if you want to rerun a STILT simulation with the same stochastic and turbulent effects on the particles, but would like to save the footprints in a different format (i.e. dense or sparse). By default, this is set to FALSE.
    - the option to --overwrite-footprints: the option to overwrite the footprint files if they already exist. By default, this is set to FALSE (and footprints are therefore skipped when they already exist, **even when STILT simulation parameters have changed**).
    - **the option to set relevant domain parameters (in the form of lat_ur, lat_ll, lon_ur, lon_ll, lon_res, lat_res, numpix_x, and numpix_y)**. If not given, it will take the default values defined in setStiltParam.r.
    - **the option to run STILT for a specific --station: takes a three-letter stationcode as input, and runs STILT only for that station**. If not given, it will run STILT for all stations in the station file. As currently implemented in <batch_scripts/stilt-multistation.sh>, STILT is ran for one core per station and loops over all stations in the station file. 
    
     **EXPLAIN MORE ABOUT BASH SCRIPT HERE**. This setup script will create the neccessstation-specific) run and output directories, and will occupy each run directory with the setting files neccessary to run STILT (LANDUSE.ASC; ROUGLEN.ASC; ASCDATA.CFG; and runhymodelc.bat), which are created in the setup.sh script (see STILT_README.md). The script copies these files from source directory, which by default is <*parent_dir*/stilt_hysplit/bdyfiles>. 

2. the setStiltParam.r script, which is used to set some neccessary STILT parameters that are not set in the STILT configuration file
   This script sets some important parameters that are not set in the submit script, namely:
   - the name of localization and footprint files created by STILT (by default '.RData' and 'footprint_<station>_<timestamp>.nc' respectively).
   - some other turbulence and mixing layer parameters that I left mostly to their default values, but can be changed if needed.
   - some options to set up a version of STILT that uses an online calculation with externally provided fluxes (e.g. VPRM fluxes for the biosphere or EDGAR fluxes for the anthroposphere). This is not used in the current version of STILT, but might be useful in the future.
  
3. the create_times.r script, which is used to define the time period and time resolution of the STILT run
   This script sets some important parameters that are not set in the submit script, namely:
   - **the start and end date of the STILT run (by default set to run for 2021)**.
   - the hours (in UTC) for which the option --filter-times filters the hours for which STILT is ran. This can be set separately for lowland and mountaineous sites, as described above.
   - the name of the file that contains the time period and time resolution of the STILT run (by default '.RDataTimes')
  