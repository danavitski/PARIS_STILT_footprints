# STILT model with added sparse footprint saving multi-core execution functionality
### Id: README.md, 29-01-2024 D. Kivits $
---
DISCLAIMER: THIS IS THE README FILE FOR THE STILT MODEL WITH ADDED SPARSE FOOTPRINT SAVING MULTI-CORE EXECUTION FUNCTIONALITY. 
For *only* installing the "regular" version of STILT follow steps 1-6 of the instructions given in the <stiltR/README.md> file. To install STILT *including* the added functionalities, follow this README file! 

This directory contains all R scripts needed to run STILT (the receptor oriented modelling package), as well as added functionality to save the footprint in a sparse format and run the model in parallel on multiple cores. The <stiltR/> directory contains all STILT-related R functions, the <batch_scripts/> directory contains all bash scripts needed to submit single- or multi-core STILT jobs on the HPC cluster, and the <stationfiles/> directory contains an exemplary file with station-specific metadata (i.e. station name, station code, latitude, longitude, STILT-corrected elevation) used in PARIS. The <merged_stilt_hysplit/> directory contains the libraries 

- Instructions to install STILT are given under the section 'INSTALLING STILT'. 
- A few examples on how to run STILT with the different functionalities are given below under the 'USING STILT + added functionalities' section.
- A quick guide into STILT's key functionalities are given in the section 'SOME KEY POINTS IN UNDERSTANDING THE STILT CODE'.
---

## INSTALLING STILT
Fetch the STILT + sparse functionality from <[this repository on Github](https://github.com/DaanKivits/STILT_PARIS.git)> into a user-defined STILT parent directory (on Snellius currenty, this is </projects/0/ctdas/PARIS/transport_models/STILT_Model/>). This can be done by using either git clone, svn checkout, or by downloading the zip file from the repository.

This version of STILT was provided by Santiago Botia and is *not* the most recent version of STILT. A tutorial on how to install the most recent version of STILT can be found on the STILT model development website, which seems to currently be offline: <[www.stilt-model.org]>(www.stilt-model.org). You will need to contact the developers to apply for a user account to access the Trac system, which is the STILT model development website.


## SOME KEY POINTS IN UNDERSTANDING THE STILT CODE
The main STILT functionality can be reduced to a few key scripts:
- <stiltR/Trajec.r>: Function to run HYSPLIT particle dispersion model and to check distribution of particles in the model domain. This function is called by the Trajecmod.r function for each timestep of the simulation.
- <stiltR/Trajecmod_sparse_dense.r>: this function calls Trajec.r for each backwards simulation timestep (as defined by create_times.r), which is later fed into the Trajecfoot.r function to save as a footprint. 
- <stiltR/Trajecfoot_sparse.r> & <stiltR/create_sparse_footprints_coords.r>: these functions are used to create the sparse footprint files. The first script filters the particles in the object from Trajec.r that are outside the model domain and returns the object as 2D dense array, and using the second script it can convert these 2D arrays to a sparse format if requested.


## USING STILT + added functionalities
To prepare for a STILT run with the added sparse and multi-core functionality, the following input data is required:

1. Meteorological driver fields, provided in <.arl> format. For now, these are prepared by Thomas Koch of the MPI BG group in Jena as the conversion from GRIB to ARL proved to be more difficult than expected. The exact format of these meteorological driver fields can vary, depending on the source. A few of the (most popular) options that are coded into STILT are the following:
   - **ECMWF (ECmetF) forecast data**: 3 hour time resolution, 72-144h forecast, 25km (~.25 degree resolution) grid
   - **Aladin (alad) forecast data**: Aladin mesoscale forecasts (MeteoFrance), 3 hour time resolution, 72h forecast, 8km grid
   - **Weather Research & Forecasting (wrf) forecast data**: WRF nested domains, variable time and horizontal resolution and forecast duration.
   - ... more options are available, have a look at the <stiltR/setStiltparam.r> script to see which options are available and how to add new ones.
2. A station file, which contains all the station-specific metadata (i.e. station name, station code, latitude, longitude, STILT-corrected elevation), in which each row represents a different station. An example of this file is located over at <stationfiles/stationfile_all.csv>.
3. STILT configuration files (i.e. <LANDUSE.ASC>; <ROUGLEN.ASC>; <ASCDATA.CFG>; and <runhymodelc.bat>), which are created in the <setup.sh> script that comes with the "default" installment of STILT (see <STILT_README.md>). The script copies these files from source directory (which by default is <stilt_hysplit/bdyfiles>) to the (station-specific) run directory.

Besides this input data, the followings files are needed (beside the default STILT R files and libraries):

1. a bash submit script, that can communicate with the SLURM scheduler on our HPC cluster (Snellius) and set most of the neccessary STILT parameters before submitting STILT to the cluster. An example of such a file can be found at <batch_scripts/stilt-multistation.sh>. 
    Using this submit script, you can set the following STILT-specific run parameters:
    - the neccessary directories: the parent directory (basedir), the run directory (rundir), the source directory that contains all the STILT-specific R functions (sourcepath), a directory that contains all template files that are later put into the run directory (bdyfiles_dir), a directory where the neccessary meteorology files are stored (metdir), and the output directory (path).
    - **the --sparse and --dense options: these options are used to save the footprint in either a dense or sparse format, the former of which means that the influence field is saved for each grid cell (even if the influence is zero); and the latter of which means that only the nonzero values of the simualted influence field are saved. This is done to save disk space, as the footprint can be quite large.**
    - **the option to set relevant domain parameters (in the form of lat_ur, lat_ll, lon_ur, lon_ll, lon_res, lat_res, numpix_x, and numpix_y)**. If not given, it will take the default values defined in setStiltParam.r.
    - **the option to run STILT for a specific --station: takes a three-letter stationcode as input, and runs STILT only for that station**. If not given, it will run STILT for all stations in the station file. 
    - **the option to --filter-times: this option is used to control whether STILT is ran for every hour of the day or only during well-mixed conditions.** For lowland sites (i.e. station with a STILT-corrected elevation below 1000m) this means STILT is ran only between 11 and 16 UTC, but for mountaineous sites (i.e. station with a STILT-corrected elevation above 1000m) STILT is ran only between 23 and 04 UTC to avoid any unwanted stratification effects in the nocturnal boundary layer (that wouldn't represent the 'average' conditions).
    - the station file name (FILENAME) that contains all the station-specific metadata (i.e. station name, station code, latitude, longitude, STILT-corrected elevation), in which each row represents a different station. An example of this file is located over at <stationfiles/stationfile_all.csv>.
    - the option to define --nhrs: this option is used to define the number of hours that STILT is ran backwards in time for. After some tests we found 240 hours to be sufficient, but the group at University of Bristol uses a hybrid scheme with 720 hours (30 days) backwards in time in their RHIME/InTEM inversions.
    - the option to define --npars: this option is used to define the number of particle sets that are released in each STILT simulation timestep. The default value is 250, as this proved to be sufficient after some tests (mostly reduced to save computational time). The group at University of Bristol uses 10000 particles per timestep in their RHIME/InTEM inversions.
    - the option to --overwrite-localization: the option to overwrite .RData localization files if they already exist. This is useful (for instance) if you want to rerun a STILT simulation with the same stochastic and turbulent effects on the particles, but would like to save the footprints in a different format (i.e. dense or sparse). By default, this is set to FALSE.
    - the option to --overwrite-footprints: the option to overwrite the footprint files if they already exist. By default, this is set to FALSE (and footprints are therefore skipped when they already exist, **even when STILT simulation parameters have changed**).
    - the option to --calc-sum: this option controls whether the influence fields are summed over all timesteps of each STILT simulation, and therefore whether single-timestep footprints are returned or not.
    - the option to define --ens-mem-num: this option controls how many times a STILT simulation is ran for (i.e. how many ensemble members are used). This is useful if you want to quantify the uncertainty in the STILT simulations caused by the stochastic and turbulent effects on the particles (i.e. the atmospheric transport itself). This also appends the footprint output name with the current ensemble member number, so that the footprints are not overwritten. By default, this option is switched off.
    
   The submit script loops over each station (row) in the user-defined stationfile and runs the <setup_multi.sh> script for each of these stations, which creates the neccessary (station-specific) run and output directories, and will occupy each run directory with the configuration files neccessary to run STILT (i.e. <LANDUSE.ASC>; <ROUGLEN.ASC>; <ASCDATA.CFG>; and <runhymodelc.bat>), which are created in the <setup.sh> script that comes with the "default" installment of STILT (see <STILT_README.md>). The script copies these files from source directory, which by default is <stilt_hysplit/bdyfiles>.
   
   After the neccessary directories have been initialized and the STILT-specific run parameters have been set, the submit script will run STILT for each station in the stationfile. Before the job is taken out of the SLURM queue, the simulations for each station should be finished first.


2. the setStiltParam.r script, which is used to set some neccessary STILT parameters that are not set in the STILT configuration file
   This script sets some important parameters that are not set in the submit script, namely:
   - the type of meteorological input data fed into STILT (by default 'ECmetF' for ECMWF forecast data)
   - the name of localization and footprint files created by STILT (by default '.RData' and 'footprint_{station}_{timestamp}.nc' respectively).
   - some other turbulence and mixing layer parameters that I left mostly to their default values, but can be changed if needed.
   - some options to set up a version of STILT that uses an online calculation with externally provided fluxes (e.g. VPRM fluxes for the biosphere or EDGAR fluxes for the anthroposphere). This is not used in the current version of STILT, but might be useful in the future.
  

3. the create_times.r script, which is used to define the time period and time resolution of the STILT run
   This script sets some important parameters that are not set in the submit script, namely:
   - **the start and end date of the STILT run (by default set to run for 2021)**.
   - the hours (in UTC) for which the option --filter-times filters the hours for which STILT is ran. This can be set separately for lowland and mountaineous sites, as described above.
   - the name of the file that contains the time period and time resolution of the STILT run (by default '.RDataTimes')


---

### **If some things related to STILT are still unclear, have a look at the stiltR/README.r document that describes the functionality of the "default" STILT installation in more detail.**
