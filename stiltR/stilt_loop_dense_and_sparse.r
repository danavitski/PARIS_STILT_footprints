# This script is based on the original stilt.r atmospheric transport model pipeline but 
# has some added functionality, such as the ability to run STILT in multi-core mode and 
# the optional creation of sparse footprints.

# Original script by C. Gerbig at MPI Jena (2010), modified by D. Kivits (2024)"""

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description='Run STILT for a single station or for all stations in a file')

# Define arguments to parse later
parser$add_argument('--station', help='Station name or file with station names')
parser$add_argument('--sparse', action='store_true', help='Save sparse footprints')
parser$add_argument('--dense', action='store_true', help='Save dense footprints')
parser$add_argument('--csv', action='store_true', help='Store output in csv file instead of netCDF4')
parser$add_argument('--filter-times', action='store_true', help='Filter times between 11 and 16 local time')
parser$add_argument('--calc-sum', action='store_true', help='Calculate sum of footprints')
parser$add_argument('--ens-members', help='Number ensemble members to run STILT with. This argument is used to run stilt.r with the same settings multiple times, and can be used to see the effect of the randomized error in STILT')
parser$add_argument('--npars', default = 250, help='Number of particles to run STILT with')
parser$add_argument('--nhrs', default = 240, help='How many hours STILT is ran backwards in time for. The model acknowledges a positive number as hours backwards in time. The default is 240 hours backwards in time, which is 10 days')
parser$add_argument('--overwrite-localization', default=FALSE, action='store_true', help='Switch to overwrite the .RData files that store the time and location of the particles released by STILT')
parser$add_argument('--overwrite-footprints', default=FALSE, action='store_true', help='Switch to overwrite footprints if they already exist')
parser$add_argument('--path', help='The path to the directory where the footprints are stored')
parser$add_argument('--rundir', help='The path to the directory where the STILT model is ran')
parser$add_argument('--altdir', help='Optional extra (project-specific) path that creates a subdir in which the STILT model is ran')
parser$add_argument('--metdir', help='The path to the directory where the meteorological data is stored')
parser$add_argument('--sourcepath', help='The path to the directory where the STILT scripts reside')
parser$add_argument('--stationfile', help='The path to the stationfile that stores all the station metadata (including 3D location)')
parser$add_argument('--lon_ur', type='double', help='The upper right longitude of the domain')
parser$add_argument('--lat_ur', type='double', help='The upper right latitude of the domain')
parser$add_argument('--lon_ll', type='double', help='The lower left longitude of the domain')
parser$add_argument('--lat_ll', type='double', help='The lower left latitude of the domain')
parser$add_argument('--lon_res', type='double', help='The longitude resolution of the domain')
parser$add_argument('--lat_res', type='double', help='The latitude resolution of the domain')
parser$add_argument('--numpix_y', type='integer', help='The amount of pixels in the y direction of the domain')
parser$add_argument('--numpix_x', type='integer', help='The amount of pixels in the x direction of the domain')

# Parse arguments
args = parser$parse_args()

# Store parsed arguments in args
sparse = args$sparse
dense = args$dense
csv = args$csv
filter_times = args$filter_times
calc_sum = args$calc_sum
ens_members = args$ens_members
npars = args$npars
nhrs = args$nhrs
path = args$path
rundir = args$rundir
metpath = args$metdir
altdir = args$altdir
sourcepath = args$sourcepath
overwrite_fp = args$overwrite_footprints
overwrite_lf = args$overwrite_localization
station = args$station

# Define grid
lon.ur = args$lon_ur
lat.ur = args$lat_ur
lon.ll = args$lon_ll
lat.ll = args$lat_ll
lon.res = args$lon_res
lat.res = args$lat_res
numpix.y = args$numpix_y
numpix.x = args$numpix_x

# Read in stationfile to get station metadata
stationpath <- args$stationfile
stationfile <- read.csv(stationpath, sep =',')

# Give path to file that sources all required R functionality
if (! exists("sourcepath")) {
   if (file.exists("stiltR"))
      sourcepath <- paste(getwd(), "/stiltR/", sep="")
   else if (file.exists("Rsc"))
      sourcepath <- paste(getwd(), "/Rsc/", sep="")
   else {
      stop('stilt.r: no stiltR or Rsc directory found.')
      quit(status=1)
   }
}

# Check if given sourcepath is valid. If yes, source all R function files
if (!file.exists(paste(sourcepath, "sourceall.r", sep=""))){
   stop('stilt.r: "sourcepath" is not a valid source path')
} else {
   cat("stilt.r: using sourcepath", sourcepath, "\n")
   source('sourceall.r')
   sourceall(sparse=sparse, dense=dense)}

# Define the runs.done.dir to store run.info and missing.footprints
runs.done.dir <- NULL
if (file.exists('./Runs.done')){runs.done.dir <- paste('./Runs.done/', station, '/', sep='')}
if (is.null(runs.done.dir) && file.exists(paste(sourcepath,'Runs.done/', sep=''))){
  runs.done.dir <- paste(sourcepath,'Runs.done/',station, '/', sep='')}
  if (!dir.exists(runs.done.dir)){dir.create(runs.done.dir)}

# Define file that stores missing footprint metadata and save a copy of it to runs.done.dir if it exists
fp_misfile <- paste(rundir,station,"/missing.footprints",sep="")
savename <- gsub(" ",".",date())
savename <- substring(savename,4)
  if (file.exists(fp_misfile)){
   
   # Copy missing.footprints file from rundir to Runs.done before starting new run!
   new_fp_misfile <- paste(runs.done.dir,"missing.footprints",savename,sep="")
   file.copy(fp_misfile, new_fp_misfile, overwrite = TRUE)
   cat("stilt.r: Saved a copy of missing.footprints file to ", runs.done.dir, "\n")}

################################
##### MULTI-CORE RUN ###########
################################
# Running STILT in this multi-core mode means that one station is ran per CPU core. The user has to define and give
# a station name, the metadata of which it will search in the user-defined stationfile. Whereas in single-core mode STILT would run all stations
# in a given stationfile on a single core, in multi-core mode STILT will divide each station over multiple cores as SLURM tasks. This has to be set up using an external 
# SLURM batch script. For both modes, the user can also define the amount of ensemble members to run STILT with. This can be useful to do a statistical 
# analysis of STILT, by running the transport model with the same settings multiple times. If the --ens-members argument is not given, STILT is ran with 
# the default amount of ensemble members, which is 1.

if (!is.null(args$station)){
cat("stilt.r: starting.. Station = ", station,"..\n")

if (!is.null(ens_members)){
for (ens_mem_num in seq(1,ens_members,1)){
cat("stilt.r: starting ensemble member = ", ens_mem_num,"..\n")

# Create a list of timestamps for which STILT is ran, depending on station location and filter_times argument
create_times(stationfile = stationfile, station = station, filter_times = filter_times, outpath = path)

# Create STILT job partitioning info to use in tracjemod() func later
partinfo <- Sys.getenv(c("STILT_PART", "STILT_TOTPART","STILT_OFFSET"), unset = NA)
if ( (is.na(partinfo[[3]])) && (!is.na(partinfo[[1]])) && (!is.na(partinfo[[2]] )) ) {
    partinfo[[3]] <-0 }

if (any(is.na(partinfo))) {
   partarg    <- NULL
   totpartarg <- NULL
   nodeoffset <- NULL
} else {
   partarg    <- as.integer(partinfo[[1]])
   totpartarg <- as.integer(partinfo[[2]])
   nodeoffset <- as.integer(partinfo[[3]])}

# Call main STILT function, store run info
run.info <- Trajecmod_sparse_dense(partarg=partarg, totpartarg=totpartarg, nodeoffset=nodeoffset, 
            csv=csv, sparse = sparse, dense = dense, calc_sum = calc_sum, ens_mem_num = ens_mem_num,
            npars = nparstilt, station=station, overwrite_fp=overwrite_fp, overwrite_lf=overwrite_lf, rundir=rundir)
            
}} else if (is.null(ens_members)){

# Create a list of timestamps for which STILT is ran, depending on station location and filter_times argument
create_times(stationfile = stationfile, station = station, filter_times = filter_times, outpath = path)

# Create STILT job partitioning info to use in tracjemod() func later
partinfo <- Sys.getenv(c("STILT_PART", "STILT_TOTPART","STILT_OFFSET"), unset = NA)
if ( (is.na(partinfo[[3]])) && (!is.na(partinfo[[1]])) && (!is.na(partinfo[[2]] )) ) {
    partinfo[[3]] <-0 }

if (any(is.na(partinfo))) {
   partarg    <- NULL
   totpartarg <- NULL
   nodeoffset <- NULL
} else {
   partarg    <- as.integer(partinfo[[1]])
   totpartarg <- as.integer(partinfo[[2]])
   nodeoffset <- as.integer(partinfo[[3]])}

# Call main STILT function, store run info
run.info <- Trajecmod_sparse_dense(partarg=partarg, totpartarg=totpartarg, nodeoffset=nodeoffset, 
            csv=csv, sparse = sparse, dense = dense, calc_sum = calc_sum, 
            npars = nparstilt, station=station, overwrite_fp=overwrite_fp, overwrite_lf=overwrite_lf, rundir=rundir)
}
} else {

##################################
##### SINGLE-CORE RUN ############
##################################
# Running STILT in single-core mode means that all stations in a given stationfile are ran on a single CPU core. This can be useful if you want to run STILT for a large
# amount of stations but don't have access to a cluster, or in case you want to test STILT for a small time period but on multiple stations. The user has to define and give
# a stationfile, which contains the metadata of all stations that the user wants to run STILT for. 

for (i in 1:nrow(stationfile)){
station <- stationfile$code[i]

cat("stilt.r: starting.. Station = ", station,"..\n")

# Create a list of timestamps for which STILT is ran, depending on station location and filter_times argument
create_times(stationfile = stationfile, station = station, filter_times = filter_times, outpath = path)

# Create STILT job partitioning info to use in tracjemod() func later
partinfo <- Sys.getenv(c("STILT_PART", "STILT_TOTPART","STILT_OFFSET"), unset = NA)
if ( (is.na(partinfo[[3]])) && (!is.na(partinfo[[1]])) && (!is.na(partinfo[[2]] )) ) {
    partinfo[[3]] <-0 }

if (any(is.na(partinfo))) {
   partarg    <- NULL
   totpartarg <- NULL
   nodeoffset <- NULL
} else {
   partarg    <- as.integer(partinfo[[1]])
   totpartarg <- as.integer(partinfo[[2]])
   nodeoffset <- as.integer(partinfo[[3]])}

# Call main STILT function, store run info
run.info <- Trajecmod_sparse_dense(partarg=partarg, totpartarg=totpartarg, nodeoffset=nodeoffset, 
            csv=csv, sparse = sparse, dense = dense, calc_sum = calc_sum, 
            npars = nparstilt, station=station, overwrite_fp=overwrite_fp, overwrite_lf=overwrite_lf, rundir=rundir)
}
}

# Save run information to Runs.done directory with timestamp as savename handle
if (!is.null(runs.done.dir)) {
  assignr(paste("run.info",savename,sep=""),run.info,runs.done.dir,printTF=T)
  cat("stilt.r: Saved run.info to ", runs.done.dir, "\n")
} else {
  cat("stilt.r: Runs.done not found in ./ or sourcepath, not saving backups of run.info and missing.footprints\n")}
