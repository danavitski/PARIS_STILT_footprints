#!/usr/bin/env Rscript

# Script to run STILT:
## edit setStiltparam.r to match the domain and folders
## edit create_times.r to match the receptor period

## set the following
library(proj4)
sourcepath <- "/home/dkivits/STILT/STILT_Model/stiltR/"
source("sourceall.r")
source("create_times.r")
source("stilt.r")
