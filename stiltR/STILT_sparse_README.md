# README file for STILT model with added sparse footprint saving multi-core execution functionality
### Id: STILT_sparse_README.md, 29-01-2024 D. Kivits $
---
DISCLAIMER: THIS IS THE README FILE FOR THE STILT MODEL WITH ADDED SPARSE FOOTPRINT SAVING MULTI-CORE EXECUTION FUNCTIONALITY.
This directory contains all R scripts needed to run STILT (the receptor oriented modelling package), as well as added functionality to save the footprint in a sparse format and run the model in parallel on multiple cores. 

- Instructions to install STILT are given under the chapter 'INSTALLING STILT'. 
- A few examples on how to run STILT with the different functionalities are given below under 'RUNNING STILT'.
- A listing of all functions used and what they do is given under 'DIRECTORY LISTING'.
---

## INSTALLING STILT
For installing the "regular" version of STILT follow steps 1-6 of the instructions given in the STILT_README.md file. That includes running the setup.sh script in the parent directory, as the setup_auto.sh script discussed in this STILT_sparse_README.md is only valid for the added functionality of sparse footprint saving and multi-core execution.

1. To prepare for a STILT run with the added sparse and multi-core functionality, run the setup_auto.sh script in the parent directory. This setup script will create the neccessary (station-specific) run and output directories, and will occupy each run directory with the setting files neccessary to run STILT (LANDUSE.ASC; ROUGLEN.ASC; ASCDATA.CFG; and runhymodelc.bat), which are created in the setup.sh script (see STILT_README.md). The script copies these files from source directory, which by default is <*parent_dir*/stilt_hysplit/bdyfiles>. 
2. <Link to Github branch where added STILT functionality is stored (everything that is added onto regular STILT)>