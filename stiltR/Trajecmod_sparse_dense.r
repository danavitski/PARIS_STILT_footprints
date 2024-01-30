#***************************************************************************************************
# Function that loops over all starting times
#***************************************************************************************************

Trajecmod_sparse_dense <- function(partarg=NULL, totpartarg=NULL, nodeoffset=NULL, csv=TRUE, 
                           sparse = TRUE, dense= FALSE, significance = 15,
                           calc_sum = FALSE, ens_mem_num=NULL, npars=NULL, station=NULL, overwrite_fp = FALSE,
                           overwrite = TRUE, rundir = NULL) {
   #---------------------------------------------------------------------------------------------------
   # Calls 'Trajec' for each starting time
   # arguments assigned from call to setStiltparam.r
   # 'csv'             controls whether footprints are stored in NetCDF (.nc) or comma-seperated (.csv) format
   # 'significance'     controls the amount of significant digits that are stored in the footprints
   # 'sparse'           controls whether sparse footprints are created
   # 'dense'            controls whether dense footprints are created
   # 'significance'     controls the significance of the stored influence, affecting the file size in the process
   # 'calc_sum'         controls whether the sum of the footprint should be automatically calculated and 
   #                    output as a separate .csv file
   # 'ens_mem_num'      controls the ensemble member number in a sensitivity run, which is used to create a unique 
   #                    name for the output files
   # 'npars'            is used solely for naming convention of the output footprints and RData files, and is read
   #                    as input from 
   # 'station'          controls the station name, which is used to create a unique name for the output files
   # 'overwrite_fp'     controls whether footprints are overwritten if they already exist
   # 'overwrite'        controls whether trajectory files are overwritten if they already exist
   # 'rundir'           controls the directory in which STILT is ran
   #---------------------------------------------------------------------------------------------------
   
   if (is.null(station)){
      stop("station is NULL, please specify a station name")
   }
   
   if(is.null(rundir)){rundir<-"~/STILT_Exe/"}
   rundir<-paste(rundir,station,"/",sep="")

   # Make sure all influence values are stored with a scientific notation
   options(scipen = 0)

   # Make sure that footprints are saved if footprintTF is set to TRUE, but both sparse and dense are set to FALSE
   if (sparse == FALSE & dense == FALSE & footprintTF == TRUE){
      print("Warning: footprintTF is set to TRUE, indicating the need to save footprints, but both sparse and dense are set to FALSE. Setting sparse to TRUE.")
      sparse = TRUE
   }
   
   # need to assign parameters; also save parameter setting in archive file with date in name
   savename <- gsub(" ", ".", date())
   savename <- substring(savename,4)
   runs.done.dir <- NULL
   if (file.exists('./Runs.done')) runs.done.dir <- './Runs.done/'
   if (is.null(runs.done.dir) && file.exists(paste(sourcepath,'Runs.done',sep='')))
     runs.done.dir <- paste(sourcepath,'/Runs.done/',sep='')
   if (is.null(runs.done.dir) && substring(path, nchar(path)-nchar("Runs.done"), nchar(path)) == "Runs.done/")
     runs.done.dir <- sourcepath
   if (!is.null(runs.done.dir)) 
   {
      file.copy("setStiltparam.r",paste(runs.done.dir, "setStiltparam", savename, ".r", sep=""), overwrite=T)
      cat("Saving copy of setStiltparam.r in ", paste(runs.done.dir, "setStiltparam", savename, ".r", sep=""), "\n",sep="")
   } 
   else 
   {
      cat("No Runs.done directory; parameter information not saved\n")
   }
   
   totpart <- 1
   if (!is.null(totpartarg)) 
   {
      cat('resetting totpart=', totpart, ' to totpartarg=', totpartarg, '\n', sep='')
      totpart <- totpartarg
   }
   part <- 1
   if (!is.null(partarg)) 
   {
      cat('Using totpart=', totpart, ' resetting part=', part, ' to partarg=', partarg, '\n', sep='')
      part <- partarg
   }
   if (!is.null(nodeoffset)) 
   {
      nummodel <- part+nodeoffset
      cat('Using nodeoffset= ', nodeoffset, ' ,results in nummmodel= ', nummodel, '\n', sep='')
   }
   else 
   {
      nummodel <- part
   }

   # get Starting info
   if (!exists("path_timesname")) {path_timesname<-path} #new tk
   if (!existsr(Timesname, path_timesname)) stop(paste("cannot find object ", Timesname, " in directory ", path_timesname, sep=""))
    cat(paste(Timesname, sep=""), path_timesname,"\n")
   StartInfo <- getr(paste(Timesname, sep=""), path_timesname) # object containing fractional julian day, lat, lon, agl for starting position and time
   # SELECTION OF A FEW Receptors for testing!
   if (Times.startrow > 0) StartInfo <- StartInfo[Times.startrow:Times.endrow,, drop=FALSE] # can be just one (Times.startrow=Times.endrow)
   
   # divide job into "totpart" parts to speed up
   if (dim(StartInfo)[1] < totpart) 
   {
      cat ('Warning: resetting totpart=', totpart, ' to dim(StartInfo)[1]=', dim(StartInfo)[1], '\n', sep='')
      totpart <- dim(StartInfo)[1]
   }
   if (part > totpart) 
   {
      stop.message <- paste('Specified part=', part, ' > totpart=', totpart, ', stopping\n')
      cat(stop.message)
      stop(stop.message)
   }
   
   if(totpart>1)
   {
      n.perpart<-rep(dim(StartInfo)[1]%/%totpart,totpart)
      if((dim(StartInfo)[1]%%totpart)>0)n.perpart[1:(dim(StartInfo)[1]%%totpart)]<-n.perpart[1:(dim(StartInfo)[1]%%totpart)]+1
      start.rows=c(1,cumsum(n.perpart)+1)
      StartInfo <- StartInfo[start.rows[part]:(start.rows[part+1]-1),, drop=FALSE]
   }
   dimnames(StartInfo)[[2]] <- toupper(dimnames(StartInfo)[[2]])
   
   if (biomassburnTF) {
      biomassburnoutmat <- matrix(nrow=dim(StartInfo)[[1]], ncol=2)
      dimnames(biomassburnoutmat) <- list(NULL, c("ident", "CO"))
   }
   
   # OVERWRITE WARNING
   if (!exists("path_stiltresult")) {path_stiltresult<-path} #new tk
   if(existsr(paste("stiltresult",part,sep=""),path=path_stiltresult)) 
   {
      warning("You are attempting to overwrite an existing stiltresult object")
      warning("Notice: If you have changed parameters and Trajecmod fails, first try to move or remove the existing stiltresult object")
   }
   
   nrows     <- length(StartInfo[,1]) # 1 row for each trajectory
   rownum    <- 1
   firsttraj <- T
   firstflux <- T
   l.remove.Trajfile <- FALSE
   if (exists('remove.Trajfile')) l.remove.Trajfile <- remove.Trajfile
   for (j in 1:nrows) 
   {
      ###############################################
      ##### run trajectories and save output ########
      ###############################################
      lat       <- StartInfo[j, "LAT"] 
      lon       <- StartInfo[j, "LON"]
      agl       <- StartInfo[j, "AGL"]
      identname <- pos2id(StartInfo[j,1], lat, lon, agl, npars=npars, ens_mem_num=ens_mem_num, sep="x")

      # check if trajectory file already exists
      if (length(list.files(path = path, pattern = identname)) == 0 | (length(list.files(path = path, pattern = identname)) > 0 & overwrite_fp==TRUE)){
      cat(paste("Trajecmod(): Trajectory file ",identname," does not exist yet \n"))
      cat("Trajecmod(): ", identname, " running at ", date(), "\n", sep="")
      
      dat       <- month.day.year(floor(StartInfo[j,1])) # from julian to mmddyy
      yr4       <- dat$year # 4 digit year
      yr        <- yr4%%100 # 2 digit year (or 1 digit...)
      mon       <- dat$month
      day       <- dat$day
      hr        <- round((StartInfo[j,1]-floor(StartInfo[j,1]))*24)

      l.ziscale <- NULL
      if (exists('ziscale')) l.ziscale <- ziscale
      l.zsg.name <- NULL
      if (exists('zsg.name')) l.zsg.name <- zsg.name
      l.create.X0 <- FALSE
      if (exists('create.X0')) l.create.X0 <- create.X0
      l.use.multi <- TRUE
      if (exists('use.multi')) l.use.multi <- use.multi
      l.hymodelc.exe <- NULL
      if (exists('hymodelc.exe')) l.hymodelc.exe <- hymodelc.exe
      if (l.use.multi) 
      {
         l.setup.list <- list()
         if (exists('setup.list')) l.setup.list <- setup.list
         info <- Trajecmulti(yr=yr, mon=mon, day=day, hr=hr, lat=lat, lon=lon, agl=agl, nhrs=nhrs,
                        numpar=nparstilt, doublefiles=T, metd=metsource, metlib=metpath,
                        conv=convect, overwrite=overwrite, outpath=path, varsout=varstrajec, rundir=rundir,
                        nummodel=nummodel, sourcepath=sourcepath, ziscale=l.ziscale, zsg.name=l.zsg.name,
                        create.X0=l.create.X0,setup.list=l.setup.list,hymodelc.exe=l.hymodelc.exe,
                        siguverr=siguverr,TLuverr=TLuverr,zcoruverr=zcoruverr,horcoruverr=horcoruverr,
                        sigzierr=sigzierr,TLzierr=TLzierr,horcorzierr=horcorzierr,outname=identname)
      }
      else 
      {
         info <- Trajec(yr=yr, mon=mon, day=day, hr=hr, lat=lat, lon=lon, agl=agl, nhrs=nhrs,
                        delt=stepsize, numpar=nparstilt, doublefiles=T, metd=metsource, metlib=metpath,
                        conv=convect, overwrite=overwrite, outpath=path, varsout=varstrajec, rundir=rundir,
                        nummodel=nummodel, sourcepath=sourcepath, ziscale=l.ziscale, zsg.name=l.zsg.name,
                        create.X0=l.create.X0,
                        siguverr=siguverr,TLuverr=TLuverr,zcoruverr=zcoruverr,horcoruverr=horcoruverr,
                        sigzierr=sigzierr,TLzierr=TLzierr,horcorzierr=horcorzierr, outname=identname)
      } # if l.use.multi
      if (firsttraj) 
      {   # set up array for run info
         run.info <- matrix(NA, nrow=nrows, ncol=length(info))
         dimnames(run.info) <- list(NULL, names(info))
         firsttraj <- F
      } 
      else 
      {
         havenames <- dimnames(run.info)[[2]]
         for (nm in names(info)) 
         {
            ndx <- which(havenames == nm)[1] # check if there are new column names
            if (is.na(ndx)) 
            { # new column name, need to add column
               run.info <- cbind(run.info, rep(NA, dim(run.info)[1]))
               dimnames(run.info)[[2]][dim(run.info)[2]] <- nm
            }
            else 
            {
            havenames[ndx] <- paste(havenames[ndx], "done", sep="")
            } # dummy
         }
      }
      run.info[j, names(info)] <- info
   
      #########################################################################
      ###### TM3-STILT ################
      if(writeBinary == T)
      {
         ####### variables for writing to binary files ######
         dlat   <-as.numeric(lat,digits=10)
         dlon   <-as.numeric(lon,digits=10)
         dagl   <-as.numeric(agl,digits=10)
         dlatres<-as.numeric(lat.res,digits=10)
         dlonres<-as.numeric(lon.res,digits=10)
         
         ####### construct filename for binary file #########
         cyr<-as.character(2000+yr)
         cmon<-as.character(mon)
         cday<-as.character(day)
         chr<-as.character(hr)
         x1<-""; x2<-""; x3<-""
         if(mon<10) x1<-paste(x1,"0",sep="")
         if(day<10) x2<-paste(x2,"0",sep="")
         if(hr<10)  x3<-paste(x3,"0",sep="")
         
         pathBinFootprintstation<-paste(pathBinFootprint,station,"/",sep="")
         if (file.access(pathBinFootprintstation,0)!=0) 
         {
            system(paste("mkdir ",pathBinFootprintstation,sep=""))
         }
         
         filename <- paste(pathBinFootprintstation,station,"_",as.character(-nhrs),"h0",as.character(ftintr),
                         "h_",cyr,x1,cmon,x2,cday,x3,chr,"_",gridtag,cendian,".d",sep="")
   #     Jan's version with height in filename
   #     filename <- paste(pathBinFootprintstation,station,"_",as.character(-nhrs),"h0",as.character(ftintr),
   #                    "h_",cyr,x1,cmon,x2,cday,x3,chr,"_",sprintf("%5.5d",as.integer(agl)),"_",gridtag,cendian,".d",sep="")
   
         if (file.exists(filename)) 
         {
            print(paste("Binary footprint file ",filename," already exists"))
            print(paste("not replaced !!"))
         }
         else
         {
            ident<-info["outname"]
            print(paste(path,".RData",ident,sep=""))
            if (file.exists(paste(path,".RData",ident,sep=""))) 
            {
               #for longer than hourly intervals for footprints, first make sure to match time intervals of flux fields
               #assume those are e.g. 0-3, 3-6 etc. UTC, or 0-24 UTC 
               #NOTE: only for hourly intervals variables "foottimes", "nfoottimes" and "nftpix" are computed in setStiltparam.r 
               if(ftintr>1)
               {
                  nfoottimes       <- -nhrs/ftintr+2                  #number of footprints computed
                  foottimes        <- rep(c(0),nfoottimes)            #vector of times (backtimes) in hours between which footprint is computed
                  nftpix           <- rep(c(0),nfoottimes)            #vector of numbers of pixels in each footprint
                  for(ft in 2:nfoottimes)
                  { 
                     foottimes[ft] <- hr+(ft-2)*ftintr 
                  }
                  foottimes[nfoottimes]<- -nhrs
                  if(hr==0)
                  {                                 #special case when starting at midnight
                     foottimes     <-foottimes[2:nfoottimes]
                     nftpix        <-nftpix[2:nfoottimes]
                     nfoottimes    <-nfoottimes-1
                  }
               }
   
               ####### call Trajecfoot ############################ 
               ident               <- info["outname"]
               foot                <- Trajecfoot(ident,      pathname=path, foottimes=foottimes, zlim=c(zbot,ztop),fluxweighting=NULL, coarse=1, vegpath=vegpath,
                                            numpix.x=numpix.x, numpix.y=numpix.y,lon.ll=lon.ll, lat.ll=lat.ll, lon.res=lon.res, lat.res=lat.res, rundir=rundir)
   #           foot                <- Trajecfoot(ident=ident,pathname=path, foottimes=foottimes, zlim=c(zbot,ztop),fluxweighting=NULL, coarse=1, vegpath=vegpath,
   #                                        numpix.x=numpix.x, numpix.y=numpix.y,lon.ll=lon.ll, lat.ll=lat.ll, lon.res=lon.res, lat.res=lat.res, rundir=rundir)
               nameslat            <- rownames(foot)
               nameslon            <- colnames(foot)
   #           print(paste("nameslat, nameslon: ",nameslat, nameslon))
               if(is.null(foot))
               {
                  print(paste("is.null(foot): ",is.null(foot),foot))
                  print(paste("No binary footprint file for TM3 written!!!!"))
               }
               else
               {
   #              #### write the output file ####
   #              cyr<-as.character(2000+yr)
   #              cmon<-as.character(mon)
   #              cday<-as.character(day)
   #              chr<-as.character(hr)
   #              x1<-""; x2<-""; x3<-""
   #              if(mon<10) x1<-paste(x1,"0",sep="")
   #              if(day<10) x2<-paste(x2,"0",sep="")
   #              if(hr<10)  x3<-paste(x3,"0",sep="")
   #             
   #              pathBinFootprintstation<-paste(pathBinFootprint,station,"/",sep="")
   #              if (file.access(pathBinFootprintstation,0)!=0) {
   #               system(paste("mkdir ",pathBinFootprintstation,sep=""))
   #              }
   #             
   ##             filename<-paste(pathBinFootprintstation,station,"_",as.character(-nhrs),"h0",as.character(ftintr),
   ##                              "h_",cyr,x1,cmon,x2,cday,x3,chr,"_",gridtag,cendian,".d",sep="")
   ##             Jan's version with height in filename
   #              filename<-paste(pathBinFootprintstation,station,"_",as.character(-nhrs),"h0",as.character(ftintr),
   #                              "h_",cyr,x1,cmon,x2,cday,x3,chr,"_",sprintf("%5.5d",as.integer(agl)),"_",gridtag,cendian,".d",sep="")
                 
                  print(paste("Writing binary footprint file: ",filename))
                  cat("binary footprint file: ", filename, "\n")
                  con    <- file(filename, "wb")
   #              write first: hours backward, hours interval, lat, lon, altitude, lat resolution, and lon resolution
                  writeBin(as.numeric(c(nhrs,ftintr,lat,lon,agl,lat.res,lon.res)),con,size=4,endian=endian) 
                  
                  nftpix <- apply(foot,c(3),function(x)return(sum(x>0)))
                  ftmax  <- which(nftpix>0)[sum(nftpix>0)]
                  for(ft in 1:(nfoottimes-1))
                  {
                     writeBin(as.numeric(c(ft,nftpix[ft])),con,size=4,endian=endian)          # index of footprint and # grids in footprint
                     #### pixel by pixel #####
                     if(nftpix[ft]>0)
                     {
                        id        <- which((foot[,,ft])>0,arr.ind=T)
                        foot_data <- cbind(as.numeric(nameslat[id[,1]]),as.numeric(nameslon[id[,2]]),foot[,,ft][id])# save coordinates of pixels of each footprint and footprint value
                        foot_data <- foot_data[order(foot_data[,1],foot_data[,2]),]
                        if(nftpix[ft]==1)foot_data<-matrix(foot_data,nrow=1) #convert back to matrix; ordering turned one-row matrices into vectors
                        writeBin(as.vector(t(foot_data)),con,size=4,endian=endian)
                        if(any(foot_data[,3]<0))browser()
                     } #if(nftpix[ft]>0)
                  } #for(ft   
                  print(paste("ftmax: ",ftmax))
           
                  close(con)
               } #if is.null(foot) 
           
            }
            else
            {
               print(paste(path,outname," does exist -> no new STILT run !!",sep=""))
            } #if (file.exists(paste(path,outname,sep=""))) 
   
         } #if (file.exists(filename))
   
   
      }  #end if(writeBinary == T)
   
      ###### end of TM3-STILT output ################
      #########################################################################
   
      #########################################################################
      ##### map trajectories to flux grids and vegetation maps ################
      ##### calculate mixing ratios at receptor points, save in result ########
      if (fluxTF) 
      {
         print(paste("Trajecmod(): rownumber j:", j))
       
       
         traj <- Trajecvprm(ident=identname, pathname=path, tracers=fluxtracers, coarse=aggregation,
                    dmassTF=T, nhrs=nhrs, vegpath=vegpath, evilswipath=evilswipath,
                    vprmconstantspath=vprmconstantspath, vprmconstantsname=vprmconstantsname, nldaspath=nldaspath,
                    nldasrad=usenldasrad, nldastemp=usenldastemp, pre2004=pre2004,
                    keepevimaps=keepevimaps, detailsTF=detailsTF, bios=fluxmod, landcov=landcov,
                    numpix.x=numpix.x, numpix.y=numpix.y, lon.ll=lon.ll, lat.ll=lat.ll,
                    lon.res=lon.res, lat.res=lat.res)
       
       
         # 'traj' is a vector
         if (existsr(paste("stiltresult", part, sep=""), path=path_stiltresult)) 
         {
            result <- getr(paste("stiltresult", part, sep=""), path=path_stiltresult)
            if (dim(result)[1] != nrows) 
            {
               if (firstflux) print("Trajecmod(): existing stiltresult has wrong dimension; creating new one.")
            } 
            else 
            {
               if (firstflux) print("Trajecmod(): found existing stiltresult, update rows in that.")
               firstflux <- FALSE
            }
         }
         if (firstflux) 
         { # at beginning create result object
            ncols <- length(traj) # all from Trajec(), + 3 from StartInfo (agl, lat, lon)
            result <- matrix(NA, nrow=nrows, ncol=ncols)
            firstflux <- F
         }
         result[rownum, ] <- traj
         dimnames(result) <- list(NULL, c(names(traj)))
         dimnames(result) <- list(NULL, dimnames(result)[[2]])
         # write the object into default database; object names are, e.g., "Crystal.1"
         assignr(paste("stiltresult", part, sep=""), result, path=path_stiltresult)
      }
      rownum <- rownum+1
     
      ##### calculate footprint, assign in object ########
      if (footprintTF)
      {
         print(" ")
         print(paste("Trajecmod(): ", identname, " running footprint at ", unix("date"), sep=""))
         #print(paste("Trajecmod(): memory in use:", memory.size()[1]))

         library(ncdf4)

         if (sparse == TRUE){
         ###############################
         ########## SPARSE #############
         ###############################
         foot <- Trajecfoot(identname, pathname=path, foottimes=foottimes, zlim=c(zbot, ztop),
                            fluxweighting=NULL, coarse=1, vegpath=vegpath,
                            numpix.x=numpix.x, numpix.y=numpix.y,
                            lon.ll=lon.ll, lat.ll=lat.ll, lon.res=lon.res, lat.res=lat.res, 
                            landcov="IGBP",wrfinput=NULL, sparse=TRUE, rundir=rundir)
         if (length(foot$infl) < 3){
            print(foot$infl)}
         if (!is.null(foot)){      
         # Check length of footprint
         #print(paste("length(foot$infl)=",length(foot$infl),sep=""))
         #print(paste("length(foot$times)=",length(foot$times),sep=""))

         #index_x = foot$index_x
         #index_y = foot$index_y
         lons     = foot$lons
         lats     = foot$lats
         infl     = foot$infl
         times    = foot$times

         # Construct longx,latx from info supplied in setStiltparam.r
         longx_1          <- seq(lon.ll,length.out=numpix.x,by=lon.res)
         latx_1           <- seq(lat.ll,length.out=numpix.y,by=lat.res)
         require(pracma) # For meshgrid
         dum            <- meshgrid(latx_1,longx_1)
         longx          <- dum$Y
         latx           <- dum$X
         nx             <- dim(latx)  [1]
         ny             <- dim(longx) [2]
         nt             <- length(foottimes)-1
         nonzerolen     <- length(foot$infl)

         if (csv == TRUE){
         # Write Output to newly created .csv file
         outfilename    <- paste(path, "footprint_", station, '_', identname, ".csv", sep="")
         
         # Construct sparse footprint
         df = data.frame(unlist(times),unlist(infl),unlist(lons),unlist(lats))
         names(df) = c("Time","Influence","Longitude","Latitude")

         # Write to output to .csv files
         if (!file.exists(outfilename)){
            write.csv(df,file=outfilename, row.names = FALSE)}
         else if (file.exists(outfilename) & overwrite_fp == TRUE){
            write.csv(df,file=outfilename, row.names = FALSE)}
         else {print(paste(outfilename, "already exists, skipping..."))}}
         
         # Define nc file global attributes for sparse format
         attrs <- list(
            summary = 'Surface influence fields (footprints) stored in sparse format (only the non-zero values)',
            model = 'STILT (v1.2), source code provided by Thomas Koch from MPI BGC',
            institution = 'Wageningen University, department of Meteorology and Air Quality, Wageningen, the Netherlands; Rijksuniversiteit Groningen, Groningen, the Netherlands; ICOS Carbon Portal, Lund, Sweden',
            contact = 'Daan Kivits; daan.kivits@wur.nl',
            conventions = 'CF-1.8',
            creation_date = format(Sys.time(), "%Y-%m-%d %H:%M"),
            crs = 'spherical earth with radius of 6370 km',
            disclaimer = 'This data belongs to the CarbonTracker project',
            history = paste('File created on', format(Sys.time(), "%Y-%m-%d %H:%M"),
                              'by dkivits, using the code on  the Subversion (SVN) repository on https://projects.bgc-jena.mpg.de/STILT/svn.
                              R version', R.version$version.string),
            creator = 'Daan Kivits, https://orcid.org/0009-0005-8856-8497',
            frequency = '1h',
            length = invisible(capture.output(paste(str(nhrs), 'h'))),
            particles = invisible(capture.output(paste(str(npars), 'particles'))),
            grid_definition = 'Sparse format of CTE-HR grid, which ranges from -15E to 35W and 33N to 72N, with a resolution of 0.1x0.2 degree',
            #longitude_pixels = paste0(str(nx)),
            #latitude_pixels = paste0(str(ny)),
            geospatial_lat_resolution = '0.1 degree',
            geospatial_lon_resolution = '0.2 degree',
            keywords = 'footprint, surface influence field, STILT',
            license = 'CC-BY-4.0',
            nominal_resolution = '0.1x0.2 degree'
         )

         # Save sparse format to .nc file
         outfilename    <- paste(path, "footprint_", station, '_', identname, ".nc", sep="")       
         dimnonzero     <- ncdim_def( "Nonzero_len", "", longname = "Amount of non-zero values", 1:nonzerolen, create_dimvar=TRUE, unlim=TRUE)
         nc_foot        <- ncvar_def( "Influence", "ppm ug ms^-2 s^-1)", longname = "Simulated influence", dimnonzero, prec="float", compression = 6)
         nc_times       <- ncvar_def( "Time", "hr", longname = "Simulation timestep", dimnonzero, compression = 6)
         #nc_index_x     <- ncvar_def( "index_x","", longname = "Longitude index of nonzero value", dimnonzero, prec="short", compression = 6)
         #nc_index_y     <- ncvar_def( "index_y","", longname = "Latitude index of nonzero value", dimnonzero, prec="short", compression = 6)
         nc_lons        <- ncvar_def( "Longitude","degrees E", longname = "Longitude of nonzero value", dimnonzero, compression = 6)
         nc_lats        <- ncvar_def( "Latitude","degrees N", longname = "Latitude of nonzero value", dimnonzero, compression = 6)
         
         if (!file.exists(outfilename) | (file.exists(outfilename) & overwrite_fp == TRUE)){
            nc <- nc_create(outfilename, list(nc_foot, nc_times, nc_lons, nc_lats))

            ncvar_put(nc = nc, varid = nc_foot, vals = infl, start=c(1), count=c(-1))
            ncvar_put(nc = nc, varid = nc_times, vals = times, start=c(1), count=c(-1))
            ncvar_put(nc = nc, varid = nc_lons, vals = lons, start=c(1), count=c(-1))
            ncvar_put(nc = nc, varid = nc_lats, vals = lats, start=c(1), count=c(-1))
            
            for (attr in names(attrs)){
               ncatt_put(nc = nc, varid = 0, attname=attr, attval=attrs[[attr]])}
            
            nc_close(nc)
            
         } else {print(paste(outfilename, "already exists, skipping..."))}
         
         if (calc_sum == TRUE){
         # Also write sum of each simulated time of the footprints to .csv file
         sumfilename      <- paste(path, 'sum_footprint_', station, '_', identname, ".csv", sep="")

         # Construct sparse footprint
         df = data.frame(unlist(times),unlist(infl),unlist(lons),unlist(lats))
         names(df) = c("Time","Influence","Longitude","Latitude")
         
         # Construct sum of each simulated footprint
         sumdf = aggregate(df$Influence ~ df$Time, data = df, FUN = 'sum')
         names(sumdf) = c("Time","Sum")
         
         # Write to output to .csv files
         if (!file.exists(sumfilename) | (file.exists(sumfilename) & overwrite_fp == TRUE)){
            write.csv(sumdf,file=sumfilename, row.names = FALSE)}
         else {print(paste(sumfilename, "already exists, skipping..."))}}
         
         } else {
            print(paste("Trajecmod(): Saving footprint info of empty influence field ...",sep=""))

            pos<-id2pos(identname)
            time<-month.day.year(floor(pos[1]))

            colnamesdf = c("ident", "yr", "mon", "day", "hr", "lat", "lon", "agl")
            fp_misfile <- paste(rundir,"missing.footprints",sep="")

            if(!file.exists(fp_misfile)){
               df <- data.frame(matrix(ncol = 8, nrow = 0, dimnames=list(NULL,colnamesdf)))
               write.csv(df, fp_misfile, row.names=FALSE)
            }

            data <- data.frame(ident = identname, yr = time$year, mon = time$month, day = time$day, hr = round((pos[1]-floor(pos[1]))*24), lat = lat, lon = lon, agl = agl)
            write.table(data, fp_misfile, append = T, col.names=FALSE, sep= ",", row.names=FALSE)

            # Create empty .nc footprint file
            outfilename    <- paste(path, "footprint_", station, '_', identname, ".nc", sep="")
            if (!file.exists(outfilename) | (file.exists(outfilename) & overwrite_fp == TRUE)){
               file.create(outfilename)
            }

            next
         }

         } else if (dense == TRUE){
         ###############################
         ########## DENSE ##############
         ###############################
         foot <- Trajecfoot(identname, pathname=path, foottimes=foottimes, zlim=c(zbot, ztop),
                            fluxweighting=NULL, coarse=1, vegpath=vegpath,
                            numpix.x=numpix.x, numpix.y=numpix.y,
                            lon.ll=lon.ll, lat.ll=lat.ll, lon.res=lon.res, lat.res=lat.res, 
                            landcov="IGBP",wrfinput=NULL, sparse=FALSE, rundir=rundir)

         if (!is.null(foot)){
         # Define nc file global attributes for sparse format
         attrs <- list(
            summary = 'Surface influence fields (footprints) stored in dense format (2D arrays over time)',
            model = 'STILT (v1.2), source code provided by Thomas Koch from MPI BGC',
            institution = 'Wageningen University, department of Meteorology and Air Quality, Wageningen, the Netherlands; Rijksuniversiteit Groningen, Groningen, the Netherlands; ICOS Carbon Portal, Lund, Sweden',
            contact = 'Daan Kivits; daan.kivits@wur.nl',
            conventions = 'CF-1.8',
            creation_date = format(Sys.time(), "%Y-%m-%d %H:%M"),
            crs = 'spherical earth with radius of 6370 km',
            disclaimer = 'This data belongs to the CarbonTracker project',
            history = paste('File created on', format(Sys.time(), "%Y-%m-%d %H:%M"),
                              'by dkivits, using the code on  the Subversion (SVN) repository on https://projects.bgc-jena.mpg.de/STILT/svn.
                              R version', R.version$version.string),
            creator = 'Daan Kivits, https://orcid.org/0009-0005-8856-8497',
            frequency = '1h',
            length = invisible(capture.output(paste(str(nhrs), 'h'))),
            particles = invisible(capture.output(paste(str(npars), 'particles'))),
            grid_definition = 'Sparse format of CTE-HR grid, which ranges from -15E to 35W and 33N to 72N, with a resolution of 0.1x0.2 degree',
            #longitude_pixels = paste0(str(nx)),
            #latitude_pixels = paste0(str(ny)),
            geospatial_lat_resolution = '0.1 degree',
            geospatial_lon_resolution = '0.2 degree',
            keywords = 'footprint, surface influence field, STILT',
            license = 'CC-BY-4.0',
            nominal_resolution = '0.1x0.2 degree'
         )   

         # Construct longx,latx from info supplied in setStiltparam.r
         longx_1          <- seq(lon.ll,length.out=numpix.x,by=lon.res)
         latx_1           <- seq(lat.ll,length.out=numpix.y,by=lat.res)
         require(pracma) # For meshgrid
         dum            <- meshgrid(latx_1,longx_1)
         longx          <- dum$Y
         latx           <- dum$X
         nx             <- dim(latx)  [1]
         ny             <- dim(longx) [2]
         nt             <- length(foottimes)-1

         if (calc_sum == TRUE){
         # Construct sum of each simulated footprint
         sumfilename      <- paste(path, 'sum_footprint_', station, '_', identname, "_dense.csv", sep="")
         summat = c()

         for (i in 1:dim(foot)[3])
         {
            sum = sum(foot[,,i], digit = significance)
            summat = append(summat,sum)
         }
            
         # Construct sum of each simulated footprint
         sumdf = data.frame(1:nt,summat)
         names(sumdf) = c("Time","Sum")
         
         # Write to output to .csv files
         if (!file.exists(sumfilename) | (file.exists(sumfilename) & overwrite_fp == TRUE)){
            write.csv(sumdf,file=sumfilename, row.names = FALSE)}
         else {print(paste(sumfilename, "already exists, skipping..."))}
         }

         # Save dense format to .nc file
         outfilename    <- paste(path, "footprint_", station, '_', identname, "_dense.nc", sep="")
         mv             <- -9999   # Fill Value
         dimt           <- ncdim_def( "Time","hr", longname = "Simulation timestep", 1:nt, create_dimvar=TRUE, unlim=TRUE)
         dimx           <- ncdim_def( "Longitude", longname = "Longitude of nonzero value", "degrees E", longx[,1])
         dimy           <- ncdim_def( "Latitude" , longname = "Latitude of nonzero value", "degrees N", latx[1,] )
         nc_foot        <- ncvar_def( "Influence","ppm/(ug/m2/s)", longname = "Influence field", list(dimx,dimy,dimt), mv, prec="float")
         
         if (!file.exists(outfilename) | (file.exists(outfilename) & overwrite_fp == TRUE)){
            nc <- nc_create(outfilename, list(nc_foot))
            ncvar_put(nc, nc_foot, aperm(foot, c(2,1,3)), start=c(1, 1, 1),count=c(-1, -1, -1))
            
            for (attr in names(attrs)){
               ncatt_put(nc = nc, varid = 0, attname=attr, attval=attrs[[attr]])}
               
            nc_close(nc)
            
         } else {print(paste(outfilename, "already exists, skipping..."))}
         } else {
            print(paste("Trajecmod(): Saving footprint info of empty influence field ...",sep=""))

            pos<-id2pos(identname)
            time<-month.day.year(floor(pos[1]))

            colnamesdf = c("ident", "yr", "mon", "day", "hr", "lat", "lon", "agl")
            fp_misfile <- paste(rundir,"missing.footprints",sep="")

            if(!file.exists(fp_misfile)){
               df <- data.frame(matrix(ncol = 8, nrow = 0, dimnames=list(NULL,colnamesdf)))
               write.csv(df, fp_misfile, row.names=FALSE)}

            data <- data.frame(ident = identname, yr = time$year, mon = time$month, day = time$day, hr = round((pos[1]-floor(pos[1]))*24), lat = lat, lon = lon, agl = agl)
            write.table(data, fp_misfile, append = T, col.names=FALSE, sep= ",", row.names=FALSE)

            # Create empty .nc footprint file
            outfilename    <- paste(path, "footprint_", station, '_', identname, ".nc", sep="")
            if (!file.exists(outfilename) | (file.exists(outfilename) & overwrite_fp == TRUE)){
               file.create(outfilename)}

            next
         }}}
         
      } else {
      print(paste("Trajectory file ",identname," already exists"))
      print(paste("Skipping trajectory!"))
      next
      }
      

      ##### plot footprint ########
      if (footplotTF) 
      {  # plot footprints
         foot <- getr(paste(identname, sep=""), path)
         footplot(foot,identname,lon.ll,lat.ll,lon.res,lat.res)
         for (foottimespos in 1:(length(foottimes)-1)) 
         {
            ############# NOT READY YET ###############
         }
      }
     
      # Specify the function parameters
      if (biomassburnTF) 
      {
         biomassburnoutmat[j, ] <- biomassburn(timesname=StartInfo, burnpath=burnpath, endpath=path, pathname=path, nhrs=nhrs, timesrow=j)
         print(paste("Biomassburning influence calculated to be ", biomassburnoutmat[j,2], " ppbv. Inserted into fireinfluence matrix row ", j, sep=""))
      }
     
      if(l.remove.Trajfile)unix(paste("rm -f ",paste(path,".RData",identname,sep=""),sep=""))
   
   } 

   # Wrap up all of the CO biomassburning calculations
   if (biomassburnTF)
     write.table(biomassburnoutmat, file=paste(path, "fireinfluencex", nhrs, "hr_", part, ".txt", sep=""), row.names=F)
   
   ##### save mixing ratios at receptor points in file, e.g. stiltresult1.csv for part=1 ########
   if (fluxTF) 
   {
      dimnames(result) <- list(NULL, dimnames(result)[[2]])
      # write the object into default database; object names are, e.g., "Crystal.1"
    
      assignr(paste("stiltresult", part, sep=""), result, path=path_stiltresult)
      print(paste("stiltresult", part, " assigned in ", path_stiltresult, sep=""))
      write.table(result, file=paste(path_stiltresult, "stiltresult", part, ".csv", sep=""), na="", row.names=F)
   }
   
   # If evi and lswi maps from vprm calculations is saved to the global environment; it should be removed here
   
   rm(list=objects(pattern="GlobalEvi"), envir=globalenv())
   rm(list=objects(pattern="GlobalLswi"), envir=globalenv())
   gc(verbose=F)
   
   if (!exists("run.info")){
      cat("All footprints are already calculated. Exiting loop. \n") 
      return(NULL)
   } else {
      return(run.info)}
}
