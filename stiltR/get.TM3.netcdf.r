#***************************************************************************************************
#  read initial values from CarbonTracker file and attach them to the "result" matrix
#***************************************************************************************************

get.TM3.netcdf <- function(yr4=NULL, mon=NULL, day=NULL, hr=NULL, tracersinifile=NULL,
                             result=NULL, result.sel=NULL,tracer="co2") {

# --------------------------------------------------------------------------------------------------
#  Interface
#  =========
#
#  yr4          YYYY of the receptor time
#  mon          month "
#  day          day   "
#  hr           hour  "
#  tracersinifile   netCDF file name (absolute path) of the TM3 analyzed fields
#  result       matrix with columns "btime", "lat", "lon" and "pres"
#               gets the column "co2ini" attached/overwritten on output
#  result.sel   selector of same length as the columns of "result": which rows should be considered
#
# --------------------------------------------------------------------------------------------------
#  $Id: get.TM3.netcdf.r,v 1.1 2009-02-24 14:27:54 gerbig Exp $
# --------------------------------------------------------------------------------------------------


   if (!is.element("package:ncdf", search())) library("ncdf") # Load ncdf library
   refdate <- month.day.year(floor(julian(mon, day, yr4) + hr/24 - result[result.sel, "btime"]/24))

   #check which files are there
   tracersinifiles<-dir(dirname(tracersinifile),full.names=T)
   years.avail<-as.numeric(substring(tracersinifiles,nchar(tracersinifiles)-6,nchar(tracersinifiles)-3))
   years.avail<-years.avail[!is.na(years.avail)]
   for(yeari in unique(refdate$year)){
      selm<-refdate$year==yeari
      resultm<-result[selm,]
      result.selm<-result.sel[selm]
      if(yeari>max(years.avail)){#limit to last available file
         print(paste("NOTICE: data in ",dirname(tracersinifile)," only until ",max(years.avail),"; replaced more recent years by ",max(years.avail),sep=""))
         yeari<-max(years.avail)
      }
      tracersinifilem<-paste(substring(tracersinifile,1,nchar(tracersinifile)-7),yeari,".nc",sep="")
      tm3file <- open.ncdf(tracersinifilem, write=F,readunlim=FALSE)
      centerlats <- get.var.ncdf(tm3file, varid="lat")           # Center
      centerlons <- get.var.ncdf(tm3file, varid="lon")           # Center
      times <- get.var.ncdf(tm3file, varid="time")           # Center (instantaneous)
      levs <- get.var.ncdf(tm3file, varid="height")           # levels
      nlevs<-length(levs)
      if(tracer=="co2")ini <- att.get.ncdf(tm3file, varid="co2mix",attname="ini")$value           # global CO2 offset
      if(tracer!="co2")ini <- 0
      # Get Ending positions
      agl4ct <- resultm[result.selm, "agl"]+resultm[result.selm, "grdht"]
      lat4ct <- resultm[result.selm, "lat"]
      lon4ct <- resultm[result.selm, "lon"]

      #-- lat lon pointer
      dlat<-unique(round(diff(centerlats),4))
      dlon<-unique(round(diff(centerlons),4))
      latpt<-round((lat4ct-centerlats[1])/dlat)+1
      lonpt<-round((lon4ct-centerlons[1])/dlon)+1
      lonpt[lonpt > as.vector(360/dlon)]<-1
      point<-cbind(lonpt,latpt) #spatial pointer (2D)

      #-- time pointer
      delt.h<-(times[2]-times[1])/3600
      timept<-round((julian(mon, day, yr4) + hr/24 - resultm[result.selm, "btime"]/24-julian(1, 1, yr4))*24/delt.h+1+0.01)#+0.01 to avoid round bug
      timept[timept<1]<-1 #default to first
      timept[timept>length(times)]<-length(times) #default to last

      tm3boundary <- rep(NA, length(timept))

      # loop over unique end times
      for (time.i in unique(timept)) {
         ctsel <- which(timept == time.i)
         if(length(ctsel)==1){
           #get height from geostrophic height
           gph<-get.var.ncdf(tm3file, varid="gph", start=c(point[ctsel,],1, time.i),count=c(1,1,nlevs,1))
           lev<-findInterval(agl4ct[ctsel], c(0,gph[-1]), rightmost.closed=TRUE,all.inside=TRUE)
           tm3boundary[ctsel] <- ini+get.var.ncdf(tm3file, varid=paste(tracer,"mix",sep=""), start=c(point[ctsel,],lev, time.i),count=c(1,1,1,1))

         } else { #more than one element to extract
           start.rd<-apply(point[ctsel,],2,min)
           count.rd<-apply(point[ctsel,],2,max)-start.rd + c(1,1)
           count.rd[count.rd<2]<-2 #read at least 2 items per dimension to have right overall dimension in object
           if(start.rd[1]+count.rd[1]-1>length(centerlons)){# stepping over dateline, need to read in last then first long. elements 
              count.rda<-count.rd;count.rdb<-count.rd;
              start.rda<-start.rd;start.rdb<-start.rd;
              count.ex<-(start.rd[1]+count.rd[1]-1-length(centerlons))
              count.rda[1]<-count.rda[1]-count.ex
              start.rdb[1]<-(start.rd[1]+1)%%length(centerlons)
              count.rdb[1]<-count.rd[1]-count.rda[1]

              gpha <- get.var.ncdf(tm3file, varid="gph", start=c(start.rd, 1, time.i),count=c(count.rda,nlevs,1))
              gphb <- get.var.ncdf(tm3file, varid="gph", start=c(start.rdb, 1, time.i),count=c(count.rdb,nlevs,1))
              if(length(dim(gpha))<3)gpha<-array(gpha,dim=c(1,dim(gpha)))
              if(length(dim(gphb))<3)gphb<-array(gphb,dim=c(1,dim(gphb)))
              require(abind)
              gph<-abind(gpha,gphb,along=1)
           } else {
              gph <- get.var.ncdf(tm3file, varid="gph", start=c(start.rd, 1, time.i),count=c(count.rd,nlevs,1))
           }
           gph[,,1]<-0 #set lowest level to zero
           hpt<-point[ctsel,1]*NA
           hpt[is.na(hpt)]<-1
           pt3d<-cbind(point[ctsel,],hpt) #turn into 3 D pointer
           point[ctsel,1]<-point[ctsel,1]-start.rd[1]+1#for use in gph as pointer
           point[ctsel,2]<-point[ctsel,2]-start.rd[2]+1
           for(l in 1:dim(gph)[3]){
             gphi<-gph[,,l][point[ctsel,]]
             pt3d[,"hpt"][agl4ct[ctsel]>gphi]<-l
           }
           start.rd<-apply(pt3d,2,min)
           count.rd<-apply(pt3d,2,max)-start.rd + c(1,1,1)
           pt3d<-pt3d+matrix(-start.rd+c(1,1,1),nrow=nrow(pt3d),ncol=ncol(pt3d),byrow=TRUE) #offset correction, to extract from small array to be read from ncdf
           if(any(count.rd>1)) tm3boundary[ctsel] <- ini+get.var.ncdf(tm3file, varid=paste(tracer,"mix",sep=""), start=c(start.rd, time.i),count=c(count.rd,1))[pt3d[,count.rd>1]]
           if(!any(count.rd>1))tm3boundary[ctsel] <- ini+get.var.ncdf(tm3file, varid=paste(tracer,"mix",sep=""), start=c(start.rd, time.i),count=c(count.rd,1)) #only one value
         }
      }

      resultm[,paste(tracer,"ini",sep="")] <- rep(0,nrow(resultm))
      resultm[result.selm,paste(tracer,"ini",sep="")] <- tm3boundary

      close.ncdf(tm3file)
      result[selm,paste(tracer,"ini",sep="")] <- resultm[,paste(tracer,"ini",sep="")]
   } #loop over months

   return(result)

}
