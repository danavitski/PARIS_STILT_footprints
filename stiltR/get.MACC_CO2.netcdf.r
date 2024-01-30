#***************************************************************************************************
#  read initial values from CarbonTracker file and attach them to the "result" matrix
#***************************************************************************************************

get.MACC_CO2.netcdf <- function(yr4=NULL, mon=NULL, day=NULL, hr=NULL, co2inifile=NULL,
                             result=NULL, result.sel=NULL, spec="co2ini") {
#ls -hl /Net/Groups/BSY/people/vbeck/WRF_runs/BC_WRF/GEMS_reanalysis
#co2inifile="/Net/Groups/BSY/people/vbeck/WRF_runs/BC_WRF/GEMS_reanalysis/co_nov_dec_08.nc"
# --------------------------------------------------------------------------------------------------
#  Interface
#  =========
#
#  yr4          YYYY of the receptor time
#  mon          month "
#  day          day   "
#  hr           hour  "
#  co2inifile   netCDF file name (absolute path) of the GEMS fields
#  result       matrix with columns "btime", "lat", "lon" and "pres"
#               gets the column "co2ini" attached/overwritten on output
#  result.sel   selector of same length as the columns of "result": which rows should be considered
#
# --------------------------------------------------------------------------------------------------
#  $Id: get.TM3.netcdf.r,v 1.1 2009-02-24 14:27:54 gerbig Exp $
# --------------------------------------------------------------------------------------------------


   if (!is.element("package:ncdf", search())) library("ncdf") # Load ncdf library
   refdate <- month.day.year(floor(julian(mon, day, yr4) + hr/24 - result[result.sel, "btime"]/24))
   yearmond<-paste(refdate$year,substring(refdate$month+100,2,3),substring(refdate$day+100,2,3),sep="")
#

   for(yearmondi in unique(yearmond)){
      selm<-yearmond==yearmondi
      resultm<-result[selm,]
      result.selm<-result.sel[selm]
      fn<-basename(co2inifile)
      co2inifilem<-paste(dirname(co2inifile),"/",substring(fn,1,nchar(fn)-11),yearmondi,".nc",sep="")
      tm3file <- open.ncdf(co2inifilem, write=F)
      centerlats <- get.var.ncdf(tm3file, varid="latitude")           # Center (assumed here, not sure!)
      centerlons <- get.var.ncdf(tm3file, varid="longitude")           # Center
      times <- get.var.ncdf(tm3file, varid="time")           # Center (instantaneous)
      levs <- get.var.ncdf(tm3file, varid="level")           # levels
      nlevs<-length(levs)
      ini <- att.get.ncdf(tm3file, varid="co2",attname="add_offset")$value           # global CO2 offset
      scl <- att.get.ncdf(tm3file, varid="co2",attname="scale_factor")$value           # global CO2 offset

      # Get Ending positions
      p4ct <- resultm[result.selm, "pres"]
      lat4ct <- resultm[result.selm, "lat"]
      lon4ct <- resultm[result.selm, "lon"]
      lon4ct[lon4ct<0] <- lon4ct[lon4ct<0]+360
      ps4ct <- 1013*exp(-resultm[result.selm,"grdht"]/8000)#get surface pressure using simple scale height

      #-- lat lon pointer
      dlat<-unique(round(diff(centerlats),4))
      dlon<-unique(round(diff(centerlons),4))
      latpt<-round((lat4ct-centerlats[1])/dlat)+1;latpt[latpt>length(centerlats)]<-length(centerlats);latpt[latpt<1]<-1
      lonpt<-round((lon4ct-centerlons[1])/dlon)+1;lonpt[lonpt>length(centerlons)]<-length(centerlons)
      point<-cbind(lonpt,latpt) #spatial pointer (2D)

      #-- time pointer
      delt.h<-(times[2]-times[1])
      timept<-round((julian(mon, day, yr4,c(1, 1, 1900)) + hr/24 - resultm[result.selm, "btime"]/24)*24/delt.h)*delt.h #hours since 1900-01-01 00:00:0.0
      timept<-(timept-times[1])/delt.h+1
      timept[timept<1]<-1 #default to first
      timept[timept>length(times)]<-length(times) #default to last

      tm3boundary <- rep(NA, length(timept))

      #prepare for height pointer
      # read heights from ECM website for 137 levels (http://old.ecmwf.int/products/data/technical/model_levels/model_def_137.html)
      plevs<-c(0.0000,0.0200,0.0310,0.0467,0.0683,0.0975,0.1361,0.1861,0.2499,0.3299,0.4288,0.5496,0.6952,0.8690,1.0742,1.3143,1.5928,1.9134,2.2797,2.6954,3.1642
        ,3.6898,4.2759,4.9262,5.6441,6.4334,7.2974,8.2397,9.2634,10.3720,11.5685,12.8561,14.2377,15.7162,17.2945,18.9752,20.7610,22.6543,24.6577,26.7735
        ,29.0039,31.3512,33.8174,36.4047,39.1149,41.9493,44.9082,47.9915,51.1990,54.5299,57.9834,61.5607,65.2695,69.1187,73.1187,77.2810,81.6182,86.1450
        ,90.8774,95.8280,101.0047,106.4153,112.0681,117.9714,124.1337,130.5637,137.2703,144.2624,151.5493,159.1403,167.0450,175.2731,183.8344,192.7389
        ,201.9969,211.6186,221.6146,231.9954,242.7719,253.9549,265.5556,277.5852,290.0548,302.9762,316.3607,330.2202,344.5663,359.4111,374.7666,390.6450
        ,407.0583,424.0190,441.5395,459.6321,478.3096,497.5845,517.4198,537.7195,558.3430,579.1926,600.1668,621.1624,642.0764,662.8084,683.2620,703.3467
        ,722.9795,742.0856,760.5996,778.4661,795.6396,812.0847,827.7756,842.6959,856.8376,870.2004,882.7910,894.6222,905.7116,916.0815,925.7571,934.7666
        ,943.1399,950.9082,958.1037,964.7584,970.9046,976.5737,981.7968,986.6036,991.0230,995.0824,998.8081,1002.2250,1005.3562,1008.2239,1010.8487,1013.2500)

      # loop over unique end times
      for (time.i in unique(timept)) {
         ctsel <- which(timept == time.i)
         if(length(ctsel)==1){
           #get vertical pointer: 
           #1) use heights from ECM website for 60 levels (in plevs)
           #2) get surface pressure using simple scale height
           psurf<-ps4ct[ctsel]
           #3) apply compression factor (simplified ...)
           pres<-p4ct[ctsel]*(1013.25/psurf);pres[pres>max(plevs)]<-max(plevs);pres[pres<=min(plevs)]<-min(plevs)+0.001
           #now get vertical pointer
           hpt<-cut(pres,breaks=plevs,labels=FALSE)
           hpt[hpt<min(levs)]<-min(levs)
           hpt[hpt>max(levs)]<-max(levs)
           tm3boundary[ctsel] <- get.var.ncdf(tm3file, varid="co2", start=c(point[ctsel,],hpt, time.i),count=c(1,1,1,1))
           #Note that scale factor and ini (global offset) has been applied in get.var.ncdf
         } else { #more than one element to extract
           start.rd<-apply(point[ctsel,],2,min)
           count.rd<-apply(point[ctsel,],2,max)-start.rd + c(1,1)
           count.rd[count.rd<2]<-2 #read at least 2 items per dimension to have right overall dimension in object
           #get vertical pointer: 
           #1) use heights from ECM website for 60 levels (in plevs)
           #2) get surface pressure using simple scale height
           psurf<-ps4ct[ctsel]
           #3) apply compression factor (simplified ...)
           pres<-p4ct[ctsel]*(1013.25/psurf);pres[pres>max(plevs)]<-max(plevs);pres[pres<=min(plevs)]<-min(plevs)+0.001
           #now get vertical pointer
           hpt<-cut(pres,breaks=plevs,labels=FALSE)
           hpt[hpt<min(levs)]<-min(levs)
           hpt[hpt>max(levs)]<-max(levs)
           pt3d<-cbind(point[ctsel,],hpt) #turn into 3 D pointer
           start.rd<-apply(pt3d,2,min)
           count.rd<-apply(pt3d,2,max)-start.rd + c(1,1,1)
           pt3d<-pt3d+matrix(-start.rd+c(1,1,1),nrow=nrow(pt3d),ncol=ncol(pt3d),byrow=TRUE) #offset correction, to extract from small array to be read from ncdf
           if(any(count.rd>1)) tm3boundary[ctsel] <- get.var.ncdf(tm3file, varid="co2", start=c(start.rd, time.i),count=c(count.rd,1))[pt3d[,count.rd>1]]
           if(!any(count.rd>1))tm3boundary[ctsel] <- get.var.ncdf(tm3file, varid="co2", start=c(start.rd, time.i),count=c(count.rd,1)) #only one value
           #Note that scale factor and ini (global offset) has been applied in get.var.ncdf
         }
      }
      resultm[,spec] <- rep(0,nrow(resultm))
      resultm[result.selm,spec] <- tm3boundary

      close.ncdf(tm3file)
      result[selm,spec] <- resultm[,spec]
   } #loop over months

   #unit conversion to ppm
   kgkg2ppb<-28.97/44.01*1E6
   result[result.sel,spec] <- result[result.sel,spec]*kgkg2ppb

   return(result)

}
