#!/usr/bin/env Rscript

create_times <- function(stationfile, station, filter_times = 'TRUE', outpath = '/scratch-shared/dkivits/STILT/footprints/')
{
   # Filter the stationfile for the station
   index <- which(stationfile$code == station)

   print(paste('Creating times for', stationfile$station[index]))
   print(paste('This is station', index, 'out of', nrow(stationfile)))

   #script to create object "Times" w/ starting times for HF
   #output is a matrix with 4 columns: 
   #-fjul (fractional julian day since 1/1/1960)
   #-lat (deg N)
   #-lon (deg E)
   #-altitude (meters above ground)
   #---------------------------------------------------------------------------------------------------
   
   ##############
   # path to store receptor information (should be the same as in setStiltparam.r)
   # path    <- "/home/dkivits/STILT/STILT_Model/Output/" 
   path        <- outpath
   outname     <- paste('.RDataTimes_',station,'.hf',sep='') #name for object with receptor information
   ##############
   # julian(mm,dd,yyyy)
   
   y1 <- 2021
   m1 <-    1
   d1 <-    1
   h1 <-    1  # keep this at 1 for now
   y2 <- 2021
   m2 <-    12
   d2 <-    31
   h2 <-    23 # keep this at 23 for now
   dt <-  1/24

  # y2 <- 2021
  # m2 <-    1
  # d2 <-    2
  # h2 <-    23

   # Get the station information
   lat = stationfile$lat[index]
   lon = stationfile$lon[index]
   agl = stationfile$alt[index]
   corrected_agl = stationfile$corrected_alt[index]
   type = stationfile$type[index]
   offset = stationfile$utc2lst_offset[index]

   if (!is.na(corrected_agl)) {
      agl = corrected_agl
   }

   print(paste('Lat, lon and agl are', lat, lon, agl))
   print(paste('The offset from UTC to local time for this station is', offset, 'hours'))

   # Create the julian dates
   fjul <- seq(julian(m1,d1,y1)+h1/24,julian(m2,d2,y2)+h2/24,dt)
   
   if (filter_times == 'TRUE'){
      #hours <- rep(seq(1,24,1), times = length(fjul)/24)
      hours <- fjul%%1 * 24
      hours[hours==0] = 24

      l <- list(jul = fjul, hours = hours)
      df <- data.frame(l)

      # Filter mountain type sites
      if (type == "2") {
         print("This is a mountain site")
         
         if (offset <= 1) {
            # Select only the hours between 23 and 03 local time for mountain sites
            dfsubset <- subset(df, hours < (4 + offset) | hours >= abs(23 + offset))
            fjul <- dfsubset$jul

         } else {
            # Select only the hours between 23 and 03 local time for mountain sites
            dfsubset <- subset(df, hours >= (23 + offset) & hours < (4 + offset))
            fjul <- dfsubset$jul
         }

      # Filter lowland type sites
      } else {
         print("This is a lowland site")

         # Select only the hours between 11 and 15 local time for lowland sites
         dfsubset <- subset(df, hours >= (11 + offset) & hours < (16 + offset))
         fjul <- dfsubset$jul
      }
   }

   fjul <- round(fjul,6)

   assignr(outname,cbind(fjul,lat,lon,agl),path=path,printTF=T)
}  