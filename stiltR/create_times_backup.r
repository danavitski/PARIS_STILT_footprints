#!/usr/bin/env Rscript

create_times <- function(station)
{
   print(paste('Creating times for ', station))
   
   #.libPaths("/home/mmolen/R/x86_64-pc-linux-gnu-library/3.4/")
   #require(base) # for julian
   #library('fields')
   #library('foreign')
   #library('methods')
   #library('methods')
   
   #script to create object "Times" w/ starting times for HF
   #output is a matrix with 4 columns: 
   #-fjul (fractional julian day since 1/1/1960)
   #-lat (deg N)
   #-lon (deg E)
   #-altitude (meters above ground)
   #
   #  $Id: create_times.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
   #---------------------------------------------------------------------------------------------------
   
   ##############
   #path to store receptor information (should be the same as in setStiltparam.r)
#   path    <- "/home/dkivits/STILT/STILT_Model/Output/" 
   path        <-"/scratch-shared/dkivits/STILT/footprints/"
   outname     <- paste('.RDataTimes_',station,'.hf',sep='') #name for object with receptor information
   ##############
   # julian(mm,dd,yyyy)
   
   y1 <- 2021
   m1 <-    7
   d1 <-    5
   h1 <-    11
   y2 <- 2021
   m2 <-    7
   d2 <-    5
   h2 <-    16
   dt <-  1/24

   fjul <- seq(julian(m1,d1,y1)+h1/24,julian(m2,d2,y2)+h2/24,dt)
#   # select only the hours between 11 and 16.
#   hours <- rep(seq(1,24,1), times = length(fjul)/24)
#   l <- list(jul = fjul, hours = hours)
#   df <- data.frame(l)
#   dfsubset <- subset(df, hours > 10 & hours < 17)
#   fjul <- dfsubset$jul
   fjul <- round(fjul,6)

   if(station == "TST") {
      # Virtual station at the boundary of the domain
      lon <- 34
      lat <- 71
      agl <- 20}
   
   if(station == "GAT") {
      #GAT
      lat <- 53.0655
      lon <- 11.4427
      agl <- 344.0
   } else if(station == "HPB") {
      #HPB
      lat <- 47.8015
      lon <- 11.0096
      agl <- 131.00
   } else if(station == "HTM") {
      #HTM
      lat <- 56.0976
      lon <- 13.4189
      agl <- 265.00
   } else if(station == "BRM") {
      #BRM
      lat <- 47.1896
      lon <- 8.1755
      agl <- 1009.00
   } else if(station == "KRE") {
      #KRE
      lat <- 49.5833
      lon <- 15.0833
      agl <- 250.00
   } else if(station == "LIN") {
      #LIN
      lat <- 52.2167
      lon <- 14.1167
      agl <- 99.00
   } else if(station == "KIT") {
      #KIT
      lat <- 49.1000
      lon <- 8.4380
      agl <- 200.00
   } else if(station == "STK") {
      #STK
      lat <- 53.0431
      lon <- 8.4588
      agl <- 200.00
   } else if(station == "CES") {
      #CES
      lat <- 51.9710
      lon <- 4.9270
      agl <- 200.00
   } else if(station == "LUT") {
      #LUT
      lat <- 53.4038
      lon <- 6.3529
      agl <-  60.00
   } else if(station == "SAC") {
      #SAC
      lat <- 48.7227
      lon <- 2.1422
      agl <- 100.00
   } else if(station == "IPR") {
      #IPR
      lat <- 45.8030
      lon <- 8.6270
      agl <- 210.00
   } else if(station == "HEI") {
      #HEI
      lat <- 49.4171
      lon <-  8.6745
      agl <- 10.00
   } else if(station == "FRE") {
      #FRE
      lat <- 49.5238
      lon <-  8.2170
      agl <- 10.00
   } else if(station == "GNS") {  
      #GNS
      lat <- 48.9850
      lon <- 2.4490
      agl <- 68.00
   } else if(station == "ROC") {
      #ROC
      lat <- 51.9327
      lon <-  3.9992
      agl <- 10.00
   } else if(station == "COU") { 
      #COU
      lat <- 48.9241
      lon <-  2.5680
      agl <- 10.00
   } else if(station == "AND") {
      #AND
      lat <- 49.0126
      lon <-  2.3017
      agl <- 10.00
   } else if(station == "OVS") {
      #OVS
      lat <- 48.7780
      lon <-  2.0482
      agl <- 10.00
   } else if(station == "ZWT") {
      #ZWT
      lat <- 51.9644
      lon <-  4.3947
      agl <- 10.00
   } else if(station == "WMS") {
      #WMS
      lat <- 51.7867
      lon <-  4.4505
      agl <- 10.00
   } else if(station == "2MV") {
      #WMS
      lat <- 51.8800
      lon <-  4.0300
      agl <- 10.00
   } else if(station == "RHI_UW") {
      #RHI_UW
      lat <- 49.50
      lon <-  8.25
      agl <- 10.00
   } else if(station == "RHI_DW") {
      #RHI_DW
      lat <- 49.42
      lon <-  8.67
      agl <- 10.00
   } else if(station == "BOR_UW") {
      #BOR_UW
      lat <- 44.53
      lon <- -0.62
      agl <- 10.00
   } else if(station == "BOR_DW") {
      #BOR_DW
      lat <- 45.07
      lon <- -0.52
      agl <- 10.00
   } else if(station == "LYO_UW") {
      #LYO_UW
      lat <- 45.48
      lon <-  4.97
      agl <- 10.00
   } else if(station == "LYO_DW") {
      #LYO_DW
      lat <- 46.04
      lon <-  4.68
      agl <- 10.00
   } else if(station == "LIL_UW") {
      #LIL_UW
      lat <- 50.12
      lon <-   3.1
      agl <- 10.00
   } else if(station == "LIL_DW") {
      #LIL_DW
      lat <- 50.98
      lon <-  3.06
      agl <- 10.00
   } else if(station == "LUX_UW") {
      #LUX_UW
      lat <- 49.46
      lon <-  5.8
      agl <- 10.00
   } else if(station == "LUX_DW") {
      #LUX_DW
      lat <- 49.85
      lon <-  6.53
      agl <- 10.00
   } else if(station == "RUR_UW") {
      #RUR_UW
      lat <- 51.01
      lon <-  6.16
      agl <- 10.00
   } else if(station == "RUR_DW") {
      #RUR_DW
      lat <- 51.85
      lon <-  7.94
      agl <- 10.00
   } else if(station == "BER_UW") {
      #BER_UW
      lat <- 52.21
      lon <- 12.92
      agl <- 10.00
   } else if(station == "BER_DW") {
      #BER_DW
      lat <- 52.7
      lon <- 13.71
      agl <- 10.00
   } else if(station == "MUN_UW") {
      #MUN_UW
      lat <- 47.95
      lon <-   11.01
      agl <- 10.00
   } else if(station == "MUN_DW") {
      #MUN_DW
      lat <- 48.34
      lon <-   12.04
      agl <- 10.00
   }
       
   assignr(outname,cbind(fjul,lat,lon,agl),path=path,printTF=T)
}   
#StationCoords    = {'GAT': ('DEU', 11.4427, 53.0655, 344.00), 'HPB': ('DEU', 11.0096, 47.8015, 131.00), 'KRE': ('CZE', 15.0833, 49.5833, 250.00), 'LIN': ('DEU', 14.1167, 52.2167,  99.00), 'KIT': ('DEU', 8.4380, 49.1000, 200.00), 
#                    'STK': ('DEU',  8.4588, 53.0431, 200.00), 'CES': ('NLD',  4.9270, 51.9710, 200.00), 'LUT': ('NLD',  6.3529, 53.4038,  60.00), 'SAC': ('FRA',  2.1422, 48.7227, 100.00), 'IPR': ('ITA', 8.6270, 45.8030, 210.00),
#                    'FRE': ('DEU',  8.2170, 49.5238,  10.00), 'HEI': ('DEU',  8.6745, 49.4171,  10.00),
#                    'ROC': ('NLD',  3.9992, 51.9327,  10.00), '2MV': ('NLD',  4.0300, 51.8800,  10.00), 'ZWT': ('NLD',  4.3947, 51.9644,  10.00), 'WMS': ('NLD', 4.4505, 51.7867,  10.00), 
#                    'OVS': ('FRA',  2.0482, 48.7780,  10.00), 'COU': ('FRA',  2.5680, 48.9241,  10.00), 'AND': ('FRA',  2.3017, 49.0126,  10.00), 'GNS': ('FRA', 2.4490, 48.9850,  68.00)}
  
