get.edgar.time<-function(result,mon,day,yr4,hr,ann.oil=NULL,ann.gas=NULL,ann.coal=NULL){
#apply time factors for different emission categories from edgar
#uses variables assigned in set.edgar.r, called from setStiltparam.r

  #get day of year, local time;1Jan=1
  doy<-floor(julian(mon,day,yr4,c(1,0,yr4))+hr/24-result[,"btime"]/24+result[, "lon"]/360) 
  doy[doy<1]<-doy[doy<1]+julian(12,31,yr4-1,c(1,1,yr4-1)) #make sure that only positive values are used as time pointers
  ly<-ifelse(julian(12,31,yr4,c(1,0,yr4))==366,TRUE,FALSE)#leap year
  sh<-ifelse(max(result[,"lat"])<0,TRUE,FALSE)
  dayvar<-get(paste("daily_var",ifelse(ly,".l",".n"),ifelse(sh,".sh",""),sep="")) #this is monthly variation, smoothly interpolated to each doy

  #get hour of day and weekday
  ltime<-weekdayhr(yr4, mon, day, hr, -result[, "btime"]*60, diffGMT=round(result[, "lon"]*24/360))# last column is weekday, sunday=0, monday=1 etc.
  ltime[ltime[,"hr"]==0,"hr"]<-24 
  ltime[ltime[,"weekd"]==0,"weekd"]<-7 #edgar time factors use 7 as sunday
  if(emisscatTF){#fluxtracers defined in set.edgar
    ipcc2time.n<-ipcc2time[tracer.info["wanted",tolower(tr.edg)]=="TRUE"]#reduce length of ipcc2time vector to only those tracers that are wanted
    t.fac<-dayvar[doy,ipcc2time.n]*hourly_var[ltime[,"hr"],ipcc2time.n]*weekly_var[ltime[,"weekd"],ipcc2time.n]
    tmp<-result[,(ncol(result)-length(fluxtracers)+1):ncol(result)][,c(grep("co2\\.",fluxtracers),grep("co\\.",fluxtracers),grep("ch4\\.",fluxtracers),grep("n2o\\.",fluxtracers))]
    result[,(ncol(result)-length(fluxtracers)+1):ncol(result)][,c(grep("co2\\.",fluxtracers),grep("co\\.",fluxtracers),grep("ch4\\.",fluxtracers),grep("n2o\\.",fluxtracers))]<-tmp*t.fac
  }
  if(!is.null(ann.oil)){
    tracer.cols<-(ncol(result)-length(fluxtracers)+1):ncol(result)
    for(fuel in c(".oil",".gas",".coal")){
      result[,tracer.cols][,grep(fuel,fluxtracers)]<-result[,tracer.cols][,grep(fuel,fluxtracers)]*get(paste("ann",fuel,sep=""))
    }
  }
  return(result)
}
