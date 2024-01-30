flttrack<-function(lon=NULL,lat=NULL,col=6,cex=0.8,cities=F,minpop=100000,pts=T,newmotifTF=T,database=NULL,xylims=NULL){
#Plots the flight tracks on top of U.S. state map
#'pts' specifies whether points are plotted or not
#If don't want to specify 'lon' & 'lat' separately, then could just give entire particle object as first argument--e.g., 'flttrack(bdat)'
#'minpop' specifies the minimum population in cities to plot
#'xylims': lower left and upper right corners (x1, x2, y1, y2)
#8/14/2000 by JCL
#
#  $Id: flttrack.r,v 1.13 2009-03-09 08:02:47 gerbig Exp $
#---------------------------------------------------------------------------------------------------

#------for debugging--------#
#lon<-c(-60,-90,-120);lat<-c(48,45,65)
#lon<-c(lon,-117.78662,-97.48695);lat<-c(lat,34.1248,24.71251)
#pts<-T;newmotifTF<-T;cities<-F
#cex<-0.8;col<-6
#------for debugging--------#

if(!exists("CountryPolySP"))load(paste(sourcepath,"CountryPolySP.rdata",sep=""))
require(sp)
require(maps)
require(mapdata)

if(!is.null(lon)&is.null(lat)){
   if(sum(dimnames(lon)[[2]]=="lat")==1)lat<-lon[,"lat"]
   if(sum(dimnames(lon)[[2]]=="LAT")==1)lat<-lon[,"LAT"]
   if(sum(dimnames(lon)[[2]]=="lon")==1)lon<-lon[,"lon"]
   if(sum(dimnames(lon)[[2]]=="LON")==1)lon<-lon[,"LON"]
}

if(newmotifTF)x11(type="dbcairo")

xlims<-range(lon,na.rm=T);ylims<-range(lat,na.rm=T)
xlims[1]<-floor(xlims[1]-0.2*diff(xlims));xlims[2]<-ceiling(xlims[2]+0.2*diff(xlims))
ylims[1]<-floor(ylims[1]-0.2*diff(ylims));ylims[2]<-ceiling(ylims[2]+0.2*diff(ylims))

#implement aspect ratio to make pretty map
yctr<-mean(lat,na.rm=T);xctr<-mean(lon,na.rm=T)
dx<-abs(diff(xlims));dy<-abs(diff(ylims))
asp<-cos(pi*yctr/180)
if(is.null(xylims)){
  if(dy<asp*dx){
     dy<-asp*dx;ylims<-c(yctr-dy/2,yctr+dy/2)
  }else{
     dx<-dy/asp;xlims<-c(xctr-dx/2,xctr+dx/2)
  }
}
#plot(0,0,xlim=xlims,ylim=ylims,type="n",axes=F,xlab="",ylab="")
#plot(0,0,xlim=xlims,ylim=ylims,type="n",xlab="",ylab="")
#par(new=T);plot(CountryPolySP,add=T,xlim=xlims,ylim=ylims)
plot(CountryPolySP,xlim=xlims,ylim=ylims,xlab="",ylab="",axes=TRUE)

if(pts)points(lon,lat,pch=16,col=col,cex=cex)
#if(cities)points(city.x,city.y,pch=16,col=3,cex=0.7,err=-1)
if(cities){
  #remove cities close to one another--otherwise letters cover one another and is ugly
  data(world.cities) #need to first load data
  tmp<-world.cities

  sel<-tmp$name%in%"Longueuil";sel<-sel|(tmp$name%in%"Laval")
  dum<-data.frame(tmp$name[!sel],tmp$lat[!sel],tmp$long[!sel],tmp$pop[!sel],
                tmp$country.etc[!sel],tmp$capital[!sel])
  names(dum)<-c("name","lat","long","pop","country.etc","capital")
  map.cities(x=dum,minpop=minpop)   #plots cities over minimum population specified by 'minpop'
 
} #if(cities){

}  #flttrack<-function(lon=NULL,lat=NULL,col=6,cex=0.8,cities=F,pts=T,newmotifTF=T){



