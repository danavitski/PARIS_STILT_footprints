#script to set edgar tracer names and flux fields for edgar emissions by fuel type and category
#for changing tracer settings make edits in section far below entitled 
#"HERE SET TRACERS FOR USE IN STILT FORWARD SIMULATIONS"

#for augmentation to multiple species (CO2, CO, CH4): editing required


#get list of all nc files for year 2010
getfs<-function(path,tr){
#function to get list of ncdf files for specific tracer
#path: the path toi search (recursively) for the files
#tr: the tracer for which files are returned
  edfs<-list.files(path,pattern="*.nc",recursive=TRUE)

  if(tr=="co2")edfs<-edfs[grepl("data_CO2\\/",edfs)|grepl("data_CO2x\\/",edfs)] #only co2 files from Greet
  if(tr=="co")edfs<-edfs[grep("data_CO/",edfs)] #only CH4 files from Greet
  if(tr=="ch4")edfs<-edfs[grep("data_CH4/",edfs)] #only CH4 files from Greet
  if(tr=="n2o")edfs<-edfs[grep("data_N2O/",edfs)] #only CH4 files from Greet

  edfs<-edfs[edfs!="data_CO2/CO2_excl_short_cycle_org_C_2010_IPCC_1A4.0.1x0.1/v42_FT2010_CO2_RCO_SOLID_nonbiowaste.0.1x0.1.nc"] #was included twice from Greet
  edfsn<-basename(edfs)
  return(list(edfs,edfsn))
}

for(tracer in c("co2","ch4","co","n2o")){
  print(paste("tracer",tracer))
  edfs<-getfs(edgar.path,tracer)
  edfsn<-edfs[[2]]
  edfs<-edfs[[1]]
#  print(length(edfsn))
############################################################################################################
#get gategories (IPCC), unify IPCC categories accross tracers
############################################################################################################
#get IPCC categories from file and path names, make necessary adjustments to have the same category names accross tracers
  ipcc<-regexpr("IPCC",edfs)
  x1<-regexpr("1x0",edfs)
  x1[x1-4-(ipcc+5)>20]<-regexpr("grmnht",edfs)[x1-4-(ipcc+5)>20]+2 #some have "1x0" close to the end of the full file name, use different filter there
  ipcc<-substring(edfs,ipcc+5,x1-4)
  if(tracer=="co2")ipcc[edfsn=="v42_FT2010_CO2_short_cycle_org_C_AWB.0.1x0.1.nc"]<-"4F" #add 4F: agricultural waste burning
  if(tracer=="co2")ipcc[edfsn=="grmnht_v43_EM_CO2_2010_TNG_CO2_2010_TNG_grmnht.0.1x0.1.nc"]<-"1A3ce"
  if(tracer=="n2o")ipcc[edfsn=="v42_FT2010_N2O_2010_IPCC_7B_7C.0.1x0.1.nc"]<-"7BC"
  #ipcc[ipcc=="7A"]<-"7" #use same cat for CH4 as for CO2
  ipcc[ipcc=="1A1"]<-"1A1a" #use same name for CH4 as for all others
  ipcc[ipcc=="1A1bcr_1B1"]<-"1A1bcr" #combine those for co2; was obviously used to separate short cycle from rest, this will seperate again via fuel types
  ipcc[ipcc=="1A2_6CD"]<-"1A2+6CD" #use same name for CO2 as for all others
  ipcc[ipcc=="1A2+6C"]<-"1A2+6CD" #use same cat for N2O as for all others (unclear meaning of 6D)
  ipcc[ipcc=="6AD"]<-"6A" #use same cat for N2O as for CH4 (contains only solid waste disposal - land fills for both tracers)
  ipcc[ipcc=="1A3a1C1"]<-"1A3a+1C1" #use same name for CH4 as for CO and N2O
  ipcc[ipcc=="1A3a"]<-"1A3a+1C1" #was only in CO2, rename as it also includes transatlantic flights 
  ipcc[ipcc=="1A3d"]<-"1A3d+1C2" #was only in CO2, rename as it probably also includes shipping (in addition to inland waterways) 
  ipcc[ipcc=="1A3d1C2"]<-"1A3d+1C2" #use same name for CH4 as for CO and N2O
  ipcc[ipcc=="1B2"]<-"1B2abc"     #for CO2, contains oil prod. + gas prod. + venting and flaring
  ipcc[ipcc=="2B"]<-"2BEFG+3"     #use same name this one (chemical industry and solvents)
  ipcc[ipcc=="2B_3"]<-"2BEFG+3"     #use same name this one (chemical industry and solvents)
  ipcc[ipcc=="2B+G+3"]<-"2BEFG+3"     #use same name this one (chemical industry and solvents)
  ipcc[ipcc=="7"]<-"7A" #use same name for CO as for CO2 and CH4 (fossil fuel fires)
  ipcc[ipcc=="7ABC"]<-"7A" #use same name for N2O as for CO (fossil fuel fires)

############################################################################################################
#exclude unwanted gategories
############################################################################################################
#which classes do we not want?
#agricultural soil emissions for CO2, IPCC 4D_AGS -> done in biosphere model
#large scale biomass burning	5+4E_BMB -> done in biomass burning model 
  ipcc.not<-c("5+4E")#excludes large scale biomass burning; needed from biomass burning model
  if(tracer=="co2")ipcc.not<-c("4C_4D", #excludes rice CH4 and agricultural soil emissions of CO2
             "5ACF_D")#excludes large scale biomass burning
  edfsn<-edfsn[!ipcc%in%ipcc.not]
  edfs<-edfs[!ipcc%in%ipcc.not]
  ipcc<-ipcc[!ipcc%in%ipcc.not]

############################################################################################################
#get fuel types
############################################################################################################
#based on file or path name
#from excel file Sanam_source_categories_14March2014f.xls
#Hard coal	Brown coal	Peat	Solid waste	heavy oil	Light oil	
#heavy+light+gas slurry production (venting & flaring emission)	
#Natural gas	Derived gas	Solid biomass	Bio liquid	Bio gas	  Other

#check ambiguities in filenames
#"solid_non_bio" or "solid_nonbiowaste" -> "solid_non"
#"solid_bio" or "solid_biomass" -> "solid_bio"
#"liquid_bio" or "bio_liquid"
#"gas_bio" or "bio_gas"
#"solid_bio" or "bio_solid"

  ftypes.long<-c("hard coal","brown coal","peat","solid waste non-bio","heavy oil","light oil",
  "venting and flaring of heavy and light oil and gas",
  "natural gas","derived gas","solid biomass","liquid biomass","biogas","other")
  ftypes<-c("coal_hard","coal_brown","coal_peat","solid_non_bio","oil_heavy","oil_light",
  "VAF",
  "gas_nat","gas_der","solid_bio","liquid_bio","gas_bio", "others")
  #define alternative names also
  ftypes2<-c("coal_hard","coal_brown","coal_peat","solid_nonbiowaste","oil_heavy","oil_light",
  "OIl_lght+hvy+gas_VAF",
  "gas_nat","gas_der","bio_solid","bio_liquid","bio_gas", "others")
  ftypes3<-c("coal_hard","coal_brown","coal_peat","solid_waste","oil_heavy","oil_light",
  "oil_lght+hvy+gas_VAF",
  "gas_nat","gas_der","bio_solid","bio_liquid","bio_gas", "cement")
  ftypes4<-c("coal_hard","coal_brown","coal_peat","solid_wst_nonbio","oil_heavy","oil_light",
  "OIl_lght+hvy+gas_VAF",
  "gas_nat","gas_der","bio_solid","bio_liquid","bio_gas", "others")
  ftypes5<-c("coal_hard","coal_brown","coal_peat","SOLIDwaste","oil_heavy","oil_light",
  "OIl_lght+hvy+gas_VAF",
  "gas_nat","gas_der","bio_solid","bio_liquid","bio_gas", "others")

  nams<-NULL
  ftypids<-rep(NA,length(edfsn)) #this is the vector of indices pointing to fuel type; picks out a fuel type for each file
  for(i in 1:length(ftypes)){
  #  print(ftypes.long[i])
    m1<-edfsn[grep(ftypes[i] ,edfsn, ignore.case=TRUE)]
    m2<-edfsn[grep(ftypes2[i],edfsn, ignore.case=TRUE)] 
    m2<-m2[!m2%in%m1]
    m3<-edfsn[grep(ftypes3[i],edfsn, ignore.case=TRUE)] 
    m3<-m3[!m3%in%m2]
    m4<-edfsn[grep(ftypes4[i],edfsn, ignore.case=TRUE)] 
    m4<-m4[!m4%in%m3]
    m5<-edfsn[grep(ftypes5[i],edfsn, ignore.case=TRUE)] 
    m5<-m5[!m5%in%m4]
  #  print(paste("1:",length(m1)))
  #  print(paste("2:",length(m2)))
    nams<-c(nams,m1,m2,m3,m4,m5)
    ftypids[edfsn%in%c(m1,m2,m3,m4,m5)]<-i
  }
  edfsn[!edfsn%in%nams]
  which(!edfsn%in%nams)

  #modify fueltype for one case, flaring emissions from light and heavy oil was included twice for CH4 and N2O:
  ftypids[edfsn=="grmnht_v43_2010_CH4__CH4_2010_PRO_OIL_light+heavy+gas_grmnht.0.1x0.1.nc"]<-which(ftypes=="gas_nat")
  ftypids[edfsn=="grmnht_v43_N2O_2010_N2O_2010_PRO_OIL_light+heavy+gas_grmnht.0.1x0.1.nc"]<-which(ftypes=="gas_nat")
  #add fueltype "others"
  #note: "other" includes Cement&Lime production (2A_NMM), chemical indurtry and solvents (2B+G+3_CHS), 
  #process emissions metall industry (2C_INF), and agricultural soil emissions (4D_AGS)
  ipcc[is.na(ftypids)] #those are the IPCC categs with not yet associated fuel type
  #cbind(ipcc[is.na(ftypids)] ,edfs[is.na(ftypids)]) #those are the IPCC categs with not yet associated fuel type
  ftypids[ipcc=="1A3a+1C1"]<-which(ftypes=="oil_light") #aviation
  ftypids[ipcc=="1A3ce"]<-which(ftypes=="oil_heavy") #non-road ground transport
  ftypids[ipcc=="1A3d+1C2"]<-which(ftypes=="oil_heavy") #inland waterways and shipping
  ftypids[ipcc=="2A"]<-which(ftypes=="others") #(process emissions of non-metallic minerals industry (cement, lime))
  ftypids[ipcc=="2BEFG+3"]<-which(ftypes=="others") #chem. industry + solvents
  ftypids[ipcc=="2C"]<-which(ftypes=="others") #process emissions metal industry
  ftypids[ipcc=="4A"]<-which(ftypes=="others") #agricultural waste burning
  ftypids[ipcc=="4B"]<-which(ftypes=="others") #(manure management in agriculture)
  ftypids[ipcc=="4C"]<-which(ftypes=="others") #(rice cultivation in agriculture)
  ftypids[ipcc=="4D"]<-which(ftypes=="others") #(agricultural soil emissions)
  ftypids[ipcc=="4F"]<-which(ftypes=="solid_bio") #agricultural waste burning
  ftypids[ipcc=="6A"]<-which(ftypes=="others") #(solid waste disposal in landfills)#agricultural waste burning
  ftypids[ipcc=="6B"]<-which(ftypes=="others") #(wastewater treatment)
  ftypids[ipcc=="7A"]<-which(ftypes=="coal_hard") #fossil fuel fires
  ftypids[ipcc=="7BC"]<-which(ftypes=="others") #fossil fuel fires

  fuelnames<-ftypes3[ftypids] #fuel names for identification of tracer
  assign(paste("fuelnames",tracer,sep="."),fuelnames)
  assign(paste("ipcc",tracer,sep="."),ipcc)
  assign(paste("edfs",tracer,sep="."),edfs)
} #loop over tracers


tr.edg<-c(paste("co2",get(paste("ipcc","co2",sep=".")),get(paste("fuelnames","co2",sep=".")),sep="."),
          paste("co" ,get(paste("ipcc","co" ,sep=".")),get(paste("fuelnames","co" ,sep=".")),sep="."),
          paste("ch4",get(paste("ipcc","ch4",sep=".")),get(paste("fuelnames","ch4",sep=".")),sep="."),
          paste("n2o",get(paste("ipcc","n2o",sep=".")),get(paste("fuelnames","n2o",sep=".")),sep="."))

############################################################################################################
############################################################################################################
#HERE SET TRACERS FOR USE IN STILT FORWARD SIMULATIONS
############################################################################################################
############################################################################################################

############################################################################################################
#which of the four species do you want to simulate based on EDGAR Categories/fueltypes? (this results in about 60 tracers per species in STILT output from Trajecvprm.r)
spec.want.catfuel<-c("co2","co","ch4","n2o") ######MODIFY HERE, any combination of "co2","co","ch4", or "n2o"
############################################################################################################


#which of the four species do you want to simulate based on EDGAR Categories/fueltypes? (this results in about 60 tracers per species in STILT output from Trajecvprm.r)
tr.edg.tmp<-NULL
for(sp in spec.want.catfuel)tr.edg.tmp<-c(tr.edg.tmp,tr.edg[grepl(paste(sp,".",sep=""),tr.edg,fixed=TRUE)])
tr.edg<-tr.edg.tmp

  tracer.names.all<-tolower(c(tr.edg[grepl("co2.",tr.edg,fixed=TRUE)],"co2",                                 #vector of names for which mixing ratios are calculated
                              tr.edg[grepl("co.",tr.edg,fixed=TRUE)],"co",
                              tr.edg[grepl("ch4.",tr.edg,fixed=TRUE)],"ch4",
                              tr.edg[grepl("n2o.",tr.edg,fixed=TRUE)],"n2o",
                              "h2","cofire","rn","ch4nat"))	
  l.traddco2<-length(tr.edg[grepl("co2.",tr.edg,fixed=TRUE)])
  l.traddco <-length(tr.edg[grepl("co.",tr.edg,fixed=TRUE)])
  l.traddch4<-length(tr.edg[grepl("ch4.",tr.edg,fixed=TRUE)])
  l.traddn2o<-length(tr.edg[grepl("n2o.",tr.edg,fixed=TRUE)])

#################################################
#Below specifies information on a) if tracer is wanted, b) if emissions are avail. as ncdf, c) the full name for emission file,
# d) the kind of lateral boundary condition (LBC), and e) the full filename for the LBC
tracer.info<-rbind(
#################################################
#### edit following part for each tracer ########
#    "co2",    "co",  "ch4",   "n2o",   "h2",   "cofire",   "rn" ,"ch4nat"
c(rep(T,l.traddco2),
      T,        							#want CO2 as tracer (using VPRM and other inventory and LBC specified below)? Select T or F
            rep(T,l.traddco),
                T,     							#want CO as tracer (using other inventory specified below)? Select T or F
                    rep(T,l.traddch4),
                        T,						#want CH4 as tracer (using other inventory specified below)? Select T or F
                            rep(T,l.traddn2o),
                                F,					#want N2O as tracer (using other inventory specified below)? Select T or F
                                         F,        F,      T,      F), #want tracer? H2 is not ready at the moment, needs proper implementation in Trajecvprm()
c(rep(T,l.traddco2),
      T,        							# CO2 fluxes from other inventory (not Edgar cat. fuel) as netCDF? Select T or F
            rep(T,l.traddco),
                T,        							# CO fluxes from other inventory (not Edgar cat. fuel) as netCDF? Select T or F
                    rep(T,l.traddch4),
                        T,        							# CH4 fluxes from other inventory (not Edgar cat. fuel) as netCDF? Select T or F
                            rep(T,l.traddn2o),
                                F,        							# N2O fluxes from other inventory (not Edgar cat. fuel) as netCDF? Select T or F
                                         F,        F,      T,      T), #surface fluxes in netCDF format (T) or as R objects (F)? (In case of CO2: same format for fossil fluxes
                                          	 #and biospheric fields (veg cover, modis indices)
c(paste(edgar.path,get(paste("edfs","co2",sep=".")),sep="")[0:l.traddco2],                               #not modify
      "/Net/Groups/BSY/people/ukarst/STILT_prepro/EDGAR/EDGAR4.1_BP2012/data/ncdf_stilt_eu2/EDGAR_4.1_0.1x0.1.CO2.2011.nc",                 #full name for emission file (CO2)
  paste(edgar.path,get(paste("edfs","co",sep=".")),sep="")[0:l.traddco],                                #not modify
              "/Net/Groups/BSY/people/ukarst/STILT_prepro/EDGAR/EDGAR4.1_BP2012/data/ncdf_stilt_eu2/EDGAR_4.1_0.1x0.1.CO.2011.nc",                 #full name for emission file (CO)
  paste(edgar.path,get(paste("edfs","ch4",sep=".")),sep="")[0:l.traddch4],                                #not modify
                    "/Net/Groups/BSY/people/ukarst/STILT_prepro/EDGAR/EDGAR4.1_BP2012/data/ncdf_stilt_eu2/EDGAR_4.1_0.1x0.1.CH4.2011.nc",              #full natme for emission file (CH4)
  paste(edgar.path,get(paste("edfs","n2o",sep=".")),sep="")[0:l.traddn2o],                                #not modify
                             "/Net/Groups/BSY/tools/STILT/fluxes_input/IER_Stuttgart/Europe/N2O.2000.nc",  #full name for emission file (N2O)
                                        "/Net/Groups/BSY/tools/STILT/fluxes_input/IER_Stuttgart/Europe/N2O.2000.nc",  #full name for H2 emissions NOT READY YET
                                                    "/Net/Groups/BSY/tools/people/cgerbig/RData/ROAM/Fluxes/BARCAfires/ncdf/CO.barcafire2009.nc", #full name for emission file (cofire)
								"/Net/Groups/BSY/people/ukarst/STILT_prepro/Radon/InGOS_Rn_map_STILT_EU2_2011_v1.0.nc", #full name for Radon emissions
									"/Net/Groups/BSY/people/svardag/STILT_modelling/CH4_other_emissions/CH4_OTHERS.0.125x0.0.083.nc"),#full name for natural CH4 emissions
#    "co2",    "co",    "ch4",   "n2o",   "h2",   "cofire",  "rn" , "ch4nat"
c(   rep("",l.traddco2), 
     "TM3", 
     rep("",l.traddco),
            "MACC",
     rep("",l.traddch4),
                      "TM3",   
     rep("",l.traddn2o),
                                 "",       "",       "",    "TM3" ,  ""),                                                    #inikind, possible values: "climat" (Gerbig et al. (2003)),
                                                                                                        #"CT" (CarbonTracker), "TM3", "LMDZ", "" (zero boundary)
c(rep("",l.traddco2),   #not modify		  	              #full names for initial & boundary cond. files; select "" for climatological LBC ("climat")
       "/User/homes/croeden/public_html/download-CO2-3D/INVERSION/OUTPUT/s96_v3.6_mix_2013.nc",		#full name for initial & boundary cond. file (CO2); functions accessing files automatically get correct year
  rep("",l.traddco),   #not modify                   							  	
          "/Net/Groups/BSY/tools/STILT/fluxes_input/MACC_CO/macc_CO_2010_01.nc",                        #full name for initial & boundary cond. file (CO); functions accessing files automatically get correct year and month
  rep("",l.traddch4),   #not modify
                    "/Net/Groups/BSY/tools/STILT/fluxes_input/TM3_CH4_nc/mu1.0_070_mix_2011.nc",        #full name for initial & boundary cond. file (CH4); functions accessing files automatically get correct year
  rep("",l.traddn2o),   #not modify
                               "",                   							#full name for initial & boundary cond. file (N2O) 
                                         "",                                                          	#full name for initial & boundary cond. file (H2)
                                                   "",							#full name for initial & boundary cond. file (COfire)
							     "/Net/Groups/BSY/tools/STILT/fluxes_input/TM3_Radon_nc/mix_2011.nc",#full name for initial & boundary cond. file (Rn)
								     "")) 				#full name for initial & boundary cond. file (ch4nat)

#################################################
#flag and path for annual emission scaling according to BP statistics
BPtf<-TRUE     # Want to use interannually varying emissions according to BP fuel statistics?
if(BPtf)BPpath<-edgar.path #now annual maps based on BP statistics written to same directory as EDGAR maps   
#NOTE: Venting and Flaring emissions are for "oil" and "gas" combined, but annual changes are taken from "oil"
############################################################################################################
############################################################################################################
#END OF SECTION TO SET TRACERS FOR USE IN STILT FORWARD SIMULATIONS
############################################################################################################
############################################################################################################


############################################################################################################
#time factors that we have used so far in COFFEE:
############################################################################################################
if(!file.exists(paste("/Net/Groups/BSY/data/EDGAR_4.3_4.2FT2010_categ/RData/edgar4.time.rdata",sep=""))){
  source(paste("/Net/Groups/BSY/data/EDGAR_4.3_4.2FT2010_categ/Rsource/edgar4.time_activity_conversion.r",sep=""))
  save(list=c("daily_var.n", "daily_var.n.sh","daily_var.l","daily_var.l.sh", "weekly_var", "hourly_var"),file=paste(edgar.path,"RData/edgar4.time.rdata",sep=""))
} else { load(file=paste("/Net/Groups/BSY/data/EDGAR_4.3_4.2FT2010_categ/RData/edgar4.time.rdata",sep=""))}

#map ipcc categories to those 12 categories used in time factors
#   1A1_ENE      Energy                                    --> 1
#   1A2_2	     Industry combustion and process emissions --> 5
#   1A3a	     aviation					     --> 9
#   1A3b_c_e     ground transport                          --> 7 
#   1A3d_SHIP    shipping					     --> 8
#   1A3d1        international shipping			     --> 8
#   1A4_5	     residential and other combustion          --> 3
#   1B           fuel production 				     --> 4
#   3            solvents					     --> 6
#   4_but_4E     agriculture without savannah burning      --> 10
#   6A_6C        solid waste                               --> 4 (equ. edgar 3 w-cats)
#   6B	     wastewater handling			     --> 4 (equ. edgar 3 w-cats)
#   7		     other(fossil fuel fires, other an.sources)--> 4 (constant) ####DIFFERENT FROM COFFEE!!! There it was 3 (but small overall contribution)

#get ipcc cats from tracer list
ipcc<-matrix(unlist(strsplit(tr.edg,".",fixed=TRUE)),ncol=3,byrow=T)[,2]
  ipcc2time<-(1:length(ipcc))*NA
  names(ipcc2time)<-ipcc
  #i<-5;par(mfrow=c(2,1));plot(daily_var.n[,i]);plot(hourly_var[,i])
  ipcc2time[names(ipcc2time)=="1A1a"]<-1   #(Power industry)
  ipcc2time[names(ipcc2time)=="1A1bcr"]<-4 #other transformation and non-energy use; assume constant
  ipcc2time[names(ipcc2time)=="1A2+6CD"]<-5 #(industrial combustion (non-power, but including waste incineration))
  ipcc2time[names(ipcc2time)=="1A3a+1C1"]<-9 #(international and domestic aviation)
  ipcc2time[names(ipcc2time)=="1A3b"]<-7 #(road transport)
  ipcc2time[names(ipcc2time)=="1A3ce"]<-8 #(non-road ground transport) (assume constant as ship transport)
  ipcc2time[names(ipcc2time)=="1A3d+1C2"]<-8 #(inland waterways and shipping)
  ipcc2time[names(ipcc2time)=="1A4"]<-3 #(buildings (residential, commercial, services, agriculture,fishing))
  ipcc2time[names(ipcc2time)=="1B1"]<-4 #(fugitive emissions from solid fuels production) assume constant
  ipcc2time[names(ipcc2time)=="1B2abc"]<-4 #assume constant
  ipcc2time[names(ipcc2time)=="1B2ac"]<-4 #(Oil production and distribution and flaring) assume constant
  ipcc2time[names(ipcc2time)=="1B2b"]<-4 #(gas production and distribution (excluding flaring)) assume constant
  ipcc2time[names(ipcc2time)=="2A"]<-4 #(process emissions of non-metallic minerals industry (cement, lime)) assume constant
  ipcc2time[names(ipcc2time)=="2BEFG+3"]<-5 #chem. industry + solvents; use industry
  ipcc2time[names(ipcc2time)=="2C"]<-5 #process emissions metal industry
  ipcc2time[names(ipcc2time)=="4A"]<-10 #(enteric fermentation in agriculture)
  ipcc2time[names(ipcc2time)=="4B"]<-10 #(manure management in agriculture)
  ipcc2time[names(ipcc2time)=="4C"]<-10 #(rice cultivation in agriculture)
  ipcc2time[names(ipcc2time)=="4D"]<-10 #(agricultural soil emissions)
  ipcc2time[names(ipcc2time)=="4F"]<-10 #(agricultural waste burning)
  ipcc2time[names(ipcc2time)=="6A"]<-4 #(solid waste disposal in landfills) assume constant
  ipcc2time[names(ipcc2time)=="6B"]<-4 #(wastewater treatment) assume constant
  ipcc2time[names(ipcc2time)=="7A"]<-4 #(foss fuel fires) assume constant
  ipcc2time[names(ipcc2time)=="7BC"]<-4 #("other" according to excel table from Greet, so unclear) assume constant
