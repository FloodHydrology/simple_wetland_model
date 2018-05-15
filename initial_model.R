##################################################################################
#Name: Simple Wetland Model
#Coder: C. Nathan Jones
#Date: 5/14/2018
#Purpose: Explore evapotranspiration effects on mixing in wetlands 
##################################################################################

##################################################################################
#Step 1: Setup Workspace ---------------------------------------------------------
##################################################################################
#Clear Memory
rm(list=ls(all=TRUE))

#Set Working Directory
wd<-"/nfs/njones-data/Research Projects/SimpleWetlandModel/Initial_Model"
setwd(paste0(wd))

#add appropriate libarary
library(raster)   #spatial analysis
library(rgdal)    #spatial analysis
library(rgeos)    #spatial analysis
library(maptools) #spatial analysis
library(plyr)     #data processing
library(dplyr)    #data processing
library(Evapotranspiration)
library(dplyr)
library(magrittr)

#Obtain Model Inputs (we'll need to redo this at some point...)
load(paste0(wd,"/Backup/inputs.RData")) #load inputs from previous model 
wd<-"/nfs/njones-data/Research Projects/SimpleWetlandModel/Initial_Model"
setwd(paste0(wd))

####################################################################################
# Step 2: Define input variables----------------------------------------------------
####################################################################################
#Define dynamic parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create vector of dynamic parameters names
dynamic_variables<-c( 
  "y_w",        #Surface water elevation(mm)
  "GW_local",   #Local ground water exchange (mm/day)
  "ET",         #ET from the water table (mm/day)
  "GW_bf",      #GW lost to "Baseflow" (mm/day)
  "runoff",     #saturation excess runoff (mm/day)
  "spill"       #surface water export from wetland (mm/day)
)

#Create individualvectors
for(i in 1:length(dynamic_variables)){
  assign(paste0(dynamic_variables[i], ".VAR"), 
         matrix(0, ncol=2, 
                nrow=length(pet.VAR), 
                dimnames=list(seq(1,length(pet.VAR), 1), c("upland","wetland"))))
}


#Define wetland parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#create input variables
giw.INFO<-c("area","invert")

#Create giw.INFO matrix
giw.INFO<-matrix(0, nrow=1, ncol=length(giw.INFO), dimnames = list(seq(1,n.wetlands,1), c(giw.INFO)))

#Assign initial Values
giw.INFO[,"area"]<-    16794*(1000^2)  #mm^2
giw.INFO[,"invert"]<-          -1000           #mm

#Define upland parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
upland.INFO<-c("area",#geometric characteristics
             "n","y_cl", "k_sat", "RD","Sy","kb" #soil characteristics
)

#Create giw.INFO matrix
upland.INFO<-matrix(0, nrow=1, ncol=length(upland.INFO), dimnames = list(c(1), c(upland.INFO)))

#Populate land.INFO matrix (lenghth units in mm)
upland.INFO[,"area"]<-  giw.INFO[,"area"]*1.5 #mm^2              #area in mm^2
upland.INFO[,"n"]<-     soils$n                             #porisity
upland.INFO[,"y_cl"]<-  -1*soils$y_cl                       #confining layer depth (mm) from SSURGO
upland.INFO[,"k_sat"]<- -soils$ksat*24                      #saturated condcuctivity (mm/day)
upland.INFO[,"RD"]<-    -soils$y_rd                         #Rooting Depth (mm)
upland.INFO[,"Sy"]<-    soils$Sy
upland.INFO[,"kb"]<-    0.046                               #Defined fromo literatture. (Going to odo with Gauge analysis eventually)


####################################################################################
# Step 3:Define model function------------------------------------------------------
####################################################################################
fun<-function(day){
  #Calculate upland water balance~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Calculate GW_local flux [positive is flux into the given control volume]
  r_ws<-(upland.INFO[,"area"]/pi)^0.5
  r_w<- (giw.INFO[,"area"]/pi)^0.5
  GW_local<- pi*upland.INFO[,"k_sat"]*((y_w.VAR[day, "wetland"]^2)-(y_w.VAR[day, "upland"]^2))/log(r_ws/r_w)
  GW_local.VAR[day,"upland"]<<- -1*GW_local/(upland.INFO[,"area"]-giw.INFO[,"area"])
  GW_local.VAR[day,"wetland"]<<- GW_local/giw.INFO[,"area"]
  
  #Flux out of the watershed (e.g. baseflow from SWAT manual)
  GW_bf.VAR[day,"land"]<<-ifelse(day==1,
                                 1,
                                 ifelse(GW_bf.VAR[day-1,"upland"]*exp(-land.INFO[,"kb"])+(-ET.VAR[day,"upland"]+GW_local.VAR[day,"upland"])*(1-exp(-land.INFO[,"kb"]))<0,
                                        GW_bf.VAR[day-1,"upland"]*exp(-land.INFO[,"kb"])+(-ET.VAR[day,"upland"]+GW_local.VAR[day,"upland"])*(1-exp(-land.INFO[,"kb"])),
                                        0))
    
  #Calculate change in upland stage for next time step
  dy<-(precip.VAR[day]-ET.VAR[day]+GW_local.VAR[day,"upland"]-GW_bf.VAR[day,"upland"])/upland.INFO[,"Sy"]
  
  #Calucate y_w for upland and runoff into wetland [n+1]
  if((y_w.VAR[day,"upland"]+dy)<=0){
    y_w.VAR[day+1,"upland"]<<-y_w.VAR[day,"upland"]+dy
  }else{
    y_w.VAR[day+1,"upland"]<<-0
    runoff.VAR[day, "upland"]<<-y_w.VAR[day,"upland"]+dy
  }
  
  #Calculate wetland water balance~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Calculate GW_local flux [positive is flux into the given control volume]
  GW_local.VAR[day,"wetland"]<<- GW_local/giw.INFO[,"area"]
  
  #Calculate change in wetland stage
  dy<-precip.VAR[day]-pet.VAR[day]+GW_local.VAR[day,"wetland"]+runoff.VAR[day, "upland"]
     
  #Calucate y_w for upland and export out of wetland [n+1]
  if((y_w.VAR[day,"wetland"]+dy)<=giw.INFO[,"invert"]){
    y_w.VAR[day+1,"wetland"]<-giw.INFO[,"invert"]
  }else{
    if((y_w.VAR[day,"wetland"]+dy)<=0){
      y_w.VAR[day+1,"wetland"]<<-y_w.VAR[day,"wetland"]+dy
    }else{
      y_w.VAR[day+1,"wetland"]<<-0
      runoff.VAR[day, "wetland"]<<-y_w.VAR[day,"wetland"]+dy
    }   
  }
}



####################################################################################
# Step 4: Run initial function------------------------------------------------------
####################################################################################
#Initial Conditions
GW_bf.VAR[1,"upland"]<-10