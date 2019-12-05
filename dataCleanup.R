#FALL 2019, DATA ANALYSIS WITH LAN NGUYEN EXAMINING SAR DATA FROM VERMILLION AREA
#GOAL: CLEAN UP / ORGANIZE DATA

rm(list=ls()) #Clear workspace

#Packages
library(tidyverse)
library(lubridate)
library(beepr)
theme_set(theme_classic())

# Load and reorganize data from Dataset 1 ---------------------------------------------------------------
#Dataset1: small dataset from random locations around Vermillion
rm(list=ls())
setwd("~/Documents/soil_moisture_analysis/Dataset1")

#Day of year
DOY <- read.table('DOY.txt') #Days of year, starting from Jan 1 = 1
#Local incident angle
LIA <- read.table('LIA.txt') #Row = location, col = observation
#Backscatter
SAR <- read.table('SAR.txt') #Same indexing as above, 0 -> NA

rownames(SAR) <- paste0('S',1:100)
rownames(LIA) <- paste0('S',1:100)

dat <- data.frame(type='SAR',DOY=DOY[,1],t(SAR)) %>% bind_rows(data.frame(type='LIA',DOY=DOY[,1],t(LIA))) %>% 
  pivot_longer(S1:S100,names_to='site',values_to='meas') %>% 
  filter(meas!=0,meas!=-9999) %>% 
  unite(siteDOY,site,DOY,sep='_') %>% 
  group_by(type,siteDOY) %>% summarize(meas=mean(meas,na.rm=TRUE)) %>% ungroup() %>% 
  pivot_wider(names_from=type,values_from=meas) %>% 
  separate(siteDOY,c('site','DOY'),sep='_') %>% mutate(site=factor(site),DOY=as.numeric(DOY)) %>% 
  arrange(site,DOY) %>% mutate(LIA=LIA/100) %>% na.omit()
save(dat,file='Rdata16.Rdata')



# Load and reorganize data from Dataset 2 ----------------------------------
#Dataset2: larger dataset from 16 individual fields, taken over 175 days
rm(list=ls())
setwd("~/Documents/soil_moisture_analysis/Dataset2")

#Problem: SAR and NDVI data are sometimes from different days

#Accessory data
sar_doy <- read.table('S1_DOY.txt') #Day of year for SAR data 
ndvi_doy <- read.table('VI_DOY.txt') #Day of year for NDVI data
dim(sar_doy); dim(ndvi_doy)
#NOTE: ndvi_doy and ndvi are from 2 years. First year should be discarded

col <- read.table('COL.txt') #Column index
row <- read.table('ROW.txt') #Row index
dim(col); dim(row)

#File numbers to append to text
numFiles <- gsub('.txt','',gsub('LIA_','',list.files(pattern='LIA')))
#Empty list
dat <- vector(mode = "list", length = length(numFiles))

for(i in 1:length(numFiles)){
  use <- numFiles[i] #Suffix of file names
  
  #Read only data from first field
  lia <- read.table(paste0('LIA_',use,'.txt'),na.strings=c('32767','0')) 
  sar <- read.table(paste0('SAR_',use,'.txt')) %>% 
    mutate_all(function(x) ifelse(x>=0,NA,x))
  ndvi <- read.table(paste0('VI_',use,'.txt'),na.strings = '-9999') #-9999 = NA
  #This is roughly twice the length of the others because 2 platforms were used to get NDVI. Should amalgamate to an average over a single day.
  
  # coords <- read.table(paste0('XYZ_',use,'.txt')) #Actual lat lon values
  dim(lia); dim(sar); dim(ndvi)
  
  #Rows are all individual days, corresponding to doy
  #Columns are individual points at a field (10201 = 101 x 101), with rows first
  
  #LIA 
  #Create matrices of matching dimensionality, just to make sure conversion is correct
  doymat <- matrix(rep(sar_doy[,1],each=ncol(lia)),ncol=ncol(sar),nrow=nrow(sar),byrow=T)
  colmat <- matrix(rep(as.vector(t(as.matrix(col))),nrow(sar)),ncol=ncol(sar),nrow=nrow(sar),byrow=T)
  rowmat <- matrix(rep(as.vector(t(as.matrix(row))),nrow(sar)),ncol=ncol(sar),nrow=nrow(sar),byrow=T)
  dim(doymat); dim(colmat); dim(rowmat)
  
  #Combine into single dataframe
  dat[[i]] <- data.frame(lia=as.vector(as.matrix(lia)),
                    sar=as.vector(as.matrix(sar)),
                    colnum=as.vector(colmat),rownum=as.vector(rowmat),
                    doy=as.vector(doymat)) %>% 
  filter(!is.na(lia)|!is.na(sar)) %>% #Strip out values with no measurements
  group_by(colnum,rownum,doy) %>% summarize(lia=mean(lia),sar=mean(sar)) %>% ungroup() %>% #Average measurements on same day
  # mutate(sarScal=scale(sar),liaScal=scale(lia),doyScal=scale(doy)) %>% #Scale sar,lia,doy
  unite(cell,colnum,rownum,remove=F) %>% mutate(cell=factor(cell)) %>% #Create discrete "cell" factor from row/col
  unite(id,cell,doy,sep='_',remove=F) %>% #Create id column for merge
  as.data.frame() 
  
  #Tests - looks OK
  # dat[[i]] %>% filter(cell=='50_50') %>% ggplot(aes(doy,sar))+geom_point()+labs(title=paste('Field',i,'at cell 50_50'))
  # chooseday <- dat[[i]] %>% select(doy) %>% distinct() %>% mutate(diff=abs(doy-190)) %>% filter(diff==min(diff)) %>% slice(1) %>% .$doy #Choose day closest to 190
  # dat[[i]] %>% filter(doy==chooseday) %>% ggplot(aes(colnum,rownum,fill=sar))+geom_raster()+labs(title=paste('Field',i,'SAR on day',chooseday))
  
  p1 <- dat[[i]] %>% ggplot(aes(colnum,rownum,fill=sar))+ geom_raster()+facet_wrap(~doy)+ labs(title=paste('Field',i,'SAR')) +scale_y_reverse()
  ggsave(paste0('./DataFigures/Field',i,'SAR.png'),p1,width=8,height=8)

  #NDVI
  #Create matrices of matching dimensionality, just to make sure conversion is correct
  doymat <- matrix(rep(ndvi_doy[,1],each=ncol(lia)),ncol=ncol(ndvi),nrow=nrow(ndvi),byrow=T)
  colmat <- matrix(rep(as.vector(t(as.matrix(col))),nrow(ndvi)),ncol=ncol(ndvi),nrow=nrow(ndvi),byrow=T)
  rowmat <- matrix(rep(as.vector(t(as.matrix(row))),nrow(ndvi)),ncol=ncol(ndvi),nrow=nrow(ndvi),byrow=T)
  dim(doymat); dim(colmat); dim(rowmat)
  
  #Combine into dataframe
  dat2 <- data.frame(ndvi=as.vector(as.matrix(ndvi)),
                         colnum=as.vector(colmat),rownum=as.vector(rowmat),
                         doy=as.vector(doymat)) %>%
    filter(!is.na(ndvi)) %>% 
    group_by(colnum,rownum,doy) %>% summarize(ndvi=mean(ndvi,na.rm=T)) %>% ungroup() %>% #Average measurements on same day
    group_by(doy) %>% mutate(percmissing=1-(n()/ncol(doymat))) %>% ungroup() %>% 
    filter(percmissing<0.20) %>% select(-percmissing) %>% #Days with more than 20% of pixels missing are excluded
    mutate(ndviScal=scale(ndvi)) %>% #Scale ndvi
    unite(cell,colnum,rownum,remove=F) %>% mutate(cell=factor(cell)) %>%  #Create discrete "cell" factor from row/col
    unite(id,cell,doy,sep='_',remove=F) %>% #Create id column for merge
    as.data.frame()
  
  #Test -looks OK
  # dat2 %>% filter(cell=='50_50') %>% ggplot(aes(doy,ndvi))+geom_point()+labs(title=paste('Field',i,'at cell 50_50'))
  # chooseday <- dat2 %>% select(doy) %>% distinct() %>% mutate(diff=abs(doy-190)) %>% filter(diff==min(diff)) %>% slice(1) %>% .$doy #Choose first day closest to 190
  # dat2 %>% filter(doy==chooseday) %>%  ggplot(aes(rownum,colnum,fill=ndvi))+geom_raster()+ facet_wrap(~doy) + labs(title=paste('Field',i,'on day',chooseday)) 
  # dat2 %>% filter(doy>=179 & doy<=211) %>% ggplot(aes(rownum,colnum,fill=ndvi))+geom_raster()+facet_wrap(~doy)
  
  p1 <- dat2 %>% ggplot(aes(colnum,rownum,fill=ndvi))+ geom_raster()+facet_wrap(~doy)+ labs(title=paste('Field',i,'NDVI'))+scale_y_reverse()
  ggsave(paste0('./DataFigures/Field',i,'NDVI.png'),p1,width=8,height=8)
  
  #Merge NDVI with SAR, LIA data - full_join keeps all row IDs
  # dat[[i]] <- 
  # temp1 <- 
  
  dat[[i]] <- dat[[i]] %>% 
    select(id,lia,sar) %>% 
    full_join(select(dat2,id,ndvi),by='id') %>% 
    separate(id,c('colnum','rownum','doy'),sep='_',convert=T) %>% 
    unite(cell,colnum,rownum,remove=F) %>%
    arrange(colnum,rownum,doy)
  
  rm(doymat,colmat,rowmat,dat2) #Cleanup extra matrices
  gc()
  print(paste0('Finished field ',i))
}
names(dat) <- paste0('site',numFiles)
beep(1)
  
save(dat,file='allFieldData.Rdata') #Save as R file



# Load and reorganize data from Dataset 3 --------------------------------

rm(list=ls())
setwd("~/Documents/soil_moisture_analysis/Dataset3/TimestampData")

#Met station data
metFiles <- list.files(pattern='acisData')

metDat <- vector(mode = "list", length = length(metFiles)) #Empty list

headers <- c('Station','time','AirTemp','AirTempFlag','AirTempComment','Precip','PrecipFlag','PrecipComment',
             'SoilTemp','SoilTempFlag','SoilTempComment','SoilMoisture','SoilMoistureFlag','SoilMoistureComment')

for(i in 1:length(metFiles)){ #Read all metdata csvs
  metDat[[i]] <- read.csv(metFiles[i],header=T,col.names=headers,stringsAsFactors=F) %>% 
    select(-contains('Comment'),-contains('Flag'),-Station)
}
metDat <- do.call('rbind',metDat) %>% #Bind together dataframes
  mutate(time=dmy_hm(time))
rm(headers,i,metFiles)

metDat %>% pivot_longer(cols=AirTemp:SoilMoisture) %>%
  ggplot(aes(time,value))+geom_line()+
  facet_wrap(~name,ncol=1,scales='free_y') #Looks ok

#LIA
lia <- read.table('localIncidenceAngle.txt',na.strings = '-9999',stringsAsFactors = F, header = F) %>% 
  setNames(.,c(paste0('cell','_',rep(1:7,each=7),'_',rep(1:7,7)))) #Read in data
liaLabs <- read.table('localIncidenceAngle_fnames.txt',header=F,stringsAsFactors=F) %>% #Get time stamp 
  mutate(V1=str_extract(V1,regex('[:digit:]{8}T[:digit:]{6}'))) %>%  #Select text with date/time info
  mutate(time=ymd_hms(V1)) %>% select(-V1) #Extract date/time
lia <- bind_cols(liaLabs,lia); rm(liaLabs)

ggplot(lia,aes(x=time,cell_1_1))+geom_point()+geom_line()

#NDVI
ndvi <- read.table('NDVI.txt',stringsAsFactors = F, header = F,na.strings='-9999') %>% 
  setNames(.,c(paste0('cell','_',rep(1:7,each=7),'_',rep(1:7,7)))) #Read in data
ndviLabs <- read.table('NDVI_fnames.txt',header=F,stringsAsFactors=F) %>% #Get time stamp
  mutate(V1=str_extract(V1,regex('[:digit:]{8}T[:digit:]{6}'))) %>%  #Select text with date/time info
  mutate(time=ymd_hms(V1)) %>% select(-V1) #Extract date/time
ndvi <- bind_cols(ndviLabs,ndvi); rm(ndviLabs)

ndvi %>% filter(complete.cases(.)) %>% 
  ggplot(aes(x=time,cell_1_1))+geom_point()+geom_line()

ndvi %>% pivot_longer(cols=contains('cell')) %>% filter(!is.na(value)) %>% 
  group_by(time) %>% summarize(value=median(value)) %>% ungroup() %>% 
  ggplot(aes(x=time,value))+geom_point()+geom_line()

#SAR (VV and VH readings)
sarVV <- read.table('SAR_VV.txt',stringsAsFactors = F, header = F,na.strings='-9999') %>% 
  setNames(.,c(paste0('cell','_',rep(1:7,each=7),'_',rep(1:7,7)))) #Read in data
sarVVLabs <- read.table('SAR_VV_fnames.txt',header=F,stringsAsFactors=F) %>% #Get time stamp
  mutate(V1=str_extract(V1,regex('[:digit:]{8}T[:digit:]{6}'))) %>%  #Select text with date/time info
  mutate(time=ymd_hms(V1)) %>% select(-V1) #Extract date/time
sarVV <- bind_cols(sarVVLabs,sarVV); rm(sarVVLabs)

sarVH <- read.table('SAR_VH.txt',stringsAsFactors = F, header = F,na.strings='-9999') %>% 
  setNames(.,c(paste0('cell','_',rep(1:7,each=7),'_',rep(1:7,7)))) #Read in data
sarVHLabs <- read.table('SAR_VH_fnames.txt',header=F,stringsAsFactors=F) %>% #Get time stamp
  mutate(V1=str_extract(V1,regex('[:digit:]{8}T[:digit:]{6}'))) %>%  #Select text with date/time info
  mutate(time=ymd_hms(V1)) %>% select(-V1) #Extract date/time
sarVH <- bind_cols(sarVHLabs,sarVH); rm(sarVHLabs)

#VH and VV highly strongly correlated
bind_rows(mutate(sarVV,type='VV'),mutate(sarVH,type='VH')) %>%
  pivot_longer(cols=contains('cell')) %>% filter(!is.na(value)) %>% 
  group_by(time,type) %>% summarize(value=median(value)) %>% ungroup() %>% 
  pivot_wider(names_from='type',values_from='value') %>% mutate(ratio=VH/VV) %>% 
  pivot_longer(cols=VH:ratio) %>% 
  ggplot(aes(x=time,value))+geom_point()+geom_line()+
  geom_smooth(method='gam',formula=y~s(x))+
  facet_wrap(~name,ncol=1,scales='free_y')

bind_rows(mutate(sarVV,type='VV'),mutate(sarVH,type='VH')) %>%
  pivot_longer(cols=contains('cell')) %>% filter(!is.na(value)) %>% 
  group_by(time,type) %>% summarize(value=median(value)) %>% ungroup() %>% 
  pivot_wider(names_from='type',values_from='value') %>% select(-time) %>% 
  
  ggplot(aes(VV,VH))+geom_point()+geom_smooth(method='lm',se=T)



# #Soil moisture and climate data from stations
# soilwater <- read.csv('ACIS_DailyData_2018.csv',stringsAsFactors = F) %>% 
#   select(-Precip..Accumulated.Source.Flag,-Precip..Source.Flag) %>% 
#   filter(grepl('Vermilion',Station.Name)) %>%  #Filter all station data except Vermillion
#   mutate(Date..Local.Standard.Time.=1:nrow(.)) %>% 
#   setNames(.,c('station','doy','precipAccum','precipAccumNote','precip','precipNote','soilwater','soilwaterFlag','soilwaterComplete')) %>% 
#   select(-station,-precipAccumNote,-precipNote,-soilwaterFlag)
# 
# soilwater %>% select(-soilwaterComplete) %>% 
#   pivot_longer(cols=precipAccum:soilwater) %>% 
#   ggplot(aes(doy,value))+geom_point()+geom_line()+
#   facet_wrap(~name,ncol=1,scales='free_y') #Looks ok
#   
# #LIA
# lia <- read.csv('Vermilion_AGDM_LIA.txt',stringsAsFactors = F,sep=' ',na.strings='32767',header=F) %>%
#   mutate(doy=read.csv('Vermilion_AGDM_S1_DOY.txt',header=F)[,]) %>% #Add day of year
#   setNames(.,c(paste0('cell','_',rep(1:7,each=7),'_',rep(1:7,7)),'doy')) %>% #Fix column names
#   group_by(doy) %>% summarize_all(~mean(.,na.rm=T)) %>% ungroup() %>% #Average over day
#   mutate_all(~ifelse(is.nan(.),NA,.)) 
# 
# ggplot(lia,aes(doy,cell_1_1))+geom_point()+geom_line() #looks ok
# 
# #NDVI 
# ndvi <- read.csv('Vermilion_AGDM_NDVI.txt',sep=' ',stringsAsFactors = F,header=F,na.strings='-9999') %>% 
#   mutate(doy=read.csv('Vermilion_AGDM_S2_DOY.txt',header=F)[,]) %>%  #Add day of year
#   setNames(.,c(paste0('cell','_',rep(1:7,each=7),'_',rep(1:7,7)),'doy')) %>% 
#   group_by(doy) %>% summarize_all(~mean(.,na.rm=T)) %>% ungroup() %>% #Average over day
#   mutate_all(~ifelse(is.nan(.),NA,.)) 
#   
#   ggplot(ndvi,aes(doy,cell_1_1))+geom_point()+geom_line()+geom_smooth(method='gam',formula=y~s(x,k=5)) #looks ok
# 
# #SAR 
# sar <- read.csv('Vermilion_AGDM_SAR.txt',stringsAsFactors = F,sep=' ',na.strings='32767',header=F) %>% 
#   mutate(doy=read.csv('Vermilion_AGDM_S1_DOY.txt',header=F)[,]) %>% #Add day of year
#   setNames(.,c(paste0('cell','_',rep(1:7,each=7),'_',rep(1:7,7)),'doy')) %>% 
#   group_by(doy) %>% summarize_all(~mean(.,na.rm=T)) %>% ungroup() %>% #Average over day
#   mutate_all(~ifelse(is.nan(.),NA,.)) 
#   
# ggplot(sar,aes(doy,cell_1_1))+geom_point() #looks ok
# 
# #Combine into single dataframe
# dat <- bind_rows(mutate(sar,meas='sar'),mutate(ndvi,meas='ndvi'),mutate(lia,meas='lia')) %>%
#   pivot_longer(cols=contains('cell')) %>% rowwise() %>% 
#   transmute(id=gsub('cell',doy,name),meas,value) %>% ungroup() %>% 
#   pivot_wider(names_from = meas, values_from = value) %>% 
#   separate(id,c('doy','row','col'),convert=T) %>%
#   left_join(select(soilwater,doy,precipAccum,precip,soilwater),by='doy')
# 
# dat %>% group_by(doy) %>% 
#   summarize(sar=mean(sar),ndvi=mean(ndvi),lia=mean(lia),precip=first(precip),soilwater=first(soilwater)) %>% 
#   pivot_longer(cols=sar:soilwater) %>% 
#   ggplot(aes(doy,value))+geom_point()+geom_line()+facet_wrap(~name,ncol=1,scales='free_y')
# 
# save(dat,file='metStationDat.Rdata') #Save as R file



  
