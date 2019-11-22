#FALL 2019, DATA ANALYSIS WITH LAN NGUYEN EXAMINING SAR DATA FROM VERMILLION AREA
#GOAL: CLEAN UP / ORGANIZE DATA

rm(list=ls()) #Clear workspace

#Packages
library(tidyr)
library(dplyr)
library(ggplot2)
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
  ggsave(paste0('./Figures/Field',i,'SAR.png'),p1,width=8,height=8)

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
  ggsave(paste0('./Figures/Field',i,'NDVI.png'),p1,width=8,height=8)
  
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

