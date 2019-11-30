#FALL 2019, DATA ANALYSIS WITH LAN NGUYEN EXAMINING SAR DATA FROM VERMILLION AREA
#GOAL: TRY TO PULL OUT TEMPORAL SIGNAL AND LIA SIGNAL FROM SAR BACKSCATTER DATA

#Packages
rm(list=ls()) #Clear workspace

library(tidyverse)
library(ggplot2)
library(mgcv)
library(lme4)
library(parallel)
library(beepr)
theme_set(theme_classic())


# Dataset 1 ----------------------------------------------------------------
#Dataset1: small dataset from random locations around Vermillion
rm(list=ls())
setwd("~/Documents/soil_moisture_analysis/Dataset1")
load('Rdata1.Rdata')

# Basic plot
dat %>% ggplot(aes(DOY,SAR))+geom_point(alpha=0.3)+
  geom_smooth(method='loess')

#Spline appears to fit well. Would be good to have random intercept/slopes for site
mod1 <- gamm(scale(SAR)~scale(LIA)+s(DOY),random=list(site=~1+LIA),data=dat)
summary(mod1$lme)
print(mod1)
plot(mod1$lme)
plot(mod1$gam,all.terms=T)


#Residuals broken up by site
dat %>% mutate(resid=resid(mod1$lme)) %>%
  ggplot(aes(site,resid))+geom_point()

#Residuals over time
dat %>% mutate(resid=resid(mod1$lme)) %>%
  ggplot(aes(DOY,resid))+geom_point()

summary(mod1$lme)


summary(mod1)
plot(mod2,all.terms=T)

plot(c(coef(mod2)[101],coef(mod2)[101]+coef(mod2)[102:200]),xlab='Site',ylab='LIA slope')
mod3 <- gam(SAR~LIA+s(DOY),data=dat,subset=site=='S1')
plot(mod3,all.terms=T)


# Dataset 2 ---------------------------------------------------------------
#Dataset2: larger dataset from 16 individual fields, taken over 175 days
#To start, only using a single field's worth of data:
rm(list=ls())
setwd("~/Documents/soil_moisture_analysis/Dataset2")
load('allFieldData.Rdata') #Load list of dataframes

library(foreach)
library(doParallel)

cl <- makeCluster(4) #8 cores
registerDoParallel(cl)

# foreach(i=1:length(dat),.packages=c('dplyr','tidyr','mgcv')) %dopar% {
for(i in 1:length(dat)){
  shortDat <- dat[[i]] #Field i
  #Step 1: create model to fill in holes in NDVI data
  
  # #Plot of SAR and NDVI as a function of DOY at 3 cells
  # shortDat %>% filter(cell=='0_0'|cell=='50_50'|cell=='100_100') %>% 
  #   select(-colnum,-rownum,-lia) %>% 
  #   pivot_longer(cols=sar:ndvi) %>% 
  #   ggplot(aes(doy,value))+geom_point()+
  #   geom_smooth(method='loess',se=F)+
  #   facet_grid(name~cell,scales='free_y')
  
  
  #Average over the entire field  
  # shortDat %>% 
  #   select(-colnum,-rownum) %>% pivot_longer(cols=lia:ndvi) %>% 
  #   group_by(doy,name) %>% summarize(mean=mean(value,na.rm=T),sd=sd(value,na.rm = T)) %>% 
  #   ungroup() %>% filter(!is.nan(mean)) %>% 
  #   ggplot(aes(doy,mean))+
  #   geom_pointrange(aes(ymax=mean+sd,ymin=mean-sd))+
  #   geom_line()+
  #   facet_wrap(~name,scales='free_y',ncol=1)+
  #   labs(x='DOY',y='Mean value',title=paste('Field ',i,'mean values'))
  
  load(paste0('./ModelResults/ndviMod1Field',i,'.Rdata'))
  # ndviMod1 <- shortDat %>% #Tensor smooth across space & time
  #   mutate(ndvi=scale(ndvi)) %>%
  #   bam(ndvi~ti(colnum,rownum,doy,bs='cr',k=12)+ti(colnum,rownum,bs='cr',k=12)+s(doy,k=12),data=.,method='REML')
  # save(ndviMod1,file=paste0('./ModelResults/ndviMod1Field',i,'.Rdata'))

  # #Save output figures
  # png(paste0('./ModelResults/ndviMod1Field',i,'.png'),height=1200,width=600)
  # par(mfrow=c(3,1)); plot(ndviMod1,se=F,residuals=T)
  # dev.off()
  # 
  # png(paste0('./ModelResults/ndviMod1Field',i,'Resid.png'),height=800,width=800)
  # par(mfrow=c(2,2)); gam.check(ndviMod1); abline(0,1,col='red')
  # dev.off()
  # summary(ndviMod1) #Not the best in terms of residual and k' checks, but probably the best we're going to get at this point. Likely driven down by cloudy days.
  
  #NOTES:   
  #Also tried tensor smooth for space, with independent spline for time: ndvi~te(colnum,rownum,bs='cr',k=12)+s(doy,k=20)
  #Results are similar, but with worse AIC and slightly lower R2. However, most of the signal taken up by time component, so leaving out some spatial variation doesn't matter much
  #All the residuals seem to be fairly heavy-tailed. Perhaps this has something to do with the day-to-day varation?
  #Next step: see what kind of smoothing is appropriate for a single day's worth of NDVI data. This may inform the spatio-temporal smoothing.
  #Result: high-res smoothing takes a very long time, but looks OK in terms of predictions. Further increasing k enhances resolution of predictions. 
  #Problems exist trying to model small "holes" in the crop (possibly areas that didn't get filled in?)
  #If just using this to "fill in" values, probably won't hurt to use ndviMod1
  
  #Gets predictions from model, fills in missing values, and truncates to temporal range of NDVI data
  shortDat <- shortDat %>%
    mutate(pred=predict(ndviMod1,newdata=shortDat)) %>%
    mutate(pred=mean(ndvi,na.rm=T)+(pred*sd(ndvi,na.rm = T))) %>% #Re-scales predicted values
    mutate(ndviNA=is.na(ndvi),ndvi=ifelse(ndviNA,pred,ndvi)) %>% #Fills in values that have NAs
    #Keep going here
    group_by(doy) %>% mutate(allNA=!any(!ndviNA)) %>% ungroup() %>% #Identify days with only NA values for NDVI
    mutate(mindoy=min(doy[!allNA]),maxdoy=max(doy[!allNA])) %>% #Get max and min values for days that had NDVI values
    filter(doy>=mindoy&doy<=maxdoy) %>% #Removes days outside of the date range
    filter(!is.na(sar)&!is.na(lia)) %>% #Remove days without LIA or SAR
    select(-pred,-allNA:-maxdoy) #Cleanup columns
  
  # #Looks OK
  # shortDat %>% group_by(doy) %>% 
  #   summarize(ndvi=mean(ndvi),imputed=!any(!ndviNA)) %>% 
  #   ggplot(aes(doy,ndvi,col=imputed))+geom_point()+geom_line(aes(group=1))
  
  #Step 2: fit SAR model; basic framework -> gam(sar ~ ndvi + lia + ti(x,y,doy) + te(x,y,doy) + s(doy))
  # load(paste0('./ModelResults/sarMod1Field',i,'.Rdata'))
  sarMod1 <- bam(sar~lia+ndvi+s(doy)+ti(colnum,rownum)+ti(colnum,rownum,doy),
                 data=mutate(shortDat,sar=scale(sar)))
  # save(sarMod1,file=paste0('./ModelResults/sarMod1Field',i,'.Rdata'))
  
  png(paste0('./ModelResults/sarMod1Field',i,'.png'),height=1200,width=600)
  par(mfrow=c(3,2)); plot(sarMod1,se=F,all.terms=T,residuals=T)
  dev.off()
  
  png(paste0('./ModelResults/sarMod1Field',i,'Resid.png'),height=800,width=800)
  par(mfrow=c(2,2)); gam.check(sarMod1); abline(0,1,col='red')
  dev.off()
  
  rm(shortDat,ndviMod1,sarMod1); gc() #Cleanup
  print(paste0('Finished field ',i))
}
beepr::beep(2)  
# stopCluster(cl) #Appears to work

#Summary: "predicted soil moisture" (i.e. spatio-temporal patterns | NDVI, LIA) seems reasonable, but we have no way of telling without looking at actual soil moisture data. Next set of analyses should deal with this.


# Dataset 3 ---------------------------------------------------------------

#Dataset 3: set of SAR, NDVI, LIA observations around climate stations, where soil moisture was recorded as well

rm(list=ls())
setwd("~/Documents/soil_moisture_analysis/Dataset3")

load('metStationDat.Rdata')

dat <- dat %>% 
  mutate(uprNDVI=max(doy[!is.na(ndvi)&!is.na(sar)]),lwrNDVI=min(doy[!is.na(ndvi)&!is.na(sar)])) %>% #Filter to range of NDVI amd SAR data
  filter(doy>=lwrNDVI & doy<=uprNDVI) %>% select(-uprNDVI,-lwrNDVI) %>% group_by(row,col) %>% 
  mutate(ndviInterp=approx(doy,ndvi,doy)$y) %>% ungroup() %>% #Interpolate NDVI values
  unite(cell,row,col,remove = F) 
  
dat %>% filter(row>4,col>4) %>% 
  ggplot(aes(x=doy,y=ndvi))+ geom_point()+
  geom_line(aes(y=ndviInterp),col='red')+
  facet_grid(row~col) #Looks OK

#Mixed effects model of SAR, using cell as a random intercept (singularity occurs if random slopes used)
#Variables linearly scaled. Log scaling provides no better fit
dat1 <- dat %>% mutate_at(vars(sar,lia,ndviInterp,soilwater),scale) %>%  #Scale variables
  filter_at(vars(sar=sar,lia,ndviInterp,soilwater),~!is.na(.)) #Get rid of NA values
# dat1 %>% select(sar,lia,soilwater,ndviInterp,precip) %>% pairs(.) #Look at distribution of data

sarMod1 <- dat1 %>% lmer(sar~ndviInterp+soilwater+lia+(1|cell),data=.,REML=F)
#Continue here 

{
  par(mfrow=c(3,1)); qqnorm(resid(sarMod1)); qqline(resid(sarMod1)); plot(fitted(sarMod1),resid(sarMod1)); plot(predict(sarMod1),predict(sarMod1)+resid(sarMod1),ylab='actual') 
  abline(0,1,col='red')
}
summary(sarMod1) #LIA


#Residuals appear to increase over time
dat1 %>% mutate(res=resid(sarMod1)) %>%
  select(doy,res) %>% 
  ggplot(aes(doy,res))+
  geom_point()+geom_smooth(method='gam',formula=y~s(x,k=8))

#Residuals appear spatially structured 
dat1 %>% mutate(res=resid(sarMod1)) %>% 
  group_by(cell) %>% summarize(res=mean(res)) %>% 
  separate(cell,c('row','col'),sep='_') %>% 
  ggplot(aes(row,col,fill=res))+geom_raster()

#Residuals mainly vary across time of the season
resMod1 <- dat1 %>% mutate(res=resid(sarMod1)) %>% 
  separate(cell,c('row','col'),sep='_',convert = T) %>% 
  select(doy,row,col,res) %>% as.data.frame() %>% 
  gam(res~s(doy)+ti(row,col),data=.)
summary(resMod1)
par(mfrow=c(2,1)); plot(resMod1,se=F,residuals=T)
par(mfrow=c(2,2)); gam.check(resMod1); abline(0,1,col='red')

#Do residuals vary with other factors?
dat1 %>% mutate(res=resid(sarMod1)) %>% 
  select(-cell,-ndvi,-row,-col,-sar,-precipAccum) %>% 
  pivot_longer(cols=doy:ndviInterp) %>% 
  ggplot(aes(value,res))+geom_point()+
  facet_wrap(~name,scales='free')

#Partial effect of LIA
dat1 %>% mutate(soilwater=0,ndviInterp=0) %>% #Marginalize across soilwater/ndviInterp
  mutate(pred=predict(sarMod1,newdata=.,re.form=~0),res=resid(sarMod1)) %>% 
  ggplot(aes(lia,pred))+ 
  geom_point(aes(y=pred+res))+geom_line(col='red',size=1)+
  labs(x='LIA',y='SAR | NDVI,Soilwater',title = 'Partial effect of LIA')

#Partial effect of NDVI
dat1 %>% mutate_at(vars(sar,lia,ndviInterp,soilwater),scale) %>% 
  mutate(soilwater=0,lia=0) %>% 
  mutate(pred=predict(sarMod1,newdata=.,re.form=~0),res=resid(sarMod1)) %>% 
  ggplot(aes(ndviInterp,pred))+geom_point(aes(y=pred+res,col=doy))+geom_line(col='red',size=1)+
  labs(x='NDVI',y='SAR | LIA,Soilwater',title = 'Partial effect of NDVI')

#Partial effect of soilwater
dat1 %>% mutate_at(vars(sar,lia,ndviInterp,soilwater),scale) %>% 
  mutate(lia=0,ndviInterp=0) %>% 
  mutate(pred=predict(sarMod1,newdata=.,re.form=~0),res=resid(sarMod1)) %>% 
  ggplot(aes(soilwater,pred))+geom_point(aes(y=pred+res,col=doy))+geom_line(col='red',size=1)+
  labs(x='Soilwater',y='SAR | LIA,NDVI',title = 'Partial effect of Soilwater')

dat2 <- dat1 %>% 


