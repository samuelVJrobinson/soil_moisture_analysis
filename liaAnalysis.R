#FALL 2019, DATA ANALYSIS WITH LAN NGUYEN EXAMINING SAR DATA FROM VERMILLION AREA
#GOAL: TRY TO PULL OUT TEMPORAL SIGNAL AND LIA SIGNAL FROM SAR BACKSCATTER DATA

#Packages
rm(list=ls()) #Clear workspace

library(tidyverse)
library(ggplot2)
library(mgcv)
library(parallel)
library(beepr)
theme_set(theme_classic())


# Dataset 1 ----------------------------------------------------------------
#Dataset1: small dataset from random locations around Vermillion
rm(list=ls())
setwd("~/Projects/LIA_analysis/Dataset1")
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
setwd("~/Projects/UofC postdoc/LIA_analysis/Dataset2")
load('allFieldData.Rdata')

shortDat <- dat[[1]] #Field 1


#Step 1: create NDVI model to fill in 


#Plot of SAR as a function of DOY at 3 cells
shortDat %>% filter(cell=='0_0'|cell=='50_50'|cell=='100_100') %>% 
  ggplot(aes(doy,sar))+geom_point()+
  geom_smooth(method='gam',formula=y~s(x,k=50),se=F)+
  facet_wrap(~cell,ncol=1)

#Raster plot of SAR
shortDat %>% 
  filter(doy==9|doy==165|doy==213|doy==261) %>% 
  ggplot(aes(colnum,rownum,fill=sar))+
  facet_wrap(~doy,ncol=2) +
  geom_raster()+scale_y_reverse()

#Dataset too large. Aggregating to lower resolution 
shorterDat <- shortDat %>% 
  # mutate(newCol=colnum%/%5,newRow=rownum%/%5) %>% #Aggregates pixels to 5x5 (50x50m)
  # mutate(newCol=colnum%/%4,newRow=rownum%/%4) %>% #Aggregates pixels to 4x4 (40x40m)
  mutate(colnum=colnum%/%3,rownum=rownum%/%3) %>% #Aggregates pixels to 3x3 (30x30m)
  group_by(colnum,rownum,doy) %>% summarize_at(vars(sar,lia),mean) %>% ungroup() %>% 
  mutate(cell=factor(paste(colnum,rownum,sep='_'))) %>% 
  mutate(sar=scale(sar),lia=(lia/100)-40)

shorterDat %>% 
  filter(cell=='3_7'|cell=='0_30'|cell=='30_30') %>% 
  ggplot(aes(doy,sar))+geom_point()+
  geom_smooth(method='gam',formula=y~s(x,k=50),se=F)+
  facet_wrap(~cell,ncol=1)

shorterDat %>% 
  filter(doy==9|doy==165|doy==213|doy==261) %>% 
  ggplot(aes(colnum,rownum,fill=sar))+
  facet_wrap(~doy,ncol=2) +
  geom_raster()+scale_y_reverse()

detectCores()
cl <- makeCluster(8) #6 CPUs

# mod1 <- gamm(sar~1+lia+s(doy,bs='cr'),random=list(cell=~1+lia),
#              data=shorterDat,control=list(msVerbose=TRUE))
mod2 <- bam(sar~1+lia+s(doy,bs='cr')+s(cell,bs='re')+s(cell,lia,bs='re'),
            data=shorterDat,cluster=cl); beep(2)
mod3 <- bam(sar~1+lia+s(doy,bs='cr')+s(cell,bs='re'),data=shorterDat,cluster=cl); beep(2)
mod4 <- bam(sar~1+lia+s(doy,bs='cr',k=30),data=shorterDat,cluster=cl); beep(2)

# summary(mod1$gam)
# summary(mod1$lme)
# plot(mod1$gam,select=1)
# plot(mod1$lme)
summary(mod2) #Between-cell variance is very low, and estimates appear iffy. Diagnosis: random effects should be removed for now.
gam.vcomp(mod2)
plot(mod2)

summary(mod3) 
gam.vcomp(mod3)
summary(mod4) 
gam.vcomp(mod4)
plot(mod4)

data.frame(doy=121:243,pred=predict(mod4,newdata=data.frame(lia=0,doy=121:243))) %>% 
  ggplot(aes(doy,pred))


# #Compare random slopes and intercepts from gamm and bam models
# data.frame(gammInt=ranef(mod1$lme)[[2]][,1],gammSlope=ranef(mod1$lme)[[2]][,2],
#            bamInt=coef(mod2)[grepl('s\\(cell\\)',names(coef(mod2)))],
#            bamSlope=coef(mod2)[grepl('s\\(cell,lia\\)',names(coef(mod2)))]) %>% 
#   pairs(.,upper.panel=function(x,y) {points(x,y); abline(0,1,col='red')},
#         lower.panel=function(x,y) text(mean(x),mean(y),round(cor(x,y),3),cex=0.1+cor(x,y)*2),
#         main='Intercepts')

# #Plot of random effects from gamm model
# ranef(mod1$lme)[[2]] %>%
#   rename('Int'='(Intercept)','liaSlope'='lia') %>% 
  # rownames_to_column() %>% mutate(rowname=gsub('1/','',rowname)) %>% 
  # separate(rowname,c('col','row'),sep='_',convert=T) %>% 
  # mutate_at(vars(Int,liaSlope),scale) %>%
  # pivot_longer(Int:liaSlope,names_to = 'type',values_to = 'value') %>% 
  # # filter(type=='liaSlope') %>% 
  # ggplot(aes(x=col,y=row,fill=value))+geom_raster()+
  # facet_wrap(~type)+
  # labs(fill='Scaled\nvalue')+
  # scale_y_reverse()


ranefMod2 <- data.frame(int=coef(mod2)[grepl('s\\(cell\\)',names(coef(mod2)))],
           liaSlope=coef(mod2)[grepl('s\\(cell,lia\\)',names(coef(mod2)))],
           cell=levels(shorterDat$cell)) %>% 
  separate(cell,c('col','row'),sep='_',convert=T)

ggplot(ranefMod2,aes(x=col,y=row,fill=int))+geom_raster()+
  scale_y_reverse()
ggplot(ranefMod2,aes(x=col,y=row,fill=liaSlope))+geom_raster()+
  scale_y_reverse()



           
#Appears to be spatial AC, both in intercepts and slopes, but need to verify. Check at 53.50933 -110.1945 (this field)

stopCluster(cl)


# Calculate NDVI --------------------------------------------------------------------
setwd("~/Projects/LIA_analysis/NDVI1")
ndvi <- read.table('NDVI.txt',stringsAsFactors=F)

ndvi <- strsplit(ndvi[c(1:nrow(ndvi)),],',')

ndvi <- lapply(1:(length(ndvi)/2),function(x){
  data.frame(site=x,doy=ndvi[[(x*2)-1]],NDVI=ndvi[[(x*2)]])
})

ndvi <- do.call('rbind',ndvi) %>% mutate(site=factor(site)) %>% 
  mutate_at(vars(doy:NDVI),function(x) as.numeric(as.character(x))) %>% 
  mutate(logNDVI=log(NDVI))

ggplot(ndvi,aes(doy,logNDVI))+geom_point()+
  facet_wrap(~site)+
  geom_smooth(method='gam',formula=y~s(x,bs='cr'))+
  labs(y='log(ndvi) value')

mod5 <- gam(logNDVI~site+s(doy,by=site,bs='cs'),data=ndvi)
summary(mod5)
gam.check(mod5)
plot(mod5,pages=1,residuals=F,se=F,rug=F,by.resids=F)

#Predictions for range of date values
ndviMod <- with(ndvi,expand.grid(site=unique(site),doy=min(doy):max(doy))) %>% 
  mutate(predNDVI=exp(predict(mod5,newdata=data.frame(site,doy)))) #Back-transform

#Check predictions - looks OK
ndviMod %>% 
  ggplot(aes(doy,predNDVI))+geom_line()+
  geom_point(data=ndvi,aes(x=doy,y=NDVI))+
  facet_wrap(~site)+labs(y='NDVI')


# Use NDVI with tensor product ----------------------------------------------------------

detectCores()
cl <- makeCluster(12)

#Range of dates to restrict analysis to
dateRange <- filter(ndviMod,site==1) %>% summarize(lwr=min(doy),upr=max(doy))

shortDat <- dat[[1]] %>%  #Field 1
  filter(doy>=dateRange$lwr,doy<=dateRange$upr) %>% #Restrict to doy range of NDVI data
  left_join(filter(ndviMod,site==1),by='doy') %>% 
  mutate(NDVIScal=scale(predNDVI)[,1]) %>%  #Scale NDVI
  mutate(sarScal=sarScal[,1],liaScal=liaScal[,1])

#More complex version with tensor product smooth  
mod6 <- bam(sarScal ~ liaScal + s(NDVIScal) + te(colnum,rownum,doy),data=shortDat,cluster=cl); beep(1)
summary(mod6)
gam.check(mod6)
plot(mod6,se=T,all.terms=T,residuals=T,pages=1)

shortDat %>% mutate(resid=resid(mod6)) %>% 
  ggplot(aes(NDVIScal,resid))+geom_point(alpha=0.5)+
  geom_smooth()

shortDat %>% mutate(pred=predict(mod6,newdata=mutate_at(shortDat,vars(liaScal,NDVIScal),mean)),resid=resid(mod6)) %>% 
  filter(cell=='50_50'|cell=='0_0'|cell=='100_100') %>% 
  ggplot(aes(x=doy))+geom_line(aes(y=pred))+geom_point(aes(y=pred+resid))+
  facet_wrap(~cell,ncol=1)
  
  

#Simpler version with single smoother
mod6a <- bam(sarScal ~ liaScal + s(NDVIScal) + s(doy),data=shortDat,cluster=cl); beep(1)
summary(mod6a)
gam.check(mod6a)
plot(mod6a,se=T,all.terms=T,residuals=T,pages=1)




# stopCluster(cl)




shortDat %>% mutate(resid=resid(mod6)) %>% 
  filter(doy==165) %>% 
  ggplot(aes(colnum,rownum,fill=resid))+geom_raster()+
  scale_y_reverse()


  
stopCluster(cl)  







