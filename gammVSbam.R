#TEST OF BAM VS GAMM PARAMETER ESTIMATION:

#Idea: y = [b0 + boi] + x[b1 + b1i] + s(z) + error, 
#   where b0i and 1i are random intercepts and slopes, and s(z) is some curvy spline

library(mgcv)

#Generate some fake data
set.seed(123)
ndat <- 900 #900 data points
ngroup <- 30 #30 groups
gInd <- sample(rep(1:30,length.out=900))#Group index
xVal <- runif(ndat,-1,1) #Range of x values
zVal <- runif(ndat,-pi,pi) #Range of z values

bo <- -1 #Intercept
b1 <- 2 #Slope of x
b0i <- rnorm(ndat/ngroup,0,0.5) #Random intercepts
b1i <- rnorm(ndat/ngroup,0,0.5) #Random slopes
s_z <- sin(zVal) #Sine curve to estimate by spline
err <- rnorm(ndat,0,0.2) #Noise
yVal <- (bo + b0i[gInd]) + (b1 + b1i[gInd])*xVal + s_z + err #Generate values
plot(xVal,yVal) #Looks OK
plot(zVal,yVal)
dat <- data.frame(y=yVal,x=xVal,z=zVal,group=factor(gInd)) #Assemble into dataframe
str(dat)

#Help file for s(x,bs='re')
# ?smooth.construct.re.smooth.spec

#GAMM MODEL
mod1 <- gamm(y~1+x+s(z,bs='cr'),random=list(group=~1+x),data=dat,
             control=list(msVerbose=TRUE))
summary(mod1$lme) #Estimates variance terms well, slope terms well
summary(mod1$gam)
plot(mod1$gam); curve(sin(x),-pi,pi,col='red',add=T) #Pretty good estimation of sine curve

#BAM MODEL

#Using "by" notation within s()
mod2a <- bam(y~1+x+s(z,bs='cr')+s(group,bs='re')+s(x,by=group,bs='re'),data=dat) 
gam.vcomp(mod2a) #Note the number of sd terms estimated -- 1 for each group
summary(mod2a)

#Using paired notation
mod2b <- bam(y~1+x+s(z,bs='cr')+s(group,bs='re')+s(x,group,bs='re'),data=dat)
gam.vcomp(mod2b) #Fewer sd terms estimated
summary(mod2b)

#Compare estimates of intercepts and slopes from :
ranInt <- data.frame(actual=b0i,gamm=ranef(mod1$lme)[[2]][,1],
                     bamA=coef(mod2a)[grepl('s\\(group\\).',names(coef(mod2a)))],
                     bamB=coef(mod2b)[grepl('s\\(group\\).',names(coef(mod2b)))])
ranSlope <- data.frame(actual=b1i,gamm=ranef(mod1$lme)[[2]][,2],
                       bamA=coef(mod2a)[grepl('s\\(x\\).',names(coef(mod2a)))],
                       bamB=coef(mod2b)[grepl('s\\(x,group\\).',names(coef(mod2b)))])

#Good estimates of random intercept terms
pairs(ranInt,upper.panel=function(x,y) {points(x,y); abline(0,1,col='red')},
      lower.panel=function(x,y) text(mean(x),mean(y),round(cor(x,y),3),cex=0.1+cor(x,y)*2),
      main='Intercepts')

#bamA slope estimates are offset by a constant value when using "by" 
pairs(ranSlope,upper.panel=function(x,y) {points(x,y); abline(0,1,col='red')},
      lower.panel=function(x,y) text(mean(x),mean(y),round(cor(x,y),3),cex=0.1+cor(x,y)*2),main='Slopes')

