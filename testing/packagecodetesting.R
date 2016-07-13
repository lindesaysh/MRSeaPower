

nsim=1000

dat<-makeToyData(200, length.panels=5, changecoef.link = log(0.2))
#newdatH0<-generateNoise(nsim, dat$mu, family='poisson', d=10)
newdat<-generateNoise(nsim, dat$mu, family='poisson', d=10)

init_glm<-glm(newdat[,1] ~ x + as.factor(evph), data=dat, family='quasipoisson')

empdistribution<-getRunsCritVals(n.sim = nsim, simData=newdat, model = init_glm, data = dat, plot=TRUE, returnDist = TRUE)

ps<-matrix(NA, nrow=nsim, ncol=2)
for(i in 1:nsim){
sim_glm<-glm(newdat[,i] ~ x + as.factor(evph), data=dat, family='quasipoisson')

ps[i,1]<-runs.test(residuals(sim_glm, type='pearson'), critvals = empdistribution)$p.value
ps[i,2]<-runs.test(residuals(sim_glm, type='pearson'))$p.value

}

length(which(ps[,1]<0.05))/nsim
length(which(ps[,2]<0.05))/nsim


rho=0.3

newdatcorr<-generateIC(dat, c(1, rho, rho^2, rho^3, rho^4), 'panelid', newdat, ncol(newdat))

init_glmc<-glm(newdatcorr[,10] ~ x + as.factor(evph), data=dat, family='quasipoisson')
acf(residuals(init_glmc, type='pearson'))
runs.test(residuals(init_glmc, type='pearson'), critvals = empdistribution)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd('C:\\Users\\lindesay\\Dropbox\\Research\\MarSco_Power_Analysis\\Data\\FoW')
dat<-read.csv('modellingoutputs/FittedData_SHCO.csv')
dat<- dat[dat$YearInt==2009,]
head(dat)
dat$newbidNum<- as.numeric(dat$newbid)

init_glm<-glm(response ~  TideState + WindStrength + SeaState + SimpPrecipitation + CloudCover, data=dat, family=quasipoisson)

# ~~~~~~~~~~~~~~~~~~
## overall change
dat2<-genOverallchangeData(log(0.5), init_glm, data = dat, panels = 'newbidNum')

require(dplyr)
plotdata<-group_by(dat2, x.pos, y.pos, eventphase)%>%
  summarise(fitsredist=mean(truth), areatime=first(areatime))
p <- ggplot(plotdata)
p + geom_tile(aes(x=x.pos/1000, y=y.pos/1000, fill=fitsredist, height=0.5*sqrt(areatime),  width=0.5*sqrt(areatime))) + facet_wrap(~eventphase) + coord_equal() + theme_bw()

group_by(dat2, eventphase)%>%
  summarise(sum(truth))

# ~~~~~~~~~~~~~~~~~~
## redist

imppoly<-list(x.pos=c(510210.0, 510944.6, 511353.8, 511166.5, 510244.5), y.pos=c(6556624, 6556997, 6555979, 6554563, 6556177))

dat3<-genRedistData(init_glm, data=dat, changecoef.link = log(0.5), panels='newbidNum', imppoly =imppoly)

plotdata<-group_by(dat3, x.pos, y.pos, eventphase)%>%
  summarise(fitsredist=mean(truth), areatime=first(areatime))
p <- ggplot(plotdata)
p + geom_tile(aes(x=x.pos/1000, y=y.pos/1000, fill=fitsredist, height=0.5*sqrt(areatime),  width=0.5*sqrt(areatime))) + facet_wrap(~eventphase) + coord_equal() + theme_bw()

group_by(dat3, eventphase)%>%
  summarise(sum(truth), mean(truth))

group_by(dat3, eventphase, impcells)%>%
  summarise(sum(truth), mean(truth))

# ~~~~~~~~~~~~~~~~~~
## redist + overall change
dat4<-genRedistData( init_glm, data=dat, c(log(0.5), log(0.75)), panels='newbidNum', imppoly=imppoly)

plotdata<-group_by(dat4, x.pos, y.pos, eventphase)%>%
  summarise(fitsredist=mean(truth), areatime=first(areatime))
p <- ggplot(plotdata)
p + geom_tile(aes(x=x.pos/1000, y=y.pos/1000, fill=fitsredist, height=0.5*sqrt(areatime),  width=0.5*sqrt(areatime))) + facet_wrap(~eventphase) + coord_equal() + theme_bw()

group_by(dat4, eventphase)%>%
  summarise(sum(truth), mean(truth))

group_by(dat4, eventphase, impcells)%>%
  summarise(sum(truth), mean(truth))

