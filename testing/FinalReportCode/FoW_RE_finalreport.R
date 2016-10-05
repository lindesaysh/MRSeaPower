devtools::load_all(pkg='C://MarineScotlandPower/MRSeaPower')
devtools::load_all(pkg='C://MarineScotlandPower/MRSea/MRSea')

require(fields)
require(splines)
require(mgcv)
require(MRSea)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
breaks<-ticks<-c(0, 0.2, 0.4, 0.6, 0.8, 1, 10, 25)
mypalette<-rev(brewer.pal(length (ticks)-1,"Spectral"))
require(sp)
require(raster)

data(fowshco)
load('testing/FoW/FoW_model.RData')
bestModel.fow<-bestModel
fow.data<-bestModel.fow$data


newdata.fow<-generateNoise(n=500, response=fitted(bestModel.fow), family='poisson', d=summary(bestModel.fow)$dispersion)

corrmat.fow<-getCorrelationMat(panel=fow.data$newbidNum, data = fow.data$response)
newdata.fow.cor<-generateIC(data = fow.data, corrs = corrmat.fow, panels = 'newbidNum', newdata=newdata.fow, nsim=500,dots = FALSE)

# ~~~~~~~~~~~~~~~~~~~
# ~~~ RE Change ~~~~
# ~~~~~~~~~~~~~~~~~~~

fow.data$index<-1:nrow(fow.data)
impcells<- data.frame(Zone=c( 'C1' , 'D1' , 'C2' , 'D2' , 'C3' ),impactcells=rep(0, 5))
mdat<-merge(fow.data, impcells, all.x=TRUE)
impactcells<-arrange(mdat, index)$impactcells
impactcells[which(is.na(impactcells))]<-1

nsim=500
truebeta<-log(0.5)
impdata.fow.re<-genRedistData(bestModel.fow, data=fow.data, changecoef.link=truebeta, panels='newbidNum', imppoly=NULL, impactcells = impactcells)

##
newdata.fow.imp.re<-generateNoise(nsim, impdata.fow.re$truth, family='poisson', d=summary(bestModel.fow)$dispersion)
corrmat.fow.imp<-rbind(corrmat.fow, corrmat.fow)
newdata.fow.cor.re<-generateIC(data = impdata.fow.re, corrs = corrmat.fow.imp, panels = 'panels', newdata=newdata.fow.imp.re, nsim=nsim)

##
fowsim_glm<-update(bestModel.fow, newdata.fow.cor.re[,1] ~. + as.factor(eventphase) -
                     LocalRadialFunction(radiusIndices, dists, radii, aR), data=impdata.fow.re)
fowsim_glm<-make.gamMRSea(fowsim_glm, panelid=1:nrow(impdata.fow.re),
                          splineParams=fowsim_glm$splineParams,
                          varshortnames=fowsim_glm$varshortnames,gamMRSea=TRUE)
require(MuMIn)
salsa2dlist<-list(fitnessMeasure = 'QAIC', knotgrid = fowsim_glm$splineParams[[1]]$knotgrid,
                  knotdim = c(100, 100),  startKnots=10, minKnots=2,
                  maxKnots=20, r_seq=fowsim_glm$splineParams[[1]]$radii, gap=0, interactionTerm='as.factor(eventphase)')

salsa2dOutputint<-runSALSA2D(fowsim_glm, salsa2dlist,
                             d2k=rbind(fowsim_glm$splineParams[[1]]$dist,
                                       fowsim_glm$splineParams[[1]]$dist),
                             k2k=fowsim_glm$splineParams[[1]]$knotDist,
                             splineParams=NULL, tol=0, chooserad=F, panels=NULL,
                             suppress.printout=TRUE)


bestModel_int<-make.gamMRSea(salsa2dOutputint$bestModel, panelid = 1:nrow(impdata.fow.re),
                             splineParams = salsa2dOutputint$splineParams,
                             varshortnames = NULL, gamMRSea=TRUE)

summary(bestModel_int)
anova(bestModel_int)


# make sure that the independent data is used to get the null distribution
empdistpower<-getRunsCritVals(n.sim = nsim, simData=newdata.fow.imp.re,
                              model = bestModel_int, data = impdata.fow.re, plot=FALSE,
                              returnDist = TRUE, dots=FALSE)

###
data("fowshco.grid")
predictdata<-rbind(data.frame(fowshco.grid, TideState=1, WindStrength=0, SeaState=1, SimpPrecipitation='NONE', CloudCover=7, MonthInt=6, areatime=fowshco.grid$Area, eventphase=0), data.frame(fowshco.grid, TideState=1, WindStrength=0, SeaState=1, SimpPrecipitation='NONE', CloudCover=7, MonthInt=6,areatime=fowshco.grid$Area, eventphase=1))

g2k<-makeDists(cbind(predictdata$x.pos, predictdata$y.pos), knotcoords = na.omit(bestModel_int$splineParams[[1]]$knotgrid), knotmat = FALSE)$dataDist



# POWER ANALYSIS
nsim=100
system.time(
  powerout.fow.re<-powerSimPll(newdata.fow.cor.re, bestModel_int, empdistpower, nsim=nsim, powercoefid=length(coef(bestModel_int)), predictionGrid=predictdata, g2k=g2k, splineParams=bestModel_int$splineParams, n.boot=500, impact.loc=c(510750, 6555700), nCores=8)

)
save(powerout.fow.re, file='testing/FinalReportCode/powerout.fow.re.RData', compress = 'bzip2')


null.output.fow.re<-pval.coverage.null(newdat.ind = newdata.fow.imp.re, newdat.corr = newdata.fow.cor.re, model = bestModel_int, nsim = nsim, powercoefid = length(coef(bestModel_int)))

save(null.output.fow.re, file='testing/FinalReportCode/null.output.fow.re.RData', compress = 'bzip2')

summary(powerout.fow.oc, null.output, truebeta=log(0.5))

powerPlot(powerout.fow.oc)

plotdata<-plot.sigdiff(powerout.fow.oc, predictdata[predictdata$eventphase==0,c('x.pos', 'y.pos')], tailed='two', error.rate = 0.05, family=FALSE)
plotdata

plotdata<-plot.sigdiff(powerout.fow.oc, predictdata[predictdata$eventphase==0,c('x.pos', 'y.pos')], tailed='two', error.rate = 0.05, family=TRUE)
plotdata

### Distance to event site plot

plot.d2imp(powerout.fow.oc)

plot.d2imp(powerout.fow.oc, pct.diff = FALSE)
