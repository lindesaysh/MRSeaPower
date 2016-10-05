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


data(fat)
fat$response<-ifelse(fat$response>0, 1, 0)

load('testing/Fat/FaT_model.RData')
bestModel.fat<-bestModel
fat.data<-bestModel.fat$data

newdata.fat<-generateNoise(n=100, response=fitted(bestModel.fat), family='binomial', size=1)



# ~~~~~~~~~~~~~~~~~~~
# ~~~ RE Change ~~~~
# ~~~~~~~~~~~~~~~~~~~
nsim=500
bnd.wf<-data.frame(x.pos=c(543289.6, 574317.5, 574892.1, 545588.0,543289.6),y.pos=c(6253767, 6253097, 6237016, 6236681,6253767))
truebeta<-log(0.5)
impdata.fat.re<-genRedistData(bestModel.fat, data=fat.data, changecoef.link=truebeta, panels='panel', imppoly=bnd.wf, impactcells = NULL)
newdata.fat.imp.re<-generateNoise(nsim, impdata.fat.re$truth, family='binomial', size=1)


fatsim_glm<-update(bestModel.fat, newdata.fat.imp.re[,1]~. + as.factor(eventphase) -
                     LocalRadialFunction(radiusIndices, dists, radii, aR), 
                   data=impdata.fat.re)
fatsim_glm<-make.gamMRSea(fatsim_glm, panelid=1:nrow(impdata.fat.re),
                          splineParams=fatsim_glm$splineParams,
                          varshortnames=fatsim_glm$varshortnames,gamMRSea=TRUE)

salsa2dlist<-list(fitnessMeasure = 'BIC', knotgrid = fatsim_glm$splineParams[[1]]$knotgrid, knotdim = c(100, 100),  startKnots=10, minKnots=2, maxKnots=20, r_seq=fatsim_glm$splineParams[[1]]$radii, gap=0, interactionTerm='as.factor(eventphase)')

salsa2dOutputint<-runSALSA2D(fatsim_glm, salsa2dlist,
                             d2k=rbind(fatsim_glm$splineParams[[1]]$dist,
                                       fatsim_glm$splineParams[[1]]$dist),
                             k2k=fatsim_glm$splineParams[[1]]$knotDist, 
                             splineParams=NULL, tol=0, chooserad=F, panels=NULL,
                             suppress.printout=TRUE)


bestModel_int<-make.gamMRSea(salsa2dOutputint$bestModel, panelid = 1:nrow(impdata.fat.re), splineParams = salsa2dOutputint$splineParams, varshortnames = bestModel.fat$varshortnames, gamMRSea=TRUE)

summary(bestModel_int)
anova(bestModel_int)


# make sure that the independent data is used to get the null distribution
empdistpower.fat<-getRunsCritVals(n.sim = nsim, simData=newdata.fat.imp.re, 
                                  model = bestModel_int, data = impdata.fat.re, plot=FALSE, 
                                  returnDist = TRUE, dots=FALSE)

grid<-expand.grid(seq(min(fat.data$x.pos), max(fat.data$x.pos), by=1000), seq(min(fat.data$y.pos), max(fat.data$y.pos), by=1000))
names(grid)<-c('x.pos', 'y.pos')

bnd<-cbind(fat.data$x.pos, fat.data$y.pos)[chull(x = fat.data$x.pos, fat.data$y.pos),]
gridsub<-grid[inout(grid, bnd),]

predictdata<-rbind(data.frame(gridsub, month=9, year=2011, eventphase=0), data.frame(gridsub, month=9, year=2011, eventphase=1))

g2k<-makeDists(cbind(predictdata$x.pos, predictdata$y.pos), knotcoords = na.omit(bestModel_int$splineParams[[1]]$knotgrid), knotmat = FALSE)$dataDist



nsim=100
system.time(powerout.fat.re<-powerSimPll(newdata.fat.imp.re, bestModel_int, empdistpower.fat, nsim=nsim, powercoefid=length(coef(bestModel_int)), predictionGrid=predictdata, g2k=g2k, splineParams=bestModel_int$splineParams, sigdif=TRUE, n.boot=500, impact.loc=c(559347.4, 6244923), nCores = 8))

save(powerout.fat.re, file='testing/FinalReportCode/powerout.fat.re.RData', compress = 'bzip2')


null.output.fat.re<-pval.coverage.null(newdat.ind = newdata.fat.imp.re, newdat.corr = NULL, model = bestModel_int, nsim = nsim, powercoefid = length(coef(bestModel_int)))

save(null.output.fat.re, file='testing/FinalReportCode/null.output.fat.re.RData', compress = 'bzip2')


summary(powerout.fat.re, null.output.fat.re, truebeta=log(0.5))

powerPlot(powerout.fat.re)

plotdata<-plot.sigdiff(powerout.fat.re, predictdata[predictdata$eventphase==0,c('x.pos', 'y.pos')], tailed='two', error.rate = 0.05, family=FALSE)
plotdata

plotdata<-plot.sigdiff(powerout.fat.re, predictdata[predictdata$eventphase==0,c('x.pos', 'y.pos')], tailed='two', error.rate = 0.05, family=TRUE)
plotdata

### Distance to event site plot

plot.d2imp(powerout.fat.re)

plot.d2imp(powerout.fat.re, pct.diff = FALSE)
