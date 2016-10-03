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
# ~~~ OC Change ~~~~
# ~~~~~~~~~~~~~~~~~~~
nsim=500
bnd.wf<-data.frame(x.pos=c(543289.6, 574317.5, 574892.1, 545588.0),y.pos=c(6253767, 6253097, 6237016, 6236681))
truebeta<-log(0.5)
impdata.fat<-genRedistData(bestModel.fat, data=fat.data, changecoef.link=truebeta, panels='panel', imppoly=bnd.wf)
newdata.fat.imp<-generateNoise(nsim, impdata.fat$truth, family='binomial', size=1)

bestModel.fat$splineParams[[1]]$dist<-rbind(bestModel.fat$splineParams[[1]]$dist, bestModel.fat$splineParams[[1]]$dist)
fatsim_glm<-update(bestModel.fat, newdata.fat.imp[,1]~. + eventphase, data=impdata.fat)
fatsim_glm$panels<-impdata.fat$panels


# make sure that the independent data is used to get the null distribution
empdistpower.fat<-getRunsCritVals(n.sim = nsim, simData=newdata.fat.imp, 
                                  model = fatsim_glm, data = impdata.fat, plot=FALSE, 
                                  returnDist = TRUE, dots=FALSE)

grid<-expand.grid(seq(min(fat.data$x.pos), max(fat.data$x.pos), by=1000), seq(min(fat.data$y.pos), max(fat.data$y.pos), by=1000))
names(grid)<-c('x.pos', 'y.pos')

bnd<-cbind(fat.data$x.pos, fat.data$y.pos)[chull(x = fat.data$x.pos, fat.data$y.pos),]
gridsub<-grid[inout(grid, bnd),]

predictdata<-rbind(data.frame(gridsub, month=9, year=2011, eventphase=0), data.frame(gridsub, month=9, year=2011, eventphase=1))

g2k<-makeDists(cbind(predictdata$x.pos, predictdata$y.pos), knotcoords = na.omit(bestModel.fat$splineParams[[1]]$knotgrid), knotmat = FALSE)$dataDist



nsim=100
system.time(powerout.fat.oc<-powerSimPll(newdata.fat.imp, fatsim_glm, empdistpower.fat, nsim=nsim, powercoefid=length(coef(fatsim_glm)), predictionGrid=predictdata, g2k=g2k, splineParams=fatsim_glm$splineParams, sigdif=TRUE, n.boot=500, impact.loc=c(581877.4, 6234143), ncores=8))

save(powerout.fat.oc, file='testing/FinalReportCode/powerout.fat.oc.RData', compress = 'bzip2')


null.output.fow.oc<-pval.coverage.null(newdat.ind = newdata.fow.imp, newdat.corr = newdata.fow.cor, model = fowsim_glm, nsim = nsim, powercoefid = length(coef(fowsim_glm)))

save(null.output.fow.oc, file='testing/FinalReportCode/null.output.fow.oc.RData', compress = 'bzip2')
# ~~~~~~~~~~~~~~~~~~~
# ~~~ RE Change ~~~~
# ~~~~~~~~~~~~~~~~~~~
nsim=500
bnd.wf<-data.frame(x.pos=c(543289.6, 574317.5, 574892.1, 545588.0,543289.6),y.pos=c(6253767, 6253097, 6237016, 6236681,6253767))
truebeta<-log(0.5)
impdata.fat.re<-genRedistData(bestModel.fat, data=fat.data, changecoef.link=truebeta, panels='panel', imppoly=bnd.wf, impactcells = NULL)