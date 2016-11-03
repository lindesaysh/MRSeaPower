devtools::load_all(pkg='C://MarineScotlandPower/MRSeaPower')
devtools::load_all(pkg='C://MarineScotlandPower/MRSea/MRSea')


setwd("C:/MarineScotlandPower/MRSeaPower")


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
# ~~~ OC Change ~~~~
# ~~~~~~~~~~~~~~~~~~~

nsim=500
impdata.fow<-genOverallchangeData(log(0.5), bestModel.fow, data = fow.data, panels = 'newbidNum')

newdata.fow.imp<-generateNoise(nsim, impdata.fow$truth, family='poisson', d=summary(bestModel.fow)$dispersion)
corrmat.fow.imp<-rbind(corrmat.fow, corrmat.fow)
newdata.fow.cor<-generateIC(data = impdata.fow, corrs = corrmat.fow.imp, panels = 'panels', newdata=newdata.fow.imp, nsim=nsim)


bestModel.fow$splineParams[[1]]$dist<-rbind(bestModel.fow$splineParams[[1]]$dist, bestModel.fow$splineParams[[1]]$dist)
fowsim_glm<-update(bestModel.fow, newdata.fow.cor[,1]~. + eventphase, data=impdata.fow)
fowsim_glm$panels<-impdata.fow$panels


# make sure that the independent data is used to get the null distribution
empdistpower<-getRunsCritVals(n.sim = nsim, simData=newdata.fow.imp,
                              model = fowsim_glm, data = impdata.fow, plot=FALSE,
                              returnDist = TRUE, dots=FALSE)

data("fowshco.grid")
predictdata<-rbind(data.frame(fowshco.grid, TideState=1, WindStrength=0, SeaState=1, SimpPrecipitation='NONE', CloudCover=7, MonthInt=6, areatime=fowshco.grid$Area, eventphase=0), data.frame(fowshco.grid, TideState=1, WindStrength=0, SeaState=1, SimpPrecipitation='NONE', CloudCover=7, MonthInt=6,areatime=fowshco.grid$Area, eventphase=1))

g2k<-makeDists(cbind(predictdata$x.pos, predictdata$y.pos), knotcoords = na.omit(fowsim_glm$splineParams[[1]]$knotgrid), knotmat = FALSE)$dataDist

# POWER ANALYSIS
nsim=500
system.time(
  powerout.fow.oc<-powerSimPll(newdata.fow.cor, fowsim_glm, empdistpower, nsim=nsim, powercoefid=length(coef(fowsim_glm)), predictionGrid=predictdata, g2k=g2k, splineParams=fowsim_glm$splineParams, n.boot=500, nCores=8)

)
save(powerout.fow.oc, file='testing/FinalReportCode/powerout.fow.oc.RData', compress = 'bzip2')


null.output.fow.oc<-pval.coverage.null(newdat.ind = newdata.fow.imp, newdat.corr = newdata.fow.cor, model = fowsim_glm, nsim = nsim, powercoefid = length(coef(fowsim_glm)))

save(null.output.fow.oc, file='testing/FinalReportCode/null.output.fow.oc.RData', compress = 'bzip2')

summary(powerout.fow.oc, null.output.fow.oc, truebeta=log(0.5))

powerPlot(powerout.fow.oc)

plotdata<-plot.sigdiff(powerout.fow.oc, predictdata[predictdata$eventphase==0,c('x.pos', 'y.pos')], tailed='two', error.rate = 0.05)
plotdata

plotdata<-plot.sigdiff(powerout.fow.oc, predictdata[predictdata$eventphase==0,c('x.pos', 'y.pos')], tailed='two', error.rate = 0.05, adjustment = 'bonferroni')
plotdata

plotdata<-plot.sigdiff(powerout.fow.oc, predictdata[predictdata$eventphase==0,c('x.pos', 'y.pos')], tailed='two', error.rate = 0.05, adjustment = 'sidak')
plotdata


pdata<-plot.preds(powerout.fow.oc, returndata=TRUE)
ggplot(pdata) + geom_tile(aes( x.pos , y.pos , fill = mean, width=500, height=500 ) ) + scale_fill_gradientn(colours=mypalette, values = breaks, space = "Lab", na.value = "grey50", guide = "colourbar", name='Mean Prediction') + theme_bw() + facet_wrap(~type+eventphase, nrow=3) + coord_equal()


ddata<-plot.diffs(powerout.fow.oc, returndata=TRUE)

ggplot(ddata) + geom_tile(aes( x.pos , y.pos , fill = mean, width=500, height=500 ) ) + scale_fill_gradientn(colours=mypalette, values = breaks, space = "Lab", na.value = "grey50", guide = "colourbar", name='Mean Prediction') + theme_bw() + facet_wrap(~type, nrow=1) + coord_equal()



