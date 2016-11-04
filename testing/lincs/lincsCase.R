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
require(splancs)


lincs<-read.csv('D:\\lindesay\\MarScoPower\\finalreportdata\\phase1Data.csv')

lincs$response<-lincs$AUK
lincs$foldid<-getCVids(lincs, 10, 'line_id')
  
initialModel<-glm(response ~ 1 + offset(log(area)), data=lincs, family=quasipoisson)

factorlist<-NULL
varlist<-c('depth', 'month')

salsa1dlist<-list(fitnessMeasure='QAIC', minKnots_1d = c(1,1), maxKnots_1d=c(5,5), startKnots_1d = c(1,1), degree=c(2,2), maxIterations=100, gaps=c(0,0))

salsa1dout<-runSALSA1D_withremoval(initialModel, salsa1dlist, varlist, factorlist, varlist_cyclicSplines = c('month'), splineParams = NULL, datain=lincs, suppress.printout = TRUE)

salsa1dout$bestModel$splineParams<-salsa1dout$splineParams

bestModel1D<-make.gamMRSea(salsa1dout$bestModel, panelid = 1:nrow(lincs), splineParams = salsa1dout$splineParams)

summary(salsa1dout$bestModel)

anova(salsa1dout$bestModel)


runPartialPlots(bestModel1D, data = lincs, factorlist.in = factorlist, varlist.in = varlist, showKnots = T)

knotgrid<-getKnotgrid(cbind(lincs$x.pos, lincs$y.pos))
distMats<-makeDists(cbind(lincs$x.pos, lincs$y.pos), na.omit(knotgrid))
# choose sequence of radii
r_seq<-getRadiiChoices(8, distMats$dataDist)

salsa2dlist<-list(fitnessMeasure = 'QAIC', knotgrid = knotgrid, knotdim = c(100, 100),  startKnots=6, minKnots=2, maxKnots=20, r_seq=r_seq, gap=0)


salsa2dOutput<-runSALSA2D(bestModel1D, salsa2dlist, d2k=distMats$dataDist, k2k=distMats$knotDist, splineParams=NULL, tol=0, chooserad=F, panels=NULL, suppress.printout=FALSE)

bestModel<-make.gamMRSea(salsa2dOutput$bestModel, panelid = 1:nrow(lincs), splineParams = salsa2dOutput$splineParams, varshortnames = varlist, gamMRSea=TRUE)

par(mfrow=c(1,2))
quilt.plot(lincs$x.pos, lincs$y.pos, lincs$response, asp=1, nrow=80, ncol=80)
quilt.plot(lincs$x.pos, lincs$y.pos, fitted(bestModel), asp=1)

save(bestModel, file='D:/lindesay/MarScoPower/finalreportdata/lincsAUKmodel.RData', compress='bzip2')

### best model - power analysis

nsim=500
impdata.auk<-genOverallchangeData(log(0.5), bestModel, data = lincs)

t<-group_by(impdata.auk, eventphase)%>%
  summarise(sum=sum(truth), mean=mean(truth), n=n())
t

bef<-impdata.auk[impdata.auk$eventphase==0,]
aft<-impdata.auk[impdata.auk$eventphase==1,]

spdf1<- SpatialPointsDataFrame( data.frame( x = bef$x.pos , y = bef$y.pos ) , data = data.frame( z = bef$truth) )
e <- extent( spdf1 )
ratio <- ( e@xmax - e@xmin ) / ( e@ymax - e@ymin )
r <- raster( nrows = 60 , ncols = floor( 60 * ratio ) , ext = extent(spdf1) )
rf <- rasterize( spdf1 , r , field = "z" , fun= mean)
rdf1 <- data.frame( rasterToPoints( rf ), evph=0 )   

spdf2<- SpatialPointsDataFrame( data.frame( x = aft$x.pos , y = aft$y.pos ) , data = data.frame( z = aft$truth) )
e <- extent( spdf2 )
ratio <- ( e@xmax - e@xmin ) / ( e@ymax - e@ymin )
r <- raster( nrows = 60 , ncols = floor( 60 * ratio ) , ext = extent(spdf2) )
rf <- rasterize( spdf2 , r , field = "z" , fun= mean)
rdf2 <- data.frame( rasterToPoints( rf ), evph=1 ) 

rdf<-rbind(rdf1, rdf2)
names(rdf)[3]<-'Mean.count'

pct.change<-((rdf$Mean.count[rdf$evph==1] - rdf$Mean.count[rdf$evph==0])/rdf$Mean.count[rdf$evph==0]
)*100

ggplot( NULL ) + geom_raster( data = rdf1 , aes( x , y , fill = pct.change ) ) + scale_fill_gradientn(colours=mypalette, values = breaks,space = "Lab", na.value = "grey50", guide = "colourbar") + theme_bw()


nsim=500
truebeta<-log(0.5)
impdata.fat<-genOverallchangeData(bestModel, data=lincs, changecoef.link=truebeta)
newdata.fat.imp<-generateNoise(nsim, impdata.fat$truth, family='poisson', d=summary(bestModel)$dispersion)

bestModel$splineParams[[1]]$dist<-rbind(bestModel$splineParams[[1]]$dist, bestModel$splineParams[[1]]$dist)
fatsim_glm<-update(bestModel, newdata.fat.imp[,1]~. + eventphase, data=impdata.fat)
fatsim_glm$panels<-impdata.fat$panels


# make sure that the independent data is used to get the null distribution
empdistpower.fat<-getRunsCritVals(n.sim = nsim, simData=newdata.fat.imp, 
                                  model = fatsim_glm, data = impdata.fat, plot=FALSE, 
                                  returnDist = TRUE, dots=FALSE)


load('D:\\lindesay\\MarScoPower\\finalreportdata\\HiDef2016_pred_grid.rdata')

