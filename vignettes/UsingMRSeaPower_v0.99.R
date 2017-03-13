## ----echo=FALSE, message=FALSE, warning=FALSE----------------------------
require(knitcitations)
cleanbib()
#biblio <- read.bibtex("newref.bib")
cite_options(citation_format = 'pandoc', cite.style = 'authoryear', max.names = 1, longnamesfirst=FALSE)

## ----setup, echo=FALSE, warning=FALSE, message=FALSE---------------------
require(MRSea)
devtools::load_all(pkg='C://MarineScotlandPower/MRSeaPower')

## ------------------------------------------------------------------------
data("nystedA_slim")

## ------------------------------------------------------------------------
data("nysted.studybnd")
data("nysted.bndwf")
data("nysted.predgrid")

## ----arrangedat, echo=TRUE-----------------------------------------------
nysted$panelid<-as.numeric(nysted$unique.transect.label)
nysted$foldid<-getCVids(nysted, 10, 'panelid')

## ----nysteddatplot, fig.cap='Figure showing the survey effort and proposed windfarm site for the Greater Wash.', comment=FALSE, message=FALSE, fig.height=6, fig.width=8----
require(splancs)
plot(nysted$x.pos, nysted$y.pos, pch=20, cex=0.8, col='grey', xlab='x-coordinate (m)',
     ylab='y-coordinate (m)', asp=1)
polymap(nysted.bndwf, add=T)
polymap(nysted.studybnd, add=T)

## ----initsetup, warning=FALSE, message=FALSE, echo=TRUE------------------
initialModel<-gamMRSea(response ~ 1 + as.factor(yearmonth)+depth + 
                    x.pos + y.pos + offset(log(area)),  data=nysted, family=quasipoisson)

## ---- warning=FALSE, message=FALSE---------------------------------------
anova(initialModel, test='F')

## ------------------------------------------------------------------------
bestModel<-initialModel

## ----echo=TRUE, comment=""-----------------------------------------------
runsTest(residuals(bestModel, type='pearson'), emp.distribution = empdistribution)

## ---- fig.height=6, fig.width=8------------------------------------------
acf(residuals(bestModel, type='pearson'))

## ------------------------------------------------------------------------
bestModel<-make.gamMRSea(bestModel, panelid=nysted$panelid)

## ----casemean, echo=TRUE, fig.cap='Histogram of the mean of the simulated data, with the red line representing the mean of the original data and the blue line is the mean of the fitted values from the data generation model. The grey lines represent 95% confidence intervals for the simulated mean.', warning=FALSE, message=FALSE, fig.height=6, fig.width=8----
require(ggplot2)
plotMean(bestModel, newdat.ic)


## ----casevar,  warning=FALSE, message=FALSE,echo=TRUE, fig.width=5, fig.height=8, fig.cap='Top: Figure showing the mean-variance relationship for the data generation model (red line) and the mean relationship for all models fitted to the simulated data (grey line), with 95 percentile based confidence intervals for the relationship (grey shading).  The blue line is the assumed mean-variance relationship in the original model used to generate the new data (here overplotting the grey line), while the red line is the original data/model combination. Bottom: Histogram of the estimated dispersion parameters from the simulated data. The blue dashed line represents the dispersion parameter from the data generation model.', fig.width=8, fig.height=4----
plotVariance(bestModel, newdat.ic, n.sim=nsim)

## ------------------------------------------------------------------------
length(badsim.id)
newdat.ic<-newdat.ic[,-badsim.id]
dim(newdat.ic)

## ----caseacf, echo=TRUE, fig.cap='ACF plots for the response data (top left), an example of the simulated data (top right), the pearsons residuals for the data generation model (bottom left) and the residuals for a model fitted to the generated data.', fig.height=6, fig.width=8----
#par(mfrow=c(2,2))
acf(nysted$response, main='Data')
acf(newdat.ic[,10], main='sim data')
acf(residuals(bestModel, type='pearson'), main='original model residuals')
acf(residuals(update(bestModel, newdat.ic[,10] ~ .), type='pearson'), main='sim model residuals')

## ----caseacf2, echo=TRUE, fig.cap='ACF plots for the response data (top left), an example of the simulated data (top right), the pearsons residuals for the data generation model (bottom left) and the residuals for a model fitted to the generated data.', fig.height=6, fig.width=8----
#par(mfrow=c(2,2))
acf(nysted$response, main='Data')
acf(newdat[,10], main='sim data')
acf(residuals(bestModel, type='pearson'), main='original model residuals')
acf(residuals(update(bestModel, newdat[,10] ~ .), type='pearson'), main='sim model residuals')

## ----casequilts, echo=TRUE, fig.height=8, fig.width= 6, fig.cap='Figure showing the spatial distribution of the original data (top left), the fitted values from the model (top right), an example of the simulated data (bottom left) and the mean of all simulated data (bottom right).', warning=FALSE, message=FALSE, fig.height=6, fig.width=8----
require(fields)
par(mfrow=c(2,2))
quilt.plot(nysted$x.pos, nysted$y.pos, nysted$response, main='Original Data', 
           ncol=20, nrow=20, asp=1)
quilt.plot(nysted$x.pos, nysted$y.pos, fitted(bestModel), main='Fitted Values', 
           ncol=20, nrow=20, asp=1)
quilt.plot(nysted$x.pos, nysted$y.pos, newdat.ic[,1], main='Simulated Data Example', 
           ncol=20, nrow=20, asp=1)
quilt.plot(nysted$x.pos, nysted$y.pos, apply(newdat.ic, 1, mean), main='Mean Simulated Data', 
           ncol=20, nrow=20, asp=1)


## ----casebias, echo=TRUE, fig.cap='Figure showing the average simulated value  minus the "true" value under the model (simulating the data) at each location.', fig.height=6, fig.width=8----
par(mfrow=c(1,1))
quilt.plot(nysted$x.pos, nysted$y.pos, apply(newdat.ic, 1, mean) - fitted(bestModel), 
           ncol=20, nrow=20)

## ----echo=TRUE-----------------------------------------------------------
truebeta<-0.80 
impdata.ny<-genChangeData(truebeta*100, model = bestModel, data = nysted, 
                           panels = "panelid")

## ----echo=TRUE, comment='', warning=FALSE, message=FALSE-----------------
require(dplyr)
t<-group_by(impdata.ny, eventphase)%>%
  summarise(sum=sum(truth), mean=mean(truth), n=n())
t

## ----echo=TRUE, warning=FALSE, message=FALSE-----------------------------
rdf1<-make.raster(ncell.y=60, 
                  xyzdata=impdata.ny[impdata.ny$eventphase==0,c("x.pos", "y.pos", "truth")], 
                  z.name = "Mean.count")

rdf2<-make.raster(ncell.y=60, 
                  xyzdata=impdata.ny[impdata.ny$eventphase==0,c("x.pos", "y.pos", "truth")],
                  z.name="Mean.count")

rdf<-rbind(data.frame(rdf1, evph=0), data.frame(rdf2, evph=1))

## ------------------------------------------------------------------------
pct.change<-((rdf$Mean.count[rdf$evph==1] - rdf$Mean.count[rdf$evph==0])/
               rdf$Mean.count[rdf$evph==0])*100

## ---- fig.height=6, fig.width=8------------------------------------------
ggplot( NULL ) + geom_raster( data = rdf1 , aes( x , y , fill = pct.change ) ) +
  scale_fill_gradientn(colours=mypalette,values=c(0,1), space = "Lab", 
                       na.value = "grey50", guide = "colourbar") + 
  theme_bw() + coord_equal()


## ----echo=TRUE, eval=TRUE------------------------------------------------
nysim_glm<-update(bestModel, newdata.ny.imp[,1]~. + eventphase, data=impdata.ny)
nysim_glm$panels<-impdata.ny$panels


## ----warning=FALSE, message=FALSE, echo=TRUE, eval=TRUE------------------
empdistpower.ny<-getEmpDistribution(n.sim = nsim, simData=newdata.ny.imp, 
                                  model = nysim_glm, data = impdata.ny, plot=FALSE,  
                                  returnDist = TRUE, dots=FALSE)

## ----echo=TRUE-----------------------------------------------------------
predictdata<-rbind(data.frame(nysted.predgrid[,2:8], yearmonth='2001/2', 
                              area=nysted.predgrid$area, eventphase=0), 
                   data.frame(nysted.predgrid[,2:8], yearmonth='2001/2', 
                              area=nysted.predgrid$area, eventphase=1))


## ------------------------------------------------------------------------
nsim=50

## ----aukpoweroc, echo=TRUE, eval=FALSE-----------------------------------
#  nysted.power.oc<-powerSimPll(newdatcor.ny.imp, nysim_glm, empdistpower.ny, nsim=nsim,
#                               powercoefid=length(coef(nysim_glm)), predictionGrid=predictdata,
#                               g2k=NULL, splineParams=NULL, n.boot=100, nCores = 2)

## ----echo=TRUE, eval=FALSE-----------------------------------------------
#  save(nysted.power.oc, file='data/powerout.nystedslim.oc.RData', compress = 'bzip2')

## ----eval=TRUE, echo=FALSE-----------------------------------------------
data("nysted.power.oc")

## ----auknull, eval=FALSE, echo=TRUE--------------------------------------
#  nsim=100
#  nysted.power.oc.null<-pval.coverage.null(newdat.ind = newdata.ny.imp,
#                                            newdat.corr = newdatcor.ny.imp, model = nysim_glm,
#                                            nsim = nsim,  powercoefid = length(coef(nysim_glm)))
#  

## ----echo=TRUE, eval=FALSE-----------------------------------------------
#  save(nysted.power.oc.null, file='data/null.output.nystedslim.oc.RData', compress = 'bzip2')

## ----eval=TRUE, echo=FALSE-----------------------------------------------
data("nysted.power.oc.null")

## ----powauksummary, comment=''-------------------------------------------
summary(nysted.power.oc, nysted.power.oc.null, truebeta=log(truebeta))

## ----powaukpowerplot, fig.cap='Figure showing how the power to detect change varies with the error rate chosen for the Forth and Tay redistribution analysis.  The first grey dashed line is at 1% and the second at 5%, traditionally values used as $p$-value cutoffs. The blue dashed lines indicate the error rate required to get a power of 80%.  The value is given in the title.', fig.height=6, fig.width=8----
powerPlot(nysted.power.oc)


## ----powaukpredplot, fig.height=12, fig.width=8, fig.cap='Figure showing the mean (middle), lower 2.5% (top) and upper 97.5% (bottom) of predicted animal counts before (left) and after (right) the event.'----
plot.preds(nysted.power.oc, cellsize=c(1000,1000))

## ----powaukdiffsplot, fig.cap='Figure showing the mean (middle), lower 2.5% (left) and upper 97.5% (right) of estimated differences between before and after the event. (difference = post - pre)', fig.height=4, fig.width=8----
plot.diffs(nysted.power.oc, cellsize=c(1000, 1000))

## ----powauksigdiff, fig.cap='Figure showing, for every grid cell, the proportion of simulations that showed a significant difference.', fig.height=6, fig.width=8----
plot.sigdiff(nysted.power.oc, 
             coordinates = predictdata[predictdata$eventphase==0,c('x.pos', 'y.pos')],
             tailed='two', error.rate = 0.05, gridcell.dim = c(1000,1000))

## ----powauksigdiffsid, fig.cap='Figure showing, for every grid cell, the proportion of simulations that showed a significant difference (with Sidak adjustment for family error rate of 0.05).', fig.height=6, fig.width=8----
plot.sigdiff(nysted.power.oc, 
             coordinates = predictdata[predictdata$eventphase==0,c('x.pos', 'y.pos')],
             tailed='two', error.rate = 0.05, adjustment='sidak', gridcell.dim = c(1000,1000))

## ---- warning=FALSE, message=FALSE---------------------------------------
require(raster)
bndwf_buf<-data.frame(buffer(SpatialPolygons(list(
                Polygons(list(Polygon(nysted.bndwf)),ID = 1))),width=2000)@
                polygons[[1]]@Polygons[[1]]@coords)
names(bndwf_buf)<-c('x.pos', 'y.pos')
detach(package:raster)

## ----warning=FALSE, message=FALSE, echo=TRUE-----------------------------
truebeta<-0.5
impdata.nyre<-genChangeData(truebeta*100, model = bestModel, data = nysted, 
                            panels = "panelid", eventsite.bnd = bndwf_buf)

## ----echo=TRUE, comment='', warning=FALSE, message=FALSE-----------------
impdata.nyre$redistid<-paste(impdata.nyre$eventphase, impdata.nyre$noneventcells)
t<-group_by(impdata.nyre, redistid)%>%
  summarise(sum=sum(truth), mean=mean(truth), n=n())
t

## ----echo=TRUE, warning=FALSE, message=FALSE, fig.height=6, fig.width=8----
rdf1<-make.raster(ncell.y = 60, 
                  xyzdata =impdata.nyre[impdata.nyre$eventphase==0,c('x.pos', 'y.pos', 'truth')], 
                  z.name = 'Mean.count')

rdf2<-make.raster(ncell.y = 60, 
                  xyzdata = impdata.nyre[impdata.nyre$eventphase==1,c('x.pos', 'y.pos', 'truth')], 
                  z.name = 'Mean.count')

rdf<-rbind(data.frame(rdf1,evph=0) , data.frame(rdf2, evph=1))

pct.change<-((rdf$Mean.count[rdf$evph==1] - rdf$Mean.count[rdf$evph==0])/
               rdf$Mean.count[rdf$evph==0])*100

ggplot( NULL ) + geom_raster( data = rdf1 , aes( x , y , fill = pct.change ) ) +
  
  scale_fill_gradientn(colours=mypalette, 
                       values = c(0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6,0.8,1, 10),
                       space = "Lab", na.value = "grey50", 
                       guide = "colourbar") + 
  theme_bw() + coord_equal()


## ----echo=TRUE, eval=TRUE------------------------------------------------
nyresim_glm<-update(bestModel, newdatcor.ny.impre[,1] ~. - as.factor(eventphase) +
                     x.pos + x.pos:as.factor(eventphase), data=impdata.nyre)
                    
nyresim_glm<-make.gamMRSea(nyresim_glm, panelid=1:nrow(impdata.nyre))


## ---- comment=''---------------------------------------------------------
anova(nyresim_glm)

## ----warning=FALSE, message=FALSE, echo=TRUE, eval=TRUE------------------
empdistpower.nyre<-getEmpDistribution(n.sim = nsim, simData=newdata.ny.impre,
                              model = nyresim_glm, data = impdata.nyre, plot=FALSE,  
                              returnDist = TRUE, dots=FALSE)


## ----echo=TRUE-----------------------------------------------------------
predictdata<-rbind(data.frame(nysted.predgrid[,2:8], yearmonth='2002/2', 
                              area=nysted.predgrid$area, eventphase=0), 
                   data.frame(nysted.predgrid[,2:8], yearmonth='2002/2', 
                              area=nysted.predgrid$area, eventphase=1))


## ----echo=TRUE, eval=FALSE-----------------------------------------------
#  save(powerout.nysted.re, file='data/powerout.nystedslim.re.RData', compress = 'bzip2')

## ----echo=TRUE, eval=FALSE-----------------------------------------------
#  save(null.output.nysted.re, file='data/null.output.nystedslim.re.RData', compress = 'bzip2')

## ----powaukresummaryre, comment=''---------------------------------------
summary(powerout.nysted.re, null.output.nysted.re)

## ----powaukrepowerplotre, fig.cap='Figure showing how the power to detect change varies with the error rate chosen for the Forth and Tay redistribution analysis.  The first grey dashed line is at 1% and the second at 5%, traditionally values used as $p$-value cutoffs. The blue dashed lines indicate the error rate required to get a power of 80%.  The value is given in the title.', fig.height=6, fig.width=8----
powerPlot(powerout.nysted.re)

## ----powaukrepredplot, fig.height=12, fig.width=8,  fig.width=8, fig.height=12, fig.cap='Figure showing the mean (middle), lower 2.5% (top) and upper 97.5% (bottom) of predicted animal counts before (left) and after (right) the event.'----
plot.preds(powerout.nysted.re, cellsize=c(1000,1000))

## ----powaukrediffsplotre, fig.cap='Figure showing the mean (middle), lower 2.5% (left) and upper 97.5% (right) of estimated differences between before and after the event. (difference = post - pre)', fig.height=6, fig.width=8----
plot.diffs(powerout.nysted.re, cellsize=c(1000,1000))

## ----powaukresigdiffre, fig.cap='Figure showing, for every grid cell, the proportion of simulations that showed a significant difference.', fig.height=6, fig.width=8----
plot.sigdiff(powerout.nysted.re, 
             coordinates = predictdata[predictdata$eventphase==0,c('x.pos', 'y.pos')],
             tailed='two', error.rate = 0.05, gridcell.dim = c(1000,1000))

## ----powaukresigdiffsidre, fig.cap='Figure showing, for every grid cell, the proportion of simulations that showed a significant difference.', fig.height=6, fig.width=8----
plot.sigdiff(powerout.nysted.re, 
             coordinates = predictdata[predictdata$eventphase==0,c('x.pos', 'y.pos')],
             tailed='two', error.rate = 0.05, adjustment='sidak', gridcell.dim = c(1000,1000))

