## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(fig=TRUE, warning=FALSE, message=FALSE,
                      eval=TRUE, cache=FALSE, echo=FALSE,
                      comment = '#>', collapse=TRUE, dev='png')

## -----------------------------------------------------------------------------
require(knitcitations)
cleanbib()
#biblio <- read.bibtex("newref.bib")
cite_options(citation_format = 'pandoc', cite.style = 'authoryear', max.names = 1, longnamesfirst=FALSE)

## ----setup--------------------------------------------------------------------
require(MRSeaPower)
require(ggplot2)
require(dplyr)

## -----------------------------------------------------------------------------
nysted <- readRDS("web/data/nystedphaseA.rds")
data("nysted.studybnd")
nysted.studybnd <- nysted.studybnd/1000
data("nysted.bndwf")
nysted.bndwf <- nysted.bndwf/1000

## ----nysteddatplot, fig.cap='Figure showing the survey effort and proposed windfarm site for the Nysted Windfarm.', comment=FALSE, message=FALSE, fig.height=4, fig.width=10, echo=FALSE----
ggplot() + 
  geom_point(data = filter(nysted, Count==0), aes(x.pos, y.pos), col="grey", size=0.7) +
  geom_point(data = filter(nysted, Count>0), aes(x.pos, y.pos, size=log(Count/Area)), colour="red", alpha=1/2) +
  geom_path(data=nysted.bndwf, aes(x.pos, y.pos)) +
  geom_path(data=nysted.studybnd, aes(x.pos, y.pos)) +
  coord_equal() + theme_bw() + facet_wrap(~YearMonth) +
  xlab("Easting") + ylab("Northing")

## -----------------------------------------------------------------------------
initialModel<-MRSea::gamMRSea(Count ~ 1 + as.factor(YearMonth) + Depth + x.pos + y.pos + offset(log(Area)), 
                      data=nysted,
                      family=quasipoisson)

## ----echo=FALSE---------------------------------------------------------------
bestModel<-make.gamMRSea(initialModel, panelid=nysted$TransectID)

impdata.ny<-genChangeData(pct.change = c(60, 90), 
                          model = bestModel, 
                          data = nysted,
                          panels = "TransectID",
                          eventsite.bnd = nysted.bndwf)

## ----warning=FALSE, message=FALSE, echo=FALSE---------------------------------
rdf1<-make.raster(ncell.y=50, 
                  xyzdata=impdata.ny %>% filter(eventphase==0) %>%
                    dplyr::select("x.pos", "y.pos", "truth"), 
                  z.name = "Mean.count")

rdf2<-make.raster(ncell.y=50, 
                  xyzdata=filter(impdata.ny, eventphase==1) %>%
                    dplyr::select("x.pos", "y.pos", "truth"),
                  z.name="Mean.count")

rdf<-rbind(data.frame(rdf1, evph=0), data.frame(rdf2, evph=1))

pct.change<-((rdf$Mean.count[rdf$evph==1] - rdf$Mean.count[rdf$evph==0])/
               rdf$Mean.count[rdf$evph==0])*100

## ----fig.height=6, fig.width=10, fig.cap="Figure showing the percentage change imposed for each grid cell."----
ggplot( NULL ) + 
  geom_raster( data = rdf1 , aes( x , y , fill = pct.change ) ) +
  scale_fill_gradientn(colours=mypalette, values=c(0,1), space = "Lab", 
                       na.value = "grey50", guide = "colourbar", name="Pct Difference") + 
  theme_bw() + coord_equal() +
  geom_path(data = nysted.bndwf, aes(x.pos, y.pos)) +
  geom_path(data = nysted.studybnd, aes(x.pos, y.pos)) + 
  xlab("Easting") + ylab("Northing")

## ----fig.height=4, fig.width=10, fig.cap="Figure showing the mean count for the baseline (0) and post impact (1)."----
ggplot() +
  geom_raster(data = rdf , aes(x , y , fill = Mean.count)) +
  facet_wrap(~evph) +
  coord_equal() +
  scale_fill_distiller(palette = "Spectral",name="Animal Counts") +
  geom_path(data = nysted.bndwf, aes(x.pos, y.pos)) +
  geom_path(data = nysted.studybnd, aes(x.pos, y.pos))+ 
  xlab("Easting") + ylab("Northing") + theme_bw()

## -----------------------------------------------------------------------------
nsim=500
# add noise
newdata.ny.imp<-generateNoise(nsim, 
                              impdata.ny$truth, family='poisson', 
                               d=summary(bestModel)$dispersion)


## -----------------------------------------------------------------------------
corrs<-getCorrelationMat(panel = impdata.ny$TransectID, data=impdata.ny$truth, dots = FALSE)

# quicker to do pre only and then bind the correlation matrices
corrs<-getCorrelationMat(panel = nysted$TransectID, data=nysted$Count, dots = FALSE)

# add correlation
newdatcor.ny.imp<-generateIC(data = impdata.ny, 
                             corrs = rbind(corrs, corrs), 
                            panels = 'panels', 
                            newdata = newdata.ny.imp, 
                            nsim = nsim, 
                              dots = FALSE)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
nysim_glm<-update(bestModel, newdatcor.ny.imp[,1]~. + eventphase, data=impdata.ny)
nysim_glm$panels<-impdata.ny$TransectID

## ----echo=FALSE---------------------------------------------------------------
nsim=300
empdistpower.ny<-getEmpDistribution(n.sim = nsim,
                                    simData=newdata.ny.imp, 
                                    model = nysim_glm, 
                                    data = impdata.ny, 
                                    plot=FALSE, returnDist = TRUE,
                                    dots=FALSE)

data("nysted.predgrid")
nysted.predgrid <- nysted.predgrid %>%
  rename(Area = area) %>%
  mutate(Depth = abs(Depth),
         x.pos = x.pos/1000,
         y.pos = y.pos/1000)

predictdata<-rbind(data.frame(nysted.predgrid[,2:8], YearMonth='2001/2', 
                              Area=nysted.predgrid$Area, eventphase=0), 
                   data.frame(nysted.predgrid[,2:8], YearMonth='2001/2', 
                              Area=nysted.predgrid$Area, eventphase=1))

powerout.nysted.re<-powerSimPll(newdat = newdatcor.ny.imp, 
                                model = nysim_glm, 
                                empdistribution = empdistpower.ny, 
                             nsim=nsim, powercoefid=length(coef(nysim_glm)), 
                             predictionGrid=predictdata, 
                             n.boot=10, nCores = 2)

## ----powaukpowerplot, fig.cap='Figure showing how the power to detect change varies with the error rate chosen.  The first grey dashed line is at 1% and the second at 5%, traditionally values used as $p$-value cutoffs. The blue dashed lines indicate the error rate required to get a power of 80%.  The value is given in the title.', fig.height=6, fig.width=8----
powerPlot(powerout.nysted.re)

## ----powaukpredplot, fig.height=8, fig.width=8, fig.cap='Figure showing the mean (middle), lower 2.5% (top) and upper 97.5% (bottom) of predicted animal counts before (left) and after (right) the event.', echo=FALSE----
plot.preds(powerout.nysted.re, cellsize=c(1,1))

## ----powaukdiffsplot, fig.cap='Figure showing the mean (middle), lower 2.5% (left) and upper 97.5% (right) of estimated differences between before and after the event. (difference = post - pre)', fig.height=4, fig.width=10, echo=FALSE----
plot.diffs(powerout.nysted.re, cellsize=c(1, 1)) 

## ----powaukresigdiffre, fig.cap='Figure showing, for every grid cell, the proportion of simulations that showed a significant difference.', fig.height=4, fig.width=8, echo=FALSE----
plot.sigdiff(powerout.nysted.re, 
             coordinates = predictdata[predictdata$eventphase==0,c('x.pos', 'y.pos')],
             tailed='two', error.rate = 0.05, gridcell.dim = c(1,1))

