devtools::load_all(pkg = 'C:/MarineScotlandPower/MRSea/MRSea')
devtools::load_all(pkg = 'C:/MarineScotlandPower/MRSeaPower')



# ~~~~~~~~~~~~~
# independent, panel size =1
# dispersion = 15
# ~~~~~~~~~~~~
dat<-makeToyData(nsim, truebeta, b0=-2,length.panels = 1)

nsim=500
nsim2<- 500
truebeta<-log(0.9)

model.power.raw1=model.power.rob1 = model.coverage1=vector(length=100)
for(i in 1:100){
  noisydat<-generateNoise(nsim, dat$mu, family='poisson', d=15)
  init_glm<-glm(noisydat[,1]~x+evph, data=dat, family='quasipoisson')
  
  #empdistribution<-getRunsCritVals(n.sim = nsim, simData=noisydat,  model = init_glm, data = dat, plot=TRUE, returnDist = TRUE, dots=FALSE)
  
  poweroutruns<-powerSimRunsTest(noisydat, init_glm, nsim=nsim2, powercoefid = 3)
  
  
  model.power.raw1[i]<-(length(which(poweroutruns$imppvals[,1]<=0.05))/nsim2)*100
  model.power.rob1[i]<-(length(which(poweroutruns$imppvals[,2]<=0.05))/nsim2)*100
  model.coverage1[i]<-impact.coverage(truebeta=truebeta, poweroutruns$betacis)*100
  
}


model.power.raw1
model.power.rob1
#model.coverage1

wilcox.test(model.power.raw1, model.power.rob1, paired = TRUE)
mean(model.power.raw1- model.power.rob1)


# ~~~~~~~~~~~~~
# independent, panel given when not necessary
# dispersion = 15
# ~~~~~~~~~~~~
dat2<-makeToyData(nsim, truebeta, b0=-2,length.panels = 10)

nsim=500
nsim2<- 500
truebeta<-log(0.9)

model.power.raw2=model.power.rob2 = model.coverage2=vector(length=100)
for(i in 1:100){
  noisydat<-generateNoise(nsim, dat2$mu, family='poisson', d=15)
  init_glm<-glm(noisydat[,1]~x+evph, data=dat2, family='quasipoisson')
  
  #empdistribution<-getRunsCritVals(n.sim = nsim, simData=noisydat,  model = init_glm, data = dat2, plot=TRUE, returnDist = TRUE, dots=FALSE)
  
  poweroutruns<-powerSimRunsTest(noisydat, init_glm, nsim=nsim2, powercoefid = 3)
  
  
  model.power.raw2[i]<-(length(which(poweroutruns$imppvals[,1]<=0.05))/nsim2)*100
  model.power.rob2[i]<-(length(which(poweroutruns$imppvals[,2]<=0.05))/nsim2)*100
  model.coverage2[i]<-impact.coverage(truebeta=truebeta, poweroutruns$betacis)*100
  
}


model.power.raw2
model.power.rob2
model.coverage2

wilcox.test(model.power.raw2, model.power.rob2, paired = T)
mean(model.power.raw2- model.power.rob2)

# ~~~~~~~~~~~~~
# independent, panel given when not necessary
# dispersion = 15
# ~~~~~~~~~~~~
data("fowshco")

nsim=500
nsim2<- 500
truebeta<-log(0.9)

# initial model
init_glm<-glm(response ~  as.factor(TideState) + WindStrength + SeaState + SimpPrecipitation + CloudCover + Depth + offset(log(areatime)), data=dat, family=quasipoisson)
# update to have panels even though data generated will be independent
init_glm<-make.gamMRSea(init_glm, panelid=dat$newbidNum)

model.power.raw3=model.power.rob3 = model.coverage3=vector(length=100)

for(i in 26:100){
  print(i)
  # induce change
  impdat<-genOverallchangeData(truebeta, model = init_glm, data = dat, panels = 'newbidNum')
  # generate noisy independent data
  noisydat<-generateNoise(n = nsim, response = impdat$truth, family = 'poisson', d=23)
  
  # update initial model to include impact
  imp_glm<-update(init_glm, noisydat[,1] ~ . + eventphase, data=impdat)
  
  # calculate empirical distribution
  #empdistribution<-getRunsCritVals(n.sim = nsim, simData=noisydat,  model = imp_glm, data = impdat, plot=TRUE, returnDist = TRUE, dots=FALSE)

  imp_glm<-make.gamMRSea(imp_glm, panelid=impdat$panels)
  
  # run power analysis to harvest p-values for eventphase
  poweroutruns<-powerSimRunsTest(noisydat, imp_glm, nsim=nsim2, powercoefid = length(coef(imp_glm)))
  
  
  # store power outputs
  model.power.raw3[i]<-(length(which(poweroutruns$imppvals[,1]<=0.05))/nsim2)*100
  model.power.rob3[i]<-(length(which(poweroutruns$imppvals[,2]<=0.05))/nsim2)*100
  model.coverage3[i]<-impact.coverage(truebeta=truebeta, poweroutruns$betacis)*100
}

model.power.raw3
model.power.rob3
model.coverage3


wilcox.test(model.power.raw3, model.power.rob3, paired = T)
mean(model.power.raw3- model.power.rob3)

