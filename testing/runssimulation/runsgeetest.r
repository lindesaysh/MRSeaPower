devtools::load_all(pkg = 'C:/MarineScotlandPower/MRSea/MRSea')
devtools::load_all(pkg = 'C:/MarineScotlandPower/MRSeaPower')

require(geepack)

# ~~~~~~~~~~~~~
# independent, panel size =10
# dispersion = 15
# ~~~~~~~~~~~~

nsim=500
nsim2<- 500
truebeta<-log(0.9)

dat<-makeToyData(nsim, truebeta, b0=-2,length.panels = 20)

model.power.raw1=model.power.rob1=model.power.gee1=model.power.geeind1 = model.coverage1=vector(length=100)

for(i in 1:100){
  print(i)
  noisydat<-generateNoise(nsim, dat$mu, family='poisson', d=15)
  init_glm<-glm(noisydat[,1]~x+evph, data=dat, family='quasipoisson')
  
  #empdistribution<-getRunsCritVals(n.sim = nsim, simData=noisydat,  model = init_glm, data = dat, plot=TRUE, returnDist = TRUE, dots=FALSE)
  
  poweroutruns<-powerSimRunsTest(noisydat, init_glm, nsim=nsim2, powercoefid = 3)
  
  
  model.power.raw1[i]<-(length(which(poweroutruns$imppvals[,1]<=0.05))/nsim2)*100
  model.power.rob1[i]<-(length(which(poweroutruns$imppvals[,2]<=0.05))/nsim2)*100
  model.power.gee1[i]<-(length(which(poweroutruns$imppvals[,3]<=0.05))/nsim2)*100
  model.power.geeind1[i]<-(length(which(poweroutruns$imppvals[,4]<=0.05))/nsim2)*100
  model.coverage1[i]<-impact.coverage(truebeta=truebeta, poweroutruns$betacis)*100
  
}


model.power.raw1
model.power.rob1
model.power.gee1
model.power.geeind1
#model.coverage1

wilcox.test(model.power.raw1, model.power.rob1, paired = TRUE)
mean(model.power.raw1- model.power.rob1)

wilcox.test(model.power.raw1, model.power.gee1, paired = TRUE)
mean(model.power.raw1- model.power.gee1)

wilcox.test(model.power.raw1, model.power.geeind1, paired = TRUE)
mean(model.power.raw1- model.power.geeind1)
