################################
# Session on model fitting
# with models of increasing complexities
###############################


############### Load libraries ############
library(rstan)
library(bayesplot)
library(brms)
library(ggplot2)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

############## Linear model ###########

## Simulate some data
dat <- data.frame(X1 = runif(100,-2,2),F1 = gl(n = 2,k = 50))
modmat <- model.matrix(~X1*F1,dat)
betas <- runif(dim(modmat)[2],-2,2)
dat$y <- rnorm(100,modmat %*% betas, 1)

## Fit the model
lin_brms <- brm(y ~ X1 * F1, dat)
#get the MCMC samples
lin_mcmc <- as.mcmc(lin_brms)
#put it into a standard matrix
lin_mcmc_m <- as.matrix(lin_mcmc) #to test

## Check the model
#convergence checks with traceplot, Rhat and effective sample size
mcmc_trace(as.mcmc(lin_brms))
rhat(lin_brms)
neff_ratio(lin_brms)
#posterior predictive check
pp_check(lin_brms,type = "dens_overlay",nsamples=100)
pp_check(lin_brms,type = "stat_2d")
pp_check(lin_brms,type = "stat")

## Do model inference
#summaries of parameters
mcmc_areas(as.mcmc(lin_brms),regex_pars="b")
#standard regression lines with 95% credible intervals
#create a new data frame where prediction should be computed
pred <- expand.grid(X1=seq(-2,2,0.1),F1=factor(1:2))
#turn this into a model matrix for easier computation
modpred <- model.matrix(~X1*F1, pred)
#MCMC linear prediction
est <- apply(lin_mcmc_m[,1:4],1,function(x) modpred %*% x)
#CrI of the linear predictor/ expected values
pred$LCI <- apply(est,1,quantile,probs=0.025)
pred$Med <- apply(est,1,quantile,probs=0.5)
pred$UCI <- apply(est,1,quantile,probs=0.975)
#a plot
ggplot(dat,aes(x=X1,y=y,color=F1))+geom_point()+
  geom_line(data=pred,aes(y=LCI),linetype="dashed")+
  geom_line(data=pred,aes(y=UCI),linetype="dashed")+
  geom_line(data=pred,aes(y=Med))

#another way of looking at this is too sample 100 MCMC values for
#all parameters and ocmpute the regression lines
rnd_mcmc <- sample(1:4000,100,replace=FALSE)
predd <- adply(lin_mcmc_m[rnd_mcmc,],1,function(x) modpred%*%x[1:4])
predd <- cbind(predd,pred) #some R magic in the background
names(predd)[1:2] <- c("MCMC_id","predict")
ggplot(predd,aes(x=X1,y=predict,color=F1))+geom_path()+
  geom_point(data=dat,aes(x=X1,y=y))

#posterior predictive distribution for each data points
predI <- apply(lin_mcmc_m[rnd_mcmc,1:5],1,function(x) rnorm(100,modmat %*% x[-5],x[5]))
predD <- melt(predI)

ggplot(predD,aes(x=Var1,y=value))+geom_line(aes(group=Var2))+
  geom_point(data=dat,aes(x=1:100,y=y),color="red",size=3)+
  labs(x="Row index",y="Y value")

#testing hypothesis
#for instance the hypothesis that the effect of F1 is stronger than the effect of X1
sum(abs(lin_mcmc_m[,"F12"]) > abs(lin_mcmc_m[,"X1"])) / dim(lin_mcmc_m)[1]
#the probability that observations in group 1 are bigger than obsevations in group 2
sum(lin_mcmc_m[,1] > (lin_mcmc_m[,1] + lin_mcmc_m[,"F12"])) / dim(lin_mcmc_m)[1]


############# Generalized linear model ############

##### Simulate some overdispersed poisson data
dat <- data.frame(X1=runif(100,-2,2),X2=runif(100,-2,2))
modmat <- model.matrix(~X1*X2,dat)
betas <- c(2,runif(dim(modmat)[2]-1,-0.75,0.75))
dat$y <- rnbinom(100,mu=exp(modmat %*% betas),size=2.5)

### Fit the model
poi_brms <- brm(y~X1*X2,data = dat,family = poisson)

#get the MCMC samples
poi_mcmc <- as.matrix(poi_brms)

## Check the model
#convergence checks with traceplot, Rhat and effective sample size
mcmc_trace(as.mcmc(poi_brms))
rhat(poi_brms)
neff_ratio(poi_brms)
#posterior predictive check
pp_check(poi_brms,type = "dens_overlay",nsamples=100)
pp_check(poi_brms,type = "stat_2d")
#mean is pretty much ok but sd is way higher than the data
#use quantile regession to get this infos
ppred <- apply(poi_mcmc[,-5],1,function(x) rpois(100,exp(modmat %*% x))) #compute post predictive distribution
qrs <- sapply(1:100,function(i) mean(ppred[i,] < dat$y[i])) #compare pp to actual data
hist(qrs,freq=FALSE,col='grey')
abline(h=1,col="blue",lty=2,lwd=2) #clear indication for underdispersion in the posterior predicted data

#clear example of overdispersion

### Use overdispersed poisson
nb_brms <- brm(y~X1*X2,data = dat,family = negbinomial)

#get the MCMC samples
nb_mcmc <- as.matrix(nb_brms)
  
## Check the model
#convergence checks with traceplot, Rhat and effective sample size
mcmc_trace(as.mcmc(nb_brms))
rhat(nb_brms)
neff_ratio(nb_brms)
#posterior predictive check
pp_check(nb_brms,type = "dens_overlay",nsamples=100)
pp_check(nb_brms,type = "stat_2d")
#again some quantile regression
ppred <- apply(nb_mcmc[,-6],1,function(x) rnbinom(100,mu=exp(modmat %*% x[1:4]),size=x[5]))
qrs <- sapply(1:100,function(i) mean(ppred[i,] < dat$y[i])) #compare pp to actual data
hist(qrs,freq=FALSE,col="grey")
abline(h=1,col="blue",lty=2,lwd=2)


## Do model inference
#summaries of parameters
mcmc_areas(as.mcmc(nb_brms),regex_pars="b")
#standard regression lines with 95% credible intervals
#create a new data frame where prediction should be computed
pred <- expand.grid(X1=seq(-2,2,0.1),X2=c(-2,0,2))
#turn this into a model matrix for easier computation
modpred <- model.matrix(~X1*X2, pred)
#MCMC linear prediction
est <- apply(nb_mcmc_m[,1:4],1,function(x) modpred %*% x)
#CrI of the linear predictor/ expected values
pred$LCI <- exp(apply(est,1,quantile,probs=0.025))
pred$Med <- exp(apply(est,1,quantile,probs=0.5))
pred$UCI <- exp(apply(est,1,quantile,probs=0.975))
#a plot
ggplot(dat,aes(x=X1,y=y,group=X2))+geom_point()+
  geom_line(data=pred,aes(y=LCI,color=factor(X2)),linetype="dashed")+
  geom_line(data=pred,aes(y=UCI,color=factor(X2)),linetype="dashed")+
  geom_line(data=pred,aes(y=Med,color=factor(X2)))


############ Generalized Linear Mixed effect model #########
### Simulate some data
dat <- data.frame(X1=runif(100,-2,2),Group=gl(n=10,k=10))
modmat <- model.matrix(~X1*Group,dat)
betas <- c(1,2,rnorm(9,0,1),rnorm(9,0,0.1))
dat$y <- rnorm(100, modmat %*% betas, 1)
#look at simulated data
ggplot(dat,aes(x=X1,y=y,color=Group))+geom_point()+stat_smooth(method="lm",se=FALSE)

### fit the model
hier_brms <- brm(y ~ X1 + (X1 | Group),dat)
#get MCMC samples
hier_mcmc <- as.matrix(hier_brms)

### Check the model
#convergence checks with traceplot, Rhat and effective sample size
mcmc_trace(as.mcmc(hier_brms))
rhat(hier_brms)
neff_ratio(hier_brms)
#posterior predictive check
pp_check(hier_brms,type = "dens_overlay",nsamples=100)
pp_check(hier_brms,type = "stat_2d")
#looks good
## Do model inference
#look at shrinkage effects
#TODO !!!!

#summaries of parameters
mcmc_areas(as.mcmc(hier_brms),regex_pars=c("b","sigma","sd"))

#plot regression line with credible intervals for an average group
#unconditional 
pred <- expand.grid(X1=seq(-2,2,0.1))
modpred <- model.matrix(~X1,pred)
predcond <- apply(hier_mcmc,1,function(x) modpred %*% x[1:2])
condint <- apply(predcond,1,quantile,probs=c(0.025,0.5,0.975))
pred <- data.frame(X1=pred$X1,LCI=condint[1,],Med=condint[2,],UCI=condint[3,],Group=factor(13))

ggplot(dat,aes(x=X1,y=y,color=Group))+geom_point()+
  geom_line(data=pred,aes(y=Med),color="black")+
  geom_ribbon(data=pred,aes(y=Med,ymin=LCI,ymax=UCI),alpha=0.2,color="grey")

#include group-level uncertainty, so conditional on groups
predgroup <- apply(hier_mcmc,1,function(x) rnorm(1,x[1],x[3]) + rnorm(1,x[2],x[4]) * pred$X1)
groupint <- apply(predgroup,1,quantile,probs=c(0.025,0.5,0.975))
pred <- cbind(pred,data.frame(LCIg=groupint[1,],Medg=groupint[2,],UCIg=groupint[3,]))

ggplot(dat,aes(x=X1,y=y,color=Group))+geom_point()+
  geom_line(data=pred,aes(y=Med),color="black")+
  geom_ribbon(data=pred,aes(y=Med,ymin=LCI,ymax=UCI),alpha=0.2,color="grey")+
  geom_ribbon(data=pred,aes(y=Medg,ymin=LCIg,ymax=UCIg),alpha=0.2,color="grey70")

