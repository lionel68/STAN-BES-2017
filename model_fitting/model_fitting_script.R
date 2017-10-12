############## An introduction to BDA using STAN ##################
#TO DO: make the relevant .stan code for each 4 models

############# Outline #########

# Four models will be explored:
# Linear model (gaussian regression)
# Generalized linear model with overdispersion 
# Generalized Linear Mixed effect model (with gaussian regression)
# Zero-inflated overdispersed generalized linear model

############### Load libraries ############
library(rstan)
library(bayesplot)
library(brms)
library(ggplot2)
library(reshape2)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) #this allows STAN to run chains on parallel cores

############## Linear model ###########
set.seed(20170927)
## Simulate some data
dat <- data.frame(X1 = runif(100,-2,2),F1 = gl(n = 2,k = 50))
modmat <- model.matrix(~X1*F1,dat)
betas <- runif(dim(modmat)[2],-2,2)
dat$y <- rnorm(100,modmat %*% betas, 1)

## Fit the model
lin_brms <- brm(y ~ X1 * F1, dat)

#in rstanarm you could use: stan_lm(y~X1*F1,dat)

#you can look at the underlying STAN code using stancode(lin_brms)
#for the pure stan code check: 
#lin_stan <- stan(file = '/media/lionel/USB_Lio/PostDoc/Workshop_BES/GitFolder/model_fitting/Models/normal_model_basic.stan',
#                 data = list(N=nrow(dat),K=ncol(modmat),X=modmat,y=dat$y))

#get the MCMC samples
lin_mcmc <- as.matrix(lin_brms) #to test

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
est <- apply(lin_mcmc[,1:4],1,function(x) modpred %*% x)
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
#all parameters and compute the regression lines
rnd_mcmc <- sample(1:4000,100,replace=FALSE)
predd <- apply(lin_mcmc[rnd_mcmc,],1,function(x) modpred%*%x[1:4])
predd <- cbind(predd,pred) #some R magic in the background
names(predd)[1:2] <- c("MCMC_id","predict")
ggplot(predd,aes(x=X1,y=predict,color=F1))+geom_path()+
  geom_point(data=dat,aes(x=X1,y=y))

#posterior predictive distribution for each data points
predI <- apply(lin_mcmc[rnd_mcmc,1:5],1,function(x) rnorm(100,modmat %*% x[-5],x[5]))
predD <- melt(predI)

#[Maxime]plot below doesn't work as there is no Var2 in predD
ggplot(predD,aes(x=Var1,y=value))+geom_line(aes(group=Var2))+
  geom_point(data=dat,aes(x=1:100,y=y),color="red",size=3)+
  labs(x="Row index",y="Y value")

#testing hypothesis
#for instance the hypothesis that the effect of F1 is stronger than the effect of X1
sum(abs(lin_mcmc[,"b_F12"]) > abs(lin_mcmc[,"b_X1"])) / dim(lin_mcmc)[1]
#the probability that observations in group 1 are bigger than observations in group 2
#[Maxime] doesn't work, no object named "lin_mcmc_m"
sum(lin_mcmc_m[,1] > (lin_mcmc_m[,1] + lin_mcmc_m[,"b_F12"])) / dim(lin_mcmc_m)[1]


############# Generalized linear model ############

##### Simulate some overdispersed poisson data
dat <- data.frame(X1=runif(100,-2,2),X2=runif(100,-2,2))
modmat <- model.matrix(~X1*X2,dat)
betas <- c(2,0.1,-0.05,-0.15)
dat$y <- rnbinom(100,mu=exp(modmat %*% betas),size=2.5)

### Fit the model, a standard poisson regression with no variation parameters
poi_brms <- brm(y~X1*X2,data = dat,family = poisson)
#in rstanarm you could do: stan_glm(y~X1*X2,dat,family=poisson)
#you can again look at the underlying STAN code using stancode(poi_brms)
#or pure stan code is available under: poisson_model_basic.stan


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
#in rstanarm: stan_glm.nb(y~X1*X2,dat)
#the pure stan code is available under: neg_binomial_basic.stan

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
#MCMC linear prediction with a log-link
est <- apply(nb_mcmc[,1:4],1,function(x) exp(modpred %*% x))
#CrI of the linear predictor/ expected values
pred$LCI <- apply(est,1,quantile,probs=0.025)
pred$Med <- apply(est,1,quantile,probs=0.5)
pred$UCI <- apply(est,1,quantile,probs=0.975)
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
#in rstanarm: stan_lmer(y~X1+(X1|Group),dat)
#pure stan code available in two versions: normal_model_allvarying_cent.stan
#hier <- stan(file = '/media/lionel/USB_Lio/PostDoc/Workshop_BES/GitFolder/model_fitting/Models/normal_model_allvarying_cent.stan',
#             data=list(N=nrow(dat),K=2,N_group=10,ID_group=as.numeric(dat$Group),X=modmat[,1:2],y=dat$y))
#hier_non <- stan(file = '/media/lionel/USB_Lio/PostDoc/Workshop_BES/GitFolder/model_fitting/Models/normal_model_allvarying_noncent.stan',
#                 data=list(N=nrow(dat),K=2,N_group=10,ID_group=as.numeric(dat$Group),X=modmat[,1:2],y=dat$y))
#the non-centered approach still does not work

#and normal_model_allvarying_noncent.stan

#get MCMC samples
hier_mcmc <- as.matrix(hier_brms)

### Check the model
#convergence checks with traceplot, Rhat and effective sample size
mcmc_trace(as.mcmc(hier_brms),regex_pars = "b|sd|sigma")
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
                  #[Maxime]geom_ribbon is giving warnings that it ignores argument y=, here and below

#include group-level uncertainty, so conditional on groups
predgroup <- apply(hier_mcmc,1,function(x) rnorm(1,x[1],x[3]) + rnorm(1,x[2],x[4]) * pred$X1)
groupint <- apply(predgroup,1,quantile,probs=c(0.025,0.5,0.975))
pred <- cbind(pred,data.frame(LCIg=groupint[1,],Medg=groupint[2,],UCIg=groupint[3,]))

ggplot(dat,aes(x=X1,y=y,color=Group))+geom_point()+
  geom_line(data=pred,aes(y=Med),color="black")+
  geom_ribbon(data=pred,aes(y=Med,ymin=LCI,ymax=UCI),alpha=0.2,color="grey")+
  geom_ribbon(data=pred,aes(y=Medg,ymin=LCIg,ymax=UCIg),alpha=0.2,color="grey70")

#some hypothesis testing
#probability that variation in intercept is higher than variation in slopes
sum(hier_mcmc[,4] > hier_mcmc[,5]) / 4000


############ Zero-inflated overdispersed poisson model #############

### simulate some data
dat <- data.frame(X1=runif(100,-2,2))
modmat <- model.matrix(~X1+I(X1^2),dat)
p_i <- rbinom(n = 100,size = 1,prob=0.7)
e_i <- rnorm(100,0,0.33)
lbd_i <- exp(modmat %*% c(2,0.5,-0.4) +e_i)
dat$N <- rpois(100,lbd_i * p_i)

### a first naive poisson model
poi_brm <- brm(N ~ X1 + I(X1^2),dat,family=poisson)
#in rstanarm: stan_glm(N~X1+I(X^2),dat,family=poisson)

#get the MCMC
poi_mcmc <- as.matrix(poi_brm)
### Check the model
#convergence checks with traceplot, Rhat and effective sample size
mcmc_trace(as.mcmc(poi_brm))
rhat(poi_brm)
neff_ratio(poi_brm)
#posterior predictive check
pp_check(poi_brm,type = "dens_overlay",nsamples=100)
pp_check(poi_brm,type = "stat_2d")
#some issue there
#count how many 0s there is vs number of zeros predicted
ppp <- apply(poi_mcmc,1,function(x) rpois(100,exp(modmat%*%x[1:3])))
obs_0 <- sum(dat$N==0)
hist(apply(ppp==0,2,sum),xlim=c(0,obs_0 + 5))
abline(v=obs_0,col="orange",lwd=2)
#look at the QRS
QRS <- sapply(1:100,function(i) mean(ppp[i,] +runif(4000,-0.5,0.5) > dat$N[i]))
hist(QRS,freq=FALSE,col="grey")
abline(h=1,col="blue",lty=2,lwd=2)
#evidence for both zero-inflation and overdispersion

### a second zero inflated overdispersed poisson model
dat$ID <- 1:nrow(dat)
zib_brm <- brm(N ~ X1 + I(X1^2) + (1 | ID), dat, family = zero_inflated_poisson)
#no option (yet) in rstanarm but relatively easy to code in STAN or JAGS
#TODO: make the STAN code for this one

zib_mcmc <- as.matrix(zib_brm)
### Check the model
#convergence checks with traceplot, Rhat and effective sample size
mcmc_trace(as.mcmc(zib_brm),regex_pars="b|sd|zi")
#posterior predictive check
pp_check(zib_brm,type = "dens_overlay",nsamples=100)
pp_check(zib_brm,type = "stat_2d")
#compute the ppp for 0s and QRS check
ppp <- apply(zib_mcmc,1,function(x){
  lbd <- exp(modmat %*% x[1:3] + rnorm(1,0,x[4]))
  ps <- rbinom(100,1,1-x[5])
  return(rpois(100,lbd*ps))
})
obs_0 <- sum(dat$N==0)
hist(apply(ppp==0,2,sum))
abline(v=obs_0,col="orange",lwd=2)
#look at the QRS
QRS <- sapply(1:100,function(i) mean(ppp[i,] +runif(4000,-0.5,0.5) > dat$N[i]))
hist(QRS,freq=FALSE,col="grey")
abline(h=1,col="blue",lty=2,lwd=2)

#look at predicted regression with credible and predicted intervals
pred <- data.frame(X1=seq(-2,2,0.1))
modpred <- model.matrix(~X1+I(X1^2),pred)
reg_pred <- apply(zib_mcmc,1,function(x) exp(modpred %*% x[1:3] +rnorm(1,0,x[4])))
reg_ci <- apply(reg_pred,1,quantile,probs=c(0.05,0.5,0.95))
pred$LCI <- reg_ci[1,]
pred$Med <- reg_ci[2,]
pred$UCI <- reg_ci[3,]
pred_pred <- sapply(1:4000,function(i) rpois(41,reg_pred[,i] * rbinom(n = 1,size = 1,prob=1-zib_mcmc[i,5])))
pred_ci <- apply(pred_pred,1,quantile,probs=c(0.05,0.5,0.95))
pred$LCI_pred <- pred_ci[1,]
pred$Med_pred <- pred_ci[2,]
pred$UCI_pred <- pred_ci[3,]

ggplot(dat,aes(x=X1,y=N))+geom_point()+
  geom_line(data=pred,aes(y=Med))+
  geom_ribbon(data=pred,aes(y=Med,ymin=LCI,ymax=UCI),color="grey30",alpha=0.5)+
  geom_ribbon(data=pred,aes(y=Med,ymin=LCI_pred,ymax=UCI_pred),color="grey10",alpha=0.2)+
  labs(x="Environmental gradient",y="Counts")
