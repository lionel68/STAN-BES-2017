############## An introduction to BDA using Stan ##################

############# Outline #########

# Four models will be explored:
# Linear model (gaussian regression) with brms
# Generalized linear model with overdispersion with rstanarm
# Generalized Linear Mixed effect model (with gaussian regression)
# Zero-inflated overdispersed generalized linear model
# Non-linear model with hierarchical structure

############### Load libraries ############
library(rstan)
library(bayesplot)
library(brms)
library(ggplot2)
library(reshape2)
library(plyr) #for adply
library(MASS) #for mvrnorm
library(rstanarm)
#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores()) #this allows Stan to run chains on parallel cores

############## Set working directory #########

setwd("~/Documents/PostDoc_Ghent/Course/Workshop_BES/GitFolder/")

############## Linear model ###########
set.seed(20170927)
## Simulate some data
dat <- data.frame(X1 = runif(100,-2,2),F1 = gl(n = 2,k = 50))
modmat <- model.matrix(~X1*F1,dat)
betas <- runif(dim(modmat)[2],-2,2)
dat$y <- rnorm(100,modmat %*% betas, 1)

## Fit the model
lin_brms <- brm(y ~ X1 * F1, dat) #where are the priors?
prior <- get_prior(y ~ X1 * F1, dat)
#flat priors for all slope parameters, only the sd get a default prior:
curve(dt(x = x,df = 3,ncp=10),0,40)
#define new priors:
priors <- c(prior(normal(0,20),class="Intercept"),prior(cauchy(0,5),class="b"))
#new model
lin_brms2 <- brm(y ~ X1 * F1, dat,prior=priors)
#in rstanarm you could use: stan_lm(y~X1*F1,dat,prior=R2(1))
#you can look at the underlying Stan code using stancode(lin_brms)
#for the pure stan code check: 
#lin_stan <- stan(file = 'model_fitting/Models/normal_model_basic.stan',
#                 data = list(N=nrow(dat),K=ncol(modmat),X=modmat,y=dat$y))

#get the MCMC samples
lin_mcmc <- as.matrix(lin_brms2) 

## Check the model
#convergence checks with traceplot, Rhat and effective sample size
mcmc_trace(as.mcmc(lin_brms2))
rhat(lin_brms2)
neff_ratio(lin_brms2)
#posterior predictive check
pp_check(lin_brms2,type = "dens_overlay",nsamples=100)
pp_check(lin_brms2,type = "stat_2d")
pp_check(lin_brms2,type = "stat")
#could also use shinystan
launch_shiny(lin_brms2)

## Do model inference
#summaries of parameters
mcmc_areas(as.mcmc(lin_brms2),regex_pars="b")
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
predd <- adply(lin_mcmc[rnd_mcmc,],1,function(x) modpred%*%x[1:4]) #note: adply is necessary here to have everything in a melted format
predd <- cbind(predd,pred) #some R magic in the background
names(predd)[1:2] <- c("MCMC_id","predict")
#in this plot with 100 realization of the regression line
ggplot(predd,aes(x=X1,y=predict,color=F1))+geom_path()+
  geom_point(data=dat,aes(x=X1,y=y))

#testing hypothesis
#for instance the hypothesis that the effect of F1 is stronger than the effect of X1
sum(abs(lin_mcmc[,"b_F12"]) > abs(lin_mcmc[,"b_X1"])) / dim(lin_mcmc)[1]
#the probability that observations in group 1 are bigger than observations in group 2
sum(lin_mcmc[,1] > (lin_mcmc[,1] + lin_mcmc[,"b_F12"])) / dim(lin_mcmc)[1]


############# Generalized linear model ############

##### Simulate some overdispersed poisson data
dat <- data.frame(X1=runif(100,-2,2),X2=runif(100,-2,2))
modmat <- model.matrix(~X1*X2,dat)
betas <- c(2,0.1,-0.05,-0.15)
dat$y <- rnbinom(100,mu=exp(modmat %*% betas),size=2.5)

### Fit the model, a standard poisson regression with no variation parameters
#this time fit using rstanarm
poi_arm <- stan_glm(y~X1*X2,dat,family=poisson)
#in brms you could do:
#poi_brms <- brm(y~X1*X2,data = dat,family = poisson)
#you can again look at the underlying Stan code using stancode(poi_brms)
#or use:
#poi_stan <- stan(file = 'model_fitting/Models/poisson_model_basic.stan',
#                 data = list(N=nrow(dat),K=ncol(modmat),X=modmat,y=dat$y))

#check the default priors
prior_summary(poi_arm)


#get the MCMC samples
poi_mcmc <- as.matrix(poi_arm)

## Check the model
#convergence checks with traceplot, Rhat and effective sample size
mcmc_trace(poi_mcmc)
rhat(poi_arm)
neff_ratio(poi_arm)
#posterior predictive check
pp_check(poi_arm,plotfun = "dens_overlay",nreps=100)
pp_check(poi_arm,plotfun = "stat_2d")
#mean is pretty much ok but sd is way lower in predictions than in the data
               

#use quantile regression to get this infos
#ppred <- apply(poi_mcmc[,-5],1,function(x) rpois(100,exp(modmat %*% x))) #compute post predictive distribution
#qrs <- sapply(1:100,function(i) mean(ppred[i,] < dat$y[i])) #compare pp to actual data
#hist(qrs,freq=FALSE,col='grey')
#abline(h=1,col="blue",lty=2,lwd=2) #clear indication for underdispersion in the posterior predicted data

#so clear example of overdispersion in data

### Use overdispersed poisson
nb_arm <- stan_glm.nb(y~X1*X2,dat)
#or use:
#nb_brms <- brm(y~X1*X2,data = dat,family = negbinomial)
#nb_stan <- stan(file = 'model_fitting/Models/neg_binomial_basic.stan',
#                 data = list(N=nrow(dat),K=ncol(modmat),X=modmat,y=dat$y))

#get the MCMC samples
nb_mcmc <- as.matrix(nb_arm)
  
## Check the model
#convergence checks with traceplot, Rhat and effective sample size
mcmc_trace(nb_mcmc)
rhat(nb_arm)
neff_ratio(nb_arm)
#posterior predictive check
pp_check(nb_arm,plotfun = "dens_overlay",nreps=100)
pp_check(nb_arm,plotfun = "stat_2d")

#again some quantile regression
#ppred <- apply(nb_mcmc[,-6],1,function(x) rnbinom(100,mu=exp(modmat %*% x[1:4]),size=x[5]))
#qrs <- sapply(1:100,function(i) mean(ppred[i,] < dat$y[i])) #compare pp to actual data
#hist(qrs,freq=FALSE,col="grey")
#abline(h=1,col="blue",lty=2,lwd=2)


## Do model inference
#summaries of parameters
mcmc_areas(nb_mcmc,pars=c("(Intercept)","X1","X2","X1:X2"))
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
mu <- c(1,2)
sigma <- matrix(c(2,1,1,3),ncol=2)
betas <- mvrnorm(n=10,mu,sigma)
dat$y <- rnorm(100, betas[dat$Group,1] + dat$X1 * betas[dat$Group,2], 1)
#look at simulated data
ggplot(dat,aes(x=X1,y=y,color=Group))+geom_point()+stat_smooth(method="lm",se=FALSE)

### fit the model
hier_brms <- brm(y ~ X1 + (X1 | Group),dat)
#in rstanarm: stan_lmer(y~X1+(X1|Group),dat)
#pure stan code available in three versions:
#a standard centered version
#hier <- stan(file = 'model_fitting/Models/normal_model_allvarying_cent.stan',
#             data=list(N=nrow(dat),K=2,N_group=10,ID_group=as.numeric(dat$Group),X=modmat[,1:2],y=dat$y))
#a version where using the cholesky matrix for the covariance between the varying effects
#hier_chol <- stan(file = 'model_fitting/Models/normal_model_allvarying_centcholesk.stan',
#             data=list(N=nrow(dat),K=2,N_group=10,ID_group=as.numeric(dat$Group),X=modmat[,1:2],y=dat$y))
#a non-centered version with a cholesky factorization
#hier_non <- stan(file = 'model_fitting/Models/normal_model_allvarying_noncent.stan',
#                 data=list(N=nrow(dat),K=2,N_group=10,ID_group=as.numeric(dat$Group),X=modmat[,1:2],y=dat$y))

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

#TODO !!!! but am wondering if this is not too much ... the session is already 300 lines long for 30min

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
#[Maxime] I have some warnings that plotting ignore some  y lines; will look into it in more detail
#some hypothesis testing
#probability that variation in intercept is higher than variation in slopes
sum(hier_mcmc[,4] > hier_mcmc[,5]) / 4000


############ Zero-inflated overdispersed poisson model #############

### simulate some data
dat <- data.frame(X1=runif(100,-2,2))
modmat <- model.matrix(~X1+I(X1^2),dat)
p_i <- rbinom(n = 100,size = 1,prob=0.8)
e_i <- rnorm(100,0,0.33)
lbd_i <- exp(modmat %*% c(2,0.6,-0.6) + e_i)
dat$N <- rpois(100,lbd_i * p_i)
dat$ID <- 1:100
#plot the effect
plot(N ~ X1,dat)

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
#QRS <- sapply(1:100,function(i) mean(ppp[i,] +runif(4000,-0.5,0.5) > dat$N[i]))
#hist(QRS,freq=FALSE,col="grey")
#abline(h=1,col="blue",lty=2,lwd=2)
#evidence for both zero-inflation and overdispersion

### a second zero inflated overdispersed poisson model
zib_brm <- brm(N ~ X1 + I(X1^2) + (1 | ID), dat, family = zero_inflated_poisson)
#no option (yet) in rstanarm but relatively easy to code in Stan or JAGS
#zib <- stan("model_fitting/Models/zero_inflated_overdispersed_poisson.stan",
#            data=list(N=nrow(dat),K=ncol(modmat),X=modmat,y=dat$N))


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
#QRS <- sapply(1:100,function(i) mean(ppp[i,] + runif(4000,-0.5,0.5) > dat$N[i]))
#hist(QRS,freq=FALSE,col="grey")
#abline(h=1,col="blue",lty=2,lwd=2)

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
#[Maxime] I have some warnings that plotting ignore some  y lines; will look into it in more detail

########### Non-linear model with hierarchical structure ############

#data simulation

#prior
prior=c(
  set_prior("uniform(0,500)", class= "b",nlpar="a",coef="Intercept"),set_prior("normal(0,5)", class= "b",nlpar="a",coef="density"),
  set_prior("normal(0,5)", class= "b",nlpar="a",coef="sex"),
  set_prior("normal(200,100)", class= "b",nlpar="b",coef="Intercept"),set_prior("normal(0,5)", class= "b",nlpar="b",coef="density"),
  set_prior("normal(0,5)", class= "b",nlpar="b",coef="sex"),
  set_prior("normal(50,20)", class= "b",nlpar="c",coef="Intercept"),set_prior("normal(0,5)", class= "b",nlpar="c",coef="density")
  set_prior("normal(0,5)", class= "b",nlpar="c",coef="sex")
)
#fitting model
mod=brm(brmsformula(Y~a/(1+exp((b-X)/c)),a~density+sex,b~density*sex,c~density+sex,nl=TRUE),data=data,prior=prior,iter=20000)


########### end of the training session script ###########



######## Code to generate some of the figures in the introduction presentation ######


#try out to make prior, likelihood and posterior distribution from
#a binomial problem
compute_dist <- function(n=10,shape1=2,shape2=2,type="1"){
  library(RColorBrewer)
  set.seed(20171023)
  
  cols <- brewer.pal(3,"Set1")
  #values of parameter p to try out
  p_vals <- seq(0,1,length=100)
  #the prior distribution
  prior <- dbeta(p_vals,shape1 = shape1, shape2 = shape2)
  #did the finger point on earth?
  is_earth <- rbinom(n=1,size=n,prob=0.8)
  #compute the likelihood of the different parameter p values
  lik <- sapply(1:100,function(x) dbinom(is_earth,size = n,prob=p_vals[x]))
  #compute the posterior distribution, making use of the conjugate prior
  post <- dbeta(p_vals,shape1 = shape1 + sum(is_earth),shape2 = shape2 + n - sum(is_earth))
  #std
  #post <- post / sum(post)
  #the output of the function
  out <- data.frame(n=n,p_vals=p_vals,prior=prior,lik=lik,post=post)
  #std for easier plotting
  out$prior <- with(out,(prior - min(prior))/(max(prior)-min(prior)))
  out$lik <- with(out,(lik - min(lik))/(max(lik)-min(lik)))
  out$post <- with(out,(post - min(post))/(max(post)-min(post)))
  par(mar=c(5,5,4,0))
  plot(prior~p_vals,out,type="l",col=cols[1],lwd=4,xlab="Values of parameter",ylab="Density (standardized)",
       main=paste0("Density distribution for a sample size of: ",n[1]),
       ylim=c(0,1.5),cex.lab=2,cex.main=2)
  if(type == "2"){
    lines(out$p_vals,out$lik,col=cols[2],lwd=4,lty=2)
  }
  else if(type == "3"){
    lines(out$p_vals,out$lik,col=cols[2],lwd=4,lty=2)
    lines(out$p_vals,out$post,col=cols[3],lwd=4,lty=3)
  }
  legend("topleft",legend=c("prior","likelihood","posterior"),lwd=4,lty=1:3,bty="n",col=cols,cex=3)
  return()
}

#low sample sizes
png("~/Documents/PostDoc_Ghent/Course/Workshop_BES/GitFolder/introduction/Figures/post1_1.png",width=800,height=800)
compute_dist(n=5,shape1=2,shape2=5)
dev.off()

png("~/Documents/PostDoc_Ghent/Course/Workshop_BES/GitFolder/introduction/Figures/post1_2.png",width=800,height=800)
compute_dist(n=5,shape1=2,shape2=5,type="2")
dev.off()

png("~/Documents/PostDoc_Ghent/Course/Workshop_BES/GitFolder/introduction/Figures/post1_3.png",width=800,height=800)
compute_dist(n=5,shape1=2,shape2=5,type="3")
dev.off()

png("~/Documents/PostDoc_Ghent/Course/Workshop_BES/GitFolder/introduction/Figures/post2.png",width=800,height=800)
compute_dist(n=10,shape1=2,shape2=5)
dev.off()

png("~/Documents/PostDoc_Ghent/Course/Workshop_BES/GitFolder/introduction/Figures/post3.png",width=800,height=800)
compute_dist(n=20,shape1=2,shape2=5,type="3")
dev.off()

png("~/Documents/PostDoc_Ghent/Course/Workshop_BES/GitFolder/introduction/Figures/post4.png",width=800,height=800)
compute_dist(n=100,shape1=2,shape2=5)
dev.off()

png("~/Documents/PostDoc_Ghent/Course/Workshop_BES/GitFolder/introduction/Figures/post5.png",width=800,height=800)
compute_dist(n=20,shape1=1,shape2=1,type="3")
dev.off()

#6 distribution to talk about the prior distribution
png("Documents/PostDoc_Ghent/Course/Workshop_BES/GitFolder/introduction/Figures/prior_choice.png",width=1500,height=900)
par(mfrow=c(2,3),mar=c(4,0.2,0.3,1),mgp=c(3,3,0))

curve(dunif(x,min=-50,max=50),-50,50,ylab="",yaxt="n",cex.axis=4)
legend("topleft",legend="A",bty="n",cex=4)

curve(dnorm(x,0,10),-50,50,ylab="",yaxt="n",cex.axis=4)
legend("topleft",legend="B",bty="n",cex=4)

curve(dnorm(x,-20,10),-50,50,ylab="",yaxt="n",cex.axis=4)
legend("topleft",legend="C",bty="n",cex=4)


curve(dnorm(x,20,10),-50,50,ylab="",yaxt="n",cex.axis=4)
legend("topleft",legend="D",bty="n",cex=4)

curve(dcauchy(x,0,20),-50,50,ylab="",yaxt="n",cex.axis=4)
legend("topleft",legend="E",bty="n",cex=4)

curve(dt(x,1),-50,50,ylab="",yaxt="n",cex.axis=4)
legend("topleft",legend="F",bty="n",cex=4)
dev.off()

#an example of classical MCMC algorithm:
#with the earth-water example
obs <- 3



png("Likelihood_start.png",width=800,height=800)
par(mar=c(4,8,2,0))
curve(dbinom(obs,10,x),0,1,xlab="Parameter p, proportion of earth on Earth",ylab="Likelihood",cex.axis=2,cex.lab=3)
dev.off()

png("Likelihood_init.png",width=800,height=800)
par(mar=c(4,8,2,0))
curve(dbinom(obs,10,x),0,1,xlab="Parameter p, proportion of earth on Earth",ylab="Likelihood",cex.axis=2,cex.lab=3)
abline(v=0.5,col="orange",lwd=3)
dev.off()

png("Likelihood_1.png",width=800,height=800)
par(mar=c(4,8,2,0))
curve(dbinom(obs,10,x),0,1,xlab="Parameter p, proportion of earth on Earth",ylab="Likelihood",cex.axis=2,cex.lab=3)
abline(v=0.5,col="orange",lwd=3)
abline(v=0.2,col="orange",lwd=2,lty=2)
dev.off()

png("Likelihood_2.png",width=800,height=800)
par(mar=c(4,8,2,0))
curve(dbinom(obs,10,x),0,1,xlab="Parameter p, proportion of earth on Earth",ylab="Likelihood",cex.axis=2,cex.lab=3)
abline(v=0.5,col="orange",lwd=3)
abline(v=0.2,col="orange",lwd=3)
dev.off()

png("Likelihood_3.png",width=800,height=800)
par(mar=c(4,8,2,0))
curve(dbinom(obs,10,x),0,1,xlab="Parameter p, proportion of earth on Earth",ylab="Likelihood",cex.axis=2,cex.lab=3)
abline(v=0.5,col="orange",lwd=3)
abline(v=0.2,col="orange",lwd=3)
abline(v=0.15,col="orange",lwd=2,lty=2)
dev.off()

#interpretation earth toss example
dat <- data.frame(IsEarth = rbinom(n = 10,size = 10,prob = 0.3))
mf <- glm(cbind(IsEarth,10-IsEarth) ~ 1, dat, family = binomial)
summary(mf)

