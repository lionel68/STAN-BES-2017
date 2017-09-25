##Note: for me, we have too many examples in this part
##I'd focus on 4 examples at most (LM, GLM poisson or binomial, mixed-effects model with varying intercepts only, 
#and one more complex where the usefulness of Bayesian is obvious

#if we show details for all 4, compare to Freq approach at low and realistic sample sizes, 
#and provide for some syntax in rstan/brms/stanarm for comparison (while focusing on only one obviously)
#it will last easily 1/2 h

# Library of working STAN models checked with simulated data
#check with rstanarm and brms
#load packages
library(rstan)
library(arm)
library(shinystan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Documents/PostDoc_Ghent/STAN_stuff/")

set.seed(20170807)
#first case simulate normal regression with 2 predictor and their interaction
N <- 100
K <- 4
dat <- data.frame(x1=runif(N,-2,2),x2=runif(N,-2,2))
X <- model.matrix(~x1*x2,dat)
beta <- rcauchy(K,0,1)
sigma <- 1
dat$y <- rnorm(N,X%*%beta,sigma)
pairs(dat)

m_lin <- stan("Models/normal_model_basic.stan",data = list(N=N,K=K,X=X,y=dat$y),chains = 1)
m_lins <- stan_lm(y~x1*x2,dat,prior = R2(0.5))
#do some standard stuff with this (see BayesianPlot)


#second case poisson regression
dat$y_p <- rpois(N,exp(X%*%beta))

m_pois <- stan("Models/poisson_model_basic.stan",data=list(N=N,K=K,X=X,y=dat$y_p),chains=1)

#third case binomial regression
W <- rpois(N,10)
y_bin <- rbinom(N,W,prob=invlogit(X%*%beta))

m_bin <- stan("Models/binomial_model_basic.stan",data=list(N=N,K=K,W=W,X=X,y=y_bin),chains=1)

#fourth case beta regression
beta <- c(1,1,-0.5,0.5)
phi <- 20
mu <- invlogit(X%*%beta)
y_bet <- rbeta(N,mu*phi,(1-mu)*phi)
dat$y_bet <- y_bet
pairs(dat[,c("x1","x2","y_bet")])

m_bet <- stan("Models/beta_model_basic.stan",data=list(N=N,K=K,X=X,y=y_bet),chains=1)

#fifth case intercept-only hierarchical model
N_group <- 10
dat$group <- gl(n = N_group,k = N/N_group)
rnd_int <- rnorm(N_group,beta[1],2)
y_hier <- rnorm(N,rnd_int[dat$group] + X[,-1] %*% beta[-1],1)
dat$y_hier <- y_hier

m_hier <- stan("Models/normal_model_interceptvarying.stan",
               data=list(N=N,K=K-1,X=X[,-1],N_group=N_group,ID_group=as.numeric(dat$group),y=y_hier),
               chains=1)

#sixth case a poisson hierarchical regression with more than one varying parameter
#including covariance between group-varying parameters

#Bram's case with explicit modelling of sd
dat <- data.frame(x1=runif(100,-2,2),ff=gl(n = 4,k = 25))
X <- model.matrix(~x1,dat)
W <- model.matrix(~ -1+ff,dat)
dat$y <- rnorm(100,X%*%c(1,2),W%*%c(1,0.5,2,1))

m_var <- stan("Models/normal_model_variationstr.stan",data=list(N=100,K=2,L=4,X=X,W=W,y=dat$y),chains=1)

#now putting plot as varyinh terms
dat <- data.frame(x1=runif(400,-2,2),ff=gl(n = 4,k = 100),plot=gl(n=100,k=4))
X <- model.matrix(~x1,dat)
W <- model.matrix(~ -1+plot,dat)

m_varp <- stan("Models/normal_model_variationstr.stan",data=list(N=100,K=2,L=4,X=X,W=W,y=dat$y),chains=1)


#now Bram's case, identical to MCMCglmm with a beta distribution
dat <- data.frame(x1=runif(400,-2,2),ff=gl(n = 4,k = 100),plot=gl(n=100,k=4))
X <- model.matrix(~ -1+x1,dat)
W <- model.matrix(~ -1+ff,dat)
rnd_int <- rnorm(100,1,1)
mu <- invlogit(rnd_int[dat$plot] + 1*dat$x1)
phi <- W%*%c(1,0.5,2,1)
shape1 <- mu*phi
shape2 <- (1-mu)*phi

dat$y <- rbeta(400,shape1,shape2)
dat$y[dat$y<1e-3] <- 1e-3 #make sure that the response stay within the bounds
dat$y[dat$y==1] <- 0.9999

m_varh <- stan("Models/beta_model_variation_interceptvarying.stan",
               data=list(N=400,K=1,L=4,X=X,W=W,y=dat$y,Nb_plots=100,ID_plots=as.numeric(dat$plot)),chains=1)

m_bet <- stan("Models/beta_model_basic.stan",data=list(N=400,K=1,X=X,y=dat$y),chains=1)
