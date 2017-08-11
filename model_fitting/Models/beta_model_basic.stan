/*
Standard beta regression for any number of predictor variables
with weakly informative priors on the betas and on the standard deviation
*/
data{
	int<lower=1> N; //number of observations
	int<lower=1> K; //number of predictor variables
	matrix[N,K] X; //the model matrix including intercept
	real<lower=0,upper=1> y[N]; //the response variable
}

parameters{
	vector[K] beta; //the regression parameters
	real<lower=0> phi;//the dispersion parameter
}

transformed parameters{
	vector<lower=0,upper=1>[N] mu; //the linear predictor
	vector<lower=0>[N] shape1;//the first shape parameter for the beta distrib
	vector<lower=0>[N] shape2;//the second shape parameter for the beta distrib
	for(n in 1:N){
		mu[n] = inv_logit(X[n] * beta);//the regression
		shape1[n] = mu[n] * phi;//re-parametrize in terms of alpha and beta
		shape2[n] = (1.0 - mu[n]) * phi;
	}
}

model{
	beta[1] ~ normal(0,10);
	beta[2:K] ~ normal(0,5);
	phi ~ gamma(0.1,0.1);

	y ~ beta(shape1,shape2); //the likelihood
}

generated quantities{
	real<lower=0,upper=1> y_gen[N];
	for(n in 1:N){
		y_gen[n] = beta_rng(shape1[n],shape2[n]);
	}
}
