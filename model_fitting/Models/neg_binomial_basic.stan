/*
Negative binomial regression for any number of predictor variables
with weakly informative priors on the betas and on the overdispersion parameter
*/
data{
	int<lower=1> N; //number of observations
	int<lower=1> K; //number of predictor variables
	matrix[N,K] X; //the model matrix including intercept
	int<lower=0> y[N]; //the response variable
}

parameters{
	vector[K] beta; //the regression parameters
	real<lower=0> phi; //the overdispersion parameter
}

model{
	vector[N] mu; //the linear predictor
	mu = X * beta; //the regression
	
	beta[1] ~ normal(0,10);
	beta[2:K] ~ normal(0,5);
	phi ~ normal(0,50);

	y ~ neg_binomial_2_log(mu,phi); //the likelihood

}

generated quantities{
	int<lower=0> y_gen[N];
	for(n in 1:N)
		y_gen[n] = neg_binomial_2_log_rng(X[n] * beta,phi);
}

