/*
Standard normal regression for any number of predictor variables
with weakly informative priors on the betas and on the standard deviation
*/
data{
	int<lower=1> N; //number of observations
	int<lower=1> K; //number of predictor variables
	matrix[N,K] X; //the model matrix including intercept
	vector[N] y; //the response variable
}

parameters{
	vector[K] beta; //the regression parameters
	real<lower=0> sigma;
}

model{
	vector[N] mu; //the linear predictor
	mu = X * beta; //the regression
	
	beta[1] ~ normal(0,10);
	beta[2:K] ~ normal(0,5);
	sigma ~ normal(0,10);

	y ~ normal(mu,sigma); //the likelihood

}

generated quantities{
	vector[N] y_gen;
	for(n in 1:N)
		y_gen[n] = normal_rng(X[n] * beta,sigma);
}

