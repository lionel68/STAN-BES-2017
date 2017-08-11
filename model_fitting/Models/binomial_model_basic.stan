/*
Standard binomial regression for any number of predictor variables
with weakly informative priors on the betas and on the standard deviation
*/
data{
	int<lower=1> N; //number of observations
	int<lower=1> K; //number of predictor variables
	matrix[N,K] X; //the model matrix including intercept
	int<lower=1> W[N]; //number of trials
	int<lower=0> y[N]; //the response variable
}

parameters{
	vector[K] beta; //the regression parameters
}

transformed parameters{
	vector[N] theta; //the linear predictor
	theta = X * beta; //the regression
}

model{
	beta[1] ~ normal(0,10);
	beta[2:K] ~ normal(0,5);

	y ~ binomial_logit(W,theta); //the likelihood
}

generated quantities{
	int<lower=0> y_gen[N];
	vector<lower=0,upper=1>[N] theta_gen;
	for(n in 1:N){
		theta_gen[n] = inv_logit(X[n] * beta);
		y_gen[n] = binomial_rng(W[n],theta_gen[n]);
	}
}


