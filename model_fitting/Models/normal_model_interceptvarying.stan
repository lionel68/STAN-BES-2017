/*
Standard intercept only hierarchical model withnormal regression for
any number of predictor variables
with weakly informative priors on the betas and on the standard deviation
using a non-centered re-parametrization
*/
data{
	int<lower=1> N; //number of observations
	int<lower=1> K; //number of predictor variables NOT including the intercept
	matrix[N,K] X; //the model matrix NOT including intercept
	int<lower=1> N_group;//number of group for the varying intercept
	int<lower=1,upper=N_group> ID_group[N];//index variable for group ID
	vector[N] y; //the response variable
}

parameters{
	real Intercept;//the population-level intercept parameter
	vector[K] beta; //the population-level regression parameters

	real<lower=0> sigma;//the population variation
	real<lower=0> group_sigma;//the deviation in intercept

	vector[N_group] error; //for the non-centered parametrization	
}
transformed parameters{
	
	vector[N_group] Intercept_group;//the group-level intercept
	vector[N] mu; //the linear predictor

	Intercept_group = Intercept + group_sigma * error;//non-centered parametrization

	for(n in 1:N)
		mu[n] = Intercept_group[ID_group[n]] + X[n] * beta;

}

model{
	Intercept ~ normal(0,10);
	error ~ normal(0,1);
	beta ~ normal(0,5);
	sigma ~ normal(0,10);
	group_sigma ~ cauchy(0,5);

	y ~ normal(mu,sigma); //the likelihood
}

generated quantities{
	vector[N] y_gen;
	for(n in 1:N)
		y_gen[n] = normal_rng(mu[n],sigma);
}
