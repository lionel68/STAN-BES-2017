/*
Normal model with regression on both the mean and the variance part
allowing explicit modelling of variation in the standard deviation
in a GLS type of way
*/

data{
	int<lower=1> N;
	int<lower=1> K;//number of predictor in the mean part
	int<lower=1> L;//number of predictor variance
	matrix[N,K] X;//model matrix of mean part
	matrix[N,L] W;//model matrix of variance part
	real y[N];
}
parameters{
	vector[K] beta;//slopes of mean part
	vector<lower=0>[L] gamma;//slopes of variance part, for now only allow positive effects
}
transformed parameters{
	vector[N] mu;//expected mean
	vector<lower=0>[N] sigma;//expected deviation

	mu = X * beta;
	sigma = W * gamma;
}
model{
	beta ~ normal(0,5);
	gamma ~ normal(0,10);

	y ~ normal(mu,sigma);
}
generated quantities{
	vector[N] y_gen;
	for(n in 1:N)
		y_gen[n] = normal_rng(mu[n],sigma[n]);
}




//mu[n] = X[n] * beta[plot_id[n]]
//sigma[n] = SpComb[n] * gamma[plot_id[n]]


//y ~ normal(mu,sigma)

