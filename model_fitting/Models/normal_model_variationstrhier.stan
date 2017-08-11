/*
Normal model with regression on both the mean and the variance part
allowing explicit modelling of variation in the standard deviation
in a GLS type of way
also with hierarchical effect on the variance part (all parameters by default)
*/

data{
	int<lower=1> N;
	int<lower=1> K;//number of predictor in the mean part
	int<lower=1> L;//number of predictor variance
	int<lower=1> n_plot;//number of plots
	int<lower=1,upper=n_plot> plot_id[N];//index of plot identity
	matrix[N,K] X;//model matrix of mean part
	matrix[N,L] W;//model matrix of variance part
	real y[N];
}
parameters{
	vector[K] beta;//slopes of mean part
	vector[L] gamma;//slopes of species combination-level variance part, for now only allow positive effects

	vector<lower=0>[L] sd_plot;//variation of plot within species combination sd

	//matrix[n_plot,L] error;//for the non-centered parametrization of the plot-level variation
}
transformed parameters{
	vector[N] mu;//expected mean
	vector<lower=0>[N] sigma;//expected deviation
	vector<lower=0>[n_plot] sigma_plotlvl;//plot-level deviation from expected values

	//sigma_plotlvl = gamma + (sd_plot * error)';//compute the plot-level variation

	mu = X * beta;
	for(n in 1:N)
		sigma[n] = sum(W[n] * sigma_plotlvl[plot_id[n]]);
}
model{
	beta ~ normal(0,5);
	gamma ~ normal(0,10);
	//to_vector(error) ~ normal(0,1);
	sd_plot ~ normal(0,1);
	sigma_plotlvl ~ normal(gamma,sd_plot);

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
