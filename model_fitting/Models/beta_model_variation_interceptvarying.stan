/*
Beta model with regression on both the mean and the variance part
allowing explicit modelling of variation 
in a GLS type of way
also with a intercept varying per plots
*/

data{
	int<lower=1> N;//number of observations
	int<lower=1> K;//number of predictor in the mean part WITHOUT the intercept
	int<lower=1> L;//number of predictor in the variance part
	int<lower=1> Nb_plots;//number of plots
	int<lower=1,upper=Nb_plots> ID_plots[N];//plot ID running form 1 to 53
	matrix[N,K] X;//model matrix of mean part WITHOUT the intercept
	matrix[N,L] W;//model matrix of variance part, best is to have no intercept
	real<lower=0,upper=1> y[N];//the beta response
}
parameters{
	vector[K] beta;//slopes of mean part
	vector<lower=0>[L] gamma;//slopes of variance part, for now only allow positive effects
	real<lower=0> plot_sigma;//the standard deviation of the varying plot effect
	real Intercept;//the population-level intercept (across plots)
	vector[Nb_plots] error;//non-centered parametrization of the random term
}
transformed parameters{
	vector<lower=0,upper=1>[N] mu;//expected mean
	vector[Nb_plots] Intercept_group;//the plot-level intercept
	vector<lower=0>[N] shape1;//first shape parameter for the beta distrib
	vector<lower=0>[N] shape2;//second shape parameter for the beta distrib
	vector<lower=0>[N] phi;//expected deviation



	Intercept_group = Intercept + plot_sigma * error;//non-centered parametrization

	phi = W * gamma;//the variance part regression

	for(n in 1:N){
		mu[n] = inv_logit(Intercept_group[ID_plots[n]] + X[n] * beta);//the mean part regression
		shape1[n] = mu[n] * phi[n];
		shape2[n] = (1 - mu[n]) * phi[n];
	}	
}
model{
	beta ~ normal(0,5);//pior on the slopes of the mean part
	gamma ~ normal(0,10);//prior on the slopes of the variance part
	plot_sigma ~ normal(0,10);//prior on the variation in plot-level intercept
	Intercept ~ normal(0,10);//prior on the intercept
	error ~ normal(0,1);//DO NOT TOUCH this prior

	y ~ beta(shape1,shape2);
	//y ~ normal(mu,phi);
}
generated quantities{
	//if needed this would generate predicted values from the model
	//vector[N] y_gen;
	//for(n in 1:N)
	//	y_gen[n] = beta_rng(shape1[n],shape2[n]);
}




//mu[n] = X[n] * beta[plot_id[n]]
//sigma[n] = SpComb[n] * gamma[plot_id[n]]


//y ~ normal(mu,sigma)

