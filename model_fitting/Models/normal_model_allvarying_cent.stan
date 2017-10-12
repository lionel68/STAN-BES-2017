/*
Hierarchical linear model with normal regression for
any number of predictor variables, all varying according to one grouping level
with weakly informative priors on the betas and on the standard deviation
using a 
*/
data{
	int<lower=1> N; //number of observations
	int<lower=1> K; //number of predictor variables NOT including the intercept
	matrix[N,K] X; //the model matrix including the intercept
	int<lower=1> N_group;//number of group for the varying intercept
	int<lower=1,upper=N_group> ID_group[N];//index variable for group ID
	vector[N] y; //the response variable
}

parameters{
	vector[K] mu; //the population-level regression parameters
	corr_matrix[K] Omega; //correlation between varying effects
	vector<lower=0>[K] tau; //variation of the varying effects
	
	vector[K] beta[N_group];//matrix of group-level effects
	
	real<lower=0> sigma_residuals;//the variation of the residuals

	//vector[N_group] error; //for the non-centered parametrization	
}
transformed parameters{
	matrix[K,K] Sigma; //the variance-covariance matrix between the varying effects

	Sigma = quad_form_diag(Omega,tau);
}

model{
	tau ~ cauchy(0,2.5);
	Omega ~ lkj_corr(2);
	mu ~ normal(0,5);

	beta ~ multi_normal(mu, Sigma);

	{
		vector[N] lin_pred; //the linear predictor
		for(n in 1:N)
			lin_pred[n] = X[n] * beta[ID_group[n]];
		y ~ normal(lin_pred,sigma_residuals);

	}
}

generated quantities{
	vector[N] y_gen;
	vector[K] sd_eff; //the standard deviation of the varying effects

	for(n in 1:N)
		y_gen[n] = normal_rng(X[n] * beta[ID_group[n]],sigma_residuals);

	
	sd_eff = sqrt(diagonal(Sigma));

}
