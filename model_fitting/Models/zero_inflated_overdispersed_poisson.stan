/*A zero-inflated overdispersed poisson model
*with a no-centered parametrization on the overdispersion part
*and assuming fixed amount of zero-inflation
*/
data{
	int<lower=0> N; //number of observation
	int<lower=0> K; //number of column in the model matrix
	matrix[N,K] X; //the model matrix
	int<lower=0> y[N]; //the counts
}

parameters{
	vector[K] beta; //the regression parameters
	real<lower=0,upper=1> p; //the probability of zeros
	real<lower=0> sigma; //the amount of overdispersion

	vector[N] z; //random deviates for the non-centered parametrization
}

transformed parameters{
	vector[N] rnd; //overdispersion values
	vector[N] lambda; //the linear predictor

	rnd = sigma * z;
	
	lambda = X * beta + rnd;
}

model{
	sigma ~ cauchy(0,2.5);
	beta ~ normal(0,5);
	p ~ beta(1, 1);
	
		
	z ~ normal(0,1);	
	
	for(n in 1:N){
		if(y[n] == 0) //if y is zero then it could either come from the zero-inflation or from the poisson part of the model
			target += log_sum_exp(bernoulli_lpmf(1 | p),
					      bernoulli_lpmf(0 | p) +
					      poisson_log_lpmf(0 | lambda[n]));
		else //if y is non-zero then it comes from the poisson part 
			target += bernoulli_lpmf(0 | p) +
					      poisson_log_lpmf(y[n] | lambda[n]);
	}
}


	
