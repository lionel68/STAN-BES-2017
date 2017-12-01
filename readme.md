# Stan Workshop for BES conference in December 2017

## Required packages

To make the most of the coding session please come to the workshop with the following packages installed:

* rstan (v2.16 or higher)
* bayesplot (v1.4 or higher)
* brms (v1.10 or higher)
* ggplot2
* reshape2
* plyr
* MASS

In addition **windows** users will most certainly need to install the RTools, please follow the instruction given on [this page](https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows).

To ensure a smooth start in the coding session check if you can run Stan models from R by running:

```
library(rstanarm)
stan_glm(mpg ~ cyl,family="gaussian",data=mtcars)
```

## Structure of the workshop

* 45min general introduction into Bayesian data analysis, 15min Lionel, 15min Maxime
* 45min interactive fitting of models using brms (or rstanarm?), using shinystan for checking model fitness, learn about convergence issues, fit models of increasing complexity:
..1. Linear model 
..2. Generalized Linear Model (Negative Binomial)
..(3. Mixed-effect models, varying intercept, varying intercept and slopes, crossed and nested random effects)
..(4. Zero-inflated overdispersed poisson model)
