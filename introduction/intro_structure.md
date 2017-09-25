# Outline of the introduction presentation

## Aim:

Give a 30min intro into Bayesian data analysis using ecological examples when possible


## Structure:

### Part 1: What is Bayesian data analysis (~10min)

Quotes from Hobbs and Hilborn: "There is a danger that questions are chosen for investigation by ecologists to fit widely sanctioned statistical methods rather than statistical methods being chosen to meet the needs of ecological questions"
What is inference? How do we gather knowledge? What is statistical inference? What information do we need from a model?

-> About inference in science: we see the data and we would like to know something about the process that generated them.
-> Posterior, likelihood and priors
-> A walkthrough example of a simple ABC computation following RB
-> Inference using Bayesian data analysis, called inverse modelling until Fischer came in in the XXth century
-> the issue of sampling the posterior
-> Example of BDA in ecology
Using and interpreting the posterior in Bayesian inference
Some definition of classical stats, likelihood, MLE, p-values ...

How does Bayes and Frequentist differ (probability vs asymptotic approximation)
How does Bayes analysis concretely work (MCMC sampling, issue of computing the integral, chains ...)
Does the Bayes/Frequentist differences matter? When? Is one better than the other?

In this part I'd like to outline:


* Why do we do stats in general, about inference, big world-small world dichotomy from StatRethinking
* The bayesian formula: Posterior info = Likelihood * prior info 
* What is the likelihood, dnorm, show some brute-force, 2D examples (mean+sd)
* What is prior infos, what kind of info do we generally have, how do we turn it into priors (in a way we always have __some__ prior info; uniform vs. weak priors)
* What is posterior, combination of likelihood and prior infos, relative weight of each = f(sample size, prior strength), P(hypothesis | data)
* It's all about uncertainty


### Part 2: Why do Bayesian data analysis (~10min)

-> Bayes (probability) vs frequentise (repeat sampling), but also Bayes and frequentist (asymptotic convergence of the posterior to the likelihood under infinite sample size)
-> Questions Bayes can answer that Freq cannot
-> Flexible model building, uncertainty in all model parameters
-> Include more information about the data

* Embracing uncertainty
* Flexibility in model building
* Getting what we are usually intersted by rather than some weird reference to null hypothesis (We tend to interpret frequentist results p(data|hypothesis) with Bayesian thinking p(hypothesis|data) already anyway, so might as well use method that actually provides what we look for)
* Fit complex models even with little data, varying weight of likelihood and prior as sample size increases
* Give one or two published examples were supplementary insights from Bayesian approach are obvious (ideally scientific benefits/insights rather than simply "technical" benefits)

### Part 3: How do we do Bayesian data analysis (~10min)

What are the options to fit Bayesian models (BUGS, JAGS, STAN), what are the differences?

-> The different options in R: BUGS, JAGS (r2jags,other packages allowing easy model fit with pre-set priors), STAN (rstan, rstanarm, brms)
-> Specificitis of model building and checking in BDA (chains, convergence, divergence)
-> The rest of the workshop

* 2 major ways: (i) coding the model in a bayesian program (WinBUGS, JAGS, STAN) and use their R interface to fit the models (R2WinBUGS, R2JAGS, rstan), or (ii) use packages that allow you to fit the models based on R formula syntax and these packages will do the translation to the bayesian program (BayesianTools, MCMCglmm, rstanarm, brms).
* About sampling (the blind man going through the likelihood landscape), what is MCMC, the Metropolis-Hasting vs Hamiltonian Monte Carlo (adv. and disadv.) (metaphors from the StatRethinking book really useful here)
* Vocabulary to know in BDA: chains, convergence, divergence, acceptance
* Model checking in BDA ...















