# Outline of the introduction presentation

## Aim:

Give a 30min intro into Bayesian data analysis using ecological examples when possible


## Structure:

### Part 1: What is Bayesian data analysis (~10min)

Quote from Hobbs and Hilborn: "There is a danger that questions are chosen for investigation by ecologists to fit widely sanctioned statistical methods rather than statistical methods being chosen to meet the needs of ecological questions"

On inference:
* What is inference? How do we gather knowledge? About inference in science: we see the data and we would like to know something about the process that generated them.
* What is statistical inference? What information do we need from a model?
* Why do we do stats in general, about inference, big world-small world dichotomy from StatRethinking

(Not sure where to place these yet
-> A walkthrough example of a simple ABC computation following RB
-> Inference using Bayesian data analysis, called inverse modelling until Fischer came in in the XXth century
-> the issue of sampling the posterior
Some definition of classical stats, likelihood, MLE, p-values ...)

* The bayesian formula: Posterior info = Likelihood * prior info 
* What is the likelihood, dnorm, show some brute-force, 2D examples (mean+sd)
* What is prior infos, what kind of info do we generally have, how do we turn it into priors (in a way we always have __some__ prior info; uniform vs. weak priors)
* What is posterior (how to use and interpret it), combination of likelihood and prior infos, relative weight of each = f(sample size, prior strength), P(hypothesis | data)
* It's all about uncertainty


### Part 2: Why do Bayesian data analysis (~10min)

Does the Bayes/Frequentist differences matter? When? Is one better than the other?
-> How does Bayes and Frequentist differ (probability vs asymptotic approximation) Bayes (probability) vs frequentist (repeat sampling), but also Bayes and frequentist (asymptotic convergence of the posterior to the likelihood under infinite sample size) Both are expected to lead to similar inference when sample size is very high and prior information negligible (how frequent is it really?)

* Embracing uncertainty in ALL model parameters
* Flexibility in model building
* Getting what we are usually interested by rather than some weird reference to null hypothesis (We tend to interpret frequentist results: p(data|hypothesis) with Bayesian thinking: p(hypothesis|data) already anyway, so might as well use method that actually provides what we look for)
* Fit complex models even with little data, varying weight of likelihood and prior as sample size increases
* Questions Bayes can answer that Freq cannot: Give one or two published examples were supplementary insights from Bayesian approach are obvious (ideally scientific benefits/insights rather than simply "technical" benefits)

### Part 3: How do we do Bayesian data analysis (~10min)

What are the options to fit Bayesian models (commonly used languages: BUGS, JAGS, STAN), what are the differences?

* 2 major ways: (i) coding the model in a bayesian program (WinBUGS, JAGS, STAN) and use their R interface to fit the models (R2WinBUGS, R2JAGS, rstan), or (ii) use packages that allow you to fit the models based on R formula syntax and these packages will do the translation to the bayesian program (BayesianTools, MCMCglmm, rstanarm, brms).
* About sampling (the blind man going through the likelihood landscape), what is MCMC, the Metropolis-Hasting vs Hamiltonian Monte Carlo (adv. and disadv.) (metaphors from the StatRethinking book really useful here). 

How does Bayes analysis concretely work (MCMC sampling, how sampling in general nicely sidesteps issue of computing the integral, chains ...)
* Vocabulary to know in BDA: chains, convergence, divergence, acceptance
* Model checking in BDA ...















