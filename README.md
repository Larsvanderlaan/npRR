# npRR

## Nonparametric inference for the conditional relative risk function using targeted machine learning
This package contains the `R` function `npRRWorkingModel()`, which implements a targeted maximum likelihood estimator (TMLE) for the coefficients of the projection of the conditional log relative risk function onto a user-specified linear working model. This TMLE can be viewed as a model-robust version of the semiparametric TMLE studied in Tuglus, Porter, van der Laan (2011).


The function requires initial estimates of the outcome regression and propensity score which can be passed directly to ``npRRWorkingModel()` 
or can be estimated internally using the `sl3` ensemble machine learning `R` package.

### Data-structure
This package is based on the data structure `(W, A, Y)` where `W` is a covariate, `A` is a binary treatment assignment, and `Y` is a binary or nonnegative outcome.
The estimand of interest is the conditional relative risk function: 
`RR(w) := E[Y | A = 1, W = w] / E[Y | A = 0, W = w]`.
This package supports observation weights, which can be used to adjust for outcome missingness/censoring.

The function `npRRWorkingModel()` provides estimate and inference for a working Poisson likelihood-based projection of the conditional relative risk
onto a user-speciifed log-linear parametric model. Coefficient estimates (`b`) are returned based on the working model `RR(w) ~ exp{b^T f(w)}` 
where `f` is a known, user-specified transformation of the covariate vector `W`.


### Inputs: 
The `npRRWorkingModel()` function requires the following arguments:
1. `formula_LRR`: A `formula` object specifying the linear working model for the log relative risk function.
2. `W`: numeric matrix containing covariate information (e.g. possible confounders).
3. `A`: binary vector of treatment assignments. 
4. `Y`: A vector of binary or nonnegative outcome values.
5. `weights`: observation weights. No weights corresponds with the argument `weights = rep(1, length(Y))`. 
6. `EY1W` : vector of estimates for `E[Y|A=1,W]` (see also sl3_Learner_EYAW argument).
7. `EY0W` : vector of estimates for `E[Y|A=0,W]` (see also sl3_Learner_EYAW argument).
8. `pA1W` : vector of estimates for `P(A=1|W)` (see also sl3_Learner_pA1W argument).


### Output of `npRRWorkingModel()`

The output of `npRRWorkingModel()` is a list object that includes:
1. Coefficient estimates  
2. Z-scores and p-values for coefficients 
3. 95% confidence intervals for coefficients

### Example code

``` r
library(npRR)
library(sl3)
# relative risk
n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rpois(n, lambda = exp( A * (1 + W + 2*W^2)  + sin(5 * W)))
formula = ~ 1 + W
# Estimate nuisance function internally using sl3 package and generalized additive models
fit <- npRRWorkingModel(formula_LRR = formula,
                        W = W, A = A, Y = Y,
                        weights = rep(1,n), 
                        sl3_Learner_EYAW = Lrnr_gam$new(),
                        sl3_Learner_pA1W = Lrnr_gam$new()
                        )
coef(fit)
```


### References

1. Tuglus, Cathy, Kristin E. Porter, and Mark J. van der Laan. "Targeted Maximum Likelihood Estimation of Conditional Relative Risk in a Semi-parametric Regression Model." (2011).
2. Van der Laan, Mark J., and Sherri Rose. Targeted learning: causal inference for observational and experimental data. Vol. 4. New York: Springer, 2011.
3. Van Der Laan, Mark J., and Daniel Rubin. "Targeted maximum likelihood learning." The international journal of biostatistics 2.1 (2006).
4. Coyle, Jeremy R and Hejazi, Nima S and Malenica, Ivana and
        Phillips, Rachael V and Sofrygin, Oleg , sl3: Modern Pipelines for Machine Learning and {Super
        Learning, 2021, R package version 1.4.2


