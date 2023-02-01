# npRR

## Nonparametric inference for the conditional relative risk function using targeted machine learning
This package contains the `R` function `npRRWorkingModel()`, which implements a targeted maximum likelihood estimator for the coefficients of the projection of the conditional log relative risk function onto a user-specified linear working model.
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
1. `W`: numeric matrix containing covariate information (e.g. possible confounders)
2. `A`: binary vector of treatment assignments 
3. `Y`: A vector of binary or nonnegative outcome values
4. `weights`: (optional) observation weights
5. `EY1W` : vector of estimates for `E[Y|A=1,W]` (see also sl3_Learner_EYAW argument)
6. `EY0W` : vector of estimates for `E[Y|A=0,W]` (see also sl3_Learner_EYAW argument)
7. `pA1W` : vector of estimates for `P(A=1|W)` (see also sl3_Learner_pA1W argument)


### Output of `npRRWorkingModel()`

The output of `npRRWorkingModel()` is a list object that includes:
1. Coefficient estimates  
2. Z-scores and p-values for coefficients 
3. 95% confidence intervals for coefficients

