% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_LRR_using_sieve.R
\name{train_learner_all_sieves}
\alias{train_learner_all_sieves}
\title{For a given base LRR learner, train the LRR learners on all sieve-based plugin and IPW empirical risk functions}
\usage{
train_learner_all_sieves(
  sl3_Learner,
  V,
  A,
  Y,
  weights,
  family_risk_function,
  outcome_function_plugin,
  weight_function_plugin,
  outcome_function_IPW,
  weight_function_IPW,
  transform_function,
  list_of_sieve_nuisances,
  Vpred = V,
  compute = FALSE
)
}
\arguments{
\item{V}{A matrix of observations of a subset of the covariates `W` for which to estimate the (possibly semi-marginalized) log relative risk (LRR).}

\item{A}{A binary vector specifying the treatment assignment. The values should be in {0,1}.}

\item{Y}{A numeric vector of binary or nonnegative observations of the outcome variable.}

\item{weights}{A numeric vector of observation weights. If no special weighting desired, supply a vector of 1's.}

\item{list_of_sieve_nuisances}{A list of sieve nuisance estimates as returned by the function \link{compute_plugin_and_IPW_sieve_nuisances}.}

\item{Vpred}{A matrix of covariates observations at which to predict the LRR. By default, \code{Vpred} equals \code{W}.}

\item{compute}{Whether to `compute` the list of delayed trained learners. See the package \code{delayed} for details.}

\item{sl3_learner_binomial}{See argument \code{sl3_learner_binomial} of \link{estimate_LRR_using_ERM} for details.}

\item{EY1W}{A numeric vector containing initial cross-fitted estimates of E[Y|A=1,W] for all observations.}

\item{EY0W}{A numeric vector containing initial cross-fitted estimates of E[Y|A=0,W] for all observations.}

\item{pA1W}{A numeric vector containing initial cross-fitted estimates of P(A=1|W) for all observations.}
}
\description{
For a given base LRR learner, train the LRR learners on all sieve-based plugin and IPW empirical risk functions
}
