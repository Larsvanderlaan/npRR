



#' Computes the data-adaptive sieve-based nuisance estimators of E[Y|A,W] and P(A=1|W)
#' that are constructed to be used with the IPW and plugin empirical risk functions.
#' @param V A matrix of observations of a subset of the covariates `W` for which to estimate the (possibly semi-marginalized) log relative risk (LRR).
#' The sieve basis will only be generated for these covariates.
#' @param A A binary vector specifying the treatment assignment. The values should be in {0,1}.
#' @param Y A numeric vector of binary or nonnegative observations of the outcome variable.
#' @param EY1W A numeric vector containing initial cross-fitted estimates of E[Y|A=1,W] for all observations.
#' @param EY0W A numeric vector containing initial cross-fitted estimates of E[Y|A=0,W] for all observations.
#' @param pA1W A numeric vector containing initial cross-fitted estimates of P(A=1|W) for all observations.
#' @param weights A numeric vector of observation weights. If no special weighting desired, supply a vector of 1's.
#' @param basis_generator A basis_generator object (see \code{fourier_basis} and \code{bspline_basis}) ...
#' that specifies a sieve basis transformation of the design matrix \code{W}.
#' @param family An R \code{family} object that specifies the link and risk function used in the data-adaptive sieve update step.
#' `family` should be a \code{binomial} object if `Y` is binary and a \code{poisson} object if `Y` is non-binary and nonnegative (e.g. a count).
#' @param debug ...
compute_plugin_and_IPW_sieve_nuisances <- function(V, A, Y, EY1W, EY0W, pA1W, weights, basis_generator, family = binomial(), debug = FALSE) {
  if(is.null(basis_generator)) {
    return(list(pA1W_star = pA1W, EY1W_star = EY1W, EY0W_star = EY0W, sieve = "no_sieve"))
  }
  # Compute sieve-transformed design matrix
  basis_generator <- basis_generator$clone()
  X <- basis_generator$set(V)$eval(V)

  # Compute data-adaptive sieve
  X_plugin <- cbind(A*X, (1-A)*X)

  EY <- ifelse(A==1, EY1W, EY0W)
  pA0 <- 1 - pA1W
  pA <- ifelse(A==1, pA1W, pA0)
  suppressWarnings(sieve_fit_plugin <- glm.fit(X_plugin, Y, weights = weights/pA, offset = family$linkfun(EY), family = family, intercept = F))
  beta_plugin <- coef(sieve_fit_plugin)
  beta_plugin[is.na(beta_plugin)] <- 0
  beta1_plugin <- beta_plugin[1:ncol(X)]
  beta0_plugin <- beta_plugin[-(1:ncol(X))]


  EY1W_star <- as.vector(family$linkinv(family$linkfun(EY1W) + X %*% beta1_plugin))
  EY0W_star <- as.vector(family$linkinv(family$linkfun(EY0W) + X %*% beta0_plugin))

  if(debug) {
    EYstar <- ifelse(A==1, EY1W_star, EY0W_star )
    print("Sieve scores plugin")
    print(colMeans(weights/pA*X_plugin*(Y - EYstar)))
  }

  X_IPW <- cbind(EY1W/pA1W * X, EY0W/pA0 * X)
  suppressWarnings(sieve_fit_IPW <- glm.fit(X_IPW, A, weights = weights, offset = qlogis(pA1W), family = binomial(), intercept = F))
  beta_IPW <- coef(sieve_fit_IPW)
  beta_IPW[is.na(beta_IPW)] <- 0
  pA1W_star <- as.vector(plogis(qlogis(pA1W) + X_IPW %*% beta_IPW))
  if(debug) {
    print("Sieve scores IPW")
    print(colMeans(weights*X_IPW*(A - pA1W_star)))
  }

  output <- list(pA1W_star = pA1W_star, EY1W_star = EY1W_star, EY0W_star = EY0W_star, sieve = basis_generator$name)
}








