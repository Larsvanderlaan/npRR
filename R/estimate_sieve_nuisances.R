



#' Computes the data-adaptive sieve-based nuisance estimators of E[Y|A,W] and P(A=1|W)
#' that are constructed to be used with the IPW and plugin empirical risk functions.
#' @param W A matrix of covariate observations
#' @param A A binary vector specifying the treatment assignment. The values should be in {0,1}.
#' @param Y A numeric vector of binary or nonnegative observations of the outcome variable.
#' @param EY1 A numeric vector containing initial cross-fitted estimates of E[Y|A=1,W] for all observations.
#' @param EY0 A numeric vector containing initial cross-fitted estimates of E[Y|A=0,W] for all observations.
#' @param pA1 A numeric vector containing initial cross-fitted estimates of P(A=1|W) for all observations.
#' @param weights A numeric vector of observation weights. If no special weighting desired, supply a vector of 1's.
#' @param basis_generator A basis_generator object (see \code{fourier_basis} and \code{bspline_basis}) ...
#' that specifies a sieve basis transformation of the design matrix \code{W}.
#' @param family An R \code{family} object that specifies the link and risk function used in the data-adaptive sieve update step.
#' `family` should be a \code{binomial} object if `Y` is binary and a \code{poisson} object if `Y` is non-binary and nonnegative (e.g. a count).
#' @param debug ...
compute_plugin_and_IPW_sieve_nuisances <- function(W, A, Y, EY1, EY0, pA1, weights, basis_generator, family = binomial(), debug = FALSE) {
  if(is.null(basis_generator)) {
    return(list(pA1_star = pA1, EY1_star = EY1, EY0_star = EY0, sieve = "no_sieve"))
  }
  # Compute sieve-transformed design matrix

  V <- basis_generator$set(W)$eval(W)

  # Compute data-adaptive sieve
  V_plugin <- cbind(A*V, (1-A)*V)

  EY <- ifelse(A==1, EY1, EY0)
  pA0 <- 1 - pA1
  pA <- ifelse(A==1, pA1, pA0)
  sieve_fit_plugin <- glm.fit(V_plugin, Y, weights = weights/pA, offset = family$linkfun(EY), family = family, intercept = F)
  beta_plugin <- coef(sieve_fit_plugin)
  beta_plugin[is.na(beta_plugin)] <- 0
  beta1_plugin <- beta_plugin[1:ncol(V)]
  beta0_plugin <- beta_plugin[-(1:ncol(V))]


  EY1_star <- as.vector(family$linkinv(family$linkfun(EY1) + V %*% beta1_plugin))
  EY0_star <- as.vector(family$linkinv(family$linkfun(EY0) + V %*% beta0_plugin))

  if(debug) {
    EYstar <- ifelse(A==1, EY1_star, EY0_star )
    print("Sieve scores plugin")
    print(colMeans(weights/pA*V_plugin*(Y - EYstar)))
  }

  V_IPW <- cbind(EY1/pA1 * V, EY0/pA0 * V)
  sieve_fit_IPW <- glm.fit(V_IPW, A, weights = weights, offset = qlogis(pA1), family = binomial(), intercept = F)
  beta_IPW <- coef(sieve_fit_IPW)
  beta_IPW[is.na(beta_IPW)] <- 0
  pA1_star <- as.vector(plogis(qlogis(pA1) + V_IPW %*% beta_IPW))
  if(debug) {
    print("Sieve scores IPW")
    print(colMeans(weights*V_IPW*(A - pA1_star)))
  }

  output <- list(pA1_star = pA1_star, EY1_star = EY1_star, EY0_star = EY0_star, sieve = basis_generator$name)
}








