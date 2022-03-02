compute_plugin_and_IPW_sieve_nuisances <- function(W, A, Y, EY1, EY0, pA1, weights, basis_generator, family = binomial(), debug = TRUE) {
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








