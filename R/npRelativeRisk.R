


#' For a given base LRR learner, train the LRR learners on all sieve-based plugin and IPW empirical risk functions
#' @param sl3_LRR_Learner_binomial See argument \code{sl3_LRR_Learner_binomial} of \link{estimate_LRR_using_ERM} for details.
#' @param W A matrix of covariate observations
#' @param A A binary vector specifying the treatment assignment. The values should be in {0,1}.
#' @param Y A numeric vector of binary or nonnegative observations of the outcome variable.
#' @param EY1 A numeric vector containing initial cross-fitted estimates of E[Y|A=1,W] for all observations.
#' @param EY0 A numeric vector containing initial cross-fitted estimates of E[Y|A=0,W] for all observations.
#' @param pA1 A numeric vector containing initial cross-fitted estimates of P(A=1|W) for all observations.
#' @param weights A numeric vector of observation weights. If no special weighting desired, supply a vector of 1's.
#' @param list_of_sieve_nuisances A list of sieve nuisance estimates as returned by the function \link{compute_plugin_and_IPW_sieve_nuisances}.
#' @param Wpred A matrix of covariates observations at which to predict the LRR. By default, \code{Wpred} equals \code{W}.
#' @param compute Whether to `compute` the list of delayed trained learners. See the package \code{delayed} for details.
train_LRR_learner_all_sieves <- function(sl3_LRR_Learner_binomial, W, A, Y, weights, list_of_sieve_nuisances, Wpred = W, compute = FALSE) {
  list_of_sieve_LRR <- lapply(list_of_sieve_nuisances, function(sieve_nuisances) {
    EY1_star <- sieve_nuisances$EY1_star
    EY0_star <- sieve_nuisances$EY0_star
    pA1_star <- sieve_nuisances$pA1_star

    delayed_plugin_LRR <- delayed_fun(estimate_LRR_using_ERM)(W, A, Y,  EY1_star, EY0_star, pA1_star, weights, sl3_LRR_Learner_binomial, learning_method = "plugin", Wpred = Wpred)

    delayed_IPW_LRR <- delayed_fun(estimate_LRR_using_ERM)(W, A, Y,  EY1_star, EY0_star, pA1_star, weights, sl3_LRR_Learner_binomial, learning_method = "IPW", Wpred = Wpred)

    output <- (list( plugin = delayed_plugin_LRR, IPW = delayed_IPW_LRR))
  })
  sieve_names <- sapply(list_of_sieve_nuisances, `[[`, "sieve")
  names(list_of_sieve_LRR) <- sieve_names
  if(compute) {
    list_of_sieve_LRR <- bundle_delayed(unlist(list_of_sieve_LRR))
    list_of_sieve_LRR <- list_of_sieve_LRR$compute()
  }
  return(list_of_sieve_LRR)
}

#' Train the LRR learners for all sieves and base learners on the full data.
#' @param W A matrix of covariate observations.
#' @param A A binary vector specifying the treatment assignment. The values should be in {0,1}.
#' @param Y A numeric vector of binary or nonnegative observations of the outcome variable.
#' @param EY1 A numeric vector containing initial cross-fitted estimates of E[Y|A=1,W] for all observations.
#' @param EY0 A numeric vector containing initial cross-fitted estimates of E[Y|A=0,W] for all observations.
#' @param pA1 A numeric vector containing initial cross-fitted estimates of P(A=1|W) for all observations.
#' @param weights A numeric vector of observation weights. If no special weighting desired, supply a vector of 1's.
#' @param list_of_sieve_nuisances A list of sieve nuisance estimates as returned by the function \link{compute_plugin_and_IPW_sieve_nuisances}.
#' @param list_of_LRR_learners A list of untrained \code{sl3_Learner} objects to be used to estimate the log relative risk LRR using the function \link{estimate_LRR_using_ERM}.
#' @param list_of_sieves A list of basis_generator objects specifying the sieve. See, for example, \code{fourier_basis} for an example and template.
#' @param Wpred A matrix of covariates observations at which to predict the LRR. By default, \code{Wpred} equals \code{W}.
#' @param compute Whether to `compute` the list of delayed trained learners. See the package \code{delayed} for details.
train_LRR_learners <- function(W, A, Y, EY1, EY0, pA1, weights, list_of_LRR_learners, list_of_sieves, Wpred = W, compute = TRUE) {

  list_of_sieve_nuisances <- lapply(list_of_sieves, function(sieve){
    compute_plugin_and_IPW_sieve_nuisances(basis_generator = sieve, W = W, A = A, Y = Y, EY1 = EY1, EY0 = EY0, pA1 = pA1, weights = weights)})
  list_of_sieve_nuisances

  all_learners_delayed <- lapply(list_of_LRR_learners, train_LRR_learner_all_sieves, list_of_sieve_nuisances = list_of_sieve_nuisances, W = W, A = A, Y = Y, weights = weights, Wpred = Wpred)
  learner_names <- lapply(list_of_LRR_learners, `[[`, "name")
  names(all_learners_delayed) <- paste0(learner_names)
  if(compute) {
    all_learners_delayed <- bundle_delayed(unlist(all_learners_delayed))
    all_learners_delayed <- all_learners_delayed$compute()
  }
  return(all_learners_delayed)
}

#' For each LRR learner, the best sieve is chosen by minimizing the double-robust one-step efficient empirical risk function over the choices of sieve-spaces.
#' @param trained_learner_list A unnested list of trained LRR learners as
#' @param learner_names Names of LRR learners. Should be `learner_names <- sapply(list_of_LRR_learners, `[[`, "name)`.
#' @param A A binary vector specifying the treatment assignment. The values should be in {0,1}.
#' This should include all observations including both the training and validation observations.
#' @param Y A numeric vector of binary or nonnegative observations of the outcome variable.
#' This should include all observations including both the training and validation observations.
#' @param EY1 A numeric vector containing initial cross-fitted estimates of E[Y|A=1,W] for all observations.
#' @param EY0 A numeric vector containing initial cross-fitted estimates of E[Y|A=0,W] for all observations.
#' @param pA1 A numeric vector containing initial cross-fitted estimates of P(A=1|W) for all observations.
#' @param weights A numeric vector of observation weights. If no special weighting desired, supply a vector of 1's.
subset_best_sieve <- function(trained_learner_list, learner_names, A, Y, EY1, EY0, pA1, weights) {
  LRR_learners <- trained_learner_list


  LRR_learners <- lapply(learner_names, function(learner_name) {
    keep <- which(stringr::str_detect(names(LRR_learners), quotemeta(learner_name)))

    sieve_learners <- LRR_learners[keep]

    # Get LRR predictions on fold-specific training set for all sieves
    # Single learner
    all_LRR_training <- lapply(sieve_learners, `[[`, "LRR_train")
    all_LRR_pred <- lapply(sieve_learners, `[[`, "LRR_pred")

    # Get DR risks on fold-specific training set for the LRR of the sieves
    all_training_risks <- unlist(lapply(all_LRR_training, function(LRR) {
      DR_risk_function_LRR(LRR, A, Y, EY1, EY0, pA1, weights)
    }))
    #print(all_training_risks)
    print("min")
    print(min(all_training_risks))
    print(which.min(all_training_risks))
    best_index <- which.min(all_training_risks)
    #names(best_index) <- names(all_training_risks)
    #print(names(all_training_risks)[best_index])
    #out <- list(sieve_learners[[best_index]])
    #names(out) <- names(all_training_risks)[best_index]
    all_LRR_training <- do.call(cbind, all_LRR_training)
    colnames(all_LRR_training) <- names(all_training_risks)
    #print(as.data.table(all_LRR_training))
    all_LRR_pred <- do.call(cbind, all_LRR_pred)
    return(list(LRR_train = all_LRR_training[,best_index], LRR_pred = all_LRR_pred[,best_index]))


  })
  names(LRR_learners) <- learner_names
  LRR_learners
}


