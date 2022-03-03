
#' Train the LRR learners for all sieves and base learners on the training set specified by a fold object.
#' @param fold \code{fold} object from the \link{origami} package specifying the fold-specific training set at which to estimate the LRR.
#' @param W A matrix of covariate observations. This should include all observations including both the training and validation observations.
#' @param A A binary vector specifying the treatment assignment. The values should be in {0,1}.
#' This should include all observations including both the training and validation observations.
#' @param Y A numeric vector of binary or nonnegative observations of the outcome variable.
#' This should include all observations including both the training and validation observations.
#' @param EY1 A numeric vector containing initial cross-fitted estimates of E[Y|A=1,W] for all observations.
#' @param EY0 A numeric vector containing initial cross-fitted estimates of E[Y|A=0,W] for all observations.
#' @param pA1 A numeric vector containing initial cross-fitted estimates of P(A=1|W) for all observations.
#' @param weights A numeric vector of observation weights. If no special weighting desired, supply a vector of 1's.
#' @param list_of_LRR_learners A list of untrained \code{sl3_Learner} objects to be used to estimate the log relative risk LRR using the function \link{estimate_LRR_using_ERM}.
#' @param list_of_sieves A list of basis_generator objects specifying the sieve. See, for example, \code{fourier_basis} for an example and template.
train_LRR_learners_using_fold <- function(fold, W, A, Y, EY1, EY0, pA1, weights, list_of_LRR_learners, list_of_sieves) {
  list_of_learners <- list_of_LRR_learners
  index <- origami::training(fold = fold)
  index_val <- origami::validation(fold = fold)
  Wfull <- W
  W <- Wfull[index,]
  Wval <- Wfull[index_val,]
  A <- A[index]
  Y <- Y[index]
  EY1 <- EY1[index]
  EY0 <- EY0[index]
  pA1 <- pA1[index]
  weights <- weights[index]

  all_learners_delayed <- train_LRR_learners(W, A, Y, EY1, EY0, pA1, weights, list_of_learners, list_of_sieves, Wpred = Wval, compute = FALSE)
  return(all_learners_delayed)
}

#' For each LRR learner and fold, the best sieve is chosen by minimizing the double-robust one-step efficient empirical risk function over the choices of sieve-spaces.
#' @param folds A \code{folds} object (list of \code{fold} objects) from the \link{origami} package specifying the fold-specific training set at which to estimate the LRR.
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
subset_best_sieve_all_folds <- function(folds, trained_LRR_learner_list, learner_names, A, Y, EY1, EY0, pA1, weights) {
  trained_learner_list <- trained_LRR_learner_list
  output <- lapply(seq_along(folds), function(fold_number) {
    fold <- folds[[fold_number]]
    training_index <- origami::training(fold = fold)
    keep <- which(stringr::str_detect(names(LRR_learners_by_fold), paste0("^", fold_number, "\\.", "+")))
    LRR_learners <- LRR_learners_by_fold[keep]
    LRR_learners <- subset_best_sieve(LRR_learners, learner_names, A[training_index], Y[training_index], EY1[training_index], EY0[training_index], pA1[training_index], weights[training_index])

    return(LRR_learners)

  })
  names(output) <- seq_along(folds)
  return(output)
}

