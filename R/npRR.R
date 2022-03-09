

npRR <- function(LRR_learners, W, A, Y, V = W, weights, Vpred = V, EY1W, EY0W, pA1W, sl3_Learner_EYAW, sl3_Learner_pA1W, list_of_sieves, cross_validate_LRR = FALSE, folds = origami::make_folds(n=length(A))) {
  if(!is.list(LRR_learners)) {
    LRR_learners <- list(LRR_learners)
  }
  list_of_LRR_learners <- LRR_learners
  learner_names <- sapply(list_of_learners, `[[`, "name")

  W <- as.matrix(W)
  A <- as.vector(A)
  Y <- as.vector(Y)
  if(missing(weights)) {
    weights <- rep(1, length(A))
  }
  weights <- as.vector(weights)

  try({
    likelihood <- estimate_initial_likelihood(W, A, Y, weights, sl3_Learner_EYAW = sl3_Learner_EYAW, sl3_Learner_pA1W = sl3_Learner_pA1W, folds = folds, outcome_type = outcome_type)
    EY1W <- likelihood$EY1W
    EY0W <- likelihood$EY0W
    pA1W <- likelihood$pA1W
  })
  EY1W <- as.vector(EY1W)
  EY0W <- as.vector(EY0W)
  pA1W <- as.vector(pA1W)

  if(missing(list_of_sieves)) {
    if(ncol(W))
      list_of_sieves <- list(
        NULL,
        fourier_basis$new(orders = c(1,0)),
        fourier_basis$new(orders = c(2,0)),
        fourier_basis$new(orders = c(3,0)),
        fourier_basis$new(orders = c(1,1)),
        fourier_basis$new(orders = c(2,1)),
        ourier_basis$new(orders = c(3,1))
      )
  }


  full_fit_LRR <- train_LRR_learners(V=V, A, Y, EY1W, EY0W, pA1W, weights, list_of_LRR_learners, list_of_sieves, Vpred = W)
  all_LRR_full_best <- subset_best_sieve(full_fit_LRR, learner_names, A, Y, EY1W, EY0W, pA1W, weights)
  all_LRR_full_predictions <- do.call(cbind, lapply(all_LRR_full_best, `[[`, "LRR_pred"))
  print(dim(all_LRR_full_predictions))
  if(cross_validate_LRR) {
    LRR_learners_all_folds <- lapply(folds, train_LRR_learners_using_fold,  W, A, Y, EY1W, EY0W, pA1W, weights, list_of_LRR_learners, list_of_sieves )
    names(LRR_learners_all_folds) <- seq_along(folds)
    LRR_learners_all_folds <- unlist(LRR_learners_all_folds)
    learner_sieve_names <- names(LRR_learners_all_folds)
    LRR_learners_all_folds_delayed <- bundle_delayed(LRR_learners_all_folds)
    LRR_learners_all_folds <- LRR_learners_all_folds_delayed$compute()
    LRR_learners_best_sieve_all_folds <- subset_best_sieve_all_folds(folds, LRR_learners_all_folds, learner_names, A, Y, EY1W, EY0W, pA1W, weights)

    cv_predictions <- cv_predict_LRR_learner(folds, LRR_learners_best_sieve_all_folds)
    best_learner_index_cv <- which.min(DR_risk_function_LRR(cv_predictions, A , Y, EY1W, EY0W, pA1W, weights))
    print(dim(cv_predictions))
    all_LRR_full_predictions <- all_LRR_full_predictions[,best_learner_index_cv]
  }
  return(all_LRR_full_predictions)






}
