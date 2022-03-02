delayed_train_LRR_fold_learners <- function(fold, W, A, Y, EY1, EY0, pA1, weights, list_of_learners, list_of_sieves) {
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

  all_learners_delayed <- delayed_train_LRR_learners(W, A, Y, EY1, EY0, pA1, weights, list_of_learners, list_of_sieves, Wpred = Wval)
  return(all_learners_delayed)
}



choose_best_sieve_LRR_all_folds <- function(folds, trained_learner_list, learner_names, A, Y, EY1, EY0, pA1, weights) {
  output <- lapply(seq_along(folds), function(fold_number) {
    fold <- folds[[fold_number]]
    training_index <- origami::training(fold = fold)
    keep <- which(stringr::str_detect(names(LRR_learners_by_fold), paste0("^", fold_number, "\\.", "+")))
    LRR_learners <- LRR_learners_by_fold[keep]
    LRR_learners <- choose_best_sieve_LRR(LRR_learners, learner_names, A[training_index], Y[training_index], EY1[training_index], EY0[training_index], pA1[training_index], weights[training_index])

    return(LRR_learners)

  })
  names(output) <- seq_along(folds)
  return(output)
}

