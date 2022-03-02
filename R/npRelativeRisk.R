






delayed_train_LRR_by_sieve <- function(learner_LRR, W, A, Y, weights, list_of_sieve_nuisances, Wpred = W) {
  list_of_sieve_LRR <- lapply(list_of_sieve_nuisances, function(sieve_nuisances) {
    EY1_star <- sieve_nuisances$EY1_star
    EY0_star <- sieve_nuisances$EY0_star
    pA1_star <- sieve_nuisances$pA1_star

    delayed_plugin_LRR <- delayed_fun(estimate_LRR_using_ERM)(W, A, Y,  EY1_star, EY0_star, pA1_star, weights, learner_LRR, learning_method = "plugin", Wpred = Wpred)

    delayed_IPW_LRR <- delayed_fun(estimate_LRR_using_ERM)(W, A, Y,  EY1_star, EY0_star, pA1_star, weights, learner_LRR, learning_method = "IPW", Wpred = Wpred)

    output <- (list( plugin = delayed_plugin_LRR, IPW = delayed_IPW_LRR))
  })
  sieve_names <- sapply(list_of_sieve_nuisances, `[[`, "sieve")
  names(list_of_sieve_LRR) <- sieve_names
  return(list_of_sieve_LRR)
}

delayed_train_LRR_learners <- function(W, A, Y, EY1, EY0, pA1, weights, list_of_learners, list_of_sieves, Wpred = W) {

  list_of_sieve_nuisances <- lapply(list_of_sieves, function(sieve){
    compute_plugin_and_IPW_sieve_nuisances(basis_generator = sieve, W = W, A = A, Y = Y, EY1 = EY1, EY0 = EY0, pA1 = pA1, weights = weights)})
  list_of_sieve_nuisances

  all_learners_delayed <- lapply(list_of_learners, delayed_train_LRR_by_sieve, list_of_sieve_nuisances = list_of_sieve_nuisances, W = W, A = A, Y = Y, weights = weights, Wpred = Wpred)
  learner_names <- lapply(list_of_learners, `[[`, "name")
  names(all_learners_delayed) <- paste0(learner_names)
  return(all_learners_delayed)
}


choose_best_sieve_LRR <- function(trained_learner_list, learner_names, A, Y, EY1, EY0, pA1, weights) {
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


