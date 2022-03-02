estimate_LRR_using_ERM <- function(W, A, Y,  EY1, EY0, pA1, weights, sl3_LRR_Learner_binomial, learning_method = c("plugin", "IPW"), Wpred = W, untransform_logit = TRUE) {
  learning_method <- match.arg(learning_method)
  data <- as.data.table(W)

  covariates <- colnames(data)
  if(learning_method == "plugin") {
    pseudo_outcome <- EY1 / (EY1 + EY0)
    pseudo_weights <- weights * (EY1 + EY0)
    data$pseudo_outcome <- pseudo_outcome
    data$pseudo_weights <- pseudo_weights
    task_LRR <- sl3_Task$new(data, covariates = covariates, outcome = "pseudo_outcome", weights = "pseudo_weights", outcome_type = "quasibinomial")
  }
  else if(learning_method == "IPW") {
    pseudo_outcome <- A
    pseudo_weights <- weights * Y / ifelse(A==1, pA1, 1 - pA1)
    data$pseudo_outcome <- pseudo_outcome
    data$pseudo_weights <- pseudo_weights
    task_LRR <- sl3_Task$new(data, covariates = covariates, outcome = "pseudo_outcome", weights = "pseudo_weights", outcome_type = "binomial")
  }
  task_LRR_pred <- sl3_Task$new(as.data.table(Wpred), covariates = covariates)

  sl3_Learner_LRR_trained <- sl3_LRR_Learner_binomial$train(task_LRR)
  LRR <- sl3_Learner_LRR_trained$predict(task_LRR)
  LRR_pred <- sl3_Learner_LRR_trained$predict(task_LRR_pred)
  if(untransform_logit) {
    LRR <- qlogis(LRR)
    LRR_pred <- qlogis(LRR_pred)
  }

  output <- list(LRR_train = as.matrix(LRR), LRR_pred = as.matrix(LRR_pred), LRR_learner = sl3_Learner_LRR_trained)
  return(output)
}



DR_risk_function_LRR <- function(LRR, A, Y, EY1, EY0, pA1, weights, debug = FALSE, return_loss = FALSE) {
  LRR <- as.matrix(LRR)
  if(!(nrow(LRR) == length(A) && nrow(LRR) == length(EY1))) {
    stop("Input lengths dont match")
  }
  EY <- ifelse(A==1, EY1, EY0)
  plugin_risk <- (EY0 + EY1) * log(1 + exp(LRR)) - EY1 * LRR
  score_comp <- (A/pA1)*(log(1 + exp(LRR)) - LRR)*(Y - EY) + ((1-A)/(1-pA1))*(log(1 + exp(LRR)) - LRR)*(Y - EY)
  if(debug){
    print(colMeans(weights * score_comp))
  }
  DR_loss <- weights * (plugin_risk + score_comp)
  if(return_loss) {
    return(DR_loss)
  } else {
    return(colMeans(DR_loss))
  }
}
