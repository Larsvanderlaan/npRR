#' Function to compute initial estimates of nuisance functions.
#' @param W A column-named matrix of baseline variables.
#' @param A A binary vector with values in {0,1} encoding the treatment assignment.
#' @param Y A numeric vector storing the outcome values.
#' @param sl3_Learner_pA1 A \code{sl3_Learner} object from the \code{tlverse/sl3} R github package that specifies the machine-learning algorithm for learning the propensity score `P(A = 1 | W)`
#' @param sl3_Learner_EY A \code{sl3_Learner} object from the \code{tlverse/sl3} R github package that specifies the machine-learning algorithm for learning the outcome conditional mean `E[Y | A, W]`. NOTE: the treatment arms are pooled in the regression. See the preprocessing sl3_Learner \code{Lrnr_stratified} if you wish to stratify the estimation by treatment.
#' @param folds A number representing the number of folds to use in cross-fitting or a fold object from the package \code{tlverse/origami}. This parameter will be passed to internal \code{sl3_Task} objects that are fed to the code{sl3_Learner}s.
estimate_initial_likelihood <- function(W, A, Y, weights = NULL, sl3_Learner_pA1, sl3_Learner_EY, folds = 10) {

  data <- data.table(W, A = A, Y = Y, weights = weights)
  covariates <- colnames(W)

  task_pA1 <- sl3_Task$new(data, covariates = covariates, outcome = "A", outcome_type = "binomial", weights = "weights", folds = folds)
  folds <- task_pA1$folds
  task_EY <- sl3_Task$new(data, covariates = c(covariates, "A"), outcome = "Y", weights = "weights", folds = folds)
  sl3_Learner_pA1 <- delayed_learner_train(sl3_Learner_pA1, task_pA1)
  sl3_Learner_EY <- delayed_learner_train(sl3_Learner_EY, task_EY)
  delayed_learner_list <- bundle_delayed(list(sl3_Learner_pA1,sl3_Learner_EY ))
  trained_learners <- delayed_learner_list$compute(progress = FALSE)

  sl3_Learner_pA1_trained <- trained_learners[[1]]
  sl3_Learner_EY_trained <- trained_learners[[2]]

  data1 <- data.table::copy(data)
  data0 <- data.table::copy(data)
  data1$A <- 1
  data0$A <- 0
  task_EY1 <- sl3_Task$new(data1, covariates = c(covariates, "A"), outcome = "Y", weights = "weights", folds = folds)
  task_EY0 <- sl3_Task$new(data0, covariates = c(covariates, "A"), outcome = "Y", weights = "weights", folds = folds)

  EY <- sl3_Learner_EY_trained$predict(task_EY)
  EY1 <- sl3_Learner_EY_trained$predict(task_EY1)
  EY0 <- sl3_Learner_EY_trained$predict(task_EY0)
  pA1 <- sl3_Learner_pA1_trained$predict(task_pA1)
  if(any(EY != ifelse(A==1, EY1, EY0))) {
    stop("EY and EY1, EY0 are inconsistent.")
  }


  internal <-  list(task_pA1 = task_pA1, task_EY = task_EY, sl3_Learner_pA1_trained = sl3_Learner_pA1_trained, sl3_Learner_EY_trained = sl3_Learner_EY_trained, folds = folds)
  output <- list(pA1 = pA1, EY1 = EY1, EY0 = EY0, internal = internal)
  return(output)
}
