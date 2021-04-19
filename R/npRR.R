npRR <- function(V, X, A, Y, RR_learners = "autoML", lrnr_A = "autoML", lrnr_Y = "autoML", lrnr_V = ifelse(ncol(V)!= ncol(X), "autoML", "NULL"), cv_nuisance = TRUE, cv_RR = TRUE, basis_list = make_fourier_basis_list(V)) {
  V <- as.matrix(V)
  X <- cbind(V,as.matrix(X))
  A <- as.vector(A)
  Y <- as.vector(Y)
  if(RR_learners == "autoML") {
    RR_learners <- default_library_RR
  }

  if(is.character(lrnr_A) && lrnr_A == "autoML") {
    lrnr_A <-  default_library_nuisance
  }
  if(is.character(lrnr_Y) && lrnr_Y == "autoML") {
    lrnr_Y <-  default_library_nuisance
  }
  if(ncol(V) > 10) {
    screen.glmnet <- TRUE
  }
  if(is.character(lrnr_V) && lrnr_V == "autoML") {
    lrnr_V <-  default_library_V
    if(screen.glmnet) {
      lrnr_V <- make_learner(Pipeline, Lrnr_LRR_glmnet.screener$new(), lrnr_V)
    }
  } else if (is.character(lrnr_V) &&lrnr_V == "NULL") {
    lrnr_V <- NULL
  }
  task <- make_task(V, X, A, Y, folds = 10)
  print("Fitting likelihoods")
  likelihood <- make_likelihood(task, lrnr_A ,lrnr_Y, lrnr_V, cv = cv_nuisance )
  print("Done Fitting likelihoods")
  genr <- make_generator(likelihood)
  task_RR <- genr(task, "validation")

  base_lrnr_lrr <- generate_plugin_learners(task_RR, basis_list, RR_learners, task, likelihood, select_sieve = TRUE, screen.glmnet = screen.glmnet, cv = cv_RR)
  lrnr_lrr <- base_lrnr_lrr$train(task_RR)
  predictions <- exp(lrnr_lrr$predict(task_RR))
  if(cv_RR) {
    reduced_lrnr <- lrnr_lrr$fit_object$full_fit$fit_object$learner_fits$Stack$fit_object$learner_fits[[which(out$learner_LRR_trained$fit_object$cv_meta_fit$coefficients==1)]]
  } else {
    reduced_lrnr <- lrnr_lrr
  }
  output <- list(RR = predictions, tmle_task = task, task_RR = task_RR, likelihood = likelihood, learner_LRR_trained = lrnr_lrr, learner_LRR_trained_squashed = reduced_lrnr)
  class(output) <- "npRR"
  return(output)
}



npRR_inference <- function(output, points = NULL) {
  V <- as.matrix(output$task_RR$X)
  output_inference <- compute_TMLE(log(output$RR), output$task_RR, output$tmle_task, output$likelihood)
  df <- as.data.frame(output_inference$estimates)
  plt <- ggplot(df, aes(x = V, y = RR)) + geom_line() + geom_ribbon(aes(ymin=lower_CI, ymax=upper_CI), alpha=0.2) + xlab("Value of V") + ylab("Relative Risk")
  output_inference$plot <- plt
  return(output_inference)
}

predict.npRR <- function(object, newV) {
  newV <- as.matrix(newV)
  Anew <- rep(1, nrow(newV))
  Ynew <- rep(1, nrow(newV))
  X <- as.matrix(object$tmle_task$get_tmle_node("X"))
  Xnew <- apply(X, 2, function(x) {
    if(length(x) < nrow(newV)) {
      x <- suppressWarnings(x + rep(0, nrow(newV)))
    } else {
      x <- x[1:nrow(newV)]
    }
    x <- as.matrix(x)
    return(x)
  })
  if(!is.matrix(Xnew)) {
    Xnew <- t(as.matrix(Xnew))
  }

  colnames(Xnew) <- colnames(X)

  print(dim(Xnew))
  print(dim(newV))
  task <- make_task(newV, Xnew, Anew, Ynew, folds = 10)
  print(task)
  likelihood <- object$likelihood
  genr <- make_generator(likelihood)
  task_RR <- genr(task, "validation")
  print(task_RR)
  return(exp(object$learner_LRR_trained_squashed$predict(task_RR)))
}



