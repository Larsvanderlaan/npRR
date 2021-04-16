
#' @import data.table
#' @import tmle3

make_task <- function(V, X = V, A, Y, weights = NULL,...) {

  if(is.vector(V)) {
    V <- as.matrix(V)
    colnames(V) <- "V"
  }
  if(is.vector(X)) {
    X <- as.matrix(X)
    colnames(X) <- "X"
  }
  colnames(V) <- paste0("V", 1:ncol(V))
  colnames(X) <- paste0("X", 1:ncol(X))

  if(is.null(weights)) {
    weights <- rep(1, length(Y))
  }
  data <- data.table(V, X, A = A, Y = Y, weights = weights)
  print(data)
  data <- data[,!duplicated(names(data)), with = F]
  print(data)
  npsem <- list(tmle3::define_node("V", colnames(V), c()),
                tmle3::define_node("X", colnames(X), c()),
                tmle3::define_node("A", "A", c("X")),
                tmle3::define_node("Y", "Y", c("A", "X")),
                tmle3::define_node("RR", c(), c("V")))

  task <- tmle3_Task$new(data, npsem, weights = "weights",...)
  return(task)
}

# Make likelihood
#' @import sl3
make_likelihood <- function(task, lrnr_A, lrnr_Y, cv = TRUE) {
  if(task$npsem$Y$variable_type$type == "binomial") {
    loss_Y <- loss_loglik_binomial
  } else {
    loss_Y <- loss_squared_error
  }
  if(cv) {
    lrnr_A <- make_learner(Pipeline, Lrnr_cv$new(lrnr_A), Lrnr_cv_selector$new(loss_loglik_binomial))
    lrnr_Y <- make_learner(Pipeline, Lrnr_cv$new(lrnr_Y), Lrnr_cv_selector$new(loss_Y))
  }
  factor_list <- list(LF_fit$new("A", lrnr_A, type = "density"),
                      LF_fit$new("Y", lrnr_Y, type = "mean"))
  likelihood <- Likelihood$new(factor_list)
  likelihood <- likelihood$train(task)
  return(likelihood)
}

# Generates revere
make_revere <- function(task, likelihood ) {
  genf <- make_generator(likelihood )
  return(sl3_revere_Task$new(genf, task))
}




make_generator <- function(likelihood) {




  gen_task_general <- function(tmle_task, fold_number) {
    task <- tmle_task$get_regression_task("RR")
    X <- task$X
    R <- tmle_task$get_tmle_node("Y")
    A <- tmle_task$get_tmle_node("A")
    g <- likelihood$get_likelihood(tmle_task, "A", fold_number)
    g <- tmle3::bound(g, 0.01)
    Q <- likelihood$get_likelihood(tmle_task, "Y", fold_number)
    cf_task1 <- tmle_task$generate_counterfactual_task(uuid::UUIDgenerate(), data.table::data.table(A= rep(1, length(g))))
    cf_task0 <- tmle_task$generate_counterfactual_task(uuid::UUIDgenerate(), data.table::data.table(A=rep(0, length(g))))
    Q1 <- bound(likelihood$get_likelihood(cf_task1, "Y", fold_number), 0.001)
    Q0 <- bound(likelihood$get_likelihood(cf_task0, "Y", fold_number),0.001)
    weightsIPW <- Y/g * task$weights
    weightsplugin <- (Q1 + Q0)
    YIPW <- A
    Yplugin <- Q1/weightsplugin
    covariates <- colnames(X)
    outcome <- c("Yplugin", "YIPW")
    weights <- c(tmle_task$nodes$weights)

    data <- cbind(data.table(g1 = ifelse(A==1, g, 1-g),  YIPW = YIPW, Yplugin = Yplugin, weightsplugin = weightsplugin, weightsIPW = weightsIPW, Q = Q, g = g, ginv = 1/g, Y = Y, A = A, Q1 = Q1, Q0 = Q0), X)
    data <- cbind(task$get_data(,weights), data)
    data$RR <- data$Q1 / data$Q0
    data$Qg1 <- data$Q1 / data$g1
    data$Qg0 <- data$Q0 / (1-data$g1)
    new_task <- sl3_Task$new(data, covariates = covariates, outcome = outcome, weights = weights, folds = task$folds)
    return(new_task)
  }


    return(gen_task_general)


}

make_eff_loss <- function(tmle_task, likelihood) {
  tmle_task <- tmle_task
  likelihood <- likelihood

  efficient_loss = function(preds, Y) {

    LRR <- preds

    R <- tmle_task$get_tmle_node("Y")
    A <- tmle_task$get_tmle_node("A")

    g <- likelihood$get_likelihood(tmle_task, "A", "validation")
    lik <- likelihood

    cf_task1 <- tmle_task$generate_counterfactual_task(uuid::UUIDgenerate(), data.table::data.table(A= rep(1, length(g))))
    cf_task0 <- tmle_task$generate_counterfactual_task(uuid::UUIDgenerate(), data.table::data.table(A=rep(0, length(g))))
    ER1 <- lik$get_likelihood(cf_task1, "Y", "validation")
    ER0 <- lik$get_likelihood(cf_task0, "Y", "validation")
    ER <- lik$get_likelihood(tmle_task, "Y", "validation")

    C1 <- A/g * (R - ER) + ER1
    C2 <- C1 + (1-A)/g * (R - ER) + ER0
    loss <- C1*-1*LRR + C2 * log(1 + exp(LRR))


    return(loss)
  }

}




make_true_loss <- function(tmle_task, likelihood) {

  n=10000
  library(simcausal)
  library(SuperLearner)
  library(data.table)
  library(sl3)
  bound <- Vectorize(tmle3::bound)
  D <- DAG.empty()
  D <- DAG.empty()
  D <- D +
    node("W1f", distr = "runif", min = -1, max = 1) +
    node("W2f", distr = "runif", min = -1, max = 1) +
    node("W3f", distr = "runif", min = -1, max = 1) +
    node("W1", distr = "rconst", const = 2*round(W1f/2,1)) +
    node("W2", distr = "rconst", const = 2*round(W2f/2,1)) +
    node("W3", distr = "rconst", const = 2*round(W3f/2,1)) +
    node("g", distr = "rconst", const = plogis(W1 + W2 +  W3 + W1*W2 + W1*W3 + W2*W3 ) )+
    node("A", distr = "rbinom", size = 1, prob = g )+
    node("gRtilde", distr = "rconst",  const =  plogis(W1 + W2 + W3 + A*(W1 + W2 + W3 ) + W1*W2 )) +
    node("gRtilde1", distr = "rconst",  const =     plogis(W1 + W2 + W3 + 1*(W1 + W2 + W3 ) + W1*W2 )) +
    node("gRtilde0", distr = "rconst",  const =   plogis(W1 + W2 + W3 + 0*(W1 + W2 + W3 ) + W1*W2 )) +
    node("gR", distr = "rconst",  const =  bound(gRtilde, 0.01) ) +
    node("R", distr = "rbinom", size = 1, prob = gR)+
    node("RR", distr = "rconst", const = gRtilde1/gRtilde0)

  setD <- set.DAG(D, vecfun = c("bound", "round"))
  data <- sim(setD, n = n)
  data <- setDT(data)
  data$id <- data$ID
  data$ID <- NULL
  data$t <- 0
  data$id <- as.factor(data$id)
  setkey(data, id ,t)
  npsem <- list(define_node("W", c("W1", "W2", "W3"), time = 0, variable_type = variable_type("continuous")), define_node("A", "A", c("W"), time = 0,  variable_type = variable_type("binomial")),

                define_node("R", "R", c("A", "W"), time = 0,  variable_type = variable_type("binomial")), define_node("RR", "R", c("W"), time = 0,  variable_type = variable_type("continuous")))
  tmle_task <- tmle3_Task$new(data, npsem, long_format = F)

  likelihood <- likelihood

  efficient_loss = function(preds, Y) {

    LRR <- preds

    R <- tmle_task$get_tmle_node("R", format = T)$R
    A <- tmle_task$get_tmle_node("A", format = T)$A

    g <- likelihood$get_likelihood(tmle_task, "A", "validation")
    lik <- likelihood

    cf_task1 <- tmle_task$generate_counterfactual_task(uuid::UUIDgenerate(), data.table::data.table(A= rep(1, length(g))))
    cf_task0 <- tmle_task$generate_counterfactual_task(uuid::UUIDgenerate(), data.table::data.table(A=rep(0, length(g))))
    ER1 <- lik$get_likelihood(cf_task1, "R", "validation")
    ER0 <- lik$get_likelihood(cf_task0, "R", "validation")
    ER <- lik$get_likelihood(tmle_task, "R", "validation")

    C1 <- A/g * (R - ER) + ER1
    C2 <- C1 + (1-A)/g * (R - ER) + ER0
    loss <- C1*-1*LRR + C2 * log(1 + exp(LRR))


    return(loss)
  }

}




