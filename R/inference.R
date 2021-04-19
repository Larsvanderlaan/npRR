compute_TMLE <- function(initial_LRR, task_RR, tmle_task, likelihood, points = quantile(unlist(task_RR$X), seq(0.1, 0.9, length.out = 10))) {
  task <- tmle_task
  if(is.null(points)) {
    points <- quantile(unlist(task_RR$X), seq(0.1, 0.9, length.out = 10))
  }
  if(ncol(task_RR$X) !=1) {
    stop("Inference only available for one dimensional `V`.")
  }
  ys <- points
  V <- unlist(task_RR$X)
  task <- tmle_task
  preds <- initial_LRR
  sigmas <- c()
  for(y in ys) {
    sigma <- get_h(V, y, preds, task_RR, task, likelihood)
    sigmas <- c(sigmas, sigma)
  }

  kernel <- function(x) {
    fun <- function(i) {
      y <- ys[i]
      sigma <- sigmas[i]
      v <- exp(-(y-x)^2/sigma^2) /  sqrt(2*3.14*sigma^2)
      v
    }
    sapply(1:length(ys), fun)
  }

  Q1V <- task_RR$get_data(,"Q1V")[[1]]
  Q0V <- task_RR$get_data(,"Q0V")[[1]]
  Q1 <- task_RR$get_data(,"Q1")[[1]]
  Q0 <- task_RR$get_data(,"Q0")[[1]]
  g1 <- task_RR$get_data(,"g1")[[1]]
  A <- task_RR$get_data(,"A")[[1]]
  Y <- task_RR$get_data(,"Y")[[1]]
  update <- preds
  RR <- exp(update)
  clever_covariates <-   (Q1V + Q0V)/(Q1V*Q0V) * kernel(V)
  EIC <- A/(g1) *(-1/(1+ RR))* (Y-Q1) + (1-A)/((1-g1))*(RR/(1+ RR))*(Y - Q0) - (1/(1+ RR))*(Q1V) + (RR/(1+ RR))*(Q0V)
  weights <- apply(as.vector(EIC)*clever_covariates,2,sd)
  eff_loss <- make_eff_loss(task, likelihood)

  max_eps <- 5e-3
  for(i in 1:20) {
    RR <- exp(update)
    clever_covariates <-   (Q1V + Q0V)/(Q1V*Q0V) * kernel(V)
    EIC <- A/(g1) *(-1/(1+ RR))* (Y-Q1) + (1-A)/((1-g1))*(RR/(1+ RR))*(Y - Q0) - (1/(1+ RR))*(Q1V) + (RR/(1+ RR))*(Q0V)
    dir <- colMeans(as.vector(EIC)*clever_covariates) / weights
    dir <- dir/sqrt(mean(dir^2) )
    clever_covariates_one_dim <- clever_covariates %*% dir / colMeans(kernel(V))

    risk0 <- mean(eff_loss(update + 0))

    risk <- function(epsilon) {
      # RR <- exp(preds + epsilon*1)
      # score <- A/(g1) *(-1/(1+ RR))* (Y-Q1) + (1-A)/((1-g1))*(RR/(1+ RR))*(Y - Q0) - (1/(1+ RR))*(Q1) + (RR/(1+ RR))*(Q0)
      # score <- abs(mean(score * clever_covariates))
      loss <- (eff_loss(update  + epsilon*clever_covariates_one_dim ))
      loss <- loss
      risk <- mean(loss)
    }



    optim_fit <- optim(
      par = list(epsilon = 0 ), fn = risk,
      lower = -max_eps, upper = max_eps,
      method = "Brent")

    epsilon <- optim_fit$par
    update <- update + epsilon* (clever_covariates%*% dir)

    RR <- exp(update)
    kern <- kernel(V)
    score1 <- A/(g1) *(-1/(1+ RR))* (Y-Q1) + (1-A)/((1-g1))*(RR/(1+ RR))*(Y - Q0) - (1/(1+ RR))*(Q1V) + (RR/(1+ RR))*(Q0V)

    score1a <- clever_covariates*as.vector(score1)

    if(sqrt(mean((colMeans(score1a)/ weights)^2)) <= 1/sqrt(n)/log(n)){
      break
    }
    if(abs(optim_fit$par) < 1e-8) {
      max_eps <- max_eps/3
    }
  }

  est_X <- kern * as.vector(update)
  psi <- colMeans(est_X)
  EIC <- score1a + t(t(est_X) - psi)
  EIC <- t(t(EIC) / colMeans(kern)) + t( -(psi/colMeans(kern)^2)*(t(kern) - colMeans(kern)))
  psi <- psi/colMeans(kern)
  radius <- 1.96*apply(EIC,2,sd)/sqrt(n)
  upper <- exp(psi + radius)
  lower <- exp(psi - radius)
  psi <- exp(psi)
  ests <- cbind(points, psi, lower, upper)
  colnames(ests) <- c("V", "RR", "lower_CI", "upper_CI")
  output <- list(estimates = ests, RR_targeted = exp(update), EIC = EIC)
  return(output)
}

get_h <- function(V, y, preds, task_RR, tmle_task, likelihood) {
  task <- tmle_task
  eff_loss <- make_eff_loss(task, likelihood)
  grid <- sd(V) * c(0.01, seq(0.05, 1, length = 25) )
  grid <- rev(grid)
  Q1V <- task_RR$get_data(,"Q1V")[[1]]
  Q0V <- task_RR$get_data(,"Q0V")[[1]]
  Q1 <- task_RR$get_data(,"Q1")[[1]]
  Q0 <- task_RR$get_data(,"Q0")[[1]]
  g1 <- task_RR$get_data(,"g1")[[1]]
  A <- task_RR$get_data(,"A")[[1]]
  Y <- task_RR$get_data(,"Y")[[1]]
  upper_CIs <- c()
  lower_CIs <- c()
  ests <- c()
  u <- max(V)
  l <- min(V)
  for(sigma in grid) {
    update <- preds

    kernel <- function(x) {
      fun <- function(y) {
        v <- exp(-(y-x)^2/sigma^2) /  sqrt(2*3.14*sigma^2)
        v
      }

      return(as.matrix(fun(y)))

    }

    for(i in 1:1) {
      clever_covariates <-   (Q1V + Q0V)/(Q1V*Q0V) * kernel(V)


      risk <- function(epsilon) {
        # RR <- exp(preds + epsilon*1)
        # score <- A/(g1) *(-1/(1+ RR))* (Y-Q1) + (1-A)/((1-g1))*(RR/(1+ RR))*(Y - Q0) - (1/(1+ RR))*(Q1) + (RR/(1+ RR))*(Q0)
        # score <- abs(mean(score * clever_covariates))
        loss <- (eff_loss(update  + epsilon*clever_covariates ))
        loss <- loss

        risk <- mean(loss)
      }


      max_eps <- 0.05
      optim_fit <- optim(
        par = list(epsilon = 0 ), fn = risk,
        lower = -max_eps, upper = max_eps,
        method = "Brent")

      epsilon <- optim_fit$par

      update <- update + epsilon* (clever_covariates)
    }
    RR <- exp(update)
    kern <- kernel(V)
    score1 <- A/(g1) *(-1/(1+ RR))* (Y-Q1) + (1-A)/((1-g1))*(RR/(1+ RR))*(Y - Q0) - (1/(1+ RR))*(Q1V) + (RR/(1+ RR))*(Q0V)
    score1a <- clever_covariates*as.vector(score1)




    est_X <- kern * as.vector(update)
    psi <- colMeans(est_X)
    EIC <- score1a + t(t(est_X) - psi)
    psi <- psi/ mean(kern)
    EIC <- EIC/mean(kern) - (psi/mean(kern)^2)*(kern - mean(kern))
    radius <- 1.96*apply(EIC,2,sd)/sqrt(n)
    upper <- psi + radius
    lower <- psi - radius
    upper_CIs <- c(upper_CIs, upper)
    lower_CIs <- c(lower_CIs, lower)
    ests <- c(ests, psi)

  }
  # Estimates increasing

  if(mean(diff(ests[5:10])) > 0) {
    best_index <- which(diff(sign(diff(lower_CIs)))==-2)+1

    best_index <- min(best_index)

  } else {

    best_index <- which(diff(sign(diff(-upper_CIs)))==-2)+1

    best_index <- min(best_index)

  }
  if(is.infinite(best_index)) {

    best_index <- length(grid)
  }

  print(best_index)

  return(grid[best_index])
}

#
#
#
#
#
#
#
# ys <- 0.2
# sigma <- 0.3
# kernel <- function(x) {
#   fun <- function(y) {exp(-(y-x)^2/sigma^2) / sqrt(2*3.14*sigma^2)}
#   sapply(ys, fun)
# }
# true <- colMeans(kernel(data_full$W1) *  log(lrnr1$predict(tsk)/ lrnr0$predict(tsk)))
#
# all_data <- task_RR$get_data()
# Q1V <- all_data$Q1V
# Q0V <- all_data$Q0V
# Q1 <- all_data$Q1
# Q0 <- all_data$Q0
# g1 <- all_data$g1
# update <- preds
# RR <- exp(update)
# clever_covariates <-   (Q1V + Q0V)/(Q1V*Q0V) * kernel(V)
# EIC <- A/(g1) *(-1/(1+ RR))* (Y-Q1) + (1-A)/((1-g1))*(RR/(1+ RR))*(Y - Q0) - (1/(1+ RR))*(Q1V) + (RR/(1+ RR))*(Q0V)
# weights <- apply(as.vector(EIC)*clever_covariates,2,sd)
# print(weights)
# for(i in 1:20) {
#   RR <- exp(update)
#   clever_covariates <-   (Q1V + Q0V)/(Q1V*Q0V) * kernel(V)
#   EIC <- A/(g1) *(-1/(1+ RR))* (Y-Q1) + (1-A)/((1-g1))*(RR/(1+ RR))*(Y - Q0) - (1/(1+ RR))*(Q1V) + (RR/(1+ RR))*(Q0V)
#   dir <- colMeans(as.vector(EIC)*clever_covariates) / weights
#   dir <- dir/sqrt(mean(dir^2) )
#   clever_covariates_one_dim <- clever_covariates %*% dir
#
#   risk0 <- mean(eff_loss(update + 0))
#
#   risk <- function(epsilon) {
#     # RR <- exp(preds + epsilon*1)
#     # score <- A/(g1) *(-1/(1+ RR))* (Y-Q1) + (1-A)/((1-g1))*(RR/(1+ RR))*(Y - Q0) - (1/(1+ RR))*(Q1) + (RR/(1+ RR))*(Q0)
#     # score <- abs(mean(score * clever_covariates))
#     risk <- mean(eff_loss(update + epsilon*clever_covariates_one_dim))
#     risk
#   }
#
#
#
#   optim_fit <- optim(
#     par = list(epsilon = 0), fn = risk,
#     lower = -0.01, upper = 0.01,
#     method = "Brent")
#   print(risk0)
#   print(optim_fit$value)
#   print(risk(optim_fit$par))
#   epsilon <- optim_fit$par
#   update <- update + epsilon* (clever_covariates%*% dir)
#   print(optim_fit$par)
#   RR <- exp(update)
#   kern <- kernel(V)
#   score1 <- A/(g1) *(-1/(1+ RR))* (Y-Q1) + (1-A)/((1-g1))*(RR/(1+ RR))*(Y - Q0) - (1/(1+ RR))*(Q1V) + (RR/(1+ RR))*(Q0V)
#   print(mean(score1))
#   score1a <- clever_covariates*as.vector(score1)
#   print(("Score"))
#   print(colMeans(score1a))
#   print(sqrt(mean((colMeans(score1a)/ weights)^2)))
#   print(1/sqrt(n)/log(n))
#   if(sqrt(mean((colMeans(score1a)/ weights)^2)) <= 1/sqrt(n)/log(n)){
#     break
#   }
#   if(abs(optim_fit$par) < 1e-9) {
#     break
#   }
# }
#
# est_X <- kern * as.vector(update)
# psi <- colMeans(est_X)
# EIC <- score1a + t(t(est_X) - psi)
# print(mean(EIC))
# radius <- 1.96*apply(EIC,2,sd)/sqrt(n)
# print(sd(EIC))
# upper <- psi + radius
# lower <- psi - radius
# ests <- cbind(psi, lower, upper)
# colnames(ests) <- c("psi", "lower", "upper")
# print(ests)
# output_list[[i]] <- ests
# pass <- (true >= lower) & (true <= upper)
# passed <- cbind(passed, pass)
# print(apply(passed,1,table))
# print(rowMeans(passed))
# }
#
# ```
#



