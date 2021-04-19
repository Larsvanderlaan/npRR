

#' @export
Lrnr_LRR_hal9001_discrete <- R6Class(
  classname = "Lrnr_LRR_hal9001", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(max_degree = 1,
                          smoothness_orders =0,
                          num_knots = c(100, 50),
                          penalize = F,

                          ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("LRR"),

    .train = function(task) {
      method <- self$params$method
      X <- task$X
      family <-  binomial()
      weights <- task$weights
      Y <- task$Y

      #smoothness <- self$params$smoothness_degree
      #smoothness <- rep(c(smoothness)[[1]], ncol(X))
      basis_list <- hal9001::enumerate_basis(as.matrix(task$X),  max_degree = self$params$max_degree, smoothness_orders = self$params$smoothness_orders, num_knots = self$params$num_knots)
      x_basis <- as.matrix(hal9001::make_design_matrix(as.matrix(task$X),basis_list))
      print("HAL dim")
      print(dim(x_basis))
      fit <- NULL
      tryCatch({
      fit <- speedglm::speedglm.wfit(Y, x_basis, family = binomial(), weights = weights, intercept = F)
      }, error = function(cond) {
        fit <- glm.fit(x_basis, Y, family = binomial(), weights = weights)
        fit <<- list(coef = coef(fit))
      })
      #fit <- glm.fit(x_basis, Y, family = binomial(), weights = weights)
      #fit <- list(coef = coef(fit))
        #fit <- hal9001::fit_hal(as.matrix(task$X), task$Y, weights = task$weights, family = binomial(), max_degree = self$params$max_degree, smoothness_orders = self$params$smoothness_orders, num_knots = self$params$num_knots)

      return(fit_object = list(fit = fit, basis_list = basis_list))
    },
    .predict = function(task = NULL) {
      fit_obj <- self$fit_object
      fit <- fit_obj$fit
      #preds <- predict(fit, new_data = as.matrix(task$X), type = "link")
       coef <- as.vector(fit$coef)
       coef[is.na(coef)] <- 0
       X <- as.matrix(task$X)
       x_basis <- as.matrix(hal9001::make_design_matrix(as.matrix(task$X),fit_obj$basis_list))

       preds <- x_basis %*% coef

      return(preds)
    },

    .required_packages = c("hal9001")
  )
)




#' @export
Lrnr_LRR_hal9001  <- R6Class(
  classname = "Lrnr_LRR_hal9001", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(max_degree = 1,
                          smoothness_orders = 1,
                          num_knots = c(100, 50),
                          ...) {
      params <- args_to_list()
      params$family <- binomial()
      lrnr_hal <- do.call(Lrnr_hal9001$new, params)
      params$lrnr_hal9001 <- lrnr_hal
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c( "LRR"),

    .train = function(task) {
      lrnr_hal <- params$lrnr_hal9001
      lrnr_hal <- lrnr_hal$train(task)
      return(fit_object = list(lrnr_hal = lrnr_hal))
    },
    .predict = function(task = NULL) {
      lrnr_hal <- fit_object$lrnr_hal
      preds <- lrnr_hal$predict(task)
      preds <- stats::qlogis(preds)
      return(preds)
    },

    .required_packages = c("hal9001")
  )
)

