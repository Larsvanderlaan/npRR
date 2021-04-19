



#' Use SuperLearner Wrappers, Screeners, and Methods, in sl3
#'
#' @description These learners provide an interface to the wrapper functions,
#'  screening algorithms, and combination methods provided by the
#'  \code{SuperLearner} package. These components add support for a range of
#'  algorithms not currently implemented natively in \code{sl3}.
#'
#' @description \code{Lrnr_pkg_SuperLearner} - Interface for \code{SuperLearner}
#'  wrapper functions. Use \code{SuperLearner::listWrappers("SL")} for a list.
#'
#' @docType class
#'
#' @rdname SuperLearner_interface
#'
#' @importFrom R6 R6Class
#'
#' @export
#'
#' @keywords data
#'
#' @return Learner object with methods for training and prediction. See
#'  \code{\link{Lrnr_base}} for documentation on learners.
#'
#' @format \code{\link{R6Class}} object.
#'
#' @family Learners
#'
#' @section Parameters:
#' \describe{
#'   \item{\code{SL_wrapper}}{The wrapper function to use.}
#'   \item{\code{...}}{Currently not used.}
#' }
#'
#
Lrnr_pkg_SuperLearner <- R6Class(
  classname = "Lrnr_pkg_SuperLearner",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(SL_wrapper, ...) {
      wrapper_fun <- get(SL_wrapper)
      params <- list(wrapper_name = SL_wrapper, wrapper_fun = wrapper_fun, ...)
      super$initialize(params = params, ...)
    }
  ),

  private = list(
    .properties = c("binomial", "continuous", "weights", "ids", "wrapper"),
    .train = function(task) {
      args <- self$params
      wrapper <- args$wrapper_fun
      # to minimize prediction costs (since we throw out predictions from here
      # anyways), newX is just a single row
      newX <- task$X[1, ]
      outcome_type <- self$get_outcome_type(task)

      if (is.null(args$family)) {
        args$family <- outcome_type$glm_family(return_object = TRUE)
      }
      args$Y <- task$Y
      args$X <- task$X
      args$newX <- newX
      args$obsWeights <- task$weights
      args$id <-  task$id
      fit_object <- do.call(wrapper, args)$fit
      # fit_object <- wrapper(
      #   task$Y, task$X, newX,
      #   family = args$family,
      #   obsWeights = task$weights, id = task$id
      # )$fit
      return(fit_object)
    },

    .predict = function(task) {
      args <- self$params
      outcome_type <- private$.training_outcome_type

      if (is.null(args$family)) {
        args$family <- outcome_type$glm_family(return_object = TRUE)
      }
      predictions <- stats::predict(
        private$.fit_object,
        newdata = task$X,
        family = args$family
      )
      return(predictions)
    },
    .required_packages = c("SuperLearner")
  )
)










glmnet_screener <- Lrnr_pkg_SuperLearner_screener$new("screen.glmnet")


valid_learners <- c("SL.xgboost", "SL.earth", "SL.gam", "SL.bayesglm", "SL.glm", "SL.glm.interaction", "SL.glmnet")

valid_SL_learners <- function() {
  print(valid_learners)
}

SL_learner_to_LRR_learner <- function(name, ...) {
  if(F && !(name %in% valid_learners)) {
    stop(paste0(name, " is not a valid SL learner. Please call valid_SL_learners() to see valid learners,"))
  }
  lrnr <- make_learner(Lrnr_pkg_SuperLearner, name, family = binomial(), outcome_type = variable_type("binomial"), ...)
  binomial_to_LRR_learner(lrnr)
}

binomial_to_LRR_learner <- function(learner) {
  make_learner(Pipeline,learner , Lrnr_chainer_link$new())
}



# Default superlearner library for the nuisance parameters associated with A and Y
default_library_nuisance <- make_learner(Pipeline, glmnet_screener, Stack$new(list(
  Lrnr_xgboost$new(max_depth = 3),
  Lrnr_xgboost$new(max_depth = 4),
  Lrnr_xgboost$new(max_depth = 5),
  Lrnr_xgboost$new(max_depth = 6),
  #Lrnr_ranger$new(),
  Lrnr_glmnet$new(),
  Lrnr_gam$new()
  #Lrnr_earth$new()
)))
# Default superlearner library for the sequential regression nuisance parameter associated with V
default_library_V <-  Stack$new(list(
  make_learner(Lrnr_pkg_SuperLearner,"SL.xgboost" , max_depth = 3, family = binomial(), outcome_type = variable_type("binomial")),
  make_learner(Lrnr_pkg_SuperLearner,"SL.xgboost" , max_depth = 4, family = binomial(), outcome_type = variable_type("binomial")),
  make_learner(Lrnr_pkg_SuperLearner,"SL.xgboost" , max_depth = 5, family = binomial(), outcome_type = variable_type("binomial")),
  make_learner(Lrnr_pkg_SuperLearner,"SL.xgboost" , max_depth = 6, family = binomial(), outcome_type = variable_type("binomial")),
  make_learner(Lrnr_pkg_SuperLearner,"SL.gam" , deg.gam = 2, family = binomial(), outcome_type = variable_type("binomial")),
  make_learner(Lrnr_pkg_SuperLearner,"SL.gam" , deg.gam = 4, family = binomial(), outcome_type = variable_type("binomial")),
  #make_learner(Lrnr_pkg_SuperLearner,"SL.earth", degree = 2, pmethod = "forward", family = binomial(), outcome_type = variable_type("binomial")),
  make_learner(Lrnr_pkg_SuperLearner,"SL.glm" , family = binomial(), outcome_type = variable_type("binomial")),
  make_learner(Lrnr_pkg_SuperLearner,"SL.glm.interaction" , family = binomial(), outcome_type = variable_type("binomial"))
))


# Default superlearner library for relative risk function.
default_library_RR <- list(
  "xgboost_3" = Lrnr_LRR_xgboost$new(max_depth = 3),
  "xgboost_4" =  Lrnr_LRR_xgboost$new( max_depth = 4),
  "xgboost_5" =  Lrnr_LRR_xgboost$new( max_depth = 5),
  "xgboost_6" =  Lrnr_LRR_xgboost$new( max_depth = 6),
  "randomForest_6" =  Lrnr_LRR_xgboost$new(random_forest = TRUE, max_depth = 6),
  "randomForest_8" =  Lrnr_LRR_xgboost$new(random_forest = TRUE, max_depth = 8),
  "gam_2" = SL_learner_to_LRR_learner("SL.gam", deg.gam = 2),
  "gam_4" = SL_learner_to_LRR_learner("SL.gam", deg.gam = 4),
  "earth_1" <- SL_learner_to_LRR_learner("SL.earth", degree = 1, pmethod = "forward"),
  "earth_2" <- SL_learner_to_LRR_learner("SL.earth", degree = 2, pmethod = "forward"),
  "bayesglm" = SL_learner_to_LRR_learner("SL.bayesglm"),
  "glm" = SL_learner_to_LRR_learner("SL.glm"),
  "glm.interaction" = SL_learner_to_LRR_learner("SL.glm.interaction"),
  #"hal9001" = Lrnr_LRR_hal9001$new(max_degree = 2, smoothness_orders = 1, num_knots = c(30, 5)),
  Lrnr_LRR_subst$new()
  #"hal9001_2" = Lrnr_LRR_hal9001$new(max_degree = 2, smoothness_orders = 1, num_knots = c(40, 10)),
)
