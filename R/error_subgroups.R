#' Divide observations into groups based on OOB predicted error of original residuals
#'
#' @description
#' `error_subgroups` runs a random forest to predict the absolute value of
#' the residuals calculated from `Y` and `preds` using all variables
#' specified in `variables`. Observations are divided into
#' subgroups based on the values of the OOB predicted error from this random
#' forest.
#'
#'
#' @param Y a vector of responses
#' @param X a matrix or data frame of covariates
#' @param preds a vector of predicted values
#' @param variables a vector of variable names in X,
#' May also be NULL, in which case all variables in X are used.
#' @param num_groups the number of subgroups to create
#'
#' @return A vector specifying subgroup assignment, of the same length as Y.
#' Smaller values indicate the predicted error was smaller.
#'
#' @export
#'
error_subgroups <- function(Y, X, preds, variables = NULL, num_groups = 10) {
  if (!is.null(variables)) {
    if (sum(variables %in% colnames(X)) != length(variables)) {
      stop("variables must be a list of column names in X")
    }
    df <- cbind(Y, X) %>%
      data.frame() %>%
      select(Y, all_of(variables))
  }
  else {
    df <- data.frame(cbind(Y, X))
  }


  df$abs_errors <- abs(df$Y - preds)
  error.model <- randomForest::randomForest(abs_errors ~ . - Y, data = df)
  df$error_preds <- predict(error.model)
  df$subgroups <- ntile(df$error_preds, num_groups)

  return(df$subgroups)
}
