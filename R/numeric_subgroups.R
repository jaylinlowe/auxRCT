
#' Divide observations into groups based on a single numeric variable
#'
#' @description `numeric_subgroups` attempts to divide the observations into
#' evenly sized groups based on `grouping col`. In some cases, there are too
#' many observations with the same value. This can happen if the variable in question
#' was 0 for a large number of observations, or if many of the observations
#' were originally missing and were replaced with mean-imputed values. In this case,
#' `numeric_subgroups` will create one subgroup of all observations with the
#' most frequent value and attempt to divide the remaining observations into evenly sized groups.
#' If the percentage of observations with the same value was extreme, observations will
#' be divided into two subgroups: those with the most frequent value and those with any other value.
#'
#' @param Y vector of responses
#' @param X a matrix or data frame of covariates
#' @param grouping_col A column in X, used to divide observations into subgroups
#' @param preds vector of predictions
#' @param max_groups number of groups to create (default 10)
#'
#' @return A vector specifying subgroup assignment, of the same length as Y.
#' Usually, the values in the vector will be intervals specifying how the numeric variable
#' was split. In the case of many observations with the exact same value, the formed subgroups
#' may instead be labeled as "mode" and "not mode".
#' @export
#'
numeric_subgroups <- function(Y, X, grouping_col, preds, max_groups = 10) {

  if(!grouping_col %in% colnames(X)) {
    stop("grouping_col must be the name of a column in X")
  }

  if (length(preds) != length(Y)) {
    stop("Y and X must be the same length")
  }
  if (length(Y) != nrow(X)) {
    stop("Incorrect dimensions for Y or X")
  }

  df <- data.frame(cbind(Y, X))

  df_num <- tryCatch({
    df$subgroups <- cut_number(df[[grouping_col]], max_groups)
    df_mod <- df},
    error = function(e) {
      sub_mode <- mode_calc(df, grouping_col)
      df$is_mode <- ifelse(df[[grouping_col]] == sub_mode, 1, 0)
      df_not_mode <- filter(df, is_mode == 0)

      df_mod <- tryCatch({
        df_not_mode$subgroups <- cut_number(df_not_mode[[grouping_col]], max_groups-1)
        df_mode <- filter(df, is_mode == 1)
        df_mode$subgroups <- rep(sub_mode, nrow(df_mode))
        df_mod <- rbind(df_mode, df_not_mode)},
        error = function(e) {
          df$subgroups <- ifelse(df[[grouping_col]] == sub_mode, as.character(sub_mode), "not mode")
          df_mod <- df
        })
    })

  #print(df_num)
  #ADJUST CODE - THIS SHOULD JUST RETURN A COLUMN WITH THE NUMERIC SPLITS, AND THEN WE RUN USER SPECIFIED GROUP
  #results <- user_specified_group(df_num[,1], df_num[,-1], grouping_col = "subgroups", preds, effect_size, alpha, beta)
  return(df_num$subgroups)
}
