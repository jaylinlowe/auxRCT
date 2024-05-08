#' Calculate sample sizes by subgroup
#'
#' @param Y a vector of responses
#' @param X matrix of covariates
#' @param grouping_col name of a column in X that defines subgroups
#' @param preds vector of predictions
#' @param effect_size desired minimal detectable difference
#' @param alpha desired Type II error rate
#' @param beta desired Type II error rate
#'
#' @return A data frame with 7 columns: `def`, `num`, `MSE`, `variance`,`resid_var`,
#' `samp_size`, and `samp_size_without`.
#' `def` is the definition of the subgroup, which may be an interval or a single
#' factor level, depending on the values in `grouping_col`. `num` is the number of
#' observations in that subgroup. `MSE` is the mean squared error of the predictions
#' for that subgroup. `variance` is the variance Y of the subgroup. `resid_var` is
#' the variance of the residuals calculated from `preds` and `Y`. `resid_var` is assumed
#' to be \eqn{\sigma^2} and is passed through the `samp_size` function. `samp_size`
#' is the resulting estimated sample size needed if the RCT population resembles this
#' subgroup. `samp_size_without` is the estimated sample size needed if the RCT
#' population resembles this subgroup and reloop is not utilized.
#' @export
#'
reloop.samp <- function(Y, X, grouping_col, preds, effect_size, alpha, beta) {

  if (!grouping_col %in% colnames(X)) {
    stop("grouping_col must be the name of a column in X")
  }

  df <- cbind(Y, X, preds)

  groups <- unique(X[[grouping_col]])
  results <- data.frame()
  indices <- list()
  for (i in 1:length(groups)) {
    subgroup_index <- which(df[[grouping_col]] == groups[i])
    indices[[i]] <- subgroup_index
    subgroup <- df[subgroup_index, ]
    r <- data.frame(def = groups[i], 'MSE' = sum((subgroup$preds - subgroup$Y)^2)/nrow(subgroup), 'num' = nrow(subgroup),
                    'variance' = var(subgroup$Y), "resid_var" = var(subgroup$preds-subgroup$Y))
    results <- rbind(results, r)
  }

  results$samp_size <- unlist(lapply(as.numeric(results$resid_var),samp_size, effect_size = effect_size, alpha = alpha, beta = beta))
  results$samp_size_without <- unlist(lapply(as.numeric(results$variance), samp_size, effect_size = effect_size, alpha = alpha, beta = beta))

  #the sorting is going to be tricky here, we have to resort indices in the same way?
  sort_order <- order(results$def)

  #got to be a more efficient way to do this? TO - DO
  indices_new <- list()
  for (j in 1:length(indices)) {
    indices_new[[j]] <- indices[[sort_order[j]]]
  }

  return(list(arrange(results, def), indices_new))
}
