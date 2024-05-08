
#' Calculate necessary sample size
#'
#' @description `samp_size` calculates the needed sample size for the entire
#' experiment to achieve the stated power, type I error, and effect size,
#'assuming equally sized groups.
#'
#' @details Sample size is calculated using the following equation:
#' \deqn{n = 4 \sigma^2 \frac{(Z_{\alpha/2} + Z_{\beta})^2}{\Delta^2}}
#' where \eqn{\Delta} is the effect size, \eqn{Z_{\alpha/2}} is \eqn{1-\alpha/2}
#' quantile value from a normal distribution, and \eqn{Z_{\beta}} is the \eqn{beta}
#' quantile value from a normal distribution. If `type = "one.sided"`, then
#' \eqn{Z_{\alpha/2}} is replaced with \eqn{Z_{\alpha}}.
#'
#' @param sigma2 An estimate of the variance of the population
#' @param effect_size desired minimal detectible difference
#' @param alpha desired Type II error rate
#' @param beta desired Type II error rate
#' @param type A character string specifying the type of test. Must be either
#' "two.sided" (default) or "one.sided"
#'
#' @return The number of observations needed for the overall experiment
#' @export
#'
samp_size <- function(sigma2, effect_size, alpha = 0.05, beta = 0.20 , type = "two.sided") {
  if (!type %in% c("one.sided", "two.sided")) {
    stop("type must be either 'one.sided' or 'two.sided'")
  }

  alpha_critical <- ifelse(type == "two.sided", qnorm(1-alpha/2), qnorm(1-alpha))
  beta_critical <- qnorm(1-beta)

  N <- (4 * sigma2 * (alpha_critical + beta_critical)^2)/(effect_size^2)
  return(ceiling(N))
}
