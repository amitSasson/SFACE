
#' @title Max lambda1 and lambda2 values
#' @description In order to perform the sensitivity analysis, the researcher must choose which values of lambda1 and lambda2 to explore. We give here an upper limit for lambda1 and lambda2, estimated from the data.
#' @param stand_formula A formula for standartization and DR, y ~ A + X, the outcome as a function of the exposure and covariates
#' @param exposure The treatment/exposure vector pf length n. Must be encoded 1 for treated and 0 for untreated.
#' @param df a data frame with columns for the outcome, exposure and covariates.
#' @param weight A numerical vector of length n, holding weights to adjust for missing subtyps, Default: 1
#' @examples
#' A <- rbinom(n = 1000, size = 1, prob = 0.5)
#' X1 <- rbinom(n = 1000, size = 1, prob = 0.5)
#' X2 <- rnorm(n = 1000, mean = 0, sd = 1)
#' y <- sample(c(0,1,2), 1000, replace=TRUE, prob=c(0.8, 0.1, 0.1) )
#' weight <- rep(1, n = 1000)
#' df <- data.frame(y, A, X1, X2, weight)
#'
#'max_lambda(stand_formula = y ~ A + X1 + X2,
#'          exposure = "A",
#'           df = df,
#'           weight = weight)
#' @seealso
#'  \code{\link[nnet]{multinom}}
#' @rdname max_lambda
#' @export
#' @importFrom nnet multinom
max_lambda <- function(stand_formula,
                       exposure,
                       df,
                       weight = 1)
{
  fit_y_by_exposure_X <- nnet::multinom(stand_formula,
                                        df,
                                        trace = FALSE,
                                        weights = weight)
  df_treat <- df_untr <- df
  df_treat[,exposure] <- 1
  df_untr[,exposure] <- 0

  pred_treat <- as.data.frame(predict(fit_y_by_exposure_X, newdata = df_treat, type = "probs"))
  colnames(pred_treat) <-c("0","1", "2")
  pred_untr <- as.data.frame(predict(fit_y_by_exposure_X, newdata = df_untr, type = "probs"))
  colnames(pred_untr) <-c("0","1", "2")

  lambda1_ratio <- sum(df[,weight]*pred_treat[,"1"])/sum(df[,weight]*pred_untr[,"2"])
  lambda2_ratio <- sum(df[,weight]*pred_treat[,"2"])/sum(df[,weight]*pred_untr[,"1"])

  return(c("lambda1" = min(1, lambda1_ratio), "lambda2" = min(1, lambda2_ratio)))
}


