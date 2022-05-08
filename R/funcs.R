
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param y The categorical outcome vector of length n.  Must be encoded o for disease-free, 1 for the first subtype and 2 for the second subtype.
#' @param A The treatment/expousre vector pf length n. Must be encoded 1 for treated and 0 for untreated.
#' @param X The n × p-matrix of covariates X, Default: NULL
#' @param subtype Indication of which subtype to estimate the SF-ACE of, Default: c(1, 2)
#' @param scale Indication of weather the SF-ACE should be estimated on the difference or risk ratio scale, Default: c("diff", "RR")
#' @param method Indication of which method to use when adjusting for covariates, Default: c("stand", "IPTW", "DR")
#' @param weights A vector of length n, holding weights to adjust for missing subtypes, Default: NULL
#' @param MultPer A numeric value indicating per how many people the effect sould be calculated on the difference scale, Default: 1
#' @return
#' @details DETAILS
#' @examples
#' A <- rbinom(n = 1000, size = 1, prob = 0.5)
#' X1 <- rbinom(n = 1000, size = 1, prob = 0.5)
#' X2 <- rnorm(n = 1000, mean = 0, sd = 1)
#' X <- matrix(c(X1,X2), nrow = 1000, byrow = FALSE)
#' y <- sample(c(0,1,2), 1000, replace=TRUE, prob=c(0.8, 0.1, 0.1) )
#' sface(y, A, X, subtype = 1, scale = "diff", method = "stand")
#'
#'
#' @seealso
#'  \code{\link[nnet]{multinom}}
#' @rdname sface
#' @export
#' @importFrom nnet multinom
sface <- function(y,
                  A,
                  X = NULL,
                  subtype = c(1,2),
                  scale = c("diff", "RR"),
                  method = c("stand", "IPTW", "DR"),
                  weights = NULL,
                  MultPer=1)
{
  model <- nnet::multinom(y~A+X,
                          trace = FALSE,
                          weights = weights)
  return(summary(model))
}
