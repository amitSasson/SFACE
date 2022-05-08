
#' @title sface
#' @description function to estimate the Subtype Free Average Causal Effect.
#' @param y The categorical outcome vector of length n.  Must be encoded o for disease-free, 1 for the first subtype and 2 for the second subtype.
#' @param A The treatment/expousre vector pf length n. Must be encoded 1 for treated and 0 for untreated.
#' @param X The n × p-matrix of covariates X, Default: NULL
#' @param subtype Indication of which subtype to estimate the SF-ACE of, Default: c(1, 2)
#' @param scale Indication of weather the SF-ACE should be estimated on the difference or risk ratio scale, Default: c("diff", "RR")
#' @param method Indication of which method to use when adjusting for covariates, Default: c("stand", "IPTW", "DR")
#' @param weight A vector of length n, holding weights to adjust for missing subtypes, Default: 1
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
                  weight = 1,
                  MultPer=1)
{
  if(method == "stand")
  {
    df <- data.frame(y, A, X, weight)
    fit_y_by_A_X <- nnet::multinom(y ~ A + X,
                                   df,
                                   trace = FALSE,
                                   weights = weight)

    df_treat <- df_untr <- df
    df_treat$A <- 1
    df_untr$A <- 0

    pred_treat <- as.data.frame(predict(fit_y_by_A_X, newdata = df_treat, type = "probs"))
    colnames(pred_treat) <-c("0","1", "2")
    pred_untr <- as.data.frame(predict(fit_y_by_A_X, newdata = df_untr, type = "probs"))
    colnames(pred_untr) <-c("0","1", "2")

    self <- ifelse(subtype == 1, "1", "2")
    other <- ifelse(subtype == 1, "2", "1")

    p_Y11_A1 <- sum(df$weight*pred_treat[,self])/sum(df$weight)
    p_Y11_A0 <- sum(df$weight*pred_untr[,self])/sum(df$weight)

    if(scale == "diff")
    {
      p_Y12_A1 <- sum(df$weight*pred_treat[,other])/sum(df$weight)
      return(MultPer*(p_Y11_A1 - p_Y11_A0)/(1-p_Y12_A1))
    }

    if(scale == "RR")
    {
      return(p_Y11_A1/p_Y11_A0)
    }
  }
  if(method == "IPTW")
  {
    df <- data.frame(y, A, X, weight)
    df$'1' <- ifelse(y == 1 ,1, 0)
    df$'2' <- ifelse(y == 2, 1, 0)

    fit_A_by_X <- glm(A ~ X,
                      df,
                      family = "binomial",
                      weights = weight)
    pred_A <-  predict(fit_A_by_X, type = "response")
    pr_A_1 <- mean(df$A)
    n <- nrow(df)

    self <- ifelse(subtype == 1, "1", "2")
    other <- ifelse(subtype == 1, "2", "1")

    df$w_A <- ifelse(df$A == 1, pr_A_1/pred_A, (1-pr_A_1)/(1-pred_A) ) #Stabilized weights
    p_Y11_A1 <- sum(df$weight*df$w_A*df$A*df[,self])/sum(df$weight*df$A)
    p_Y11_A0 <- sum(df$weight*df$w_A*(1-df$A)*df[,self])/sum(df$weight*(1-df$A))

    if(scale == "diff")
    {
      p_Y12_A1 <- sum(df$weight*df$w_A*df$A*df[,other])/sum(df$weight*df$A)
      return(MultPer*(p_Y11_A1 - p_Y11_A0)/(1-p_Y12_A1))
    }

    if(scale == "RR")
    {
      return(p_Y11_A1/p_Y11_A0)
    }
  }

}
