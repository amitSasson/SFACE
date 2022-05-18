
#' @title Subtype Free Average Causal Effect
#' @description A function to estimate the Subtype Free Average Causal Effect.
#' @param y The categorical outcome vector of length n.  Must be encoded 0 for disease-free, 1 for the first subtype and 2 for the second subtype.
#' @param A The treatment/expousre vector pf length n. Must be encoded 1 for treated and 0 for untreated.
#' @param X The n × p-matrix of covariates X, Default: NULL
#' @param subtype Should the SF-ACE be estimated for subtype 1 or subtype 2
#' @param scale Should the SF-ACE be estimated on the difference or risk ratio scale.
#' @param method Which method to use when adjusting for covariates, possibilities include standardization ("stand"), Inverse Probability Treatment Weighting ("IPTW"), and doubly robust estimation ("DR")
#' @param lambda1 sensitivity parameter for subtype 1. Can range between 0 (S-Monotonicity for subtype 1) and 1 (D-Monotonicity for subtype 1), Default: 0
#' @param lambda2 sensitivity parameter for subtype 2. Can range between 0 (S-Monotonicity for subtype 2) and 1 (D-Monotonicity for subtype 2), Default: 0
#' @param weight A numerical vector of length n, holding weights to adjust for missing subtypes, Default: 1
#' @param MultPer A numeric value indicating per how many people the effect should be calculated on the difference scale, Default: 1
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
sface(stand_formula = y ~ A + X1 + X2,
      iptw_formula = A ~ X1 + X2,
      exposure = "A",
      outcome = "y",
      df = df,
      subtype = c(1),
      scale = c("diff","RR"),
      method = c("stand", "IPTW"),
      weight = "weight")
sface(stand_formula = y ~ A + X1 + X2,
      iptw_formula = A ~ X1 + X2,
      exposure = "A",
      outcome = "y",
      df = df,
      subtype = 1,
      weight = "weight")
sface <- function(stand_formula,
                  iptw_formula,
                  exposure,
                  outcome,
                  df,
                  subtype = c(1,2),
                  scale = c("diff", "RR"),
                  method = c("stand", "IPTW", "DR"),
                  lambda1 = 0,
                  lambda2 = 0,
                  weight = 1,
                  MultPer=1)
{
  #checking the input:
  #if(length(exposure) != length(y)) {stop("exposure and y are of different lenghts")}
  #if(!all(exposure %in% c(0,1))) {stop("exposure can only contain 0 and 1 values")}
  #if(!all(y %in% c(0,1,2))) {stop("y can only contain 0, 1 and 2 values")}
  #if(!(subtype %in% c(1,2))) {stop("The subtype should be 1 or 2")}
  #if(!all(scale %in% c("diff", "RR"))) {stop("The scale should be 'diff' or 'RR' ")}
  #if(!all(method %in% c("stand", "IPTW", "DR"))) {stop("The scale should be 'stand','IPTW', or 'DR' ")}
  #if(any(weight <= 0)) {stop("weights can't be negative")}
  #if(lambda1 < 0 | lambda1 > 1) {stop("Lambda 1 should be between 0 and 1")}
  #if(lambda2 < 0 | lambda2 > 1) {stop("Lambda 2 should be between 0 and 1")}
  #if(MultPer <= 0 ) {stop("MultPer can't be negative")}

  sface_list <- list()
  subtype <- as.character(subtype)

  if (any(c("stand","DR") %in% method))
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

    n_w <- sum(df[,weight])
  }

  if (any(c("IPTW","DR") %in% method))
  {
    fit_exposure_by_X <- glm(iptw_formula,
                             df,
                             family = "binomial",
                             weights = weight)
    df$'1' <- ifelse(df[,outcome] == 1 ,1, 0)
    df$'2' <- ifelse(df[,outcome] == 2, 1, 0)

    pred_exposure <-  predict(fit_exposure_by_X, type = "response")

  }
  #calculate the expectations needed using stand
  if("stand" %in% method)
  {
    for (current_subtype in subtype)
    {
      self <- ifelse(current_subtype == 1, "1", "2")
      other <- ifelse(current_subtype == 1, "2", "1")

      p_Y11_exposure1 <- sum(df[,weight]*pred_treat[,self])/n_w
      p_Y21_exposure0 <- sum(df[,weight]*pred_untr[,other])/n_w
      p_Y11_exposure0 <- sum(df[,weight]*pred_untr[,self])/n_w

      if("diff" %in% scale )
      {
        p_Y20_exposure1 <- 1-sum(df[,weight]*pred_treat[,other])/n_w
        sface_list[["diff"]][["stand"]][[current_subtype]] <- MultPer*(p_Y11_exposure1-lambda2*p_Y21_exposure0+(lambda1-1)*p_Y11_exposure0)/(p_Y20_exposure1-lambda2*p_Y21_exposure0)
      if(length(sface_list[["diff"]][["stand"]]) == 2)
      {
        sface_list[["diff"]][["stand"]][["theta"]] <- sface_list[["diff"]][["stand"]][["1"]]-sface_list[["diff"]][["stand"]][["2"]]
      }
      }

      if("RR" %in% scale )
      {
        sface_list[["RR"]][["stand"]][[current_subtype]] <- (p_Y11_exposure1-lambda2*p_Y21_exposure0)/((1-lambda1)*p_Y11_exposure0)
        if(length(sface_list[["RR"]][["stand"]]) == 2)
        {
          sface_list[["RR"]][["stand"]][["theta"]] <- sface_list[["RR"]][["stand"]][["1"]]-sface_list[["RR"]][["stand"]][["2"]]
        }
      }
    }
  }

  #calculate the expectations needed using IPTW
  if("IPTW" %in% method)
  {
    pr_exposure_1 <- mean(df[,exposure])
    n <- nrow(df)

    df$w_exposure <- ifelse(df[,exposure] == 1, pr_exposure_1/pred_exposure, (1-pr_exposure_1)/(1-pred_exposure) ) #Stabilized weights

    #q99 <- quantile(df$w_exposure, .99)
    #df$w_exposure <- ifelse(df$w_exposure > q99, q99, df$w_exposure)

    n_w <- sum(df[,weight])

    for (current_subtype in subtype)
    {
      self <- ifelse(current_subtype == 1, "1", "2")
      other <- ifelse(current_subtype == 1, "2", "1")
      p_Y11_exposure1 <- sum(df[,weight]*df$w_exposure*df[,exposure]*df[,self])/sum(df[,weight]*df[,exposure])
      p_Y11_exposure0 <- sum(df[,weight]*df$w_exposure*(1-df[,exposure])*df[,self])/sum(df[,weight]*(1-df[,exposure]))
      p_Y21_exposure0 <- sum(df[,weight]*df$w_exposure*(1-df[,exposure])*df[,other])/sum(df[,weight]*(1-df[,exposure]))

      if("diff" %in% scale)
      {
        p_Y20_exposure1 <- 1- sum(df[,weight]*df$w_exposure*df[,exposure]*df[,other])/sum(df[,weight]*df[,exposure])
        sface_list[["diff"]][["IPTW"]][[current_subtype]] <- MultPer*(p_Y11_exposure1-lambda2*p_Y21_exposure0+(lambda1-1)*p_Y11_exposure0)/(p_Y20_exposure1-lambda2*p_Y21_exposure0)
        if(length(sface_list[["diff"]][["IPTW"]]) == 2)
        {
          sface_list[["diff"]][["IPTW"]][["theta"]] <- sface_list[["diff"]][["IPTW"]][["1"]]-sface_list[["diff"]][["IPTW"]][["2"]]
        }
      }

      if("RR" %in% scale )
      {
        sface_list[["RR"]][["IPTW"]][[current_subtype]] <- (p_Y11_exposure1-lambda2*p_Y21_exposure0)/((1-lambda1)*p_Y11_exposure0)
        if(length(sface_list[["RR"]][["IPTW"]]) == 2)
        {
          sface_list[["RR"]][["IPTW"]][["theta"]] <- sface_list[["RR"]][["IPTW"]][["1"]]-sface_list[["RR"]][["IPTW"]][["2"]]
        }
      }
    }
  }

  #calculate the expectations needed using DR
  if("DR" %in% method)
  {
    for (current_subtype in subtype)
    {
      self <- ifelse(current_subtype == 1, "1", "2")
      other <- ifelse(current_subtype == 1, "2", "1")
      p_Y11_exposure1 <- sum(df[,weight]*(df[,exposure]*df[,self]/(pred_exposure) - ((df[,exposure]-pred_exposure)*pred_treat[,self])/pred_exposure))/n_w
      p_Y11_exposure0 <- sum(df[,weight]*((1-df[,exposure])*df[,self]/(1-pred_exposure) + ((df[,exposure]-pred_exposure)*pred_untr[,self])/(1-pred_exposure)))/n_w
      p_Y21_exposure0 <- sum(df[,weight]*((1-df[,exposure])*df[,other]/(1-pred_exposure) + ((df[,exposure]-pred_exposure)*pred_untr[,other])/(1-pred_exposure)))/n_w

      if("diff" %in% scale )
      {
        p_Y20_exposure1 <- 1 - sum(df[,weight]*(df[,exposure]*df[,other]/(pred_exposure) - ((df[,exposure]-pred_exposure)*pred_treat[,other])/pred_exposure))/n_w
        sface_list[["diff"]][["DR"]][[current_subtype]] <- MultPer*(p_Y11_exposure1-lambda2*p_Y21_exposure0+(lambda1-1)*p_Y11_exposure0)/(p_Y20_exposure1-lambda2*p_Y21_exposure0)
        if(length(sface_list[["diff"]][["DR"]]) == 2)
        {
          sface_list[["diff"]][["DR"]][["theta"]] <- sface_list[["diff"]][["DR"]][["1"]]-sface_list[["diff"]][["DR"]][["2"]]
        }
      }

      if("RR" %in% scale )
      {
        sface_list[["RR"]][["DR"]][[current_subtype]] <- (p_Y11_exposure1-lambda2*p_Y21_exposure0)/((1-lambda1)*p_Y11_exposure0)
        if(length(sface_list[["RR"]][["DR"]]) == 2)
        {
          sface_list[["RR"]][["DR"]][["theta"]] <- sface_list[["RR"]][["DR"]][["1"]]-sface_list[["RR"]][["DR"]][["2"]]
        }
      }
    }
  }

  subtype <- paste("subtype ", subtype)
  if(length(subtype) == 2)
  {
    subtype[3] <- "theta"
  }

  cat("The estimates SF-ACEs are:","\n","\n")
  if ("diff" %in% scale)
  {
    cat("On the difference scale:","\n")
    diff_table <- do.call(rbind.data.frame, sface_list[["diff"]])
    colnames(diff_table) <- subtype
    print(diff_table)
    cat("\n")
  }
  if ("RR" %in% scale)
  {
    cat("On the Risk Ratio scale:","\n")
    RR_table <- do.call(rbind.data.frame, sface_list[["RR"]])
    colnames(RR_table) <- subtype
    print(RR_table)
  }
  return(invisible(sface_list))
}


