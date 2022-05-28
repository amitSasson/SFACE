#' @title Subtype Free Average Causal Effect
#' @description A function to estimate the Subtype Free Average Causal Effect.
#' @param stand_formula A formula for standartization and DR, y ~ A + X, the outcome as a function of the exposure and covariates
#' @param iptw_formula  A formula for IPTW and DR, A ~ X, the exposure as a function of the covariates.
#' @param exposure The treatment/exposure vector pf length n. Must be encoded 1 for treated and 0 for untreated.
#' @param outcome The categorical outcome vector of length n.  Must be encoded 0 for disease-free, 1 for the first subtype and 2 for the second subtype.
#' @param df a data frame with columns for the outcome, expousre and covariates.
#' @param subtype Should the SF-ACE be estimated for subtype 1 or subtype 2
#' @param scale Should the SF-ACE be estimated on the difference or risk ratio scale.
#' @param method Which method to use when adjusting for covariates, possibilities include standardization ("stand"), Inverse Probability Treatment Weighting ("IPTW"), and doubly robust estimation ("DR")
#' @param lambda1 sensitivity parameter for subtype 1. Can range between 0 (S-Monotonicity for subtype 1) and 1 (D-Monotonicity for subtype 1), Default: 0
#' @param lambda2 sensitivity parameter for subtype 2. Can range between 0 (S-Monotonicity for subtype 2) and 1 (D-Monotonicity for subtype 2), Default: 0
#' @param weight A numerical vector of length n, holding weights to adjust for missing subtypes, Default: 1
#' @param MultPer A numeric value indicating per how many people the effect should be calculated on the difference scale, Default: 1
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' A <- rbinom(n = 1000, size = 1, prob = 0.5)
#' X1 <- rbinom(n = 1000, size = 1, prob = 0.5)
#' X2 <- rnorm(n = 1000, mean = 0, sd = 1)
#' X <- matrix(c(X1,X2), nrow = 1000, byrow = FALSE)
#' y <- sample(c(0,1,2), 1000, replace=TRUE, prob=c(0.8, 0.1, 0.1) )
#' weight <- runif(n = 1000, 0,1)
#' df <- data.frame(y, A, X1, X2, weight)
#'
#' sface(stand_formula = y ~ A + X1 + X2,
#' iptw_formula = A ~ X1 + X2,
#' exposure = "A",
#' outcome = "y",
#' df = df,
#' scale = c("diff","RR"),
#' method = c("stand", "IPTW"),
#' weight = "weight")
#' @seealso
#'  \code{\link[nnet]{multinom}}
#' @rdname sface
#' @export
#' @importFrom nnet multinom
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
  if(!all(df$exposure %in% c(0,1))) {stop("exposure can only contain 0 and 1 values")}
  if(!all(df$y %in% c(0,1,2))) {stop("y can only contain 0, 1 and 2 values")}
  if(!all(subtype %in% c(1,2))) {stop("The subtype should be 1 or 2")}
  if(!all(scale %in% c("diff", "RR"))) {stop("The scale should be 'diff' or 'RR' ")}
  if(!all(method %in% c("stand", "IPTW", "DR"))) {stop("The scale should be 'stand','IPTW', or 'DR' ")}
  if(any(df$weight <= 0)) {stop("weights can't be negative")}
  if(any(lambda1 < 0) | any(lambda1 > 1)) {stop("Lambda 1 should be between 0 and 1")}
  if(any(lambda2 < 0) | any(lambda2 > 1)) {stop("Lambda 2 should be between 0 and 1")}
  if(MultPer <= 0 ) {stop("MultPer can't be negative")}

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
        estiamtors <- outer(X = lambda1,
                            Y = lambda2,
                            FUN = diff_calc,
                            p_Y11_exposure1,
                            p_Y21_exposure0,
                            p_Y11_exposure0,
                            p_Y20_exposure1,
                            MultPer)
        if(length(lambda1) > 1) {rownames(estiamtors) <- lambda1}
        if(length(lambda2) > 1) {colnames(estiamtors) <- lambda2}
        if(length(lambda1) == 1 & length(lambda2) == 1) {estiamtors <- c(estiamtors)}
        sface_list[["sface"]][["diff"]][["stand"]][[current_subtype]] <- estiamtors
      if(length(sface_list[["sface"]][["diff"]][["stand"]]) == 2)
      {
        sface_list[["sface"]][["diff"]][["stand"]][["theta"]] <- sface_list[["sface"]][["diff"]][["stand"]][["1"]]-sface_list[["sface"]][["diff"]][["stand"]][["2"]]
      }
      }

      if("RR" %in% scale )
      {
        estiamtors <- outer(X = lambda1,
                            Y = lambda2,
                            FUN = RR_calc,
                            p_Y11_exposure1,
                            p_Y21_exposure0,
                            p_Y11_exposure0)
        if(length(lambda1) > 1) {rownames(estiamtors) <- lambda1}
        if(length(lambda2) > 1) {colnames(estiamtors) <- lambda2}
        if(length(lambda1) == 1 & length(lambda2) == 1) {estiamtors <- c(estiamtors)}

        sface_list[["sface"]][["RR"]][["stand"]][[current_subtype]] <- estiamtors

        if(length(sface_list[["sface"]][["RR"]][["stand"]]) == 2)
        {
          sface_list[["sface"]][["RR"]][["stand"]][["theta"]] <- sface_list[["sface"]][["RR"]][["stand"]][["1"]]-sface_list[["sface"]][["RR"]][["stand"]][["2"]]
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
        estiamtors <- outer(X = lambda1,
                            Y = lambda2,
                            FUN = diff_calc,
                            p_Y11_exposure1,
                            p_Y21_exposure0,
                            p_Y11_exposure0,
                            p_Y20_exposure1,
                            MultPer)
        if(length(lambda1) > 1) {rownames(estiamtors) <- lambda1}
        if(length(lambda2) > 1) {colnames(estiamtors) <- lambda2}
        if(length(lambda1) == 1 & length(lambda2) == 1) {estiamtors <- c(estiamtors)}

        sface_list[["sface"]][["diff"]][["IPTW"]][[current_subtype]] <- estiamtors
        if(length(sface_list[["sface"]][["diff"]][["IPTW"]]) == 2)
        {
          sface_list[["sface"]][["diff"]][["IPTW"]][["theta"]] <- sface_list[["sface"]][["diff"]][["IPTW"]][["1"]]-sface_list[["sface"]][["diff"]][["IPTW"]][["2"]]
        }
      }

      if("RR" %in% scale )
      {
        estiamtors <- outer(X = lambda1,
                            Y = lambda2,
                            FUN = RR_calc,
                            p_Y11_exposure1,
                            p_Y21_exposure0,
                            p_Y11_exposure0)
        if(length(lambda1) > 1) {rownames(estiamtors) <- lambda1}
        if(length(lambda2) > 1) {colnames(estiamtors) <- lambda2}
        if(length(lambda1) == 1 & length(lambda2) == 1) {estiamtors <- c(estiamtors)}

        sface_list[["sface"]][["RR"]][["IPTW"]][[current_subtype]] <- estiamtors
        if(length(sface_list[["sface"]][["RR"]][["IPTW"]]) == 2)
        {
          sface_list[["sface"]][["RR"]][["IPTW"]][["theta"]] <- sface_list[["sface"]][["RR"]][["IPTW"]][["1"]]-sface_list[["sface"]][["RR"]][["IPTW"]][["2"]]
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
        estiamtors <- outer(X = lambda1,
                            Y = lambda2,
                            FUN = diff_calc,
                            p_Y11_exposure1,
                            p_Y21_exposure0,
                            p_Y11_exposure0,
                            p_Y20_exposure1,
                            MultPer)
        if(length(lambda1) > 1) {rownames(estiamtors) <- lambda1}
        if(length(lambda2) > 1) {colnames(estiamtors) <- lambda2}
        if(length(lambda1) == 1 & length(lambda2) == 1) {estiamtors <- c(estiamtors)}

        sface_list[["sface"]][["diff"]][["DR"]][[current_subtype]] <- estiamtors
        if(length(sface_list[["sface"]][["diff"]][["DR"]]) == 2)
        {
          sface_list[["sface"]][["diff"]][["DR"]][["theta"]] <- sface_list[["sface"]][["diff"]][["DR"]][["1"]]-sface_list[["sface"]][["diff"]][["DR"]][["2"]]
        }
      }

      if("RR" %in% scale )
      {
        estiamtors <- outer(X = lambda1,
                            Y = lambda2,
                            FUN = RR_calc,
                            p_Y11_exposure1,
                            p_Y21_exposure0,
                            p_Y11_exposure0)
        if(length(lambda1) > 1) {rownames(estiamtors) <- lambda1}
        if(length(lambda2) > 1) {colnames(estiamtors) <- lambda2}
        if(length(lambda1) == 1 & length(lambda2) == 1) {estiamtors <- c(estiamtors)}
        sface_list[["sface"]][["RR"]][["DR"]][[current_subtype]] <- estiamtors
        if(length(sface_list[["sface"]][["RR"]][["DR"]]) == 2)
        {
          sface_list[["sface"]][["RR"]][["DR"]][["theta"]] <- sface_list[["sface"]][["RR"]][["DR"]][["1"]]-sface_list[["sface"]][["RR"]][["DR"]][["2"]]
        }
      }
    }
  }
  sface_list[["additional_info"]][["lambda1"]] <- lambda1
  sface_list[["additional_info"]][["lambda2"]] <- lambda2
  sface_list[["additional_info"]][["subtype"]] <- subtype


  class(sface_list) <- "sface"
  return(sface_list)
}

diff_calc <- function(lambda1, lambda2, p_Y11_exposure1,p_Y21_exposure0, p_Y11_exposure0, p_Y20_exposure1, MultPer)
{
  MultPer*(p_Y11_exposure1-lambda2*p_Y21_exposure0+(lambda1-1)*p_Y11_exposure0)/(p_Y20_exposure1-lambda2*p_Y21_exposure0)
}

RR_calc <- function(lambda1, lambda2, p_Y11_exposure1,p_Y21_exposure0, p_Y11_exposure0)
{
  (p_Y11_exposure1-lambda2*p_Y21_exposure0)/((1-lambda1)*p_Y11_exposure0)
}


#' @title Print SF-ACE
#' @description this function prints lists of the class "sface"
#' @param x a list of class "sface", usually the output of the function sface
#' @param \dots not used
#' @return
#' @details
#' @examples
#' A <- rbinom(n = 1000, size = 1, prob = 0.5)
#' X1 <- rbinom(n = 1000, size = 1, prob = 0.5)
#' X2 <- rnorm(n = 1000, mean = 0, sd = 1)
#' X <- matrix(c(X1,X2), nrow = 1000, byrow = FALSE)
#' y <- sample(c(0,1,2), 1000, replace=TRUE, prob=c(0.8, 0.1, 0.1) )
#' weight <- runif(n = 1000, 0,1)
#' df <- data.frame(y, A, X1, X2, weight)
#'
#' lst <- sface(stand_formula = y ~ A + X1 + X2,
#' iptw_formula = A ~ X1 + X2,
#' exposure = "A",
#' outcome = "y",
#' df = df,
#' weight = "weight",
#' lambda1 = c(0.3, 0.5, 0.7))
#'
#' print(lst)
#' @rdname print.sface
#' @export print.sface
#' @export
print.sface <- function(x, ...)
{
  sface_list <- x
  lambda1 <- sface_list[["additional_info"]][["lambda1"]]
  lambda2 <- sface_list[["additional_info"]][["lambda2"]]
  subtype <- paste0("subtype", sface_list[["additional_info"]][["subtype"]])
  if (length(subtype) == 2) {subtype[3] <- "theta"}

  cat("The estimates SF-ACEs are:","\n","\n")
  for (sc in names(sface_list[["sface"]]))
  {
    if(sc == "diff") {cat("On the difference scale:","\n")}
    if(sc == "RR") {cat("On the RR scale:","\n")}

    for(m in names(sface_list[["sface"]][[1]]))
    {
      cat("Using ", as.character(m),",", "\n")
      if(length(lambda1) == 1 & length(lambda2) == 1)
      {
        ans <- do.call(cbind.data.frame, sface_list[["sface"]][[sc]][[m]])
        colnames(ans) <- subtype
        print(ans)
        cat("\n")
      }
      else
      {
        for(su in names(sface_list[["sface"]][[1]][[1]]))
        {
          if(su %in% 1:2) {cat("For Subtype", as.character(su),",", "\n")}
          if(su == "theta") {cat("For", as.character(su),",", "\n")}
        ans <- sface_list[["sface"]][[sc]][[m]][[su]]
        rownames(ans) <- paste0("lambda1=",as.character(lambda1))
        colnames(ans) <- paste0("lambda2=",as.character(lambda2))
        print(ans)
        cat("\n")
        }
      }
      cat("\n")
    }
    cat("\n")
  }
}

#' @title Plot  SF-ACE
#' @description this function plots lists of the class "sface"
#' @param x a list of class "sface", usually the output of the function sface
#' @param \dots not used
#' @details
#' @examples
#' A <- rbinom(n = 1000, size = 1, prob = 0.5)
#' X1 <- rbinom(n = 1000, size = 1, prob = 0.5)
#' X2 <- rnorm(n = 1000, mean = 0, sd = 1)
#' X <- matrix(c(X1,X2), nrow = 1000, byrow = FALSE)
#' y <- sample(c(0,1,2), 1000, replace=TRUE, prob=c(0.8, 0.1, 0.1) )
#' weight <- runif(n = 1000, 0,1)
#' df <- data.frame(y, A, X1, X2, weight)
#'
#' lst <- sface(stand_formula = y ~ A + X1 + X2,
#' iptw_formula = A ~ X1 + X2,
#' exposure = "A",
#' outcome = "y",
#' df = df,
#' weight = "weight",
#' lambda1 = c(0.3, 0.5, 0.7))
#'
#' plot(lst)
#' @rdname plot.sface
#' @importFrom gridExtra grid.arrange
#' @importFrom dplyr group_by
#' @importFrom tibble rownames_to_column
#' @importFrom purrr map2
#' @importFrom tidyr nest
#' @import ggplot2
#' @export plot.sface
#' @export
plot.sface <- function(x, ...)
{
  sface_list <- x
  lambda1 <- sface_list[["additional_info"]][["lambda1"]]
  lambda2 <- sface_list[["additional_info"]][["lambda2"]]
  subtype <- paste0("subtype", sface_list[["additional_info"]][["subtype"]])
  if (length(subtype) == 2) {subtype[3] <- "theta"}

  if(length(lambda1) > 1 & length(lambda2) == 1)
  {
    plot_one_lambda(sface_list, lambda = "lambda1", lambda_vals = lambda1)
  }

  if(length(lambda1) == 1 & length(lambda2) > 1)
  {
    plot_one_lambda(sface_list, "lambda2", lambda2)
  }

  if(length(lambda1) > 1 & length(lambda2) > 1)
  {
    plot_two_lambdas(sface_list, lambda1, lambda2)
  }

  if(length(lambda1) == 1 & length(lambda2) == 1)
  {
    print("plotting method requieres multiple lambda1 or lambda2 values")
  }
}

plot_one_lambda <- function(sface_list, lambda, lambda_vals)
{
  full_ans <-list()
  i <- 1

  for (sc in names(sface_list[["sface"]]))
  {
    for(m in names(sface_list[["sface"]][[1]]))
    {
      for(su in names(sface_list[["sface"]][[1]][[1]]))
      {
        ans <- sface_list[["sface"]][[sc]][[m]][[su]]
        if(lambda == "lambda2") {ans <- t(ans)}
        rownames(ans) <- lambda_vals
        colnames(ans) <- "value"
        ans <- tibble::rownames_to_column(as.data.frame(ans), "lambda")
        ans$method <- m
        ans$scale <- sc
        ans$subtype <- su
        full_ans[[i]] <- ans
        i <- i + 1
      }
    }
  }
  full_ans <- do.call(rbind, full_ans)
  p <- ggplot(full_ans, aes(x = lambda, y = value, color = subtype, group = subtype)) +
    facet_wrap(scale ~ method, scales = "free") +
    geom_point(alpha=0.6) +
    geom_line() +
    ylab("SF-ACE") +
    theme_bw()

  print(p)
}

plot_two_lambdas <- function(sface_list, lambda1_vals, lambda2_vals)
{
  full_ans <-list()
  i <- 1

  for (sc in names(sface_list[["sface"]]))
  {
    for(m in names(sface_list[["sface"]][[1]]))
    {
      for(su in names(sface_list[["sface"]][[1]][[1]]))
      {
        ans <- sface_list[["sface"]][[sc]][[m]][[su]]
        rownames(ans) <- lambda1_vals
        colnames(ans) <- lambda2_vals
        ans <- as.data.frame(as.table(ans))
        colnames(ans) <- c("lambda1", "lambda2", "value")
        ans$method <- m
        ans$scale <- sc
        ans$subtype <- su
        full_ans[[i]] <- ans
        i <- i + 1
      }
    }
  }
  full_ans <- do.call(rbind, full_ans)

    plot_func <- function(df, name)
  {
   ggplot(data = df, aes(x = lambda1, y = lambda2, fill = value)) +
      geom_tile() +
      scale_fill_continuous(name = name) +
      facet_wrap(.~ method)
  }

  nested_tmp <- mutate(tidyr::nest(dplyr::group_by(full_ans, scale)), plots = purrr::map2(data, scale, plot_func))
  gridExtra::grid.arrange(grobs = nested_tmp$plots)
}



