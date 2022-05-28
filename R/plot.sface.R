#' @title Plot  SF-ACE
#' @description this function plots lists of the class "sface"
#' @param x a list of class "sface", usually the output of the function sface
#' @param \dots not used
#' @details
#' @examples
#' A <- rbinom(n = 1000, size = 1, prob = 0.5)
#' X1 <- rbinom(n = 1000, size = 1, prob = 0.5)
#' X2 <- rnorm(n = 1000, mean = 0, sd = 1)
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
