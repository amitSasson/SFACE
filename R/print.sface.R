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
print.sface <- function(x, digits = 4, ...)
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
        print(ans, digits = digits)
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
          print(ans, digits = digits)
          cat("\n")
        }
      }
      cat("\n")
    }
    cat("\n")
  }
}
