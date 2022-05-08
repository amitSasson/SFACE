# this function returns X
#' @title self
#' @description this function returns X
#' @param x PARAM_DESCRIPTION
#' @return the param x
#' @details DETAILS
#' @examples self(2)
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname self
#' @export

self <- function(x)
{
  #and that's it
  return(x)
}

#' @title sface
#' @description FUNCTION_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[nnet]{multinom}}
#' @rdname SFACE
#' @export
#' @importFrom nnet multinom
sface <- function(df)
{
  model <- nnet::multinom(y~X, df)
  return(summary(model))
}
