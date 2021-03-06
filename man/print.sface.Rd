% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.sface.R
\name{print.sface}
\alias{print.sface}
\title{Print SF-ACE}
\usage{
\method{print}{sface}(x, digits = 4, ...)
}
\arguments{
\item{x}{a list of class "sface", usually the output of the function sface}

\item{digits}{number of digits chosen for the print. Default: 4 .}

\item{\dots}{not used}
}
\value{

}
\description{
this function prints lists of the class "sface"
}
\examples{
A <- rbinom(n = 1000, size = 1, prob = 0.5)
X1 <- rbinom(n = 1000, size = 1, prob = 0.5)
X2 <- rnorm(n = 1000, mean = 0, sd = 1)
y <- sample(c(0,1,2), 1000, replace=TRUE, prob=c(0.8, 0.1, 0.1) )
weight <- runif(n = 1000, 0,1)
df <- data.frame(y, A, X1, X2, weight)

lst <- sface(stand_formula = y ~ A + X1 + X2,
iptw_formula = A ~ X1 + X2,
exposure = "A",
outcome = "y",
df = df,
weight = "weight",
lambda1 = c(0.3, 0.5, 0.7))

print(lst)
}
