% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.sface.R
\name{plot.sface}
\alias{plot.sface}
\title{Plot  SF-ACE}
\usage{
\method{plot}{sface}(x, ...)
}
\arguments{
\item{x}{a list of class "sface", usually the output of the function sface}

\item{\dots}{not used}
}
\description{
this function plots lists of the class "sface"
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

plot(lst)
}
