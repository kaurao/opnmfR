% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opnmfR.R
\name{opnmfR_test}
\alias{opnmfR_test}
\title{opnmf test}
\usage{
opnmfR_test(X = NULL, y = NULL, r = 2, W0 = "nndsvd", ...)
}
\arguments{
\item{X}{A matrix, if NULL the "iris" data is used (default NULL)}

\item{y}{A vector with class labels same length as nrow of X, if X=NULL then iris$Species is used (default NULL)}

\item{r}{A number, rank to use (default 2)}

\item{W0}{A string or matrix for initialization (default "nndsvd")}

\item{...}{additional parameters passed to \code{opnmfR} and \code{opnmfRcpp}}
}
\value{
A list with factorization results using R and Rcpp function calls
}
\description{
opnmf test
}
\examples{
result <- opnmfR_test()
}
