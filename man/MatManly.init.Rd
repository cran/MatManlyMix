\name{MatManly.init}
\alias{MatManly.init}
% Also NEED an '\alias' for EACH other topic documented here.
\title{ Initialization for the EM algorithm for matrix clustering}
\description{
Runs the initialization for the EM algorithm for matrix clustering}
\usage{
MatManly.init(Y, X = NULL, K, la = NULL, nu = NULL, Mu.type = 0, 
Psi.type = 0, n.start = 10, tol = 1e-05)

}
% maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y }{ dataset of random matrices (p x T x n), n random matrices of dimensionality (p x T)}
  \item{X }{ dataset of explanatory variables (T x q x n), q explanatory variables for modeling Y}
  \item{K }{ number of clusters}
  \item{la }{ initial transformation parameters (K x p)}
  \item{nu }{ initial transformation parameters (K x T)}
  \item{Mu.type }{ mean structure: 0-unrestricted, 1-additive}
  \item{Psi.type }{ covariance structure of the Psi matrix: 0-unrestricted, 1-AR1}
  \item{n.start }{ initial random starts}
  \item{tol }{ tolerance level}
}
\details{
Random starts are used to obtain different starting values. The number of clusters, the skewness parameters, and number of random starts need to be specified. In the case when transformation parameters are not provided, the function runs the EM algorithm without any transformations, i.e., it is equivalent to the EM algorithm for a matrix Gaussian mixture.
Notation: n - sample size, p x T - dimensionality of the random matrices, K - number of mixture components. 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Initialization }% __ONLY ONE__ keyword per line
