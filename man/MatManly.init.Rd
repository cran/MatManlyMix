\name{MatManly.init}
\alias{MatManly.init}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Initialization for the EM algorithm for matrix clustering}
\description{
Runs the initialization for the EM algorithm for matrix clustering}
\usage{
MatManly.init(Y, X = NULL, K, la = NULL, Psi.type = 0, n.start = 10, tol = 1e-05)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y }{ dataset matrix (p x T x n) predictor in linear model}
  \item{X }{ dataset matrix (q x p x n) independent variable in linear model}
  \item{K }{ number of clusters}
  \item{la }{ initial transformation parameters (K x p)}
  \item{Psi.type }{ covariance structure of the Psi matrix: 0-unrestricted, 1-diagonal, 2-AR1}
  \item{n.start }{ initial random starts}
  \item{tol }{ tolerance level}
}
\details{
Random starts are used to obtain different starting values. The number of clusters, the skewness parameters, and number of random starts need to be specified. In the case when transformation parameters are not provided, the function runs the EM algorithm without any transformations, i.e., it is equivalent to the EM algorithm for a Gaussian mixture model.
Notation: n - sample size, p x T - dimensionality of the dataset Y, K - number of mixture components. 
}

\examples{

set.seed(123)

data(crime)

Y <- crime[c(2,7),,] / 1000

p <- dim(Y)[1]
T <- dim(Y)[2]
n <- dim(Y)[3]
K <-  2

init1 <- MatManly.init(Y, K = K, la = matrix(0.1, K, p), n.start = 1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Initialization }% __ONLY ONE__ keyword per line
