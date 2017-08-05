\name{MatManly.EM}
\alias{MatManly.EM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ EM algorithm for matrix clustering}
\description{
Runs the EM algorithm for matrix clustering}
\usage{
MatManly.EM(Y, X = NULL, initial = NULL, id = NULL, la = NULL, tau = NULL, 
Mu = NULL, beta = NULL, Sigma = NULL, Psi = NULL, Psi.type = 0, 
tol = 1e-05, max.iter = 1000, size.control = 0, silent = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y }{ dataset matrix (p x T x n) predictor in linear model}
  \item{X }{ dataset matrix (q x p x n) independent variable in linear model}
  \item{initial }{ initialization parameters provided by function MatManly.init}
  \item{id }{ initial membership vector}
  \item{la }{ initial transformation parameters (K x p)}
  \item{tau }{ initial mixing proportions (length K)}
  \item{Mu }{ initial mean matrices (p x T x K)}
  \item{beta }{ initial coefficient matrices (q x p x K)}
  \item{Sigma }{ initial array of sigma (p x p x K)}
  \item{Psi }{ initial array of psi (T x T x K)}
  \item{Psi.type }{ covariance structure of the Psi matrix: 0-unrestricted, 1-diagonal, 2-AR1}
  \item{tol }{ tolerance level}
  \item{max.iter }{ maximum number of iterations}
  \item{size.control }{ minimum size of clusters allowed for controlling spurious solutions}
  \item{silent }{ whether to produce output of steps or not}
}
\details{
Runs the EM algorithm for clustering matrices for a provided dataset. Both Gaussian mixture model and Manly mixture models are employed.  
A user has three options to initialize the EM algorithm. The user can use the MatManly.init() function to get initial parameters and input them as 'initial'. The second  choice is to specify either initial id vector 'id' and transformation parameters 'la'. The third option is to input initial mode parameters 'la', 'tau', 'Mu', and 'Sigma' and 'Psi'. In the case when transformation parameters are not provided, the function runs the EM algorithm without any transformations, i.e., it is equivalent to the EM algorithm for a Gaussian mixture model. If some transformation parameters have to be excluded from the consideration, in the corresponding positions of matrix 'la', the user has to specify value 0.
Notation: n - sample size, p x T - dimensionality of the dataset Y, K - number of mixture components.

}
\value{
  \item{tau }{ vector of mixing proportions (length K)}
  \item{la }{ matrix of the skewness parameters (K x p)}
  \item{Mu }{ matrix of the estimated mean matrices (p x T x K)}
  \item{beta }{ matrix of the coefficient parameters (q x p x K)}
  \item{Sigma }{ array of the estimated sigma (p x p x K)}
  \item{Psi }{ array of the estimated psi (T x T x K)}
  \item{Psi.type }{ covariance structure of the Psi matrix: 0-unrestricted, 1-diagonal, 2-AR1}
  \item{gamma }{ matrix of the estimated posterior probabilities (n x K)}
  \item{id }{ estimated membership vector (length n)}
  \item{ll }{ log likelihood value}
  \item{bic }{ Bayesian Information Criterion}
}

\examples{
set.seed(123)

data(crime)


Y <- crime[c(2,7),,] / 1000

p <- dim(Y)[1]
T <- dim(Y)[2]
n <- dim(Y)[3]
K <-  2

init <- MatManly.init(Y, K = K, la = matrix(0.1, K, p), Psi.type = 2, n.start = 1)
M1 <- MatManly.EM(Y, initial = init, max.iter = 1000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EM algorithm }% __ONLY ONE__ keyword per line
