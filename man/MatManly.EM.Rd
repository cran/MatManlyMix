\name{MatManly.EM}
\alias{MatManly.EM}
% Also NEED an '\alias' for EACH other topic documented here.
\title{ EM algorithm for matrix clustering}
\description{
Runs the EM algorithm for matrix clustering}
\usage{
MatManly.EM(Y, X = NULL, initial = NULL, id = NULL, la = NULL, nu = NULL, tau = NULL, 
Mu = NULL, beta = NULL, Sigma = NULL, Psi = NULL, Mu.type = 0, Psi.type = 0, 
tol = 1e-05, max.iter = 1000, size.control = 0, silent = TRUE)
}
% maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y }{ dataset of random matrices (p x T x n), n random matrices of dimensionality (p x T)}
  \item{X }{ dataset of explanatory variables (T x q x n), q explanatory variables for modeling Y}
  \item{initial }{ initialization parameters provided by function MatManly.init()}
  \item{id }{ initial membership vector}
  \item{la }{ initial transformation parameters (K x p)}
  \item{nu }{ initial transformation parameters (K x T)}
  \item{tau }{ initial mixing proportions (length K)}
  \item{Mu }{ initial mean matrices (p x T x K)}
  \item{beta }{ initial coefficient matrices (q x p x K)}
  \item{Sigma }{ initial array of sigma (p x p x K)}
  \item{Psi }{ initial array of psi (T x T x K)}
  \item{Mu.type }{ mean structure: 0-unrestricted, 1-additive}
  \item{Psi.type }{ covariance structure of the Psi matrix: 0-unrestricted, 1-AR1}
  \item{tol }{ tolerance level}
  \item{max.iter }{ maximum number of iterations}
  \item{size.control }{ minimum size of clusters allowed for controlling spurious solutions}
  \item{silent }{ whether to produce output of steps or not}
}
\details{
Runs the EM algorithm for modeling and clustering matrices for a provided dataset. Both matrix Gaussian mixture and matrix Manly mixture with given explanatory variables (data matrix X) or without explanatory variables (X is null) can be employed.  
A user has three options to initialize the EM algorithm. The user can use the MatManly.init() function to get initial parameters and input them as 'initial'. The second  choice is to specify either initial id vector 'id' and transformation parameters 'la'. The third option is to input initial mode parameters 'la', 'tau', 'Mu', and 'Sigma' and 'Psi'. In the case when transformation parameters are not provided, the function runs the EM algorithm without any transformations, i.e., it is equivalent to the EM algorithm for a Gaussian mixture model. If some transformation parameters have to be excluded from the consideration, in the corresponding positions of matrix 'la', the user has to specify value 0.
A user also has three options to specify the covariance structure of the 'Psi' matrices, including unrestricted case, spherical matrices and autoregressive structured matrices. 
Notation: n - sample size, p x T - dimensionality of the random matrices, K - number of mixture components.
}
\value{
  \item{la }{ matrix of the skewness parameters (K x p)}
  \item{nu }{ matrix of the skewness parameters (K x T)}
  \item{tau }{ vector of mixing proportions (length K)}
  \item{Mu }{ matrix of the estimated mean matrices (p x T x K)}
  \item{beta }{ matrix of the coefficient parameters (q x p x K)}
  \item{Sigma }{ array of the estimated sigma (p x p x K)}
  \item{Psi }{ array of the estimated psi (T x T x K)}
  \item{Mu.type }{ mean structure: 0-unrestricted, 1-additive}
  \item{Psi.type }{ covariance structure of the Psi matrix: 0-unrestricted, 1-AR1}
  \item{gamma }{ matrix of the estimated posterior probabilities (n x K)}
  \item{id }{ estimated membership vector (length n)}
  \item{ll }{ log likelihood value}
  \item{bic }{ Bayesian Information Criterion}
 \item{iter }{ number of EM iterations run}
 \item{flag }{ convergence flag (0 - success, 1 - failure)}
}

\examples{
set.seed(123)

data(crime)


Y <- crime$Y[c(2,7),,] / 1000

p <- dim(Y)[1]
T <- dim(Y)[2]
n <- dim(Y)[3]
K <-  2

#init <- MatManly.init(Y, K = K, la = matrix(0.1, K, p), nu = matrix(0.1, K, T))
#M1 <- MatManly.EM(Y, initial = init, max.iter = 1000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EM algorithm }% __ONLY ONE__ keyword per line
