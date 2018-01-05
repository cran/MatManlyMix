%  file MatManlyMix/man/IMDb.Rd  
%  This file is a component of the package 'MatManlyMix' for R
%---------------------
\name{Satellite}
\alias{Satellite}
\docType{data}
\encoding{UTF-8}
\title{Satellite data}
\description{Data publicly available at the University of California - Irvine machine learning repository (http://archive.ics.uci.edu/ml), was originally obtained by NASA.
}
\usage{data(IMDb)}
\format{ 
  A list of 2 objects: Y and id, where Y represents the data array of spectral values and id represents the true id of three classes: Soil with vegetation stubble, damp grey soil, and grey soil. Y is the of dimensionality 4 x 9 x 845. 
}
\details{The data are publicly available on http://archive.ics.uci.edu/ml.}

\examples{

data(Satellite)

}
\keyword{datasets}
