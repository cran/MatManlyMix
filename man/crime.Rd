%  file MatManlyMix/man/cime.Rd  
%  This file is a component of the package 'MatManlyMix' for R
%---------------------
\name{crime}
\alias{crime}
\docType{data}
\encoding{UTF-8}
\title{Crime data}
\description{Data collected by FBI's Uniform Crime on the violent and property crimes of 236 cities.
}
\usage{data(crime)}
\format{ 
  A data frame with 236 observations on the following 10 variables from year 2000 through 2012. The data has been cleaned with no missing values. 
  \describe{
   \item{Population}{Population of each city;}
   \item{Violent Crime rate}{Total number of violent crimes;}
   \item{Murder and non-negligent manslaughter rate}{Number of murders;}
   \item{Forcible rape rate}{Number of rape crimes;}
   \item{Robbery rate}{Number of robberies;}
   \item{Aggravated assault rate}{Number of assaults;}
   \item{Property crime rate}{Total number of property crimes;}
   \item{Burglary rate}{Number of burglary crimes;}
   \item{Larceny-theft rate}{Number of theft crimes;}
   \item{Motor vehicle theft rate}{Number of vehicle theft crimes;}
  }
}
\details{The data have been made publicly available by FBI's Uniform Crime Reports.}

\examples{

data(crime)

}
\keyword{datasets}
