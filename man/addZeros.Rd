\name{addZeros}
\alias{addZeros}

\title{
Add zeros to a numeric vector when transforming to character. 
}
\description{
Makes simple addition of zeros in front of numbers to ensure proper ordering 
of labels based on numbers. 
}
\usage{
addZeros(x)
}

\arguments{
  \item{x}{
Numeric vector.
}

}
\details{
Makes simple addition of zeros in front of numbers to ensure proper ordering 
of labels based on numbers. 
}
\value{

A new vector of type character.

\describe{
\item{$xNew}{the new vector with added zeros in the front}
}
}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
}

\examples{
x <- c(1,10,100)
addZeros(x)
}
\seealso{\code{\link{obs}} -- the core package }
