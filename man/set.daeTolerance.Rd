\name{set.daeTolerance}
\alias{set.daeTolerance}
\title{Sets the value of daeTolerance for the package dae}
\description{A function that sets the value such that, in \pkg{dae} functions, 
values less than it are considered to be zero. Initially, \code{daeTolerance} 
is set to 1e-10.}
\usage{set.daeTolerance(tolerance = daeTolerance)}
\arguments{
 \item{tolerance}{The value to which \code{daeTolerance} is to be set.}
}
\value{The value of {daeTolerance} is returned invisibly.}
\author{Chris Brien}
\examples{
## set daeTolerance.
set.daeTolerance(.Machine$double.eps ^ 0.5)
}
\keyword{manip}
\keyword{projector}