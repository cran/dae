\name{is.allzero}
\alias{is.allzero}
\title{Tests whether all elements are approximately zero}
\description{A single-line \code{function} that tests whether all elements are zero 
             (approximately).}
\usage{is.allzero(x)}
\arguments{
 \item{x}{An \code{object} whose elements are to be tested.}
}
\details{
All the elements of \code{x} are tested as being less than \code{daeTolerance}, 
which is initially set to 1e-10. The function \code{\link{set.daeTolerance}} can 
be used to change \code{daeTolerance}.}
\value{A \code{logical}.}
\author{Chris Brien}
\examples{
## create a vector of 9 zeroes and a one
y <- c(rep(0,9), 1)

## check that vector is only zeroes is FALSE 
is.allzero(y)
}
\keyword{manip}