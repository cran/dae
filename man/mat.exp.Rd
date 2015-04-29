\name{mat.exp}
\alias{mat.exp}
\title{Forms an exponential correlation matrix}
\description{Form the correlation \code{\link{matrix}} of order equal to the length of 
             \code{coordinates}. The \code{\link{matrix}} has diagonal 
             elements equal to one and the off-diagonal element in the ith row 
             and jth column equal to \eqn{\rho^k}{rho^(k)} where 
             \eqn{k = |coordinate[i]- coordinate[j]|}{k = abs(coordinate[i]- coordinate[j])}.}
\usage{mat.exp(coordinates, rho)}
\arguments{
 \item{coordinates}{The coordinates of points whose correlation \code{\link{matrix}} 
                    is to be formed.}
 \item{rho}{The correlation for points a distance of one apart.}
}
\value{A correlation \code{\link{matrix}} whose elements depend on the power of the 
       absolute distance apart.}
\seealso{\code{\link{mat.I}}, \code{\link{mat.J}}, \code{\link{mat.ar1}}, \code{\link{mat.banded}}}
\examples{
    corr <- mat.exp(coordinates=c(3:6, 9:12, 15:18), rho=0.1)
}
\keyword{array}
