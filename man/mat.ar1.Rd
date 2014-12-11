\name{mat.ar1}
\alias{mat.ar1}
\title{Forms an ar1 correlation matrix}
\description{Form the correlation \code{\link{matrix}} of order \code{order} whose 
             correlations follow the ar1 pattern. The \code{\link{matrix}} has diagonal 
             elements equal to one and the off-diagonal element in the ith row 
             and jth column equal to \eqn{\rho^k}{rho^(k)} where 
             \eqn{k = |i- j|}{k = abs(i- j)}.}
\usage{mat.ar1(order, rho)}
\arguments{
 \item{order}{The order of the \code{\link{matrix}} to be formed.}
 \item{rho}{The correlation on the first off-diagonal.}
}
\value{A correlation \code{\link{matrix}} whose elements follow an ar1 pattern.}
\seealso{\code{\link{mat.I}}, \code{\link{mat.J}}}
\examples{
    corr <- mat.ar1(order=4, rho=0.4)
}
\keyword{array}
