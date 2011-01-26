\name{rmvnorm}
\alias{rmvnorm}
\title{generates a vector of random values from a multivariate normal distribution}
\description{Generates a vector of random values from an n-dimensional 
             multivariate normal distribution whose mean is given by the 
             n-vector \code{mean} and variance by the 
             n xn symmetric matrix \code{V}.}
\usage{rmvnorm(mean, V)}
\arguments{
 \item{mean}{The mean vector of the multivariate normal distribution from which 
               the random values are to be generated.}
 \item{V}{The variance matrix of the multivariate normal distribution from which 
               the random values are to be generated..}
}
\value{An \code{\link{vector}} of length n, equal to the length of \code{mean}.}
\details{
The method is:
 a) obtain the upper-triangular factor of the Cholesky decompostion of \code{V};
 b) generate a vector of length equal to \code{mean} of standard normal values;
 c) premultiply the vector of standard normal values by the transpose of the 
    upper triangular factor and, to the result, add \code{mean}.
}
\author{Chris Brien}
\seealso{\code{\link{fac.ar1mat}},  \code{\link{fac.vcmat}}, 
         in package \pkg{dae}, and \code{\link{rnorm}}.}
\examples{
## set up a two-level factor and a three-level factor, both of length 12
A <- factor(rep(1:2, each=6))
B <- factor(rep(1:3, each=2, times=2))

## generate random values from a multivariate normal for which 
#the mean is 20 for all variables and 
#the variance matrix has random effects for factor A, ar1 pattern for B and 
#residual random variation
mean <- rep(20, 12)
V <- fac.vcmat(A, 5) + fac.ar1mat(B, 0.6) + 2*mat.I(12)
y <- rmvnorm(mean, V)
}
\keyword{datagen}
