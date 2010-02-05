\name{as.numfac}
\alias{as.numfac}
\title{Convert a factor to a numeric vector}
\description{Converts a \code{factor} to a numeric \code{vector} with integers 1, 2, ...
             for the levels that correspond approximately to its original 
             numeric values. It uses the method described in 
             \code{\link{factor}}.}
\usage{as.numfac(factor)}
\arguments{
 \item{factor}{The \code{factor} to be converted.}
}
\value{A numeric \code{vector}.}
\section{Warning}{Labels with non-numeric characters will convert to \code{NA}. 
                  \code{\link{fac.recode}} can be used to recode labels to 
                  numeric characters.}
\author{Chris Brien}
\seealso{\code{\link{fac.recode}} in package \pkg{dae}, \code{\link{factor}}.}
\examples{
## set up a factor and convert it to a numeric vector
a <- factor(rep(1:3, 4))
x <- as.numfac(a)
}
\keyword{factor}
\keyword{manip}