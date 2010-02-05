\name{fac.divide}
\alias{fac.divide}
\title{Divides a factor into several individual factors}
\description{Takes a \code{factor} and divides it into several individual 
 \code{factors} as if the \code{levels} in the combined \code{factor} correspond 
 to the numbering of the \code{levels} combinations of the individual 
 \code{factors} when these are arranged in standard or Yates order.}
\usage{fac.divide(combined.factor, factor.names, order="standard")}
\arguments{
 \item{combined.factor}{A \code{factor} that is to be divided into the individual 
       \code{factors} listed in \code{factor.names}.}
 \item{factor.names}{A \code{list} of \code{factors} to be formed. The names in 
       the \code{list} are the names of the \code{factors} and the component of 
       a name is either a) a single \code{numeric} value that is the number of 
       \code{levels}, b) a numeric \code{vector} that contains the \code{levels} 
       of the \code{factor}, or c) a character \code{vector} that contains the 
       \code{labels} of the \code{levels} of the \code{factor}.}
  \item{order}{Either \code{standard} or \code{yates}. The order in which the 
        \code{levels} combinations of the 
        \code{factors} in \code{factor.names} are to be considered as numbered; 
        \code{standard} numbers them as if they are arranged in standard order, 
        that is with the first factor moving slowest and the last factor 
        moving fastest; \code{yates} numbers them 
        as if they are arranged in Yates order, that is with the first factor 
        moving fastest and last factor moving slowest.}
}
\value{A \code{data.frame} whose columns consist of the \code{factors} listed in 
  \code{factor.names} and whose values have been computed from the combined 
  \code{factor}. All the \code{factors} will be of the same length.}
\author{Chris Brien}
\seealso{\code{\link{fac.combine}} in package \pkg{dae}.}
\note{A single \code{factor} name may be supplied in the \code{list} in which case 
      a \code{data.frame} is produced that contains the single \code{factor} 
      computed from the numeric \code{vector}. This may be useful when calling 
      this \code{function} from others.}
\examples{
## generate a small completely randomized design for 6 treatments 
n <- 12
CRD.unit <- list(Unit = n)
treat <- factor(rep(1:4, each = 3))
CRD.lay <- fac.layout(unrandomized=CRD.unit, randomized=treat, seed=956)

## divide the treatments into two two-level factor A nd B
CRD.facs <- fac.divide(CRD.lay$treat, factor.names = list(A = 2, B = 2))
}
\keyword{factor}
\keyword{manip}