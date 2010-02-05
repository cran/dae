\name{fac.recode}
\alias{fac.recode}
\title{Recodes the values of a factor using each value in a supplied vector to 
        replace the corresponding level of the factor}
\description{Recodes factor values using values in a vector.}
\usage{fac.recode(factor, newlevels)}
\arguments{
 \item{factor}{The \code{factor} to be recoded.}
 \item{newlevels}{A \code{vector }of length \code{levels(factor)} containing values to use in 
                  the recoding.}
}
\value{A \code{factor}.}
\author{Chris Brien}
\seealso{\code{\link{as.numfac}} and \code{\link{mpone}} in package \pkg{dae}, 
\code{\link{factor}}, \code{\link{relevel}}.}
\examples{
## set up a factor with labels
a <- factor(rep(1:4, 4), labels=c("A","B","C","D"))
 
## recode "A" and "D" to 1 and "B" and "C" to 2
b <- fac.recode(a, c(1,2,2,1))
}
\keyword{factor}
\keyword{manip}