\name{print.summary.pcanon}
\alias{print.summary.pcanon}
\title{Prints the values in an \code{\link{summary.pcanon}} object}
\description{Prints a \code{summary.pcanon} object, which is also a 
             \code{data.frame}, in a pretty format.}
\usage{\method{print}{summary.pcanon}(x, ...)}
\arguments{
 \item{x}{A \code{summary.pcanon} object.}
 \item{...}{further arguments passed to \code{print}.}
}
\value{No value is returned.}


\seealso{\code{\link{summary.pcanon}}}
\examples{
## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Block = 6, Unit = 4)
PBIBD2.nest <- list(Unit = "Block")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- fac.layout(unrandomized = PBIBD2.unit, 
                         nested.factors=PBIBD2.nest, 
                         randomized = trt)

##obtain combined decomposition and summarize
unit.trt.canon <- projs.canon(list(unit=~ Block/Unit, trt=~ trt), data = PBIBD2.lay)
summ <- summary(unit.trt.canon, which = c("aeff","eeff","order"))
print(summ)
}
\keyword{design}
\keyword{projector}