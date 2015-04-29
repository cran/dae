\name{summary.pcanon}
\alias{summary.pcanon}
\alias{summary,pcanon-method}
\title{A summary of the results of a canonical analysis of the relationships between two or more sets of projectors}
\description{Produces a summary of the efficiency criteria computed from the 
             canonical efficiency factors for the joint decomposition of two or 
             more sets of projectors (Brien and Bailey, 2009) obtained using 
             \code{\link{projs.canon}}.}
\usage{\method{summary}{pcanon}(object, which.criteria = c("aefficiency", "eefficiency"), ...)}
\arguments{
 \item{object}{A \code{list} of class \code{pcanon} produced by \code{projs.canon}.}
 \item{which.criteria}{A character \code{vector} nominating the efficiency criteria to be 
                       included in the summary. It can be \code{none}, \code{all} or some 
                       combination of \code{aefficiency}, \code{mefficiency}, 
                       \code{eefficiency}, \code{sefficiency} and \code{order}.}
 \item{...}{further arguments affecting the summary produced.}
}
\value{An object of classes \code{summary.pcanon} and \code{data.frame}, whose rows correspond 
       to subspaces in the decomposition; it is of the nature of a skeleton analysis of variance 
       table. It has an \code{attribute} named \code{ntiers} that is equal to the 
       number of tiers.}
\references{Brien, C. J. and R. A. Bailey (2009). Decomposition tables for 
            multitiered experiments. I. A chain of randomizations.  
            \emph{The Annals of Statistics}, \bold{36}, 4184 - 4213.}
\author{Chris Brien}
\seealso{\code{\link{projs.canon}}, \code{\link{proj2.efficiency}}, 
         \code{\link{efficiency.criteria}}, \code{\link{proj2.combine}},  
         \code{\link{proj2.eigen}}, \code{\link{projs.structure}}, \cr 
         \code{\link{print.summary.pcanon}}in package \pkg{dae}, 
         \code{\link{eigen}}.

         \code{\link[=projector-class]{projector}} for further information about this class.}
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
summary(unit.trt.canon, which = c("aeff","eeff","order"))
}
\keyword{array}
\keyword{design}
\keyword{projector}
