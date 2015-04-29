\name{efficiencies.pcanon}
\alias{efficiencies.pcanon}
\title{Extracts the canonical efficiency factors from a list of class \code{pcanon}.}
\description{Produces a \code{list} containing the canonical efficiency factors 
             for the joint decomposition of two or more sets of projectors 
             (Brien and Bailey, 2009) obtained using \code{\link{projs.canon}}.}
\usage{efficiencies.pcanon(object, which = "adjusted")}
\arguments{
 \item{object}{A \code{list} of class \code{pcanon} produced by \code{\link{projs.canon}}.}
 \item{which}{A character string, either \code{adjusted} or \code{pairwise}. For 
              \code{adjusted}, the canonical efficiency factor are adjusted for other 
              projectors from from the same set. For \code{pairwise}, they are the unadjusted 
              canonical efficiency factors between pairs of projectors consisting of one 
              projector from each of two sets.}
}
\value{A \code{list} with a component for each component of \code{object}, except for the 
       last component -- \code{\link{projs.canon}} for more information about the components 
       of a \code{pcanon} object.}
\references{Brien, C. J. and R. A. Bailey (2009). Decomposition tables for 
            multitiered experiments. I. A chain of randomizations.  
            \emph{The Annals of Statistics}, \bold{36}, 4184 - 4213.}
\author{Chris Brien}
\seealso{\code{\link{projs.canon}}, \code{\link{summary.pcanon}}, \code{\link{proj2.efficiency}}, \code{\link{proj2.combine}}, \code{\link{proj2.eigen}},  \cr 
         \code{\link{projs.structure}} in package \pkg{dae}, 
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
efficiencies.pcanon(unit.trt.canon)
}
\keyword{array}
\keyword{design}
\keyword{projector}
