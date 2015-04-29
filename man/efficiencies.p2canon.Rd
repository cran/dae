\name{efficiencies.p2canon}
\alias{efficiencies.p2canon}
\title{Extracts the canonical efficiency factors from a list of class \code{p2canon}.}
\description{Produces a \code{list} containing the canonical efficiency factors 
             for the joint decomposition of two sets of projectors 
             (Brien and Bailey, 2009) obtained using \code{\link{projs.2canon}}.}
\usage{efficiencies.p2canon(object, which = "adjusted")}
\arguments{
 \item{object}{A \code{list} of class \code{p2canon} produced by \code{\link{projs.2canon}}.}
 \item{which}{A character string, either \code{adjusted} or \code{pairwise}. For 
              \code{adjusted}, the canonical efficiency factor are adjusted for other 
              projectors from \code{Q2}. For \code{pairwise}, they are the unadjusted 
              canonical efficiency factors between pairs of projectors consisting of one 
              projector from each of two sets.}
}
\value{A \code{list} with a component for each element of the \code{Q1} argument from 
       \code{\link{projs.2canon}}. Each component is \code{list}, each its components 
       corresponding to an element of the \code{Q2} argument from \code{\link{projs.2canon}}}
\references{Brien, C. J. and R. A. Bailey (2009). Decomposition tables for 
            multitiered experiments. I. A chain of randomizations.  
            \emph{The Annals of Statistics}, \bold{36}, 4184 - 4213.}
\author{Chris Brien}
\seealso{\code{\link{projs.2canon}}, \code{\link{summary.p2canon}}, \code{\link{proj2.efficiency}}, \code{\link{proj2.combine}}, \code{\link{proj2.eigen}},  \cr 
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

##obtain projectors using projs.structure
Q.unit <- projs.structure(~ Block/Unit, data = PBIBD2.lay)
Q.trt <- projs.structure(~ trt, data = PBIBD2.lay)

##obtain combined decomposition and summarize
unit.trt.p2canon <- projs.2canon(Q.unit, Q.trt)
efficiencies.p2canon(unit.trt.p2canon)
}
\keyword{array}
\keyword{design}
\keyword{projector}
