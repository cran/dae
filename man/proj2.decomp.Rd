\name{proj2.decomp}
\alias{proj2.decomp}
\title{Canonical efficiency factors and eigenvectors in joint decomposition of two projectors}
\description{Computes the canonical efficiency factors for the joint decomposition of two 
       projection matrices and the eigenvectors corresponding to the first projector 
       (James and Wilkinson, 1971).}
\usage{proj2.decomp(Q1, Q2)}
\arguments{
 \item{Q1}{An object of class "\code{\link[=projector-class]{projector}}".}
 \item{Q2}{An object of class "\code{\link[=projector-class]{projector}}".}
}
\value{A \code{list} with components efficiencies and eigenvectors. }
\details{
The component efficiencies is a \code{vector} containing the nonzero canonical 
efficiency factors for the joint decomposition of the two projectors. 
The nonzero canonical efficiency factors are the nonzero eigenvalues of 
Q1 \%*\% Q2 \%*\% Q1 (James and Wilkinson, 1971). An eigenvalue is regarded 
as zero if it is less than \code{daeTolerance}, which is initially set to 1e-10. 
The function \code{\link{set.daeTolerance}} can be used to change \code{daeTolerance}.

The component eigenvectors is an n x r \code{\link{matrix}}, where n is the order of the 
projectors and r is the number of nonzero canonical efficiency factors; 
it contains the eigenvectors of Q1 corresponding to the nonzero canonical 
efficiency factors. The eigenvectors for Q2 can be obtained by premultiplying 
those for Q1 by Q2.}
\references{James, A. T. and Wilkinson, G. N. (1971) Factorization of the 
            residual operator and canonical decomposition of nonorthogonal 
            factors in the analysis of variance. \emph{Biometrika}, \bold{58}, 279-294.}
\author{Chris Brien}
\seealso{\code{\link{proj2.efficiency}}, \code{\link{proj2.ops}} in package \pkg{dae}, 
\code{\link{eigen}}.

         \code{\link[=projector-class]{projector}} for further information about this class.}
\examples{
## PBIBD(2) from p. 379 of Cochran and Cox (1957) Experimental Designs. 
## 2nd edn Wiley, New York
PBIBD2.unit <- list(Blocks = 6, Units = 4)
PBIBD2.nest <- list(Units = "Blocks")
trt <- factor(c(1,4,2,5, 2,5,3,6, 3,6,1,4, 4,1,5,2, 5,2,6,3, 6,3,4,1))
PBIBD2.lay <- fac.layout(unrandomized = PBIBD2.unit, 
                         nested.factors=PBIBD2.nest, 
                         randomized = trt)

## obtain projectors for units
Q.G <- projector(matrix(1, nrow=24, ncol=24)/24)                         
Q.B <- projector(fac.meanop(PBIBD2.lay$Blocks) - Q.G)
Q.BP <- projector(diag(1, nrow=24) - Q.B - Q.G)

## obtain projector for trt
Q.T <- projector(fac.meanop(PBIBD2.lay$trt) - Q.G)

## obtain intra- and inter-block decompositions
decomp.inter <- proj2.decomp(Q.B, Q.T)
decomp.intra <- proj2.decomp(Q.BP, Q.T)

#extract intrablock efficiencies
decomp.intra$efficiencies
}
\keyword{array}
\keyword{design}
\keyword{projector}
