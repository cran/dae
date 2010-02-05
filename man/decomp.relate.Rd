\name{decomp.relate}
\alias{decomp.relate}
\title{Examines the relationship between the eigenvectors for two decompositions}
\description{Two decompositions produced by \code{\link{proj2.decomp}} are compared 
             by computing all pairs of crossproduct sums of eigenvectors from the 
             two decompositions. It is most useful when the calls to 
             \code{\link{proj2.decomp}} have the same Q1.}
\usage{decomp.relate(decomp1, decomp2)}
\arguments{
 \item{decomp1}{A \code{list} containing components efficiencies and eigenvectors such as 
                is produced by \code{\link{proj2.decomp}}.}
 \item{decomp2}{Another \code{list} containing components efficiencies and eigenvectors such as 
                is produced by \code{\link{proj2.decomp}}.}
}
\value{A \code{matrix} that is r1 x r2 where r1 and r2 are the numbers of efficiencies 
       of \code{decomp1} and \code{decomp2}, respectively. The \code{rownames} 
       and \code{columnnames} of the \code{matrix} are the values of the 
       efficiency factors from \code{decomp1} and \code{decomp2}, respectively.}

\details{
Each element of the r1 x r2 \code{matrix} is the sum of crossproducts of a pair of 
eigenvectors, one from each of the two decompositions. A sum is regarded as zero 
if it is less than \code{daeTolerance}, which is initially set to 1e-10. 
The function \code{\link{set.daeTolerance}} can be used to change \code{daeTolerance}.}
\author{Chris Brien}
\seealso{\code{\link{proj2.decomp}}, \code{\link{proj2.ops}} in package \pkg{dae}, 
\code{\link{eigen}}.}
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
Q.B <- projector(meanop(PBIBD2.lay$Blocks) - Q.G)
Q.BP <- projector(diag(1, nrow=24) - Q.B - Q.G)

## obtain projector for trt
Q.T <- projector(meanop(PBIBD2.lay$trt) - Q.G)

## obtain intra- and inter-block decompositions
decomp.inter <- proj2.decomp(Q.B, Q.T)
decomp.intra <- proj2.decomp(Q.BP, Q.T)

## check that intra- and inter-block decompositions are orthogonal
decomp.relate(decomp.intra, decomp.inter) 
}
\keyword{array}
\keyword{design}
\keyword{projector}
