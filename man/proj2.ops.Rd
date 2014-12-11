\name{proj2.ops}
\alias{proj2.ops}
\title{Compute the projection and Residual operators for two, possibly nonorthogonal, 
       projection matrices}
\description{A procedure that computes the projection operators that decompose 
             the range of Q1 into a part that pertains to Q2 and a part that is 
             orthogonal to Q2. It also produces the nonzero canonical 
             efficiency factors for the joint decomposition of Q1 and Q and the 
             corresponding eigenvectors of Q1 (James and Wilkinson, 1971). 
             Q1 and Q2 may be nonorthogonal.}
\usage{proj2.ops(Q1, Q2)}
\arguments{
 \item{Q1}{A symmetric \code{projector} whose range is to be decomposed.}
 \item{Q2}{A symmetric \code{projector} whose range in Q1 is required.}
}
\value{A \code{list} with the following components:
\enumerate{
\item{\bold{efficiencies:} a \code{vector} containing the nonzero canonical efficiency factors;}
\item{\bold{eigenvectors:} an n x r \code{\link{matrix}}, where n is the order of the projectors and 
            r is the number of nonzero canonical efficiency factors; it contains 
            the eigenvectors of Q1 corresponding to the nonzero canonical 
            efficiency factors.} 
\item{\bold{Qconf:} a \code{projector} onto the part of the range of Q1 with 
            which Q2 is confounded;}
\item{\bold{Qres:} a \code{projector} onto the part of the range of Q1 that is 
            orthogonal to the range of Q2.}}}
\details{
The nonzero canonical efficiency factors are the nonzero eigenvalues of 
Q1 \%*\% Q2 \%*\% Q1 (James and Wilkinson, 1971). An eigenvalue is regarded 
as zero if it is less than \code{daeTolerance}, which is initially set to 1e-10. 
The function \code{\link{set.daeTolerance}} can be used to change \code{daeTolerance}.

The eigenvectors are the eigenvectors of Q1 corresponding to the nonzero canonical 
efficiency factors. The eigenvectors for Q2 can be obtained by premultiplying 
those for Q1 by Q2.

Qres is computed using equation 4.10 from James and Wilkinson (1971) and Qconf 
is obtained by subtracting Qres from Q1.}
\references{James, A. T. and Wilkinson, G. N. (1971) Factorization of the 
            residual operator and canonical decomposition of nonorthogonal 
            factors in the analysis of variance. \emph{Biometrika}, \bold{58}, 279-294.}
\author{Chris Brien}
\seealso{\code{\link{proj2.decomp}}, \code{\link{proj2.efficiency}}, \code{\link{decomp.relate}} 
in package \pkg{dae}.

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

## obtain the projection operators for the interblock analysis
PBIBD2.Bops <- proj2.ops(Q.B, Q.T)
Q.B.T <- PBIBD2.Bops$Qconf
Q.B.res <- PBIBD2.Bops$Qres

## demonstrate their orthogonality
is.allzero(Q.B.T \%*\% Q.B.res)
}
\keyword{array}
\keyword{design}
\keyword{projector}
