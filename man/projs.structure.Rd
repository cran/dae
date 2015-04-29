\name{projs.structure}
\alias{projs.structure}
\title{orthogonalised projectors for the terms in a formula}
\description{Produces a set of mutually orthogonal projectors, one for each 
             term in the formula. These specify a structure, or an orthogonal
             decomposition of the data space. Firstly, the primary projector 
             \eqn{\mathbf{X(X'X)^-X'}}{X(X'X)^-X'}, where X is the design 
              matrix for the term, is calculated for each term. Then this 
              projector is made orthogonal to terms aliased with it. 
              Two methods of orthogonalization are available: 
              \code{differencingfirst} and \code{eigenmethods}.}
\usage{projs.structure(formula, orthogonalize = "differencingfirst", data = NULL, ...)}
\arguments{
 \item{formula}{An object of class \code{formula} from which the terms will be 
                obtained.}
 \item{orthogonalize}{A string indicating the method for orthogonalizing a 
                      projector to those for terms that occurred previously 
                      in the formula. Two options are available:  
                      \code{differencingfirst}; \code{eigenmethods}.}
 \item{data}{A data frame contains the values of the factors and variables 
             that occur in \code{formula}.}
 \item{...}{further arguments passed to \code{terms}.}
}
\value{A \code{list} with a component for each term in the \code{formula},  
       the component containing the orthogonalized projector for the term. }
\details{In orthogonalizing a projector to others in the set, the 
         \code{differencingfirst} method subtracts from a primary projector 
         each orthogonalized projector for a term whose factors/variables are a 
         subset of its factors/variables. This relies on ensuring that all 
         projectors whose factors/variables are a subset of the current 
         projector occur before it in the expanded formula. It is checked that 
         the set of matrices are mutually orthogonal. If they are not then 
         eigenmenthods are rmployed. 

         The \code{eigenmethods} uses 
         equation 4.10 of James and Wilkinson (1971), which involves 
         calculating the canonical efficiency factors for pairs of primary 
         projectors. The latter method is needed when differencing 
         does not result in a set of mutually orthogonal projectors. 
         It can be invoked directly, especially when it is known that 
         differencing will not produce an orthogonal set of projectors.}
\references{James, A. T. and Wilkinson, G. N. (1971) Factorization of the 
            residual operator and canonical decomposition of nonorthogonal 
            factors in the analysis of variance. \emph{Biometrika}, \bold{58}, 279-294.}
\author{Chris Brien}
\seealso{\code{\link{proj2.efficiency}}, \code{\link{proj2.combine}}, \code{\link{proj2.eigen}}, 
           \code{\link{projs.2canon}} in package \pkg{dae}, \cr 
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

## manually obtain projectors for units
Q.G <- projector(matrix(1, nrow=24, ncol=24)/24)                         
Q.B <- projector(fac.meanop(PBIBD2.lay$Block) - Q.G)
Q.BP <- projector(diag(1, nrow=24) - Q.B - Q.G)

## manually obtain projector for trt
Q.T <- projector(fac.meanop(PBIBD2.lay$trt) - Q.G)

##compute intrablock efficiency criteria
effic <- proj2.efficiency(Q.BP, Q.T)
effic
efficiency.criteria(effic)

##obtain projectors using projs.structure
Q.unit <- projs.structure(~ Block/Unit, data = PBIBD2.lay)
Q.trt <- projs.structure(~ trt, data = PBIBD2.lay)

##obtain combined decomposition and summarize
unit.trt.p2canon <- projs.2canon(Q.unit, Q.trt)
summary(unit.trt.p2canon, which = c("aeff","eeff","order"))
}
\keyword{array}
\keyword{design}
\keyword{projector}
