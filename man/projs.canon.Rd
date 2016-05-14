\name{projs.canon}
\alias{projs.canon}
\title{A canonical analysis of the relationships between sets of projectors}
\description{Computes the canonical efficiency factors for the joint 
       decomposition of two or more structures or sets of mutually orthogonally 
       projectors (Brien and Bailey, 2009), orthogonalizing projectors in a set 
       to those earlier in the set of projectors with 
       which they are partially aliased. The results can be summarized in the 
       form of a skeleton ANOVA table.}
\usage{
projs.canon(formulae, orthogonalize = "differencing", 
            which.criteria = c("aefficiency","eefficiency","order"), 
            omit.projectors = c("p2canon", "combined"), data = NULL, ...)}
\arguments{
 \item{formulae}{An object of class \code{list} whose compoents are of class \code{formula}.}
 \item{orthogonalize}{A string indicating the method for orthogonalizing a 
                      projector to those for terms that occurred previously 
                      in the formula. Two options are available:  
                      \code{differencing}; \code{eigenmethods}. 
                      If a single value is given, it is used for all formulae.}
 \item{which.criteria}{A character \code{vector} nominating the efficiency 
                       criteria to be included in the summary of aliasing 
                       between terms within a structure. It can be 
                       \code{none}, \code{all} or some combination of 
                       \code{aefficiency}, \code{mefficiency},
                       \code{sefficiency}, \code{eefficiency}, \code{xefficiency}, 
                       \code{order} and \code{dforthog} -- for details see
                       \code{\link{efficiency.criteria}}. 
                       If \code{none}, no summary is printed.}
 \item{omit.projectors}{A \code{character vector} of the types of  projectors to 
                      omit from the returned \code{pcanon} object. If 
                      \code{pcanon} is included in the vector then the 
                      projectors in these objects will be replaced with a 
                      \code{numeric} containing their degrees of freedom. 
                      If \code{combined} is included in the vector 
                      then the projectors for the combined decomposition will be 
                      replaced with a \code{numeric} containing their degrees of 
                      freedom. If \code{none} is included in the vector then no 
                      projectors will be omitted.}
 \item{data}{A data frame contains the values of the factors and variables 
             that occur in \code{formulae}.}
 \item{...}{further arguments passed to \code{terms}.}
}
\value{A \code{list} of class \code{pcanon}. It has as many components as there 
       are  formulae. The first component is the joint decomposition of two 
       structures derived from the first two formulae, being the \code{p2canon}
       object produced by \code{\link{projs.2canon}}. Then there is a component 
       for each further formulae that contains the \code{p2canon} object obtained 
       by applying \code{\link{projs.2canon}} to the structure for a formula and 
       the already established joint decomposition of the structures for the 
       previous formulae in the \code{formulae}. The last component contains the 
       the \code{list} of the projectors that give the combined canonical 
       decomposition derived from the \code{formulae}.}
\details{For each formula supplied in \code{formulae}, the set of projectors is 
         obtained using \code{\link{projs.structure}}; there is one projector 
         for each term in a formula. Then \code{\link{projs.2canon}} is used 
         to perform an analysis of the canonical relationships between two sets 
         of projectors for the first two formulae. If there are further formulae, 
         the relationships between its projectors and the already established 
         decomposition is obtained using \code{\link{projs.2canon}}. The core 
         of the analysis is the determination of eigenvalues of the product of 
         pairs of projectors using the results of James and Wilkinson (1971).}
\references{Brien, C. J. and R. A. Bailey (2009). Decomposition tables for 
            multitiered experiments. I. A chain of randomizations.  
            \emph{The Annals of Statistics}, \bold{36}, 4184 - 4213.

           James, A. T. and Wilkinson, G. N. (1971) Factorization of the residual 
           operator and canonical decomposition of nonorthogonal factors in the 
           analysis of variance. \emph{Biometrika}, \bold{58}, 279-294.}
\author{Chris Brien}
\seealso{\code{\link{summary.pcanon}}, \code{\link{efficiencies.pcanon}}, 
         \code{\link{projs.2canon}}, \code{\link{projs.structure}} , \cr
         \code{\link{proj2.efficiency}}, \code{\link{proj2.combine}}, 
         \code{\link{proj2.eigen}},  \code{\link{efficiency.criteria}}, 
         in package \pkg{dae}, \code{\link{eigen}}.

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

## Three-phase sensory example from Brien and Payne (1999)
\dontrun{data(Sensory3Phase.dat)
         Eval.Field.Treat.canon <- projs.canon(list(
                                     eval=~ ((Occasions/Intervals/Sittings)*Judges)/Positions, 
                                     field=~ (Rows*(Squares/Columns))/Halfplots,
                                     treats=~ Trellis*Method), data=Sensory3Phase.dat)
         summary(Eval.Field.Treat.canon, which.criteria =c("aefficiency", "order"))
}
}
\keyword{array}
\keyword{design}
\keyword{projector}
