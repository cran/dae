\name{mat.Vpred}
\alias{mat.Vpred}
\title{Calculates the variances of a set of predicted effects from a mixed model}
\description{Calculates the variances of a set of predicted effects, given the 
             incidence matrix for the effects to be predicted and, optionally, 
             a variance matrix of the effects, a incidence matrix for the 
             nuisance fixed factors and covariates, the variance matrix of the nuisance 
             random effects in the mixed model and the residual variance matrix.}
\usage{
mat.Vpred(W, Gg = 0, X = matrix(1, nrow = nrow(W), ncol = 1), Vu = 0, R)}
\arguments{
 \item{W}{The incidence \code{\link{matrix}} for the effects to be predicted.}
 \item{Gg}{The variance \code{\link{matrix}} of the effects to be predicted. 
           If the effects to be predicted are fixed, set to 0.}
 \item{X}{The incidence \code{\link{matrix}} for the nuisance fixed factors 
          and covariates. The default is a column vector of ones.}
 \item{Vu}{The variance \code{\link{matrix}} of the nuisance random effects. 
           If there are none, set to zero.}
 \item{R}{The residual variance \code{\link{matrix}}.}
}
\value{A \code{\link{matrix}} containing the variances and covariances of the 
       predictions.}
\details{The matrix is calculated as \cr
         \code{Vpred = t(W) \%*\% Vinv \%*\% W + ginv(Gg) - A\%*\%ginv(t(X)\%*\%Vinv\%*\%X)\%*\%t(A)}, 
         where \code{Vinv = ginv(Vu + R)}, \code{A = t(W) \%*\% Vinv \%*\% X} and ginv(B) is the unique Moore-Penrose inverse of B formed using the eigendecomposition of B.}
\author{Chris Brien}
\seealso{\code{\link{designAmeasures}}.}
\references{Smith, A. B., D. G. Butler, C. R. Cavanagh and B. R. Cullis (2015). 
            Multi-phase variety trials using both composite and individual replicate 
            samples: a model-based design approach. 
            \emph{Journal of Agricultural Science}, \bold{153}, 1017-1029.}
\examples{
## Reduced example from Smith et al. (2015)
## Generate two-phase design
mill.fac <- fac.gen(list(Mrep = 2, Mday = 2, Mord = 3))
field.lay <- fac.gen(list(Frep = 2, Fplot = 4))
field.lay$Variety <- factor(c("D","E","Y","W","G","D","E","M"), 
                            levels = c("Y","W","G","M","D","E"))
start.design <- cbind(mill.fac, field.lay[c(3,4,5,8,1,7,3,4,5,8,6,2),])
rownames(start.design) <- NULL

## Set up matrices
W <- model.matrix(~ -1+ Variety, start.design)
ng <- ncol(W)
Gg<- diag(1, ng)
Vu <- with(start.design, fac.vcmat(Mrep, 0.3) + 
                         fac.vcmat(fac.combine(list(Mrep, Mday)), 0.2) + 
                         fac.vcmat(Frep, 0.1) + 
                         fac.vcmat(fac.combine(list(Frep, Fplot)), 0.2))
R <- diag(1, nrow(start.design))
  
## Caluclate information matrix
Vp <- mat.Vpred(W = W, Gg = Gg, Vu = Vu, R = R)
}
\keyword{array}
\keyword{design}
