\name{designAmeasures}
\alias{designAmeasures}
\title{Calculates the average variance of pairwise differences from the variance matrix 
       for predictions}
\description{Calculates the average variance of pairwise differences between,  or of 
             elementary contrasts of, predictions using the variance matrix for the 
             predictions. The weighted average variance of pairwise differences can be 
             computed from a vector of \code{replications}, as described by Williams and 
             Piepho (2015). It is possible to compute either 
             A-optimality measure for different subgroups of the predictions. If groups 
             are specified then the A-optimality measures are calculated for the differences 
             between predictions within each group and for those between predictions from 
             different groups. If groupsizes are specified, but groups are not, the
             predictions will be sequentially broken into groups of the size specified by 
             the elements of groupsizes. The groups can be named.
}
\usage{designAmeasures(Vpred, replications = NULL, groupsizes = NULL, groups = NULL)}
\arguments{
 \item{Vpred}{The variance \code{\link{matrix}} of the predictions. It can be obtained 
              using \code{\link{mat.Vpredicts}}.}
 \item{replications}{A \code{\link{numeric}}  vector whose length is equal to the  
              number of rows (columns) in \code{Vpred} and whose elements are to be 
              used to calculate weights for each pair of differences.}
 \item{groupsizes}{A \code{\link{numeric}} containing group sizes. The sum of the 
                   elements of \code{groupsizes} must be less than or equal to the 
                   order of \code{Vpred}. If \code{groupsizes} is a named 
                   \code{vector}, the names are used to label the groups. 
                   If \code{NULL}, either \code{groups} is used or 
                   the average for all pairwise differences is obtained.}
 \item{groups}{A \code{\link{list}}, each element of which is a \code{\link{numeric}}, 
               vector with integers that specify the subroup of the predictions 
               over whose pairwise differences the variances are 
               to be averaged. If there is more than one group, the variances of 
               all between and within group pairwise differences are averaged. 
               If the elements of \code{groups} are named , the names are used to 
               label the groups. If \code{groups} is \code{NULL}, either 
               \code{groupsizes} is used or the average for all pairwise differences 
               is obtained.}
}
\value{A \code{\link{matrix}} containing the within and between group A-optimality measures.}
\details{The variance matrix of pairwise differences is calculated as 
         \eqn{v_{ii} + v_{jj} - 2 v_{ij}}{v_{ii} + v_{jj} - 2 v_{ij}}, 
         where \eqn{v_{ij}}{v_{ij}} is the element from the ith row and jth column of 
         \code{Vpred}. if \code{replication} is not \code{NULL} then weights are computed as 
         \eqn{r_{i} * r_{j} / \mathrm{mean}(\mathbf{r})}{r_{i} * r_{j} / mean(r)}, 
         where \eqn{\mathbf{r}}{r} is the \code{replication} vector and \eqn{r_{i}}{r_{i}} 
         and \eqn{r_{j}}{r_{j}} are elements of \eqn{\mathbf{r}}{r}. The \eqn{(i,j)}{(i,j)} 
         element of the variance matrix of pairwise differences is multiplied by the  
         \eqn{(i,j)}{(i,j)}th weight. Then the mean of the variances of the pairwise 
         differences is computed for the nominated \code{groups}.}
\author{Chris Brien}
\seealso{\code{\link{mat.Vpred}}, \code{\link{designAnatomy}}.}
\references{Smith, A. B., D. G. Butler, C. R. Cavanagh and B. R. Cullis (2015). 
            Multi-phase variety trials using both composite and individual replicate 
            samples: a model-based design approach. 
            \emph{Journal of Agricultural Science}, \bold{153}, 1017-1029.
            
            Williams, E. R., and Piepho, H.-P. (2015). Optimality and contrasts in 
            block designs with unequal treatment replication. 
            \emph{Australian & New Zealand Journal of Statistics}, \bold{57}, 203-209.}
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
n <- nrow(start.design)
W <- model.matrix(~ -1+ Variety, start.design)
ng <- ncol(W)
Gg<- diag(1, ng)
Vu <- with(start.design, fac.vcmat(Mrep, 0.3) + 
                         fac.vcmat(fac.combine(list(Mrep, Mday)), 0.2) + 
                         fac.vcmat(Frep, 0.1) + 
                         fac.vcmat(fac.combine(list(Frep, Fplot)), 0.2))
R <- diag(1, n)
  
## Calculate the variance matrix of the predicted random Variety effects
Vp <- mat.Vpred(W = W, Gg = Gg, Vu = Vu, R = R)
  
## Calculate A-optimality measure
designAmeasures(Vp)
designAmeasures(Vp, groups=list(fldUndup = c(1:4), fldDup = c(5,6)))
grpsizes <- c(4,2)
names(grpsizes) <- c("fldUndup", "fldDup")
designAmeasures(Vp, groupsizes = grpsizes)
designAmeasures(Vp, groupsizes = c(4))
designAmeasures(Vp, groups=list(c(1,4),c(5,6)))

## Calculate the variance matrix of the predicted fixed Variety effects, elminating the grand mean
Vp.reduc <- mat.Vpred(W = W, Gg = 0, Vu = Vu, R = R, 
                      eliminate = projector(matrix(1, nrow = n, ncol = n)/n))
## Calculate A-optimality measure
designAmeasures(Vp.reduc)


}
\keyword{design}
