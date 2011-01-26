\name{mat.dirprod}
\alias{mat.dirprod}
\title{Forms the direct product of two matrices}
\description{Form the direct product of the \eqn{m \times n}{m x n} \code{matrix} 
             \bold{A} and the \eqn{p \times q}{p x q} \code{matrix} \bold{B}. 
             It is also called the Kroneker product and the right direct product. 
             It is defined to be the result of replacing each element of 
             \bold{A}, \eqn{a_{ij}}{a_{ij}}, with \eqn{a_{ij}\bold{B}}{a_{ij}\bold{B}}. 
             The result \code{matrix} is  \eqn{mp \times nq}{mp \times nq}.
             
             The method employed uses the \code{rep} function to form two 
             \eqn{mp \times nq}{mp \times nq} matrices: (i) the direct 
             product of \bold{A} and \bold{J}, and (ii) the direct product of 
             \bold{J} and \bold{B}, where each \bold{J} is a matrix of ones 
             whose dimensions are those required to produce an 
             \eqn{mp \times nq}{mp \times nq} matrix. Then the 
             elementwise product of these two matrices is taken to yield the result.}
\usage{mat.dirprod(A, B)}
\arguments{
 \item{A}{The left-hand \code{matrix} in the product.}
 \item{B}{The right-hand \code{matrix} in the product.}
}
\value{An \eqn{mp \times nq}{mp \times nq} \code{matrix}.}
\seealso{\code{matmult}}
\examples{
    col.I <- mat.I(order=4)
    row.I <- mat.I(order=28)
    V <- mat.dirprod(col.I, row.I)
}
\keyword{array}
