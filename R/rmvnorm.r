"rmvnorm" <- function(mean, V)
#function to generate a vector of values from an n-dimensional normal 
#distribution whose expectation is given by the n-vector mean and variance by 
#the n x n symmetric matrix V.
#The method:
# a) obtaines the upper-triangular factor of the Cholesky decompostion of V
# b) generates a vector of length n of standard normal values
# c) premultiplies the vector of standard normal values by the transpose of the
#    upper triangular factor and, to the result, adds mean.
{ n <- length(mean)
  if (n!=nrow(V) | n!=ncol(V))
    stop("length of mean must equal the number of rows and the number of columns in V")
  if (!isSymmetric(V))
    stop("V must be a symmetric matrix")
  R <- chol(V)
  z <- rnorm(n)
  u <- t(R) %*% z
  y <- mean + u
  y
}