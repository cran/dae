"elements" <- function(x, subscripts)
#a function that returns, in a vector, the elements of x specified by the subscripts.
#x is an array
#subscripts is a two dimensional object interpreted as elements by dimensions
{ if (!(class(subscripts) == "matrix" | class(subscripts) == "data.frame"))
    stop("subscripts must be in a matrix or data.frame")
  n <- nrow(subscripts)
  nsub <- ncol(subscripts)
  if (nsub == 2)
  { sapply(1:n, function(i)x[subscripts[i,1],subscripts[i,2]])
  }
  else
    sapply(1:n, function(i)eval(parse(text=paste("x[", paste(subscripts[i,], collapse=","), "]"))))
}

"mat.I" <- function(order)
{ diag(rep(1, order))
}

"mat.J" <- function(order)
{ n <- order*order
  matrix(rep(1, n), nrow=order, ncol=order)
}

"mat.ar1" <- function(order, rho)
#function to form the correlation matrix of size order with an ar1 pattern for
#correlation parameter rho
{ n <- order*order
  ar1 <- matrix(rep(rho, n), nrow=order, ncol=order)
  row.no <- matrix(rep(1:order, times=order), nrow=order, ncol=order)
  col.no <- matrix(rep(1:order, each=order), nrow=order, ncol=order)
  power <- abs(row.no - col.no)
  ar1 <- ar1^power
}

"mat.dirprod" <- function(A, B)
#function create the direct product of A and B
{ rA <- nrow(A); cA <- ncol(A)
  rB <- nrow(B); cB <- ncol(B)
  Aexp <- A[rep(1:rA, each=rB), rep(1:cA, each=cB)]
  Bexp <- eval(parse(text=paste("cbind(", paste(rep("B", cA), collapse=","), ")")))
  Bexp <- eval(parse(text=paste("rbind(", paste(rep("Bexp", rA), collapse=","), ")")))
  Aexp*Bexp
}
