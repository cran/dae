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

"mat.banded" <- function(x, nrow, ncol)
  #Function to form a banded matrix from a vector of values
  #- band 1 is the diagonal, band 2 the first subdiagonal and so on
{ nband <- length(x)
  if (nband > min(nrow, ncol))
    stop("Have supplied values for more than ",min(nrow, ncol) ," bands")
  matrix <- matrix(0, nrow=nrow, ncol=ncol)
  for (i in 1:nband)
    matrix[row(matrix)==col(matrix)+i-1 | row(matrix)+i-1 == col(matrix)] <- x[i]
  return(matrix)
}

"mat.ar1" <- function(rho, order)
#function to form the correlation matrix of size order with an ar1 pattern for
#correlation parameter rho
{ if (abs(rho) > 1)
    stop(paste("abs(rho) must be less than 1 \n",
                "Note order of mat.ar1 arguments changed to (rho, order)",sep=""))
  n <- order*order
  if (order < 2)
    stop("The order must be 2 or more")
  ar1 <- matrix(rep(rho, n), nrow=order, ncol=order)
  power <- abs(outer(1:order, 1:order, "-"))
  ar1 <- ar1^power
  return(ar1)
}

"mat.ar2" <- function(ARparameters, order)
{ #check ARparameters values
  if (length(ARparameters) <2)
    stop("Must supply two values in ARparameters")
  if (abs(ARparameters[1]) >= 1)
    stop("ARparameters[1] must be less than 1")
  if (abs(ARparameters[2]) >= 1)
    stop("ARparameters[2] must be less than 1")
  if (order < 3)
    stop("The order must be 3 or more")
  
  #Calculate autocorrelations
  corrs <- vector(mode = "numeric", length = order)
  corrs[1] <- 1
  corrs[2] <- ARparameters[1]/(1-ARparameters[2])
  for (k in 3:order)
    corrs[k] <- ARparameters[1]*corrs[k-1] + ARparameters[2]*corrs[k-2]
  
  #Form correlation matrix
  ar2 <- mat.banded(corrs, order, order)  
  return(ar2)
}

"mat.ar3" <- function(ARparameters, order)
{ #check ARparameters values
  if (abs(ARparameters[1]) >= abs(1-ARparameters[2]))
    stop("ARparameters[1] must be less than 1 - ARparameters[2]")
  if (abs(ARparameters[2]) >= 1)
    stop("ARparameters[2] must be less than 1")
  if (abs(ARparameters[3]) >= 1)
    stop("ARparameters[3] must be less than 1")
  if (order < 4)
    stop("The order must be 4 or more")
  
  #Calculate autocorrelations
  corrs <- vector(mode = "numeric", length = order)
  omega <- 1 - ARparameters[2] - ARparameters[3] * (ARparameters[1] + ARparameters[3])
  corrs[1] <- 1
  corrs[2] <- (ARparameters[1]+ARparameters[2]*ARparameters[3])/omega
  corrs[3] <- (ARparameters[1] * (ARparameters[1] + ARparameters[3]) + ARparameters[2] * (1 - ARparameters[2])) / omega
  for (k in 4:order)
    corrs[k] <- ARparameters[1]*corrs[k-1] + ARparameters[2]*corrs[k-2] + ARparameters[3]*corrs[k-3]
  
  #Form correlation matrix
  ar3 <- mat.banded(corrs, order, order)  
  return(ar3)
}

"mat.exp" <- function(rho, coordinates) 
  { if (abs(rho) > 1)
      stop(paste("abs(rho) must be less than 1 \n",
               "Note order of mat.exp arguments changed to (rho, coordinates)",sep=""))
    order <- length(coordinates)
    n <- order * order
    mat <- matrix(rep(rho, n), nrow = order, ncol = order)
    rownames(mat) <- colnames(mat) <- coordinates
    power  <- abs(outer(coordinates,coordinates,'-'))
    mat <- mat^power
    return(mat)
  }

"mat.gau" <- function(rho, coordinates) 
{ if (abs(rho) > 1)
  stop(paste("abs(rho) must be less than 1 \n",
             "Note order of mat.exp arguments changed to (rho, coordinates)",sep=""))
  order <- length(coordinates)
  n <- order * order
  mat <- matrix(rep(rho, n), nrow = order, ncol = order)
  rownames(mat) <- colnames(mat) <- coordinates
  power  <- (outer(coordinates,coordinates,'-'))^2
  mat <- mat^power
  return(mat)
}

"mat.sar" <- function(SARparameter, order)
{ #check SARparameter values
  if (abs(SARparameter) >= 1)
    stop("SARparameter must be less than 1")

  #Calculate autocorrelations
  corrs <- vector(mode = "numeric", length = order)
  corrs[1] <- 1
  corrs[2] <- SARparameter/(1 + (SARparameter*SARparameter/4))
  for (k in 3:order)
    corrs[k] <- SARparameter*corrs[k-1] - (SARparameter*SARparameter/4)*corrs[k-2]
  
  #Form correlation matrix
  sar <- mat.banded(corrs, order, order)  
  return(sar)
}

"mat.sar2" <- function(gamma, order, print = NULL)
{ options <- c("ar3parameters")
  if (is.null(print))
    opt <- "none"
  else  
    opt <- options[unlist(lapply(print, check.arg.values, options=options))]
  
  #Calculate phi values
  phi <- vector(mode = "numeric", length = 3)
  phi[1] <- gamma[1] + 2*gamma[2]
  phi[2] <- -gamma[2] * (2*gamma[1] + gamma[2])
  phi[3] <- gamma[1] * gamma[2] * gamma[2]
  if (opt == "ar3parameters")
    cat(paste("Calculated  parameter values for the associated AR3 process\n", 
              paste(phi, collapse=", "),"\n"))
  
  #Calculate autocorrelations
  sar2 <- mat.ar3(phi, order)
  return(sar2)
}

"mat.ma1" <- function(MAparameter, order)
  #function to form the correlation matrix of size order with an ma1 pattern for
  #MAparameter
{ #check MAparameters value
  if (abs(MAparameter) > 1)
  stop("abs(MAparameter) must be less than 1 \n")
  n <- order*order
  if (order < 2)
    stop("The order must be 2 or more")

  #Calculate autocorrelations
  rho <- -MAparameter / (1 + MAparameter*MAparameter)
  
  #Form correlation matrix
  ma1 <- mat.banded(x = c(1, rho), nrow = order, ncol = order)
  return(ma1)
}

"mat.ma2" <- function(MAparameters, order)
  #function to form the correlation matrix of size order with an ma1 pattern for
  #MAparameter
{ #check MAparameters values
  if (length(MAparameters) <2)
    stop("Must supply two values in MAparameters")
  if (abs(MAparameters[1]) >= 1)
    stop("MAparameters[1] must be less than 1")
  if (abs(MAparameters[2]) >= 1)
    stop("MAparameters[2] must be less than 1")
  if ((MAparameters[1] + MAparameters[2] >= 1) | 
      (MAparameters[1] - MAparameters[2] >= 1))
    stop("The sum and difference of MAparameters must be less than 1")
  if (order < 3)
    stop("The order must be 3 or more")
  
  #Calculate autocorrelations
  div <- 1 + MAparameters[1]*MAparameters[1] + MAparameters[2]*MAparameters[2]
  rho1 <- -MAparameters[1]*(1 - MAparameters[2]) / div
  rho2 <- -MAparameters[2] / div

  #Form correlation matrix
  ma2 <- mat.banded(x = c(1, rho1, rho2), nrow = order, ncol = order)
  return(ma2)
}

"mat.arma" <- function(ARparameter, MAparameter, order)
  #function to form the correlation matrix of size order with an arma(1,1) pattern for
  #given ARparameter and MAparameter
{ #check ARparameter and MAparameter values
  if (abs(MAparameter) > 1)
    stop("abs(MAparameter) must be less than 1 \n")
  if (abs(ARparameter) > 1)
    stop("abs(ARparameter) must be less than 1 \n")
  n <- order*order
  if (order < 3)
    stop("The order must be 3 or more")
  
  corrs <- vector(mode = "numeric", length = order)
  corrs[1] <- 1
  corrs[2] <- (MAparameter - ARparameter)*(1 - MAparameter*ARparameter) / 
                 (1 + MAparameter*MAparameter - 2*MAparameter*ARparameter)
  for (k in 3:order)
    corrs[k] <- ARparameter * corrs[k-1]

  #Form correlation matrix
  arma <- mat.banded(corrs, order, order)  
  return(arma)
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


"mat.dirsum" <- function(matrices)
  #Function to form the direct sum of a list of matrices
{ if (!is.list(matrices))
    stop("Must supply a list for matrices")
  if (!all(unlist(lapply(matrices, is.matrix))))
    stop("All elements of the matrices list must be matrices")
  nr <- lapply(matrices, nrow)
  nc <- lapply(matrices, ncol)
  m <- sum(unlist(nr))
  n <- sum(unlist(nc))
  dsum <- matrix(0, nrow=m, ncol=n)
  r1 <- r2 <- c1 <- c2 <- 0
  for (i in 1:length(matrices))
  { r1 <- r2 + 1
    c1 <- c2 + 1
    r2 <- r2 + nr[[i]]
    c2 <- c2 + nc[[i]]
    dsum[r1:r2, c1:c2] <- matrices[[i]]
  }
  return(dsum)  
}

#Function form variance matrix for random effects for a natural cubic smoothing spline
mat.ncssvar <- function(sigma2s = 1, knot.points, print = FALSE)
{
  #Check knot.points  
  if (is.matrix(knot.points))
    if (ncol(knot.points) != 1)
      stop("knot.points should be a single column matrix")
  if (is.unsorted(knot.points))
    stop("the knot.points should be in increasing order")
  
  h <- as.vector(knot.points)
  r <- length(h)
  krange <- (h[r] - h[1]) / (r - 1)
  h <- diff(h, lag=1) / krange
  
  TD <- matrix(0, nrow = r, ncol = r)
  diag(TD)[2:(r-1)] <- (h[2:(r-1)] + h[1:(r-2)])/3
  TD[row(TD)==col(TD)-1] <- c(0,h[2:(r-2)], 0)/6
  TD[row(TD)==col(TD)+1] <- c(0,h[2:(r-2)], 0)/6
  TD <- TD[2:(r-1), 2:(r-1)]
  TD <- sigma2s * TD
  if (print)
  {
    cat("\n\n#### TD \n\n")
    print(TD)
  }
  
  return(TD)
}

#design matrix for a natural cubic smoothing splines
Zncsspline <- function(knot.points, Gpower = 0, print = FALSE)
{
  #Check knot.points  
  if (is.matrix(knot.points))
    if (ncol(knot.points) != 1)
      stop("knot.points should be a single column matrix")
  if (is.unsorted(knot.points))
    stop("the knot.points should be in increasing order")
  
  h <- as.vector(knot.points)
  r <- length(h)
  krange <- (h[r] - h[1]) / (r - 1)
  h <- diff(h, lag=1) / krange
  
  delta <- diag(as.vector(1/h), nrow = r, ncol = (r-2))
  delta[row(delta)==col(delta)+1] <- -(1/h[1:(r-2)] + 1/h[2:(r-1)])
  delta[row(delta)==col(delta)+2] <- 1/h[2:(r-1)]
  Z <- delta %*% ginv(t(delta) %*% delta)
  if (print)
  {
    cat("\n\n#### delta\n\n")
    print(delta)
  }
  
  TD <- matrix(0, nrow = r, ncol = r)
  diag(TD)[2:(r-1)] <- (h[2:(r-1)] + h[1:(r-2)])/3
  TD[row(TD)==col(TD)-1] <- c(0,h[2:(r-2)], 0)/6
  TD[row(TD)==col(TD)+1] <- c(0,h[2:(r-2)], 0)/6
  TD <- TD[2:(r-1), 2:(r-1)]
  if (print)
  {
    cat("\n\n#### TD \n\n")
    print(TD)
  }
  
  #operator to get power of a matrix
  "%^%" <- function(x, n) 
    with(eigen(x), vectors %*% (values^n * t(vectors))) 
  
  if (Gpower != 0)  
    Z <- Z %*% (TD %^% Gpower)
  return(Z)
}

### Function to calculate the variance of predictions for Genotypes
"mat.Vpred" <- function(W, Gg = 0, X = matrix(1, nrow = nrow(W), ncol = 1), Vu = 0, R)
{ #set W or X to a column vector of 0s and Gg or Vu to a matrix of 0s if not effects for them not random
  if (all(Gg < 1e-08))
    Gginv <- Gg
  else
    Gginv <- solve(Gg)
  Vinv <- ginv(Vu + R)
  A <- t(W)%*%Vinv 
  Vpred <- A%*%W + Gginv
  if (!all(X < 1e-08))
  {
    AX <- A%*%X
    Vpred <- Vpred - AX%*%ginv(t(X)%*%Vinv%*%X)%*%t(AX)
  }
  Vpred <- ginv(Vpred)
  return(Vpred)
}



### Function to calculate 
"designAmeasures" <- function(Vpred, groupsizes = NULL, groups = NULL)
{ 
  #determine any groupings of the variances
  n <- nrow(Vpred)
  if (n != ncol(Vpred))
    stop("Vpred must be square")  
  if (is.null(groupsizes) & is.null(groups))
  {
    groupsizes <- n
    groups <- list(all = 1:n)
    ngrp <- 1
  } else
  {
    if (!is.null(groups)) #working on groups
    {
      if (!is.list(groups))
        stop("the groups argument should be a list")
      if (!is.null(groupsizes))
        warning("Both groups and grouspsizes are set - using groups")
      ngrp <- length(groups)
      groupsizes <- unlist(lapply(groups, length))
    } else #working on group sizes
    {
      ngrp <- length(groupsizes)
      end <- cumsum(groupsizes)
      if (end[ngrp] > n)
        stop("The sum of the group sizes must less than or equal to the size of Vpred")
      if(ngrp == 1)
        start <- 1
      else
        start <- c(0, end[-ngrp]) +1
      groups <- mapply(function(start,end)
                             {start:end},
                       start, end, SIMPLIFY = FALSE)
      if (!is.null(names(groupsizes)))
        names(groups) <- names(groupsizes)
    }
  }
  
  # Calculate the variance matrix for differences between predictions
  varDiff <- matrix(rep(diag(Vpred), each = n), nrow = n) + 
             matrix(rep(diag(Vpred), times = n), nrow = n) - 2 * Vpred
  
  # Calculate all within and between group A-measure values
  A <- matrix(0, nrow = ngrp, ncol = ngrp)
  for (i in 1:ngrp)
  {
    for(j in i:ngrp)
    {
      A[i, j] <- sum(varDiff[groups[[i]], groups[[j]]])
      if (i==j)
        A[i, i] <- A[i, i]/(groupsizes[i]*(groupsizes[i]-1))
      else
      {
        A[i, j] <- A[i, j]/(groupsizes[i]*groupsizes[j])
        A[j, i] <- A[i, j]
      }
    }
  }
  if (!is.null(names(groups)))
  {
    rownames(A) <- colnames(A) <- names(groups)
  }
  return(A)
}
