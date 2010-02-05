daeTolerance = 1e-10

set.daeTolerance <- function(tolerance = daeTolerance)
{ proj <- environment(set.daeTolerance)
  unlockBinding("daeTolerance", proj)
  daeTolerance <<- tolerance
  lockBinding("daeTolerance", proj)
  invisible(daeTolerance)
}

#function to test whether all elements are zero
"is.allzero" <- function(x)
{ all(x < daeTolerance)
}

correct.degfree <- function(object)
{
#check that degrees of freedom are correct
  check.df <- FALSE
  if (!is.projector(object))
    stop("Must supply a valid object of class projector")
  if (is.na(object@degfree))
    stop("Degrees of freedom are missing. Can use degfree<- to set.")
  e <- eigen(object@.Data, symmetric=T, only.values=T)
  nonzero.e <- e$values[e$values > daeTolerance]
	dflen <- length(nonzero.e)
  if ((object@degfree - dflen) > daeTolerance)
      stop("Degrees of freedom of projector are not correct")   
  check.df <- TRUE
  check.df
 }

validProjector <- function(object)
{
#checks that a matrix is a square projection operator
  Q <- object@.Data
#  tol <- 1e-10
  isproj <- TRUE
  if (nrow(Q) != ncol(Q))
     isproj <- stop("Matrix is not square")
     else if (!is.allzero(t(Q)-Q))
              isproj <- "Matrix is not symmetric"
          else 
          { if (!is.allzero(Q%*%Q - Q))
               isproj <- "Matrix is not idempotent"
          }
   isproj
}

projector <- function(Q)
{ p <- new("projector", .Data=Q)
  validity <- validProjector(p)
  if (class(validity) == "character")
     stop(validity)
  e <- eigen(Q, symmetric=T, only.values=T)
	dfsum <- as.integer(sum(e$values))
  nonzero.e <- e$values[e$values > daeTolerance]
	dflen <- length(nonzero.e)
  if ((dfsum - dflen) > daeTolerance)
    stop("Inconsistency in calculating degrees of freedom of projector")   
  p@degfree=dflen
  p
}

setClass("projector", representation("matrix", degfree = "integer"), prototype(degfree=as.integer(NA)))
setAs(from="projector", to="matrix", 
      def=function(from){m <- from@.Data; m})
setValidity("projector", validProjector, where=".GlobalEnv")

is.projector <- function(object)
{ inherits(object, "projector") & validObject(object)
}

degfree <- function(object)
{ if (!inherits(object, "projector"))
    stop("Must supply an object of class projector")
  object@degfree
}

"degfree<-" <- function(object, value)
#A function to replace supplied or computed the degrees of freedom of the objectector object
{ if (!is.projector(object))
    stop("Must assign to a valid object of class projector")
  if (length(value) == 1)
    object@degfree <- as.integer(value)
  else
  { e <- eigen(object@.Data, symmetric=T, only.values=T)
	  nonzero.e <- e$values[e$values > daeTolerance]
	  object@degfree <- length(nonzero.e)
  }
  object
}

print.projector <- function(x, ...)
{ if (!inherits(x, "projector"))
    stop("Must supply an object of class projector")
  print(as(x, "matrix"))
  cat("degfree: ",x@degfree,"\n")
  invisible(x)
}
#setMethod("print", signature=(x = "projector"), print.projector)
setMethod("show", "projector", function(object) print.projector(object))

meanop <- function(factor)
{
#computes the projection matrix that produces the means corresponding to a (generalized) factor
  if (!is.factor(factor))
    stop("Must supply a single factor as the argument")
  n <- length(factor)
	repl <- table(factor)
  lev.fac <- levels(factor)
  if (0 %in% repl)
  { lev.fac <- names(repl)[!(repl == 0)]
  }
  fact.new <- factor(factor, labels=lev.fac)
	l.fac <- length(lev.fac)
	X.fac <- matrix(rep(fact.new, times=l.fac), ncol=l.fac, dimnames=list(1:n, lev.fac))
	X.fac <- sapply(1:l.fac, function(i, X.fac) as.numeric(X.fac[,i] == lev.fac[i]), X.fac=X.fac)
	repl.new <- table(fact.new)
	div <- diag(1/repl.new, nrow=length(repl.new))
	M.fac <- X.fac %*% div %*% t(X.fac)
	M.fac <- projector(M.fac)
	if (degfree(M.fac) != l.fac)
	   stop("Degrees of freedom of projector not equal to observed number of levels of ",deparse(substitute(factor)),"\n")
	M.fac
}

proj2.decomp <- function(Q1, Q2)
{ #A procedure to compute the eigenvalues and eigenvectors for the decomposition 
  #of Q1 pertaining to Q2  i.e. the common eigenvalues of Q1Q2Q1 and Q2Q1Q2 and the 
  #eigenvectors of Q1 in thier joint dcomposition.
  #They are stored in a list with elements named efficiencies and eigenvectors
  if (!is.projector(Q1) | !is.projector(Q2))
    stop("Must supply valid objects of class projector")
  if (nrow(Q1) != nrow(Q2))
    stop("Matrices not conformable.")
	Q121 <- Q1 %*% Q2 %*% Q1
	eff <- eigen(Q121, symmetric=T)
	nonzero.eff <- eff$values[eff$values > daeTolerance]
	r <- length(nonzero.eff)
	nonzero.eigen <- eff$vectors[,1:r]
	nonzero.eigen <- (abs(nonzero.eigen) > daeTolerance)* nonzero.eigen
	list(efficiencies = nonzero.eff, eigenvectors = nonzero.eigen)
}

proj2.efficiency <- function(Q1, Q2)
{ #A procedure to compute the canonical efficiency factors (eigenvalues) in the joint 
  #decomposition of Q1 and Q2 
  proj.Q1Q2 <- proj2.decomp(Q1, Q2)
  nonzero.eff <- proj.Q1Q2$efficiencies
	nonzero.eff
}

decomp.relate <- function(decomp1, decomp2)
{ #A procedure to examine the the relationship between the eigenvectors in 
  #decomp1 and decomp2
  #decomp is a list produced by proj2.decomp or proj2.ops
  relmat <- crossprod(decomp1$eigenvectors, decomp2$eigenvectors)
  relmat <- (abs(relmat) > daeTolerance)* relmat
  dimnames(relmat) <- list(as.character(round(decomp1$efficiencies,4)), as.character(round(decomp2$efficiencies,4)))
  relmat
}

"proj2.ops" <- function(Q1, Q2)
{ #A procedure to compute the Residual operator for P remove Q when P and Q are nonorthogonal.
  #  Corresponding projection operator for Q in P is also obtained.
  n <- nrow(Q1)
  if (n != nrow(Q2))
    stop("Matrices not conformable.")
  isproj <- is.projector(Q1) & is.projector(Q2)
  #compute efficiencies
  decomp <- proj2.decomp(Q1, Q2)
  Eff.Q1.Q2 <- decomp$efficiencies 
  if (length(Eff.Q1.Q2) == 0) #check matrices are orthogonal
  { stop("Matrices are orthogonal.")
  }
  EffUnique.Q1.Q2 <- unique(Eff.Q1.Q2)
  K <- length(EffUnique.Q1.Q2)
  if (K == 1 & EffUnique.Q1.Q2[1] == 1) #check for just confounded (i.e. eff = 1)
  { Qconf <- projector(Q2)
    Qres <- projector(Q1 - Q2)
  }
  else      #compute projection operators for partially confounded case
  { I <- diag(1, nrow = n, ncol = n)
    Qres <- Q1
    for(i in 1:K)
      Qres <- (Q1 %*% (I - Q2/EffUnique.Q1.Q2[i]) %*% Qres)
    Qres <- projector(Qres)
    Qconf <- projector(Q1 - Qres)
  }
  list(efficiencies = Eff.Q1.Q2, eigenvectors=decomp$eigenvectors, Qconf = Qconf, Qres = Qres)
}

