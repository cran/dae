\name{fac.match}
\alias{fac.match}
\title{Match, for each combination of a set of columns in \code{x}, the row that has the same combination in \code{table}}
\description{Match, for each combination of a set of columns in \code{x}, 
             the row that has the same combination in \code{table}. 
             That is, there must be only one row in \code{table} for 
             each combination of the specified set of columns in \code{x}. 
             It can be viewed as a generalization of the match function 
             from a single vector to multiple vectors.}
\usage{fac.match(x, table, col.names)}
\arguments{
 \item{x}{an R object, normally a \code{data.frame}, possibly a 
          \code{matrix}.}
 \item{table}{an R object, normally a \code{data.frame}, possibly a 
       \code{matrix}.}
 \item{col.names}{A \code{character vector} giving the columns in \code{x} 
                  and \code{table} that are to be matched.}
}
\value{A \code{\link{vector}} of length equal to \code{x} that gives the 
       rows in \code{table} that match the combinations of 
       \code{col.names} in \code{x}.}
\author{Chris Brien}
\seealso{\code{\link{match}}}
\examples{\dontrun{
index <- fac.match(adj.water.dat, longi.dat,  c("Cart","Days"))
longi.dat[index, c("Weight.Before", "Weight.After")] <-  
                    adj.water.dat[,c("Weight.Before", "Weight.After")]
}}
\keyword{design}
\keyword{factor}
