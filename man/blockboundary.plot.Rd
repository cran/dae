\name{blockboundary.plot}
\alias{blockboundary.plot}
\title{This function plots a block boundary on a plot produced by \code{\link{design.plot}}.}
\description{This function plots a block boundary on a plot produced by \code{\link{design.plot}}.
             It allows control of the starting unit, through rstart and cstart, 
             and the number of rows (nr) and columns (nc) from the starting unit 
             that the blocks to be plotted are to cover.}
\usage{blockboundary.plot(bdef = NULL, bseq = FALSE, rstart= 0, cstart = 0, 
                                nr, nc, bcol = 1, bwd = 2)}
\arguments{
 \item{bdef}{A \code{\link{matrix}} of block sizes:
    \itemize{
       \item if there is only one row, then the first element is interpreted as the no. rows in 
         each block and blocks with this number of rows are to be repeated  across the rows of the design.
       \item if there is more than one row, then each row of the matrix specifies a block, 
         with the sequence of rows in the matrix specifying a corresponding 
         sequence of blocks down the rows of the design.}
    Similarly, a single value for a column specifies a repetition of blocks of that size 
         across the columns of the design, while several column values specifies a 
         sequence of blocks across the columns of the size specified.}
 \item{bseq}{A \code{logical} thatt determines whether block numbers are repetitions 
              or sequences of block numbers.}
 \item{rstart}{A \code{numeric} speccifying the row after which the plotting 
               of block boundaries is to start.}
 \item{cstart}{A \code{numeric} speccifying the column after which the plotting 
               of block boundaries is to start.}
 \item{nr}{A \code{numeric} the number of rows (nr), from the starting unit, 
             that the blocks to be plotted are to cover.}
 \item{nc}{A \code{numeric} the number of columns (nc), from the starting unit, 
             that the blocks to be plotted are to cover.}
 \item{bcol}{A \code{character string} specifying the colour of the block boundary.
  
             See \code{Colour specification} under the \code{\link{par}} function.}
 \item{bwd}{A \code{numeric} giving the width of the block boundary to be plotted.}
}
\value{no values are returned, but modifications are made to the currently active plot.}


\seealso{\code{\link{design.plot}}, \code{\link{par}}, \code{DiGGer}}
\examples{\dontrun{
    SPL.Lines.mat <- matrix(as.numfac(Lines), ncol=16, byrow=T)
    colnames(SPL.Lines.mat) <- 1:16
    rownames(SPL.Lines.mat) <- 1:10
    SPL.Lines.mat <- SPL.Lines.mat[10:1, 1:16]
    windows()
    design.plot(SPL.Lines.mat,trts=1:10,new=TRUE,
                rstr="Rows",cstr="Columns", chtdiv=3, rprop = 1,cprop=1,
                plotbndry = TRUE)
    #Plot Mainplot boundaries
    blockboundary.plot(bdef = cbind(4,16), rstart = 1, bwd = 3, bcol = "green", 
                       nr = 9, nc = 16)
    blockboundary.plot(bdef = cbind(1,4), bwd = 3, bcol = "green", nr = 1, nc = 16)
    blockboundary.plot(bdef = cbind(1,4), rstart= 9, bwd = 3, bcol = "green", 
                       nr = 10, nc = 16)
    #Plot all 4 block boundaries            
    blockboundary.plot(bdef = cbind(8,5,5,4), bseq=T, cstart = 1, rstart= 1, 
                       bwd = 3,bcol = "blue", nr = 9, nc = 15)
    blockboundary.plot(bdef = cbind(10,16), bwd=3,bcol="blue", nr=10, nc=16)
    #Plot border and internal block boundaries only
    blockboundary.plot(bdef = cbind(8,14), cstart = 1, rstart= 1,
                       bwd = 3, bcol = "blue", nr = 9, nc =  15)
    blockboundary.plot(bdef = cbind(10,16), bwd = 3, bcol = "blue", 
                       nr = 10, nc = 16)}}

\keyword{design}
\keyword{plot}