\name{design.plot}
\alias{design.plot}
\title{This function plots treatments within a matrix. It is being superseded by \code{\link{design.plot}}.}
\description{This function plots treatments within a matrix and may be used to 
             build a graphical representation of a matrix, highlighting the position 
             of certain treatments and the blocking factors used in the design. 
             It is a modified version of the function supplied with DiGGer. 
             It includes more control over the labelling of the rows and columns 
             of the design and allows for more flexible plotting of designs with 
             unequal block size. However, it is being supersed by \code{\link{design.plot}},
             which has more descriptive argument names.}
\usage{design.plot(dsgn, trts = NULL, rprop = 1, cprop = 1, label = TRUE,
                 plotchar = NULL, plotbndry = TRUE, 
                 chtdiv = 2, rchtdiv = 1, cchtdiv = 1, 
                 bseq = FALSE, bdef = NULL, bcol = 1, bwd = 2, 
                 rotate = FALSE, new = TRUE, 
                 cstr = "Range",rstr = "Row", rlab = TRUE, clab = TRUE,  
                 font = 1, rdecrease = FALSE, cdecrease = FALSE, ...)
}
\arguments{
 \item{dsgn}{a \code{\link{matrix}} containing a set of integers or characters 
             representing the treatments.}
 \item{trts}{A \code{integer} or \code{character vector} giving specific treatment 
             labels to be plotted.}
 \item{rprop}{a value giving the proportion of the row boundary of cell to plot.}
 \item{cprop}{a value giving the proportion of the column boundary of cell to plot.}
 \item{label}{a \code{logical} to indicate whether treatment labels are to be 
              plotted in the cells. If TRUE, print a label for all treatments or 
              specific treatments listed in \code{trts}. If FALSE, no labels are 
              not printed in the cells.}
 \item{plotchar}{Either a \code{character} vector containing labels for the 
              whole set of treatments or a single \code{integer} specifying a symbol to 
              be used in plotting treatments.}
 \item{plotbndry}{A \code{logical} indicting whether a boundary is to plotted 
              around a cell.}
 \item{chtdiv}{A \code{numeric} specifying the amount by which 
              plotting text and symbols in the cells should be magnified/reduced relative to the default.}
 \item{rchtdiv}{A \code{numeric} specifying the amount by which the labels of the rows of the
              design should be magnified/reduced relative to the default.}
 \item{cchtdiv}{A \code{numeric} specifying the amount by which the labels of the columns of the
              design should be magnified/reduced relative to the default.}
 \item{bseq}{A \code{logical} that determines whether block numbers are repetitions 
              or sequences of block numbers.}
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
 \item{bcol}{A \code{character string} specifying the colour of the block boundary.
  
             See \code{Colour specification} under the \code{\link{par}} function.}
 \item{bwd}{A \code{numeric} giving the width of the block boundary to be plotted.}
 \item{rotate}{A \code{logical} which, if \code{TRUE}, results in the matrix being rotated 
              90 degrees for plotting.}
 \item{new}{A \code{logical} indicating if a new plot is to be produced or the current 
              plot is added to.}
 \item{cstr}{A \code{character string} to use as a label for columns of the \code{\link{matrix}}.}
 \item{rstr}{A \code{character string} to use as a label for rows of the \code{\link{matrix}}.}
 \item{rlab}{A \code{logical} indicating each row of the design is labelled. If the 
              rows of the matrix are labelled, these are used; otherwise 1:nrow is used.}
 \item{clab}{A \code{logical} indicating each column of the design is labelled. If the 
              columns of the matrix are labelled, these are used; otherwise 1:ncol is used.}
 \item{font}{An \code{integer} specifying the font family to be used for row and column labelling.}
 \item{rdecrease}{A \code{logical} indicating whether to reverse the row labels.}
 \item{cdecrease}{A \code{logical} indicating whether to reverse the column labels.}
 \item{...}{further arguments passed to \code{\link{polygon}} in plotting the cell.}
}
\value{no values are returned, but a plot is produced.}

\references{Coombes, N. E. (2009). \emph{DiGGer design search tool in R}. 
            \url{http://www.austatgen.org/files/software/downloads/}}


\seealso{\code{\link{blockboundary.plot}}, \code{\link{par}}, \code{\link{polygon}}, 
         \code{DiGGer}}
\examples{\dontrun{
   design.plot(des.mat, trts=1:4, col="lightblue", new=TRUE, 
               rstr="Lanes", cstr="Positions", chtdiv=3, rprop = 1,cprop=1,
               plotbndry = TRUE)
   design.plot(des.mat, trts=5:87, label=T, col="grey", chtdiv=3, new=FALSE,
               plotbndry = TRUE)
   design.plot(des.mat, trts=88:434, label=T, col="lightgreen", chtdiv=3, 
               new=FALSE, plotbndry = TRUE,
               bseq=TRUE, bdef=cbind(4,10,12), bwd=3, bcol="blue")}}
\keyword{design}
\keyword{plot}