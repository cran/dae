\name{interaction.ABC.plot}
\alias{interaction.ABC.plot}
\title{Plots an interaction plot for three factors}
\description{Plots a \code{function} (the mean by default) of the \code{response} for
       the combinations of the three \code{\link{factor}s} specified as the \code{x.factor}
       (plotted on the x axis of each plot), the \code{groups.factor} (plotted
       as separate lines in each plot) and the \code{trace.factor} (its \code{levels}
       are plotted in different plots). Interaction plots for more than three
       \code{\link{factor}s} can be produced by using \code{\link{fac.combine}} to combine all but
       two of them into a single \code{\link{factor}} that is specified as the
       \code{trace.factor}.}
\usage{interaction.ABC.plot(response, x.factor, groups.factor, 
       trace.factor,data, fun="mean", title="A:B:C Interaction Plot", 
       xlab, ylab, key.title, lwd=4, columns=2, ...)}
\arguments{
 \item{response}{A numeric \code{vector} containing the response variable from 
                 which a function (the mean by
                 default) is computed for plotting on the y-axis.}
 \item{x.factor}{The \code{\link{factor}} to be plotted on the x-axis of each plot.}
 \item{groups.factor}{The \code{\link{factor}} plotted as separate lines in each plot.}
 \item{trace.factor}{The \code{\link{factor}} for whose \code{levels} there are separate plots.}
 \item{data}{A \code{\link{data.frame}} containing the three \code{factors} and the 
             \code{response}.}
 \item{fun}{The \code{function} to be computed from the \code{response} for each
            combination of the three factors \code{x.factor},
            \code{groups.factor} and \code{trace.factor}. By default, the mean
            is computed for each combination.}
 \item{title}{Title for plot window. By default it is "A:B:C Interaction Plot".}
 \item{xlab}{\code{Label} for the x-axis. By default it is the name of the
             \code{x.factor}.}
 \item{ylab}{\code{Label} for the y-axis. By default it is the name of the
             \code{response}.}
 \item{key.title}{\code{Label} for the xkey to the lines in each plot. By default it
                  is the name of the \code{groups.factor}.}
 \item{lwd}{The width of the \code{lines}. By default it is 4.}
 \item{columns}{The number of columns for arranging the several plots for the
                 levels of the \code{groups.factor}. By default it is 2.}
 \item{\dots}{Other arguments that are passed down to the function \code{\link{xyplot}}.}
}
\value{An object of class "\code{\link[=Lattice]{trellis}}", which is 
       automatically plotted by \code{print.trellis}.}
\note{A \code{\link{data.frame}} called \code{data.means} is created, attached and detached
      during execution of this function.
}
\author{Chris Brien}
\seealso{\code{\link{fac.combine}} in package \pkg{dae}, \code{\link{interaction.plot}}.}
\examples{
## plot for generated data
## use ?ABC.Interact.dat for data set details
data(ABC.Interact.dat)
interaction.ABC.plot(MOE, A, B, C, data=ABC.Interact.dat)

## plot for Example 14.1 from Mead, R. (1990). The Design of Experiments: 
## Statistical Principles for Practical Application. Cambridge, 
## Cambridge University Press.  
## use ?SPLGrass.dat for details
data(SPLGrass.dat)
interaction.ABC.plot(Main.Grass, x.factor=Period,
                     groups.factor=Spring, trace.factor=Summer,
                     data=SPLGrass.dat,
                     title="Effect of Period, Spring and Summer on Main Grass")
}
\keyword{aplot}
\keyword{hplot}
\keyword{design}
