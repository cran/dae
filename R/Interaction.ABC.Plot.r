"interaction.ABC.plot" <- function(response, x.factor, groups.factor, trace.factor, data, 
                                 fun="mean", title="A:B:C Interaction Plot", 
                                 xlab, ylab, key.title, lwd=4, columns=2, ...)
{
# form data.frame containing the means of the response variable for the three factors
  name.r <- deparse(substitute(response))
  name.x <- deparse(substitute(x.factor))
  name.g <- deparse(substitute(groups.factor))
  name.t <- deparse(substitute(trace.factor))
  no.x <- length(levels(as.factor(data[[match(name.x, names(data))]])))
  no.g <- length(levels(as.factor(data[[match(name.g, names(data))]])))
  no.t <- length(levels(as.factor(data[[match(name.t, names(data))]])))
  fnames <- list(x.factor = 1:no.x, groups.factor = 1:no.g, trace.factor = 1:no.t)
  data.means <- fac.gen(generate = fnames, order="yates")
  data.means <- data.frame(data.means, as.vector(tapply(data[[match(name.r, names(data))]], 
                           list(data[[match(name.x, names(data))]],
                                data[[match(name.g, names(data))]], 
                                data[[match(name.t, names(data))]]), FUN=mean, simplify=T)))
  dimnames(data.means)[[2]] <- c(name.x, name.g, name.t, name.r)
  attach(data.means)
  on.exit(detach("data.means"))
# set up arguments for plot
  if (missing(xlab)) xlab <- deparse(substitute(x.factor))
  if (missing(ylab)) ylab <- deparse(substitute(response))
  if (missing(key.title)) key.title <- deparse(substitute(groups.factor))
  formula.plot <- formula(paste(deparse(substitute(response)), " ~ as.numeric(", 
      deparse(substitute(x.factor)), ") | ", deparse(substitute(trace.factor))))
# initiate plot and set options for lines
  trellis.device()
  superpose.line <- trellis.par.get("superpose.line")
  superpose.line$col <- rep(1, 7)
  superpose.line$lty <- 1:7
  superpose.line$lwd <- rep(lwd, 7)
  trellis.par.set("superpose.line", superpose.line)
#  guiModify(class="GraphSheet", GraphSheetColor="Transparent")
# do the plot
  xyplot(formula.plot, groups=groups.factor, main=title, xlab = xlab, ylab=ylab, 
    panel = function(x,y,subscripts,groups) {
                  panel.superpose(x,y,subscripts,groups, type="l")},
    key=list(title=key.title, cex.title=1, transparent=T, columns=columns, 
             text=list(levels(as.factor(groups.factor))),
             lines=Rows(trellis.par.get("superpose.line"), 1:no.g)),
    strip = function(...) strip.default(..., strip.names=c(T,T), style=1),
    as.table=T)
}
