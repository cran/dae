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
  levels.x <- levels(data.means[[name.x]])
  if (any(is.na(levels.x)))
    stop("The x.factor has NA as a level")
  if (any(is.na(as.numeric(levels.x))))
    data.means[[name.x]] <- fac.recode(data.means[[name.x]], 1:length(levels.x))
  data.means[[name.x]] <- as.numeric(data.means[[name.x]])
# set up arguments for plot
  if (missing(xlab)) xlab <- deparse(substitute(x.factor))
  if (missing(ylab)) ylab <- deparse(substitute(response))
  if (missing(key.title)) key.title <- deparse(substitute(groups.factor))
  formula.plot <- formula(paste(deparse(substitute(response)), " ~ as.numeric(", 
      deparse(substitute(x.factor)), ") | ", deparse(substitute(trace.factor))))
# do the plot
  lbl.fn <- function(variable, value) {
  value <- paste(name.t,": ", as.character(value), sep="")
  return(value)}
  int.plot <- ggplot(data=data.means, aes_string(x = name.x, y = name.r, linetype=name.g, colour=name.g)) +
                     geom_line() + geom_point() +
                     facet_grid(formula(paste("~ ", name.t, sep="")),, labeller = lbl.fn)
  print(int.plot)
}
