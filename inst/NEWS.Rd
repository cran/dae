\name{news}
\title{News for Package dae}
\section{Changes in dae version 2.7-16 (2016-05-16)}{
  \itemize{
    \item Fixed a bug in fac.gen when there is a factor with a single character value.
  }
}
\section{Changes in dae version 2.7-15 (2016-05-14)}{
  \itemize{
    \item Revised fac.recode and fac.divide to deal with factors with both numeric and character levels and to improve efficiency.
  }
}
\section{Changes in dae version 2.7-14 (2016-05-10)}{
  \itemize{
    \item Added designLatinSqrSys. 
    \item Revised some documentation.
  }
}
\section{Changes in dae version 2.7-13 (2016-04-25)}{
  \itemize{
    \item Added functions mat.Vpred and Ameasures. 
    \item Fixed bugs in fac.vcmat. 
  }
}
\section{Changes in dae version 2.7-12 (2016-03-10)}{
  \itemize{
    \item Reversed order of parameters in mat.ar1 and mat.exp. 
    \item Added functions mat.ar2, mat.ar3, mat.ma1, mat.ma2, mat.arma and mat.sar2. 
  }
}
\section{Changes in dae version 2.7-11 (2016-02-11)}{
  \itemize{
    \item Added \code{xefficiency}, the maximum efficiency factor projector-based functions. 
  }
}
\section{Changes in dae version 2.7-10 (2016-01-19)}{
  \itemize{
    \item Improved the code for analyzing the relationships between pairs of projectors. 
  }
}
\section{Changes in dae version 2.7-9 (2015-01-03)}{
  \itemize{
    \item Include new functions designPlot and blockboundaryPlot, that are the same as design.plot and blockboundary.plot, except that they have less cryptic argument names and some arguments have been modified. It is envisaged that the new functions will supersede the old functions with the old functions being deprecated in due course.
    \item Fixed bug in projs.structure to ensure that arguments, including \code{keep.order}, are passed to the terms function. 
  }
}
\section{Changes in dae version 2.7-8 (2015-12-16)}{
  \itemize{
    \item Add arguments rchtdiv and cchtdiv to design.plot so that can size row and column labels independetly of the cell labels. 
  }
}
\section{Changes in dae version 2.7-7 (2015-11-15)}{
  \itemize{
    \item Improved the documentation of efficiency criteria in a number of functions. 
  }
}
\section{Changes in dae version 2.7-6 (2015-11-02)}{
  \itemize{
    \item Allow a single formula to be specified in formula of projs.canon.
    \item Improved proj2.combine to equate Qconf to Q1 when all of Q1 is confounded by Q2. 
  }
}
\section{Changes in dae version 2.7-5 (2015-10-29)}{
  \itemize{
    \item Fix bug in mat.banded in which the elements were being placed starting with the first subdiagonal instead of the diagonal. 
  }
}
\section{Changes in dae version 2.7-4 (2015-10-01)}{
  \itemize{
    \item Add arrgument \code{ggplotFuncs} to interaction.ABC.plot to allow \code{ggplot} functions to be used to 
    modify the plot. 
    \item Add detect.diff to allow the calculation of the detectable difference for given replication, power, delta and alpha.
  }
}
\section{Changes in dae version 2.7-3 (2015-08-23)}{
  \itemize{
    \item Change option \code{differencingfirst} to \code{differencing} in projs.structure and projs.canon. 
          Code modified so that differencing only orthogonalizes by subtracting from a projector those 
          projectors whose factors are a subset of its factors. 
    \item Added new criteria \code{dforth} to efficiency.criteria.
    \item Added \code{order} to the default setting of \code{which.criteria} in all functions that have this argument.
    \item NEWS added.
  }
}
\section{Changes in dae version 2.7-2 (2015-06-24)}{
  \itemize{
    \item minor reprogramming.
    \item documentation of deprecated functions added.
    \item NEWS added.
  }
}
\section{Changes in dae Version 2.7-1 (2015-05-26)}{
  \itemize{
    \item An except argument added to fac.layout so that nominated unrandomized factors are not
          permuted.
  }
}
\section{Changes in dae Version 2.7-0 (2015-04-29)}{
  \itemize{
    \item Added functions to perform an eigenalysis of the projection matrices for an experimental  
          design to produce the decomposition or skeleton-anova table. These functions include 
          projs.combine.p2canon, projs.structure, projs.2canon, projs.canon, summary.p2canon, 
          summary.pcanon, print.summary.p2canon, print.summary.pcanon. The top-level functions 
          are those with the suffix 'pcanon'.
    \item This version uploaded to CRAN.
  }
}
