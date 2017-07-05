\name{news}
\title{News for Package dae}
\section{Changes in dae version 3.0-04 (2017-07-05)}{
  \itemize{
    \item{Apply DAvid Stern's fix to fac.nested so that it deals with missing values in the nesting factor.}
}}
\section{Changes in dae version 3.0-02 (2017-06-30)}{
  \itemize{
    \item{Added links to documentation.}
}}
\section{Changes in dae version 3.0-01 (2017-06-29)}{
  \itemize{
    \item{Added a component named \code{aliasing} to \code{pstructure.object} that stores the information about (partial) aliasing between sources derived a single \code{formula}. Added a component named \code{aliasing} to  \code{pcanon.object} that collects all the aliasing that has occured in the decomposition.}
    \item{Added the function print.aliasing to print a data.frame of class \code{aliasing} that contains aliasing information.}
    \item{Renamed Ameasures to designAmeasures}
    \item{Replaced the vignette for randomzing designs with a vignette that discusses the facilities for designing experiments that includes the use of designRandomize and designAnatomy.}
}}
\section{Changes in dae version 3.0-0 (2017-06-11)}{
  \itemize{
    \item{Added the new function designRandomize with the intention that it will eventaully replace fac.layout. The new function replaces the arguments \code{unrandomized} and \code{randomized} with \code{recipients} and \code{allocated}, respectively, and \code{nested.factors} with \code{nested.recipients}. Also, it changes the default for \code{unit.permutation} to \code{FALSE}.}
    \item{Deprecated projs.structure, replacing it with pstructure.formula, creating a pstructure.object that stores the projectors, the lists of their term names and source names, and the marginality matrix for the terms. Added arguments \code{marginality} and \code{check.marginality} to allow the marginality matrix to be supplied, overwriting any that is produced by pstructure.formula. The argument \code{check.marginality} allows one to decide if the supplied and calculated marginality matrices are the same. Methods for as.data.frame, marginality and print were added for objects of class \code{pstructure}.}
    \item{Deprecated projs.canon, it being replaced by designAnatomy.}
    \item{Modified the function designAnatomy to accomodate the changes to pstructure. Added new argument \code{labels} to designAnatomy, which has option \code{sources} that will result in sources indicating whether they represent interactions and/or nested effects. Added argument \code{marginality} that allows a list of marginality matrices to be supplied. Added the argument \code{check.marginality} allows one to decide if the supplied and calculated marginality matrices are the same. }
    \item{Added argument \code{labelswap} to summary.pcanon so that can swap between terms and sources in the output.}
    \item{Added new method marginality.pcannon to extract the list of marginality matrices corresponding to the \code{formulae} in the call to designAnatomy. Added the argument \code{labelswap} to the summary.pcanon method so that either terms or sources can be output.}
    \item{Added the oats split-plot design from Yates (1937).}
  }
}
\section{Changes in dae version 2.7-26 (2017-05-04)}{
  \itemize{
    \item Converted projs.canon to designAnatomy, changing the position of \code{data} to be second in the argument list and changing \code{meanTerm} to \code{grandMean}. Retained projs.canon as a synonym with the same arguments as before. Also changed \code{meanTerm} to \code{grandMean} in projs.structure.
    \item Added new functions mat.ncssvar and Zspline that calculate the varaince matrix of the random spline effects and the design matrix for fitting random spline effects, respectively.
  }
}
\section{Changes in dae version 2.7-25 (2017-04-12)}{
  \itemize{
    \item In designPlot, have the whole \code{designMatrix} plotted if \code{labels} is \code{NULL}.
  }
}
\section{Changes in dae version 2.7-24 (2017-02-10)}{
  \itemize{
    \item Add the more general \code{hybrid} option to projs.structure and  projs.canon and made it the default option.
    \item Changed fac.combine so that a single factor would pass through without generating an error message.
    \item Revise interaction.ABC.plot to account for ggplot deprecation.
    \item Add functions mat.sar, to form the correlation matrix from an SAR parameter, and mat.gau, to form a correlation function based on the square of the distance between coordinates.
  }
}
\section{Changes in dae version 2.7-23 (2017-01-05)}{
  \itemize{
    \item Fix a bug in tukey.1df when there is no Error function in the model.
    \item Add tips.
  }
}
\section{Changes in dae version 2.7-22 (2016-10-03)}{
  \itemize{
    \item Fix bug in no.reps when the number of reps to achieve the desired power is less than 2.
  }
}
\section{Changes in dae version 2.7-21 (2016-09-19)}{
  \itemize{
    \item Fix bug in fac.combine that can introduce a level equal to \code{NA} when none of the factors has such a level.
  }
}
\section{Changes in dae version 2.7-20 (2016-09-16)}{
  \itemize{
    \item Rebuild under R 3.4.0.
  }
}
\section{Changes in dae version 2.7-19 (2016-09-07)}{
  \itemize{
    \item Include the argument \code{meanTerm} in projs.structure and projs.canon to allow the option to include the grand mean term.
  }
}
\section{Changes in dae version 2.7-18 (2016-08-07)}{
  \itemize{
    \item Fixed a bug for aliased terms in projs.jandw.
  }
}
\section{Changes in dae version 2.7-17 (2016-05-18)}{
  \itemize{
    \item Extended Ameasures to include a \code{groups} argument that allows explicit specification of the elements of groups for which A values are required.
  }
}
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
