\name{dae-package}
\alias{dae-package}
\alias{dae}
\docType{package}
\title{
Functions Useful in the Design and ANOVA of Experiments
}
\description{
The content falls into the following groupings: (i) Data, (ii) Factor manipulation 
    functions, (iii) Design functions, (iv) ANOVA functions, (v) Matrix functions, 
    (vi) Projector and canonical efficiency functions, and (vii) Miscellaneous functions. 
    A document 'daeDesignRandomization.pdf', available in the doc subdirectory of the
    installation directory for 'dae', describes the use of the package for generating randomized 
    layouts for experiments. The ANOVA functions facilitate the extraction of information when 
    the 'Error' function has been used in the call to 'aov'.
}
\details{
\tabular{ll}{
Package: \tab dae\cr
Version: \tab 2.7-2\cr
Date: \tab 2015-06-24\cr
Depends: \tab 
R (>= 2.10.0),
ggplot2,
methods\cr
License: \tab GPL (>=2)\cr
URL: \tab http://chris.brien.name\cr
Built: \tab R 3.2.0; 2015-05-26 UTC; windows\cr
}

Index:
\preformatted{
(i) Data

ABC.Interact.dat        Randomly generated set of values indexed by
                        three factors
Fac4Proc.dat            Data for a 2^4 factorial experiment
Sensory3Phase.dat       Data for the three-pahse sensory evaluation 
                        experiment in Brien, C.J. and Payne, R.W. (1999)
SPLGrass.dat            Data for an experiment to investigate the
                        effects of grazing patterns on pasture
                        composition

(ii) Factor manipulation functions

as.numfac               Convert a factor to a numeric vector
fac.combine             Combines several factors into one
fac.divide              Divides a factor into several individual
                        factors
fac.gen                 Generate all combinations of several factors
fac.match               Match, for each combination of a set of columns
                        in 'x', the row that has the same combination
                        in 'table'
fac.nested              creates a factor whose values are generated
                        within those of the factor nesting.fac
fac.recode              Recodes the 'levels' and values of a factor
                        using the value in position i of the
                        'newlevels' vector to replace the ith 'level'
                        of the 'factor'.
mpone                   Converts the first two levels of a factor into
                        the numeric values -1 and +1

(iii) Design functions

blockboundary.plot      This function plots a block boundary on a plot
                        produced by 'design.plot'.
design.plot             This function plots treatments within a matrix.
fac.layout              Generate a randomized layout for an experiment
no.reps                 Computes the number of replicates for an
                        experiment
power.exp               Computes the power for an experiment

(iv) ANOVA functions 

fitted.aovlist          Extract the fitted values for a fitted model
                        from an aovlist object
fitted.errors           Extract the fitted values for a fitted model
interaction.ABC.plot    Plots an interaction plot for three factors
qqyeffects              Half or full normal plot of Yates effects
resid.errors            Extract the residuals for a fitted model
residuals.aovlist       Extract the residuals from an aovlist object
strength                Generate paper strength values
tukey.1df               Performs Tukey's
                        one-degree-of-freedom-test-for-nonadditivity
yates.effects           Extract Yates effects

(v) Matrix functions

elements                Extract the elements of an array specified by
                        the subscripts
fac.ar1mat              forms the ar1 correlation matrix for a
                        (generalized) factor
fac.sumop               computes the summation matrix that produces
                        sums corresponding to a factor
fac.vcmat               forms the variance matrix for the variance
                        component of a (generalized) factor
mat.I                   Forms a unit matrix
mat.J                   Forms a square matrix of ones
mat.ar1                 Forms an ar1 correlation matrix
mat.dirprod             Forms the direct product of two matrices
mat.dirsum              Forms the direct sum of a list of matrices


(vi) Projector and canonical efficiency functions

Projector class:

projector               Create projectors
projector-class         Class projector
is.projector            Tests whether an object is a valid object of
                        class projector
print.projector         Print projectors
correct.degfree         Check the degrees of freedom in an object of
                        class projector
degfree                 Degrees of freedom extraction and replacement

Accepts two or more formulae:

projs.canon             A canonical analysis of the relationships between
                        two or more sets of projectors
summary.pcanon          A summary of the results of an analysis of
                        the relationships between two or more sets of 
                        projectors
print.summary.pcanon    Prints the values in an 'summary.pcanon' object
efficiencies.pcanon     Extracts the canonical efficiency factors from a 
                        list of class 'pcanon'

Accepts exactly two formulae:

projs.2canon            A canonical analysis of the relationships between
                        two sets of projectors
projs.combine.p2canon   Extract, from a p2canon object, the projectors
summary.p2canon         A summary of the results of an analysis of
                        the relationships between two sets of projectors
print.summary.p2canon   Prints the values in an 'summary.p2canon' object
                        that give the combined decomposition
efficiencies.p2canon    Extracts the canonical efficiency factors from
                        a list of class 'p2canon'.

Others:

decomp.relate           Examines the relationship between the
                        eigenvectors for two decompositions
efficiency.criteria     Computes efficiency criteria from a set of
                        efficiency factors
fac.meanop              computes the projection matrix that produces
                        means
proj2.eigen             Canonical efficiency factors and eigenvectors
                        in joint decomposition of two projectors
proj2.efficiency        Computes the canonical efficiency factors for
                        the joint decomposition of two projectors
proj2.combine           Compute the projection and Residual operators
                        for two, possibly nonorthogonal, projectors
projs.structure         Orthogonalised projectors for the terms in a formula
show-methods            Methods for Function 'show' in Package dae

(vii) Miscellaneous functions

extab                   Expands the values in table to a vector
get.daeTolerance        Gets the value of daeTolerance for the package
                        dae
harmonic.mean           Calcuates the harmonic mean.
is.allzero              Tests whether all elements are approximately
                        zero
rmvnorm                 generates a vector of random values from a
                        multivariate normal distribution
set.daeTolerance        Sets the value of daeTolerance for the package
                        dae
}
}
\author{
Chris Brien <Chris.Brien@unisa.edu.au>.

Maintainer: Chris Brien <Chris.Brien@unisa.edu.au>
}
\keyword{ package }
