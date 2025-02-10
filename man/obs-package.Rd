\name{obs-package}
\alias{obs}
\docType{package}
\author{
Giovanny Covarrubias-Pazaran
}
\title{

\strong{O}ptimization  of \strong{B}reeding  \strong{S}chemes
\cr

}
\description{

The obs package is nice wrapper of the AlphaSimR package that enables the 
quick evaluation of different breeding scheme optimization questions. The current
questions that can be optimized are:

\strong{Sparse testing}: with the \code{\link{sparset}} function which creates
sparse testing scenarios under different target population of environments for 
a given population structure and returns the expected accuracies and expected
genetic gains (influenced by selection intensities) under the different designs.

}

\section{Keeping obs updated}{

The obs package is updated on CRAN every 4-months due to CRAN policies but you can find the latest source at https://github.com/covaruber/obs. This can be easily installed typing the following in the R console:

library(devtools)

install_github("covaruber/obs")

This is recommended if you reported a bug, was fixed and was immediately pushed to GitHub but not in CRAN until the next update.
}

\section{Tutorials}{

For tutorials on how to perform different analysis with obs please look at the vignettes by typing in the terminal:

\strong{vignette("obs.intro")}

}

\section{Models Enabled}{

The machinery behind the scenes is AlphaSimR and the Markovian Coalescent Simulator.

}

\section{Bug report and contact}{

If you have any questions or suggestions please post it in https://stackoverflow.com or https://stats.stackexchange.com

I'll be glad to help or answer any question. I have spent a valuable amount of time developing this package. Please cite this package in your publication. Type 'citation("obs")' to know how to cite it.


}
\references{

Giovanny Covarrubias-Pazaran, Dorcus Gemenet, Christian Werner, Marlee Labroo  (2025).  obs: a simple evolutionary algorithm for complex problems. To be submitted to Bioinformatics.

Gaynor, R. Chris, Gregor Gorjanc, and John M. Hickey. 2021. AlphaSimR: an R package for breeding program simulations. G3 Gene|Genomes|Genetics 11(2):jkaa017. https://doi.org/10.1093/g3journal/jkaa017.

Chen GK, Marjoram P, Wall JD (2009). Fast and Flexible Simulation of DNA Sequence Data. Genome Research, 19, 136-142. http://genome.cshlp.org/content/19/1/136.

}
\keyword{ R package }
