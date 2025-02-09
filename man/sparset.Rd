\name{sparset}
\alias{sparset}

\title{
Fits a grid of treatments to optimize the the sparse testing design.
}
\description{
Using the AlphaSimR machinery it creates a population structure which is then subject to
different sparse testing designs to tacke a specified target population of environments.
}
\usage{

sparset(gridEnvs, gridIndsPerEnv, lb=10, ub=Inf, 
        nCrosses=30, nProgeny=10, 
        G=NULL, h2=0.5, p=0.1,
        nItersMacs=10, nItersTrt=5,
        verbose=TRUE)

}

\arguments{

  \item{gridEnvs}{
A numeric vector specifying the grid of values to explore in terms of number of environments where we should test.
}
  \item{gridIndsPerEnv}{
A numeric vector specifying the grid of values to explore in terms of number of individuals grown per environment.
}
  \item{lb}{
A scalar specifying the minimum number of plots that a treatment defined by the grid should have to go into the simulation.
}
  \item{ub}{
A scalar specifying the maximum number of plots that a treatment defined by the grid should have to go into the simulation.
}
  \item{nCrosses}{
A numeric value indicating how many crosses our population consists of. See details.
}
  \item{nProgeny}{
A numeric value indicating in average how many progeny per cross our population consists of. See details.
}
  \item{G}{
A square symmetric matrix for the genetic correlations expected in the target population of environments (TPE).
}
  \item{h2}{
A scalar value defining the heritability assumed per environment.
}
  \item{p}{
A scalar value defining the assumed proportion to be selected to calculate the expected genetic gain.
}
  \item{nItersMacs}{
A scalar value defining the number of coalescent simulations to use. We recommend to set this value to your required level of certainty. Normally a value of 10 to 20 is good enough.
}
  \item{nItersTrt}{
A scalar value defining the number of different randomizations to use when running the sparse testing. We recommend to set this to a value between 5 and 10.
}
  \item{verbose}{
A logical value indicating if messages should be printed.
}

}
\details{

Using the \code{AlphaSimR} machinery (runMacs) it creates a population structure given by the arguments
\code{nCrosses} and \code{nProgeny} which is important to the accuracy of the genomic selection models.
Then this populations is grown in a target population of environments defined by the G matrix.

The individuals are grown in different sparse designs defined by a grid formed by the arguments \code{gridEnvs}
and \code{gridIndsPerEnv} as:

\code{expand.grid(gridEnvs,gridIndsPerEnv)}

and the resulting grid of treatments are randomized and grown \code{nItersTrt} times. The phenotypes harvesed 
from those trials is assumed to be harvested at a certain \code{h2} value. This process is repeated \code{nItersMacs} times.

}
\value{

\describe{

\item{accuracy}{the accuracy value measured as the Pearson correlation between true and estimated breeding values for a given treatment. }
\item{$nFarms}{the number of farms where the individuals were grown}
\item{$nInds}{the number of individuals that were grown in each environment.}
\item{$repTrt}{the replicate for a given sparse testing randomization.}
\item{$repMacs}{the replicate for a given coalescent simulation.}
\item{$nPlots}{the total number of plots used.}
\item{$propMaxPlot}{the total number of plots used expressed as a proportion with respect to the maximum number of plots used in the biggest scenario of the simulation.}
\item{$trt}{a definition of the treament expressed in terms of percentage of plots and total number of plots used.}
\item{$overlapInds}{the average percentage of individuals overlapping between environments.}
\item{$sparsityInds}{the average percentage of sparsity between environments.}

}

}

\references{

Giovanny Covarrubias-Pazaran (2025).  obs: a library of optimized breeding schemes questions. To be submitted to Bioinformatics.

Gaynor, R. Chris, Gregor Gorjanc, and John M. Hickey. 2021. AlphaSimR: an R package for breeding program simulations. G3 Gene|Genomes|Genetics 11(2):jkaa017. https://doi.org/10.1093/g3journal/jkaa017.

Chen GK, Marjoram P, Wall JD (2009). Fast and Flexible Simulation of DNA Sequence Data. Genome Research, 19, 136-142. http://genome.cshlp.org/content/19/1/136.

}

\examples{

# define the TPE complexity
G = simGECorMat0(nEnv =300,nMegaEnv = 1,mu=0.2, v=0.2)
hist(G[lower.tri(G)])

 \donttest{

# fit the simulations
res <- sparset( 
  # define grids to explore
  gridEnvs =  seq(5,30,5), gridIndsPerEnv = seq(100,300,50), 
  # define constraints in plots
  lb=0, ub=Inf, 
  # define population structure parameters
  nCrosses = 30, nProgeny = 10,
  # define TPE structure parameters
  G=NULL, h2=0.15, p=0.1,
  # define replications for the simulations
  nItersMacs = 5, nItersTrt = 5 
)

plot.sparsetMod(res) # see expected gain from each scenario
plot.sparsetMod(res, y="accuracy") # what is the accuracy of each scenario

 }

}
\seealso{\code{\link{obs}} -- the information of the package }
