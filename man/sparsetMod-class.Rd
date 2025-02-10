\name{plot.sparsetMod}
\title{Diagnostic Plots for 'sparsetMod' Fits}
\alias{plot.sparsetMod}
\usage{
\method{plot}{sparsetMod}(object, y="expGain", 
             lb=0, ub=Inf, box=FALSE, 
             color="nIndsAcrossFarms", \dots)
}
\arguments{

  \item{object}{a fitted sparset model}
  
  \item{y}{a character value to indicate which column should be plotted as
  y-axis variable. By default is the expected genetic gain.}

  \item{lb}{A scalar specifying the minimum number of plots that a 
  treatment defined by the grid should have to go into the simulation.}
  
    \item{ub}{
  A scalar specifying the maximum number of plots that a treatment
  defined by the grid should have to go into the simulation.
  }
  
  \item{box}{a logical value to indicate if the plot should be of
  type boxplot or violin plot}

  \item{color}{a character value indicating which column should be 
  used to color the treatments.}

  \item{\dots}{optional arguments passed to the lattice
  plot function.}
}
\description{
  diagnostic plots for sparsetMod fits
}
\details{
  Diagnostic plots for the sparse testing simulations aim to show the 
  trade off between different strategies in terms of accuracy and 
  selection intensity.
}
\author{
  original version in \CRANpkg{obs} package by Giovanny Covarrubias-Pazaran.
}
\examples{


# define the TPE complexity
G = simGECorMat0(nEnv =100,nMegaEnv = 1,mu=0.2, v=0.2)
hist(G[lower.tri(G)])

# fit the simulations
res <- sparset( 
  # define grids to explore
  gridEnvs =  seq(5,15,5), gridIndsPerEnv = seq(100,200,50), 
  # define constraints in plots
  lb=0, ub=Inf, 
  # define population structure parameters
  nCrosses = 30, nProgeny = 10,
  # define TPE structure parameters
  G=G, h2=0.15,
  # define replications for the simulations
  nItersMacs = 5, nItersTrt = 5
)

plot.sparsetMod(res) # see expected gain from each scenario
plot.sparsetMod(res, y="accuracy") # what is the accuracy of each scenario



}
