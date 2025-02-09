# obs: Optimized Breeding Schemes

The obs package is nice wrapper of the AlphaSimR package that enables the 
quick evaluation of different breeding scheme optimization questions. The current
questions that can be optimized are:

Sparse testing: with the 'sparset' function which creates
sparse testing scenarios under different target population of environments for 
a given population structure and returns the expected accuracies and expected
genetic gains (influenced by selection intensities) under the different designs.

More to come.

## Installation

You can install the development version of `obs` from GitHub:

``` r
devtools::install_github('covaruber/obs')
```

## Vignettes

 - [Quick start for the obs package](https://cran.r-project.org/package=obs/vignettes/obs_intro.pdf)
 
## Development

The obs package is under active development. If you have ideas of breeding scheme optimization
questions that can be automated please share those with me. I'll make sure you are credited with
the idea.
