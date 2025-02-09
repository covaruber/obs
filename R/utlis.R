
##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "3.5.0"))
    stop("This package requires R 3.5.0 or later")
  if(interactive()) {
    packageStartupMessage(blue(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(blue(paste("[] Optimized Breeding Schemes  (obs) 1.0.1 (2025-03)                []",sep="")),appendLF=TRUE)
    packageStartupMessage(paste0(blue("[] Author: Giovanny Covarrubias-Pazaran",paste0(bgGreen(white(" ")), bgWhite(magenta("*")), bgRed(white(" "))),"                        []")),appendLF=TRUE)
    packageStartupMessage(blue("[] Dedicated to the CGIAR, BMGF and USAID.                           []"),appendLF=TRUE)
    packageStartupMessage(blue("[] Type 'vignette('lobs.intro')' for a short tutorial                []"),appendLF=TRUE)
    packageStartupMessage(blue(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(blue("lobs is updated on CRAN every 4-months due to CRAN policies"),appendLF=TRUE)
    packageStartupMessage(blue("Source code is available at https://github.com/covaruber/lobs"),appendLF=TRUE)
  }
  invisible()
}

.onLoad <- function(library, pkg){
  Sys.setenv("OMP_THREAD_LIMIT"=2)
}
##################################################################################################
##################################################################################################

simGECorMat0 <- function (nEnv, nMegaEnv, mu = 0.7, v = 0.2, mu2 = 0, v2 = 0.3) 
{
  ff <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  G = matrix(NA, nEnv, nEnv)
  (nEnv2 <- nEnv/nMegaEnv)
  G
  starts <- seq(1, nEnv, nEnv/nMegaEnv)
  ends <- c((starts - 1)[-1], nEnv)
  for (i in 1:nMegaEnv) {
    corsprov <- rnorm((nEnv2 * (nEnv2 - 1))/2, mu, v)
    counter = 1
    for (j in starts[i]:ends[i]) {
      for (k in j:ends[i]) {
        if (j == k) {
          G[j, k] <- 1
        }
        else {
          G[j, k] <- corsprov[counter]
          counter <- counter + 1
        }
      }
    }
  }
  G <- ff(G)
  tofill <- which(is.na(G), arr.ind = TRUE)
  G[tofill] <- rnorm(nrow(tofill), mu2, v2)
  G[which((G) > 1)] <- 0.98
  G[which((G) < -1)] <- -0.98
  G <- ff(G)
  return(G)
}

addZeros <- function(x){
  nz <- nchar(max(x))
  add <- abs(nchar(x)-nz)
  toAdd <- gsub( "1", "", ifelse(add > 0 , 10^(add) , "") )
  newX <- paste0(toAdd,as.character(x))
  return(newX)
}

