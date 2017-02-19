# Circular Drift-diffusion Model 

This package implements circular drift-diffusion model, using Armadillo C++. A
parallle MATLAB toolbox is tentatively stored in inst/matlab.

## Getting Started

The following examples are extracted from CircularDDM help pages

```
require(CircularDDM)

###################
## dddm example  ##
###################
x <- cbind(
RT= c(1.2595272, 0.8693937, 0.8009044, 1.0018933, 2.3640007, 1.0521304),
R = c(1.9217430, 1.7844653, 0.2662521, 2.1569724, 1.7277440, 0.8607271)
)
pVec <- c(a=2.45, vx=1.5, vy=1.25, t0=.1, s=1)
dddm(x, pVec)

###################
## rddm example  ##
###################
pVec <- c(a=2, vx=1.5, vy=1.25, t0=.25, s=1)
den  <- rddm(1e3, pVec);
hist(den[,1], breaks = "fd", xlab="Response Time", main="Density")
hist(den[,3], breaks = "fd", xlab="Response Angle", main="Density")

```

## Installation 

```
## From github
devtools::install_github("TasCL/CircularDDM")
## From source: 
install.packages("CircularDDM_0.0.1.tar.gz", repos = NULL, type="source")
```

## Prerequisites
 - R (>= 3.0.2)
 - Rtools
 - Rcpp (>= 0.12.3)
 - RcppArmadillo (>= 0.6.700.6.0)
 
## References
* Smith, P. L. (2016). Diffusion Theory of Decision Making in Continuous Report.
Psychological Review, 123(4), 425--451, doi:  http://dx.doi.org/10.1037/rev0000023.

