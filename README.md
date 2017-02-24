# Circular Drift-diffusion Model 

This package implements circular drift-diffusion model, using Armadillo C++. 
MATLAB callable functions reside in inst/matlab folder.

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

## Installation and Prerequisites


### R package 
```
## From github
devtools::install_github("TasCL/CircularDDM")
## From source: 
install.packages("CircularDDM_0.0.2.tar.gz", repos = NULL, type="source")
```

R package requires:
- R (>= 3.0.2)
- Rtools
- Rcpp (>= 0.12.3)
- RcppArmadillo (>= 0.6.700.6.0)

### MATLAB toolbox 
1. Rename inst/matlab as inst/CircularDDM
2. Copy the "CircularDDM" folder to MATLAB root folder. You can find it our 
by enter "matlabroot" in MATLAB prompt. In Linux with a R2016b version copy, it 
is usually at /usr/local/MATLAB/R2016a/. 
3. Use a text editor to open /usr/local/MATLAB/R2016a/toolbox/local/pathdef.m.
Note that pathdef.m might have not set a write flag. "chmod +w pathdef.m" will 
add w flag on the file.
4. Add "matlabroot,'/toolbox/CircularDDM:', ..." before "%%% END ENTRIES %%%"
5. Save pathdef.m and leave.
6. Enter the following three command in MATLAB prompt to test if CircularDDM
package is installed properly.

```
help('rcircularddm')
help('dcircularddm')
help('rvonmises')
```

MATLAB package requires:
- GNU gsl 1.16 (gsl require CBLAS)
- Armadillo ( >= 0.6.700.6.0)
- Armadillo requires LAPACK and BLAS or OpenBLAS (alread including LAPACK), or
Intel Math Kernel Library (MKL)
- See [Armadillo C++](http://arma.sourceforge.net/download.html) for details.


## References
* Smith, P. L. (2016). Diffusion Theory of Decision Making in Continuous Report.
Psychological Review, 123(4), 425--451, doi:  http://dx.doi.org/10.1037/rev0000023.

