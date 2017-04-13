# Circular Drift-Diffusion Model 

CircularDDM implements circular drift-diffusion/decision-diffusion model, using 
Armadillo C++.  

MATLAB callable functions reside in inst/matlab folder.  Currently, we provide 
only Linux binary (mexa64 files compiled by g++ 5.4.0). If the libraries, 
GNU gsl and armadillo C++ (requiring LAPACK and BLAS), are installed,
one should be able to compile the source codes in other systems, such as 
Windows (mexw64 files) and OS X (mexmaci64 files).  

## Getting Started

The following examples are from CircularDDM help pages

```
require(CircularDDM)

###################
## dcddm example ##
###################
x <- cbind(
RT= c(1.2595272, 0.8693937, 0.8009044, 1.0018933, 2.3640007, 1.0521304),
R = c(1.9217430, 1.7844653, 0.2662521, 2.1569724, 1.7277440, 0.8607271)
)
pVec <- c(a=2.45, vx=1.5, vy=1.25, t0=.1, s=1)
dcddm(x, pVec)

###################
## rcddm example ##
###################
pVec <- c(a=2, vx=1.5, vy=1.25, t0=.25, s=1)
den  <- rcddm1(1e3, pVec);
hist(den[,1], breaks = "fd", xlab="Response Time",  main="Density")
hist(den[,2], breaks = "fd", xlab="Response Angle", main="Density")

t0    <- .25
thre  <- runif(10, 0, 2)
angle <- rvm(1e3, 0, 10)
sp    <- cbind(runif(10, 0, 1), runif(10, 0,1))
dat   <- rcddm2(1e3, thre, angle, sp, t0)


## Fix starting point at a constant
t0    <- .25
thre  <- runif(1, 0, 2); threshold
angle <- rvm(1e3, 0, 10)
sp    <- cbind(runif(1, 0, 1),  runif(1, 0, 1)) 
dat   <- rcddm2(1e3, threshold, angle, startpoint, t0)

```

## Installation and Prerequisites

### R package 
```
## From github
devtools::install_github("TasCL/CircularDDM")
## From source: 
install.packages("CircularDDM_0.1.4.tar.gz", repos = NULL, type="source")
```

R package requires:

- R (>= 3.0.2)
- Rtools
- Rcpp (>= 0.12.3)
- RcppArmadillo (>= 0.6.700.6.0)

### MATLAB toolbox 
The installation for MATLAB is simply copying the files in 'matlab' folder to 
MATLAB toolbox folder. For example, in Linux default path, MATLAB toolbox is at
'/usr/local/MATLAB/R2016a/toolbox' (version number may vary, e.g., 'R2015', etc.). 
The following is a step-by-step instruction.

1. Create a folder (e.g., 'mkdir') called 'CircularDDM' under
'/usr/local/MATLAB/R2016a/toolbox'
2. Copy all files in the 'inst/matlab' in the source tar ball to
'/usr/local/MATLAB/R2016a/toolbox/CircularDDM'. Entering 'matlabroot' in MATLAB
prompt will show MATLAB root folder.  
3. Use a text editor to open /usr/local/MATLAB/R2016a/toolbox/local/pathdef.m.
Note that pathdef.m might have not set a write flag. "chmod +w pathdef.m" will 
add write ('w') flag on the file.
4. Add "matlabroot,'/toolbox/CircularDDM:', ..." before "%%% END ENTRIES %%%"
5. Save and close the file.
6. Enter the following three command in MATLAB prompt to test if CircularDDM
package is installed properly.

```
help('rcircularddm')
help('dcircularddm')
help('rvonmises')
```

MATLAB package requires:

- GNU gsl 1.16 (gsl requires CBLAS)
- Armadillo ( >= 0.6.700.6.0)
- Armadillo requires either (1) LAPACK and BLAS, (2) OpenBLAS or (3) Intel Math
Kernel Library (MKL). I have not tested the third option. However Aramadillo's 
author recommends it as a perhaps even faster option. 
- See [Armadillo C++](http://arma.sourceforge.net/download.html) for details.

## Contributors

- Yi-Shin Lin <yishin.lin@utas.edu.au> 
- Andrew Heathcote <andrew.heathcote@utas.edu.au> 
- Peter Kvam <kvam.peter@gmail.com> 

## References
* Smith, P. L. (2016). Diffusion Theory of Decision Making in Continuous Report.
Psychological Review, 123(4), 425--451, http://dx.doi.org/10.1037/rev0000023.

