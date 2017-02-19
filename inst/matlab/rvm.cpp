#include "mex.h"
#include <iostream>
#include <armadillo>

void rvm(int n, double mu, double k, double y[]) {
  double U, U1, U2, a, b, r, z, f, c;

  a = 1 + std::sqrt(1+4 * k * k);
  b = (a - std::sqrt(2*a))/(2* k);
  r = (1 + b*b)/(2*b);

  arma::vec out(n);
  arma::vec::iterator i = out.begin() ;
  do {
    int idx = std::distance(out.begin(), i);
    z = std::cos(M_PI * arma::as_scalar(arma::randu<arma::vec>(1)));
    f  = (1. + r * z)/(r + z);
    c  = k * (r - f);

    U = arma::as_scalar(arma::randu<arma::vec>(1));
    // U = R::runif(0, 1);
    if(c * (2 - c) > U) {
      U1 = arma::as_scalar(arma::randu<arma::vec>(1));
      *i = (U1 > .50) ? std::acos(f) + mu : -std::acos(f) + mu;
      y[idx] = *i;
      if(k == 0) {
        *i = 2*M_PI * arma::as_scalar(arma::randu<arma::vec>(1));
        y[idx] = *i;
      }
      i++;
    } else {
      if(std::log(c/U) + 1 >= c) {
        U2 = arma::as_scalar(arma::randu<arma::vec>(1));
        *i = (U2 > .50) ? std::acos(f) + mu : -std::acos(f) + mu;
        y[idx] = *i;
        if(k == 0) {
          *i = 2*M_PI * arma::as_scalar(arma::randu<arma::vec>(1));
          y[idx] = *i;
        }
        i++;
      }
    }
  } while(i < out.end());

  return;
}

void mexFunction(int nlhs, mxArray *plhs[], /*output*/
                 int nrhs, const mxArray *prhs[]) /* Input variables */
{
    double *n, *mu, *k, *out;
    if (nrhs != 3) {      /* Check for proper number of arguments */
	    mexErrMsgIdAndTxt( "MATLAB:rvm:invalidNumInputs",
                "Three input arguments required.");
    }

    n  = mxGetPr(prhs[0]);    // Assign pointers three inputs */
    mu = mxGetPr(prhs[1]);
    k  = mxGetPr(prhs[2]);

    /* Create an m-by-n mxArray; you will copy existing data into it */
    plhs[0] = mxCreateNumericMatrix(*k, 1, mxDOUBLE_CLASS, mxREAL);
    out     = mxGetPr(plhs[0]);
    rvm(*n, *mu, *k, out);  // call rvm
    return;
}


