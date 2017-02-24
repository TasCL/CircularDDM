#include <iostream>
#include <random>
#include "mex.h"

void rvm(int n, double mu, double k, double y[]) {
  double U, a, b, r, z, f, c;
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution1(0.0, 1.0);
  std::uniform_real_distribution<double> distribution2(0.0, 2.0*M_PI);

  a = 1 + std::sqrt(1 + 4 * k * k);
  b = (a - std::sqrt(2*a))/(2*k);
  r = (1 + b*b)/(2*b);

  int i = 0;
  do {
    z = std::cos(M_PI * distribution1(generator));
    f = (1. + r * z)/(r + z);
    c = k * (r - f);
    U = distribution1(generator);

    if(c * (2 - c) > U) {
      y[i] = (distribution1(generator) > .50) ? std::acos(f) + mu : -std::acos(f) + mu;
      if(k == 0) {y[i] = distribution2(generator);}
      i++;
    } else {
      if(std::log(c/U) + 1 >= c) {
        y[i] = (distribution1(generator) > .50) ? std::acos(f) + mu : -std::acos(f) + mu;
        if(k == 0) {y[i] = distribution2(generator);}
        i++;
      }
    }
  } while(i < n);
}

void mexFunction(int nlhs, mxArray *plhs[], /*output*/
                 int nrhs, const mxArray *prhs[]) /* Input variables */
{
    double *n, *mu, *k, *out;
    if (nrhs < 3) {      /* Check for proper number of arguments */
	    mexErrMsgIdAndTxt( "MATLAB:rvm:invalidNumInputs",
                "Three input arguments required.");
    }

    n  = mxGetPr(prhs[0]);    // Assign pointers three inputs */
    mu = mxGetPr(prhs[1]);
    k  = mxGetPr(prhs[2]);

    /* Create an m-by-n mxArray; you will copy existing data into it */
    plhs[0] = mxCreateNumericMatrix(*n, 1, mxDOUBLE_CLASS, mxREAL);
    out     = mxGetPr(plhs[0]);
    rvm(*n, *mu, *k, out);  // call rvm
    return;
}


