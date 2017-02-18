#include "mex.h"
#include <iostream>
#include <armadillo>
#include <gsl/gsl_sf_bessel.h>

arma::vec besselJ(arma::vec x, double nu) {
  arma::vec out(x.n_elem);
  for(arma::vec::iterator i=x.begin(); i!=x.end(); ++i)
  {
    int idx  = std::distance(x.begin(), i);
    out[idx] = gsl_sf_bessel_Jnu(nu, *i); // R's besselJ
  }
  return out;
}

inline double findzero(double n, double x0, int kind, double tol=1e-12,
  int MAXIT=100, double err=1) {
  // Tolerance; Maximum number of times to iterate; Initial error

  double a, b, x, n1=n+1;
  int iter=0;

  do {
    switch (kind) {
    case 0 :
      a = gsl_sf_bessel_Inu(n, x0); // Regular Modified Bessel Functions—Fractional Order
      b = gsl_sf_bessel_Inu(n1, x0);
      break;
    case 1:
      a = gsl_sf_bessel_Jnu(n, x0); // Regular Bessel Function—Fractional Order
      b = gsl_sf_bessel_Jnu(n1, x0);
      break;
    case 2:
      a = gsl_sf_bessel_Ynu(n, x0); // Irregular Bessel Functions—Fractional Order
      b = gsl_sf_bessel_Ynu(n1, x0);
      break;
    default:
      a = 0;
      b = 0;
    }
    err  = (2*a*x0*(n*a - b*x0) ) /
      ( 2*b*b*x0*x0 - a*b*x0*(4*n1) + (n*n1+x0*x0)*a*a );
    x    = x0 - err;
    x0   = x;
    iter = iter + 1;
  } while ( (std::abs(err) > tol) & (iter < MAXIT) );

  if (iter > (MAXIT - 1)) {
    std::cout << "Failed to converge to within tolerance.\n" <<
      "Try a different initial guess";
    x=INFINITY ;
  }

  return x;

}

void besselzero(double nu, int k, int kind, double y[]) {
  double x0;
  arma::vec x = arma::zeros<arma::vec>(3*k);
  for (int j=1; j<=3*k; j++) {     // Initial guess of zeros
    x0     = 1 + std::sqrt(2) + (j-1) * M_PI + nu + std::pow(nu, 0.4);
    x(j-1) = findzero(nu, x0, kind);     // Halley's method
    if (x(j-1) == INFINITY) { std::cout << "Bad guess.\n";}
  }
  if(!x.is_sorted()) { x = sort(x); };
  arma::vec onevec = arma::ones<arma::vec>(1);
  arma::vec dx     = arma::join_vert(onevec, arma::abs(arma::diff(x)));
  arma::vec out    = x(arma::find(dx > 1e-8));
  if( out.has_nan() ) {std::cout << "NA found.\n";}

  arma::vec tmp = out.rows(0, k-1);
  for(arma::vec::iterator i=tmp.begin(); i!=tmp.end(); ++i)
  {
    int idx = std::distance(tmp.begin(), i);
    y[idx] = *i;
  }
  return;
}

void mexFunction(int nlhs, mxArray *plhs[], /*output*/
                 int nrhs, const mxArray *prhs[]) /* Input variables */
{
    double *nu, *k, *kind, *out;
    if (nrhs != 3) {      /* Check for proper number of arguments */
	    mexErrMsgIdAndTxt( "MATLAB:besselzero:invalidNumInputs",
                "Three input arguments required.");
    }

    nu   = mxGetPr(prhs[0]);    // Assign pointers three inputs */
    k    = mxGetPr(prhs[1]);
    kind = mxGetPr(prhs[2]);

    /* Create an k-by-1 mxArray */
    plhs[0] = mxCreateNumericMatrix(*k, 1, mxDOUBLE_CLASS, mxREAL);
    out     = mxGetPr(plhs[0]);
    besselzero(*nu, *k, *kind, out);  // call besselzero
    return;
}
