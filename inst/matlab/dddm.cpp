//   Copyright (C) <2017>  <Yi-Shin Lin>
//
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; version 2
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License along
//   with this program; if not, write to the Free Software Foundation, Inc.,
//   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#include "mex.h"
#include <iostream>
#include <armadillo>
#include <gsl/gsl_sf_bessel.h>

arma::vec getVec(double *x, size_t nx) {
  arma::vec out(nx);
  for(int i=0; i<nx; i++) { out[i]=*(x+i); }
  return out;
}

arma::mat getMat(double *x, size_t nrow, size_t ncol) {
  int nx = nrow * ncol, i=0, j=0, k=0;
  arma::mat out(nrow, ncol);

  do {
    for(j=0; j<ncol; j++) {
      for(k=0; k<nrow; k++) {
        out.col(j).row(k) = *(x+i);
        i++;
      }
    }
  } while (i < nx);

  return out;
}

arma::vec besselJ(arma::vec x, double nu) {
    arma::vec out(x.n_elem);
    for(arma::vec::iterator i=x.begin(); i!=x.end(); ++i)
    {
        int idx  = std::distance(x.begin(), i);
        out[idx] = gsl_sf_bessel_Jn(nu, *i); // R's besselJ
        /* if(std::isnan(out[idx])) {
          std::cout << "In besselJ: " << out[idx] << "\t nu is " << nu <<
            "i is " << *i << std::endl;
        } */
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
      a = gsl_sf_bessel_In(n, x0); // Regular Modified Bessel Functions—Fractional Order
      b = gsl_sf_bessel_In(n1, x0);
      break;
    case 1:
      a = gsl_sf_bessel_Jn(n, x0); // Regular Bessel Function—Fractional Order
      b = gsl_sf_bessel_Jn(n1, x0);
      break;
    case 2:
      a = gsl_sf_bessel_Yn(n, x0); // Irregular Bessel Functions—Fractional Order
      b = gsl_sf_bessel_Yn(n1, x0);
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

arma::vec besselzero(double nu, int k, int kind) {
    double x0;
    arma::vec x = arma::zeros<arma::vec>(3*k);
    for (int j=1; j<=3*k; j++) {         // Initial guess of zeros
        x0     = 1 + std::sqrt(2) + (j-1) * M_PI + nu + std::pow(nu, 0.4);
        x(j-1) = findzero(nu, x0, kind);     // Halley's method
        if (x(j-1) == INFINITY) {std::cout << "Bad guess.";}
    }
    if(!x.is_sorted()) { x = sort(x); };
    arma::vec onevec = arma::ones<arma::vec>(1);
    arma::vec dx     = arma::join_vert(onevec, arma::abs(arma::diff(x)));
    arma::vec out    = x(arma::find(dx > 1e-8));
    if( out.has_nan() ) {std::cout << "NA found.";}
    return out.rows(0, k-1);
}

arma::vec logLik_resp(arma::mat x, arma::vec pVec) {
  // This is the first part of equation (23) with exp in Smith (2016)
  // pVec [a, vx, vy, t0, s] == [thresh, mu1, mu2, ndt, sigmasq]
  arma::vec rts = x.col(0);
  arma::vec choices = x.col(1);
  int n = choices.n_elem;

  // Each row of pMat is a replicates of eg pVec[0] of the number of choices
  arma::mat pMat = arma::repmat(pVec, 1, n);
  arma::vec term0, term1, term2, term3, term4, t0_vec, term0_vec, term3_vec;
  term0 = pVec[0] / pVec[4];
  term1 = trans(pMat.row(1)) % arma::cos(choices);
  term2 = trans(pMat.row(2)) % arma::sin(choices);
  term3 = (0.5 * pVec[4]) * ( std::pow(pVec[1], 2) + std::pow(pVec[2], 2) );

  term0_vec = arma::repmat(term0, n, 1);
  term3_vec = arma::repmat(term3, n, 1);
  t0_vec    = arma::repmat(pVec.row(3), n, 1);
  term4     = rts - t0_vec;
  arma::vec out = term0_vec % (term1 + term2) - (term3_vec % term4) ;
  return out;
}
arma::vec logLik_dt(arma::mat x, arma::vec pVec, int k=141) {
  // This is the first part of equation (23) with exp in Smith (2016)
  // pVec [a, vx, vy, t0, s] == [thresh, mu1, mu2, ndt, sigmasq]
  int n, idx; double dt, tmp;
  arma::vec rts, j0k, j0k2, J1, scalar, scalar_vec;
  rts     = x.col(0);
  n       = rts.n_elem;
  j0k     = besselzero(0, k, 1);
  j0k2    = arma::pow(j0k, 2); // squared j0k
  J1      = besselJ(j0k, 1);
  J1(k-1) = J1(k-1) / 2; // replace the last element

  arma::vec out(n);
  for(arma::vec::iterator i=rts.begin(); i!=rts.end(); ++i)
  {
    idx        = std::distance(rts.begin(), i);
    dt         = *i - pVec[3];
    tmp        = -0.5 * pVec[4] * dt / std::pow(pVec[0] , 2);
    scalar     = tmp;
    scalar_vec = arma::repmat(scalar, k, 1);
    out[idx]   = pVec[4] / std::pow(pVec[0], 2) *
      arma::accu(j0k / J1 % arma::exp(scalar_vec % j0k2));
    // When RT and t0 is almost identical, we deem it unlikely.
    if ((*i - pVec[3]) < 0.01) {out[idx] = 1e-10;} //
  }

  return arma::log(out);

}

void dddm(double *X, size_t nrow, size_t ncol, double *pvec, size_t npvec,
  int k, double y[])
{
  arma::mat x    = getMat(X, nrow, ncol);
  arma::vec pVec = getVec(pvec, npvec);

  arma::vec LL_dt   = logLik_dt(x, pVec, k);
  arma::vec LL_resp = logLik_resp(x, pVec);
  arma::vec tmp = LL_dt + LL_resp;
  for(arma::vec::iterator i=tmp.begin(); i!=tmp.end(); ++i)
  {
    int idx = std::distance(tmp.begin(), i);
    y[idx] = *i;
  }
}

void mexFunction(int nlhs, mxArray *plhs[], /*output*/
                 int nrhs, const mxArray *prhs[]) /* Input variables */
{
    double *X, *pvec, *k, *out; /* pointers to input matrices and output*/
    size_t nrow, ncol, npvec;     /* matrix dimensions */

    /* Check for proper number of arguments */
    if (nrhs != 3) {
       mexErrMsgIdAndTxt("MATLAB:dddm:rhs",
                      "This function requires 3 input arguments.");
    }

  X     = mxGetPr(prhs[0]); /* pointer to first input matrix  */
  pvec  = mxGetPr(prhs[1]); /* pointer to second input matrix */
  k     = mxGetPr(prhs[2]); /* precision */
  nrow  = mxGetM(prhs[0]);  /* dimensions of input matrices   */
  ncol  = mxGetN(prhs[0]);
  npvec = mxGetN(prhs[1]);  /* number of parameter */

  /* Validate input arguments */
  if (ncol != 2) {
    mexErrMsgIdAndTxt("MATLAB:dddm:invalidInputType",
      "Input matrix must be 2 columns.");
  }
  if (npvec != 5) {
    mexErrMsgIdAndTxt("MATLAB:dddm:invalidInputType",
      "pVec must be a 5-parameter vector.");
  }

  /* Create an m x 1 mxArray */
  plhs[0] = mxCreateNumericMatrix( mwSize(nrow), 1, mxDOUBLE_CLASS, mxREAL);
  out     = mxGetPr(plhs[0]);
  dddm(X, nrow, ncol, pvec, npvec, *k, out);  // call logLik_dt
  return;
}


