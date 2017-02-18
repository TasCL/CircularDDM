#include "mex.h"
#include <iostream>
#include <armadillo>

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

void logLik_resp(double *X, size_t nrow, size_t ncol, double *pvec,
  size_t npvec, double y[])
{
  // This is the first part of equation (23) with exponentiation in Smith (2016)
  // pVec [a, vx, vy, t0, s] == [thresh, mu1, mu2, ndt, sigmasq]
  arma::mat x    = getMat(X, nrow, ncol);
  arma::vec pVec = getVec(pvec, npvec);

  arma::vec rts     = x.col(0);
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

  for(arma::vec::iterator i=out.begin(); i!=out.end(); ++i)
  {
    int idx = std::distance(out.begin(), i);
    y[idx] = *i;
  }
  return;
}

void mexFunction(int nlhs, mxArray *plhs[], /*output*/
                 int nrhs, const mxArray *prhs[]) /* Input variables */
{

  double *X, *pvec, *out; /* pointers to input matrices and output*/
  size_t nrow, ncol, npvec;     /* matrix dimensions */

  if (nrhs != 2) {      /* Check for proper number of arguments */
    mexErrMsgIdAndTxt("MATLAB:logLik_resp:rhs",
     "This function requires 2 input matrices.");
  }

  X     = mxGetPr(prhs[0]); /* pointer to first input matrix */
  pvec  = mxGetPr(prhs[1]); /* pointer to second input matrix */
  nrow  = mxGetM(prhs[0]); /* dimensions of input matrices */
  ncol  = mxGetN(prhs[0]);
  npvec = mxGetN(prhs[1]); /* number of parameter */

  /* Validate input arguments */
  if (ncol != 2) {
    mexErrMsgIdAndTxt("MATLAB:logLik_resp:invalidInputType",
      "Input matrix must be 2 columns.");
  }
  if (npvec != 5) {
    mexErrMsgIdAndTxt("MATLAB:logLik_resp:invalidInputType",
      "pVec must be a 5-parameter vector.");
  }

  X    = mxGetPr(prhs[0]); /* pointer to first input matrix */
  pvec = mxGetPr(prhs[1]);

  /* Create an m x 1 mxArray */
  plhs[0] = mxCreateNumericMatrix( mwSize(nrow), 1, mxDOUBLE_CLASS, mxREAL);
  out     = mxGetPr(plhs[0]);
  logLik_resp(X, nrow, ncol, pvec, npvec, out);  // call logLik_resp
  return;
}


