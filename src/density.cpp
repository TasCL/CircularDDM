#include <CircularDDM.hpp>

using namespace Rcpp;

// #include <cstdlib>      // std::rand, std::srand
// #include <ctime>        // std::time
// #include <algorithm>    // std::random_shuffle

void set_seed(unsigned int seed) {
  Environment base_env("package:base");
  Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

arma::vec besselJ(arma::vec x, double nu = 1.0) {
  arma::vec out(x.n_elem);   // Get J1k
  for (size_t i = 0; i < x.n_elem; i++) out(i) = R::bessel_j(x(i), nu);
  return out;
}

inline double findzero(double n, double x0, unsigned int kind,
  double tol = 1e-12, unsigned int MAXIT = 100, double err = 1.0) {
  // Halley's method: tol= tolerance; MAXIT = maximum iterations; err = initial error
  double a, b, x, n1 = n + 1.0;
  unsigned int iter = 0;

  do {
    switch (kind) {
    case 0 :
      a = R::bessel_i(x0, n,  0);
      b = R::bessel_i(x0, n1, 0);
      break;
    case 1:
      a = R::bessel_j(x0, n);
      b = R::bessel_j(x0, n1);
      break;
    case 2:
      a = R::bessel_y(x0, n);
      b = R::bessel_y(x0, n1);
      break;
    default:
      a = 0;
      b = 0;
    }

    err  = (2.0*a*x0*(n*a - b*x0)) / ( 2.0*b*b*x0*x0 - a*b*x0*(4.0 * n1) + (n*n1+x0*x0)*a*a );
    x    = x0 - err;
    x0   = x;
    iter++;

  } while ( (std::abs(err) > tol) & (iter < MAXIT) );

  if (iter > (MAXIT - 1)) {
    Rcpp::Rcout << "Failed to converge. Try a different initial guess\n";
    x = INFINITY;
  }
  return x;
}

// [[Rcpp::export]]
arma::vec besselzero(double nu, unsigned int k, unsigned int kind) {
  double x0;
  unsigned int k3 = 3 * k;
  arma::vec x(k3); x.zeros();
  for (size_t j = 0; j < k3; j++) {
    // Initial guess of zeros
    x0   = 1.0 + std::sqrt(2.0) + (double)j * M_PI + nu + std::pow(nu, 0.4);
    x(j) = findzero(nu, x0, kind);     // Halley's method
    if (x(j) == INFINITY) {Rcpp::stop("Bad guess.");}
  }

  if(!x.is_sorted()) { x = sort(x); };
  arma::vec onevec = arma::ones<arma::vec>(1);
  arma::vec dx     = arma::join_vert(onevec, arma::abs(arma::diff(x)));
  arma::vec out    = x(arma::find(dx > 1e-8));
  if( out.has_nan() ) {Rcpp::stop("NA found.");}
  return out.rows(0, k-1);
}

// [[Rcpp::export]]
arma::vec dhamana(arma::vec RT, double a, double t0, double sigma,
  unsigned int k) {

  unsigned int n = RT.n_elem;
  double dt, sigma2 = sigma*sigma, a2 = a*a;
  double n5siga2 = -0.5*sigma2/a2;

  arma::vec out(n), tmp(k);
  arma::vec j0k  = besselzero(0.0, k, 1); // the zeros of a zero-order Bessel function, J_0(x)
  arma::vec j0k2 = j0k % j0k;
  arma::vec J1k  = besselJ(j0k, 1);

  for (size_t i = 0; i < n; i++) {
    dt = RT(i) - t0;
    tmp.fill(dt * n5siga2);
    out(i) = (dt < .02) ? 1e-100 : sigma2/a2 * arma::accu((j0k/J1k) % arma::exp(tmp % j0k2));
  }
  return out;
}

// [[Rcpp::export]]
arma::vec dcircle(arma::mat x, double a, double mu1, double mu2,
  double t0, double sigma, unsigned int k) {

  arma::vec RT   = x.col(0);
  arma::vec R    = x.col(1);
  unsigned int n = x.n_rows;
  double scaling = 0.5 * (mu1*mu1 + mu2*mu2) / (sigma*sigma);

  arma::vec mu1_vec(n); mu1_vec.fill(mu1);
  arma::vec mu2_vec(n); mu2_vec.fill(mu2);
  arma::vec t0_vec(n);  t0_vec.fill(t0);
  arma::vec term0(n);   term0.fill(a / (sigma*sigma));

  arma::vec term1 = mu1_vec % arma::cos(R) + mu2_vec % arma::sin(R);
  arma::vec term2(n); term2.fill(scaling);
  arma::vec dpta = dhamana(RT, a, t0, sigma, k);
  arma::vec out  = dpta % arma::exp(term0 % term1 - (term2 % (RT - t0_vec)));
  return out;
}

// [[Rcpp::export]]
arma::mat rcircle(unsigned int n, double a, double mu1, double mu2,
  double t0, double sigma, unsigned int k, double lower, double upper,
  unsigned int nth) {

  unsigned int nrt  = 1 + (upper - lower) / lower;
  arma::vec rtRange = arma::linspace(lower, upper, nrt);
  arma::vec thRange = arma::linspace(-M_PI, M_PI, nth);
  arma::mat mat1(1, 2);
  arma::mat x(nrt*nth, 2);

  for (size_t i = 0; i < nth; i++) {
    for (size_t j = 0; j < nrt; j++) {
      unsigned int idx = j + nrt * i;
      mat1.col(0) = rtRange(j);
      mat1.col(1) = thRange(i);
      x.row(idx)  = mat1;
    }
  }

  arma::vec den    = dcircle(x, a, mu1, mu2, t0, sigma, k);
  arma::vec genCDF = arma::cumsum(den / arma::accu(den));
  arma::vec RT(n), R(n);

  for (size_t i = 0; i < n; i++) {
    arma::uword idx = arma::index_min(arma::abs(genCDF - R::runif(0, 1)));
    unsigned int ridx = (idx % nrt);
    unsigned int cidx = std::floor(idx / nrt);
    RT(i) = rtRange(ridx);
    R(i)  = thRange(cidx);
  }

  return arma::join_horiz(RT, R);
}

// [[Rcpp::export]]
arma::vec convDrift(arma::vec mean_v, bool mag2xy) {

  arma::vec out(2);

  if (mag2xy) {
    out.row(0) = mean_v(0)*std::cos(mean_v(1));
    out.row(1) = mean_v(0)*std::sin(mean_v(1));
  } else {
    out.row(0) = std::sqrt(arma::accu(arma::square(mean_v)));
    out.row(1) = std::atan2(mean_v(1), mean_v(0));
  }
  return out;
}

struct generateIntVec {
    int x ;
    generateIntVec() {x = 0;}
    int operator()() {return x++;}
} ; // integer generator:


/*
std::vector<int> genInt (int n) {
    // generate an integer vector based on C index from 0 to n-1
    std::vector<int> out(n); // a empty vector with n integer mem space.
    generateIntVec generator; //
    std::generate(out.begin(), out.end(), generator);
    return out;
}

std::vector<int> genInt2 (int n) {
    // generate an integer vector based on C index from 0 to n-1
    std::vector<int> out(n); // a empty vector with n integer mem space.
    generateIntVec generator; //
    std::generate(out.begin(), out.end(), generator);

    std::srand (unsigned (std::time(0)));
    std::random_shuffle(out.begin(), out.end());
    return out;
}

std::vector<int> genInt3 (int n) {
    std::vector<int> out(n); // a empty vector with n integer mem space.
    for(std::vector<int>::iterator i=out.begin(); i!=out.end(); ++i)
    {
        *i = std::ceil((double)n * R::runif(0.0, 1.0)); // idx == 0
    }
    return out;
}

*/

// #include <RcppArmadilloExtensions/sample.h>
// arma::vec test_sample(int a) {
//   arma::vec x   = arma::randn<arma::vec>(5);
//   int size = 1;
//   bool replace = 1;
//   arma::vec prob = arma::randu<arma::vec>(5);
//   arma::vec out = Rcpp::RcppArmadillo::sample(x, size, replace, prob);
//   return out;
// }

