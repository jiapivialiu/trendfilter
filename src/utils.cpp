#include <stdexcept>
#include <Eigen/Sparse>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <dspline.h>
#include "utils.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(dspline)]]

typedef Eigen::COLAMDOrdering<int> Ord;

using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::ArrayXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/* General utilities */
Eigen::SparseMatrix<double> row_scale(Eigen::SparseMatrix<double> A,
    Eigen::ArrayXd v) {
  return v.matrix().asDiagonal() * A;
}
Eigen::MatrixXd row_scale(Eigen::MatrixXd A,
    Eigen::ArrayXd v) {
  return v.matrix().asDiagonal() * A;
}
Eigen::SparseMatrix<double> col_scale(Eigen::SparseMatrix<double> A,
    Eigen::ArrayXd v) {
  return A * v.matrix().asDiagonal();
}
Eigen::MatrixXd col_scale(Eigen::MatrixXd A,
    Eigen::ArrayXd v) {
  return A * v.matrix().asDiagonal();
}

/* Matrix construction */
Eigen::SparseMatrix<double> identity(int n) {
  SparseMatrix<double> Id(n, n);
  Id.setIdentity();
  return Id;
}

Eigen::SparseMatrix<double> diagonal(Eigen::ArrayXd diag) {
  Rcpp::Rcout << "getting size" << std::endl;
  int n = diag.size();
  Rcpp::Rcout << "initializating sparse matrix" << std::endl;
  SparseMatrix<double> D(n, n);
  Rcpp::Rcout << "setting diagonal entries" << std::endl;
  D.diagonal() = diag;
  return D;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> get_dk_mat(int k, NumericVector xd,
    bool tf_weighting) {
  int n = xd.size();
  return dspline::rcpp_b_mat(k, xd, tf_weighting, Rcpp::seq(0, n-k-1), true);
}

Eigen::SparseMatrix<double> get_penalty_mat(int k, NumericVector xd) {
  int n = xd.size();
  return dspline::rcpp_b_mat(k, xd, true, Rcpp::seq(0, n-k-1), true);
}

/* Polynomial subspace projection */
Eigen::VectorXd legendre_polynomial(Eigen::VectorXd x, int k, double a,
    double b) {
  ArrayXd xa = 2*(x.array()-a)/(b-a) - 1;
  if (k==0) {
    return VectorXd::Ones(x.size());
  } else if (k==1) {
    return xa.matrix();
  } else if (k==2) {
    return (1.5*xa.pow(2) - 0.5).matrix();
  } else if (k==3) {
    return (2.5*xa.pow(3) - 1.5*xa).matrix();
  } else {
    throw std::invalid_argument("`k` must be 0, 1, 2, or 3.");
  }
}

Eigen::MatrixXd polynomial_basis(const Eigen::VectorXd& x, int k, double a=0,
    double b=1) {
  int n = x.size();
  MatrixXd basis_mat(n, k+1);
  for (int j = 0; j < k+1; j++) {
    basis_mat.col(j) = legendre_polynomial(x, j, a, b);
  }
  return basis_mat;
}

Eigen::VectorXd project_polynomials(const NumericVector& x, const VectorXd& y,
    const ArrayXd& weights, int k) {
  Eigen::ColPivHouseholderQR<MatrixXd> qr;
  ArrayXd sqrt_weights = weights.sqrt();
  VectorXd x_vec(Rcpp::as<Eigen::VectorXd>(x));
  MatrixXd basis_mat = polynomial_basis(x_vec, k, x_vec.minCoeff(),
      x_vec.maxCoeff());
  // If this isn't accurate enough, can also use SVD.
  qr.compute(row_scale(basis_mat, sqrt_weights));
  VectorXd beta = qr.solve((y.array()*sqrt_weights).matrix());
  VectorXd projection = basis_mat*beta;
  if (qr.info() > 0) {
    Rcpp::Rcout << "Eigen QR solve returned nonzero exit status.\n" << std::endl;
  }
  return projection;
}

/* Tridiagonal matrix solve */
Eigen::VectorXd tridiag_forward(const Eigen::VectorXd& a,
    const Eigen::VectorXd& b, const Eigen::VectorXd& c) {
  int n = a.size();
  Eigen::VectorXd cp(n-1);

  // Forward sweep part 1
  cp[0] = c[0]/b[0];
  for (int i=1; i<n-1; i++) {
    cp[i] = c[i]/(b[i] - a[i]*cp[i-1]);
  }
  return cp;
}

// Technically, construction of dp is also part of the forward sweep,
// but it makes sense to include it in "backsolve" since d can change
// over iterations.
Eigen::VectorXd tridiag_backsolve(
    const Eigen::VectorXd& a, const::VectorXd& b,
    const Eigen::VectorXd& cp, const Eigen::VectorXd& d) {
  int n = d.size();
  Eigen::VectorXd dp(n);
  Eigen::VectorXd x(n);

  // Forward sweep part 2
  dp[0] = d[0]/b[0];
  for (int i=1; i<n; i++) {
    dp[i] = (d[i] - a[i]*dp[i-1])/(b[i] - a[i]*cp[i-1]);
  }

  // Backsolve
  x[n-1] = dp[n-1];
  for (int i=n-2; i>=0; i--) {
    x[i] = dp[i] - cp[i]*x[i+1];
  }
  return x;
}

std::tuple<Eigen::VectorXd,Eigen::VectorXd,Eigen::VectorXd> extract_tridiag(
    Eigen::SparseMatrix<double> A) {
  int n = A.cols();
  VectorXd a(n);
  VectorXd b(n);
  VectorXd c(n);
  // Extract diagonal into b
  b = A.diagonal();
  // Extract (-1)-diagonal into a
  a[0] = 0;
  for (int i=1; i<n; i++) {
    a[i] = A.coeff(i, i-1);
  }
  // Extract (+1)-diagonal into c
  for (int i=0; i<n-1; i++) {
    c[i] = A.coeff(i, i+1);
  }
  c[n-1] = 0;
  return std::make_tuple(a, b, c);
}

/* Miscellaneous */
// [[Rcpp::export]]
double get_lambda_max(const NumericVector& x, const Eigen::VectorXd& y,
    const Eigen::ArrayXd& weights, int k) {
  ArrayXd sqrt_weights = weights.sqrt();
  SparseMatrix<double> ck1_mat = get_dk_mat(k+1, x, true);
  SparseQR<SparseMatrix<double>, Ord> qr;
  qr.compute(col_scale(ck1_mat, sqrt_weights.inverse()).transpose());
  VectorXd u_infty = qr.solve((y.array() * sqrt_weights).matrix());
  return u_infty.lpNorm<Eigen::Infinity>();
}

// [[Rcpp::export]]
Eigen::ArrayXd get_lambda_seq(const NumericVector& x, const Eigen::VectorXd& y,
    const Eigen::ArrayXd& weights, int k,
    double lambda_min_ratio=1e-5, int n_lambda=50) {
  if (lambda_min_ratio >= 1) {
    throw std::invalid_argument("lambda_min_ratio must be less than 1.");
  }
  double lambda_max = get_lambda_max(x, y, weights, k);
  double lambda_min = lambda_max * lambda_min_ratio;
  return ArrayXd::LinSpaced(n_lambda,
      std::log(lambda_max), std::log(lambda_min)).exp();
}


