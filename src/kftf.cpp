#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <dspline.h>
#include "utils.h"
#include "kftf.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(dspline)]]

using Eigen::ArrayXd;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::VectorXd;
using Rcpp::_;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::Named;
using Rcpp::NumericVector;


// Always prefix the types in function signatures
void f1step(double y,
            const Eigen::RowVectorXd& Z,
            const Eigen::MatrixXd& ZZ,
            double H,
            const Eigen::MatrixXd& A,
            const Eigen::MatrixXd& RQR,
            Eigen::VectorXd& a,
            Eigen::MatrixXd& P,
            double& vt,
            double& Ft,
            Eigen::VectorXd& Kt) {

  MatrixXd Ptemp = P;
  vt = y - Z * a;
  Ft = Z * P * Z.transpose() + H;
  Kt = A * P * Z.transpose() / Ft;
  a = A * a + Kt * vt;
  Ptemp -= P * ZZ * P / Ft;
  P = A * Ptemp * A.transpose() + RQR;

  // symmeterize
  Ptemp = P;
  Ptemp += P.transpose();
  P = Ptemp / 2;

  // Some entries of P can be _really_ small in magnitude, set them to 0
  // but maintain symmetric / posdef
  for (int i = 0; i < P.rows(); i++) {
   for (int j = 0; j <= i; j++) {
     if (abs(P(i,j)) < 1e-30) {
       if (i == j) {
         P.row(i).setZero();
         P.col(i).setZero();
       } else {
         P(i,j) = 0;
         P(j,i) = 0;
      }
     }
   }
  }
}

// [[Rcpp::export]]
List f1step_test(double y,
                 Eigen::RowVectorXd Z,
                 double H,
                 Eigen::MatrixXd A,
                 Eigen::MatrixXd RQR,
                 Eigen::VectorXd a,
                 Eigen::MatrixXd P) {
  MatrixXd ZZ = Z.transpose() * Z;
  double vt;
  double Ft;
  VectorXd Kt(a);
  f1step(y, Z, ZZ, H, A, RQR, a, P, vt, Ft, Kt);

  List list = List::create(Named("a") = a, _["P"] = P, _["vt"] = vt,
                           _["Ft"] = Ft, _["Kt"] = Kt);
  return list;
}

void df1step(double y,
             const Eigen::RowVectorXd& Z,
             double H,
             const Eigen::MatrixXd& A,
             const Eigen::MatrixXd& RQR,
             Eigen::VectorXd& a,
             Eigen::MatrixXd& P,
             Eigen::MatrixXd& Pinf,
             int& rankp,
             double& vt,
             double& Ft,
             double& Finf) {
  double tol = 1.490116e-08;  // sqrt(.Machine$double.eps);
  int k = a.size();
  MatrixXd Ptemp(k, k);

  VectorXd Mt = P * Z.transpose();
  Ft = Z * P * Z.transpose() + H;
  VectorXd Minf = Pinf * Z.transpose();
  Finf = Z * Minf;
  vt = y - Z * a;

  double finv;
  if (Finf > tol) {  // should always happen
    finv = 1 / Finf;
    a += vt * finv * Minf;
    P += Ft * pow(finv, 2) * Minf * Minf.transpose();
    P -= finv * Mt * Minf.transpose() + finv * Minf * Mt.transpose();
    Pinf -= finv * Minf * Minf.transpose();
    rankp--;
  } else {  // should never happen
    Finf = 0;
    if (Ft > tol) {
      finv = 1 / Ft;
      a += vt * finv * Mt;
      P -= finv * Mt * Mt.transpose();
    }
  }
  if (Ft < tol)
    Ft = 0;

  a = A * a;
  P = A * P * A.transpose() + RQR;
  Pinf = A * Pinf * A.transpose();

  // symmeterize
  Ptemp = P;
  Ptemp += P.transpose();
  P = Ptemp / 2;
  Ptemp = Pinf;
  Ptemp += Pinf.transpose();
  Pinf = Ptemp / 2;


  // Fix possible negative definiteness, should never happen
  for (int i = 0; i < k; i++) {
    if (P(i,i) < 0) {
      P.row(i).setZero();
      P.col(i).setZero();
    }
    if (Pinf(i,i) < 0) {
      Pinf.row(i).setZero();
      Pinf.col(i).setZero();
    }
  }
}


// [[Rcpp::export]]
List df1step_test(double y,
                  Eigen::RowVectorXd Z,
                  double H,
                  Eigen::MatrixXd A,
                  Eigen::MatrixXd RQR,
                  Eigen::VectorXd a,
                  Eigen::MatrixXd P,
                  Eigen::MatrixXd Pinf,
                  int rankp) {
  double vt;
  double Ft;
  double Finf;
  df1step(y, Z, H, A, RQR, a, P, Pinf, rankp, vt, Ft, Finf);

  List out = List::create(
    Named("vt") = vt, _["Ft"] = Ft, _["Finf"] = Finf, _["a"] = a,
                      _["P"] = P, _["Pinf"] = Pinf, _["rankp"] = rankp);
  return out;
}

void kftfcpp(const Eigen::VectorXd& y,
             int k,
             double lambda,
             Eigen::VectorXd& theta) {

  // define transition matrix A
  MatrixXd A = MatrixXd::Zero(k, k);
  // For uneven spacing, we'll have to extract or rebuild d_mat/b_mat
  // The below is not ideal
  // IntegerVector xi = Rcpp::seq_len(k + 1);
  // NumericVector xd = Rcpp::as<NumericVector>(xi);
  // Eigen::SparseMatrix Dseq = dspline::rcpp_b_mat(k, xd, false, k, false);

  for (int i = 0; i < k - 1; i++) {
    A(i + 1, i) = 1;
  }
  double s = -1;
  for (int i = 0; i < k; i++) {
    s *= -1;
    A(0, i) = Rf_choose(k, i + 1) * s;
  }


  int n = y.size();
  RowVectorXd Z = RowVectorXd::Zero(k);
  Z(0) = 1;
  MatrixXd ZZ = Z.transpose() * Z;
  VectorXd H = VectorXd::Ones(n);
  VectorXd R = VectorXd::Zero(k);
  R(0) = 1;
  double Q = 1 / lambda;
  MatrixXd RQR = R * R.transpose() * Q;
  // Rcpp::Rcout << "RQR = " << RQR << std::endl;
  VectorXd a1 = VectorXd::Zero(k);
  MatrixXd at = MatrixXd::Zero(k, n + 1);
  MatrixXd P1 = MatrixXd::Zero(k, k);
  MatrixXd Pt = MatrixXd::Zero(k * k, n + 1);
  MatrixXd P1inf = MatrixXd::Identity(k, k);
  MatrixXd Pinf = MatrixXd::Zero(k * k, n + 1);
  // matrix.reshaped() is in Eigen 3.4.0 (2+ yrs old), but not in current
  // RcppEigen (see https://github.com/RcppCore/RcppEigen/issues/103)
  // but we can seemingly convert to array and assign back to matrix
  Pinf.col(0) = Map<VectorXd>(P1inf.data(), k*k);
  // save first row of each P1 & P1inf for easier computation of theta:
  // MatrixXd Pt_res = MatrixXd::Zero(k, n + 1);    // for P1
  // MatrixXd Pinf_res = MatrixXd::Zero(k, n + 1);  // for P1inf
  // Pinf_res.col(0) = P1inf.row(0);


  // forward
  int d = 0;
  int rankp = k;
  double vt_b = 0.0;
  double Ft_b = 0.0;
  double Finf_b = 0.0;
  VectorXd vt = VectorXd::Zero(n);
  VectorXd Ft = VectorXd::Zero(n);
  VectorXd Finf = VectorXd::Zero(n);
  MatrixXd Kt = MatrixXd::Zero(k, n);
  MatrixXd Kinf = MatrixXd::Zero(k, n);
  VectorXd Kt_b = VectorXd::Zero(k);


  while (rankp > 0 && d < n) {
    df1step(y(d), Z, H(d), A, RQR, a1, P1, P1inf, rankp, vt_b, Ft_b, Finf_b);
    at.col(d + 1) = a1;
    vt(d) = vt_b;
    Ft(d) = Ft_b;
    Finf(d) = Finf_b;
    Kt.col(d) = A * P1 * Z.transpose();
    Kinf.col(d) = A * P1inf * Z.transpose();
    Pt.col(d + 1) = Map<VectorXd>(P1.data(), k*k);
    Pinf.col(d + 1) = Map<VectorXd>(P1inf.data(), k*k);
    // Rcpp::Rcout << "d = " << d << std::endl;
    d++;
  }

  // Rcpp::Rcout << "P1 = " << P1 << std::endl;
  // Rcpp::Rcout << "a1 = " << a1 << std::endl;


  for (int i = d; i < n; i++) {
    f1step(y(i), Z, ZZ, H(i), A, RQR, a1, P1, vt_b, Ft_b, Kt_b);
    vt(i) = vt_b;
    Ft(i) = Ft_b;
    Kt.col(i) = Kt_b;
    at.col(i + 1) = a1;
    Pt.col(i + 1) = Map<VectorXd>(P1.data(), k*k);
  }


  // backward
  VectorXd r = VectorXd::Zero(k);
  VectorXd r1 = VectorXd::Zero(k);
  VectorXd rtmp = VectorXd::Zero(k);
  MatrixXd L0 = MatrixXd::Zero(k, k);
  MatrixXd L1 = MatrixXd::Zero(k, k);

  for (int i = n - 1; i >= d; i--) {
    // Rcpp::Rcout << "i = " << i << std::endl;
    L0 = A - Kt.col(i) * Z;
    r1 = Z.transpose() * vt(i) / Ft(i) + L0.transpose() * r;
    r = r1;
    P1 = Map<MatrixXd>(Pt.col(i).data(), k, k);
    theta(i) = at.col(i)(0) + P1.row(0) * r;
  }


  for (int i = d - 1; i >= 0; i--) {
    // Rcpp::Rcout << "i = " << i << std::endl;
    P1 = Map<MatrixXd>(Pt.col(i).data(), k, k);
    P1inf = Map<MatrixXd>(Pinf.col(i).data(), k, k);
    L0 = A - A * P1inf * ZZ / Finf(i);
    L1 = A * P1inf * ZZ * Ft(i) / pow(Finf(i), 2);
    L1 -= A * P1 * ZZ / Finf(i);
    r1 = Z.transpose() * vt(i) / Finf(i);
    rtmp = L0.transpose() * r1 + L1.transpose() * r;
    r1 += rtmp;
    rtmp = L0.transpose() * r;
    r = rtmp;
    theta(i) = at.col(i)(0) + P1.row(0) * r + P1inf.row(0) * r1;
  }
}

// [[Rcpp::export]]
Eigen::VectorXd kftf_cpp_test(Eigen::VectorXd y, int k, double lambda) {
  int n = y.size();
  VectorXd theta(n);
  kftfcpp(y, k, lambda, theta);

  return theta;
}