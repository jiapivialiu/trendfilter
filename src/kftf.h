#ifndef __KFTF_H
#define __KFTF_H


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
            Eigen::VectorXd& Kt);

Rcpp::List f1step_test(double y,
                       Eigen::RowVectorXd Z,
                       double H,
                       Eigen::MatrixXd A,
                       Eigen::MatrixXd RQR,
                       Eigen::VectorXd a,
                       Eigen::MatrixXd P);
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
             double& Finf);

Rcpp::List df1step_test(double y,
                        Eigen::RowVectorXd Z,
                        double H,
                        Eigen::MatrixXd A,
                        Eigen::MatrixXd RQR,
                        Eigen::VectorXd a,
                        Eigen::MatrixXd P,
                        Eigen::MatrixXd Pinf,
                        int rankp);

void kftfcpp(const Eigen::VectorXd& y,
             int k,
             double lambda,
             Eigen::VectorXd& theta);

Eigen::VectorXd kftf_cpp_test(Eigen::VectorXd y,
                              int k,
                              double lambda);

#endif
