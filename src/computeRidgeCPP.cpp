#include <RcppEigen.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]  
size_t computeRidgeCPP(Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> w, 
                     Eigen::Map<Eigen::MatrixXd> At, Eigen::Map<Eigen::MatrixXd> bt, 
                     Eigen::Map<Eigen::VectorXd> y){
  const size_t T=experts.rows();
  const size_t N=experts.cols();
  

  auto tw=w.transpose();
  auto texperts=experts.transpose();
  
  //Pre-allocate LU matrix storage
  Eigen::FullPivLU<Eigen::MatrixXd> lu(At);
  //Eigen::PartialPivLU<Eigen::MatrixXd> lu(At);

  for (size_t t=0 ; t<T ; t++){ 
    lu.compute(At); //compute factorization (no allocation)
    tw.col(t)=lu.solve(bt);
    At+=texperts.col(t) * experts.row(t);
    bt+=y[t]*texperts.col(t);
  }
  
  return 0;
}


  

