#include "progress.h"
#include <RcppEigen.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]  
size_t computeRidgeCPP(Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> w, 
                     Eigen::Map<Eigen::MatrixXd> At, Eigen::Map<Eigen::MatrixXd> bt, 
                     Eigen::Map<Eigen::VectorXd> y, bool quiet){
  const size_t T=experts.rows();
  const size_t N=experts.cols();
  

  auto tw=w.transpose();
  auto texperts=experts.transpose();

  Eigen::VectorXd a(N);
    
  // init progress
  IntegerVector steps;
  if (! quiet) steps = init_progress_cpp(T);
  
  for (size_t t=0 ; t<T ; t++){ 
    // update progress
    if (! quiet) update_progress_cpp(t+1, steps);
    
    tw.col(t) = At * bt;
    a =  At * texperts.col(t);
    At+= - a * a.transpose() / (1 + experts.row(t) * a);
    bt+= y[t] * texperts.col(t);
    
  }
  // end progress
  if (! quiet) end_progress_cpp();
  
  return 0;
}


  

