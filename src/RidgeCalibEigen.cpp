#include <RcppEigen.h>
using namespace Rcpp;
#include <limits>
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
size_t RidgeCalibStep1(size_t tp1,double dbestlambda,
    Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights,
    Eigen::Map<Eigen::MatrixXd> wlambda, Eigen::Map<Eigen::MatrixXd> w0,
    Eigen::Map<Eigen::MatrixXd> bt, 
    Eigen::Map<Eigen::MatrixXd> gridlambda,
    Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> lambda,
    Eigen::Map<Eigen::VectorXd> cumulativeloss, Eigen::Map<Eigen::VectorXd> prediction){
  
  const size_t t=tp1-1;
  size_t bestlambda=std::floor(dbestlambda)-1;
  
  const size_t nlambda=wlambda.cols();
  
  weights.row(t).array() = wlambda.col(bestlambda);
  prediction[t] = experts.row(t) * weights.row(t).transpose();
  lambda[t] = gridlambda(bestlambda,0);
    
  cumulativeloss.array() += ((experts.row(t) * wlambda).array() - y[t]).square();

  bt+=y[t]*experts.transpose().col(t);

  bestlambda=0;
  double mincl=std::numeric_limits<double>::max();
  for (size_t l=0 ; l<nlambda ; l++){
    if (mincl>cumulativeloss[l]){
      bestlambda=l;
      mincl=cumulativeloss[l];
    }
  }
  return bestlambda+1;
}

  


  

