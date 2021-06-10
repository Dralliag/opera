#include <RcppEigen.h>
using namespace Rcpp;
#include <limits>
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
size_t RidgeCalibStep1(size_t tp1,double dbestlambda,
    Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights,
    Eigen::Map<Eigen::MatrixXd> wlambda, Eigen::Map<Eigen::MatrixXd> w0,
    Eigen::Map<Eigen::MatrixXd> At, Eigen::Map<Eigen::MatrixXd> bt, 
    Eigen::Map<Eigen::MatrixXd> gridlambda, Eigen::Map<Eigen::MatrixXd> predlambda, 
    Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> lambda,
    Eigen::Map<Eigen::VectorXd> cumulativeloss, Eigen::Map<Eigen::VectorXd> prediction){
  
  const size_t t=tp1-1;
  size_t bestlambda=std::floor(dbestlambda)-1;
  
  const size_t T=experts.rows();
  const size_t N=experts.cols();
  const size_t nlambda=wlambda.cols();
  
  weights.row(t).array() = wlambda.col(bestlambda);
  prediction[t] = experts.row(t) * weights.row(t).transpose();
  lambda[t] = gridlambda(bestlambda,0);
    
  predlambda.row(t).array() = experts.row(t) * wlambda;
  cumulativeloss.array() += (predlambda.row(t).array()-y[t]).square();

  At+=experts.transpose().col(t) * experts.row(t);
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

// [[Rcpp::export]]
void RidgeCalibStep2(
    Eigen::Map<Eigen::MatrixXd> wlambda, Eigen::Map<Eigen::VectorXd> w0,
    Eigen::Map<Eigen::MatrixXd> At, Eigen::Map<Eigen::MatrixXd> bt, 
    Eigen::Map<Eigen::VectorXd> gridlambda){
  
  const size_t N=wlambda.rows();
  const size_t nlambda=wlambda.cols();
  
  //LP : this is probably costly...
  Eigen::MatrixXd cAt(At);
  Eigen::MatrixXd cbt(bt);
  
  //LP : this is probably costly...
  Eigen::FullPivLU<Eigen::MatrixXd> lu(At);
 //Eigen::PartialPivLU<Eigen::MatrixXd> lu(At);
 
  for (size_t k=0 ; k<nlambda ; k++){
    const double diagA=gridlambda[k];
    for (size_t i=0 ; i<N ; i++){
        cAt(i,i)=At(i,i)+diagA;
        cbt(i,0)=bt(i,0)+diagA*w0[i];
    }
    
    lu.compute(cAt);
    if (lu.isInvertible()){
      wlambda.col(k)=lu.solve(cbt);
    }
    else{
      wlambda.col(k).array()=std::numeric_limits<double>::quiet_NaN();
    }
  }
  
  return ;
} 
  


  

