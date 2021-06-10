#include "closspred.h"
#include <RcppEigen.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]



template <class LT, bool G>
double EWAEigen( Eigen::Map<Eigen::MatrixXd> awake, 
               Eigen::Map<Eigen::MatrixXd> experts, 
               Eigen::Map<Eigen::MatrixXd> weights,
               Eigen::Map<Eigen::VectorXd> y, 
               Eigen::Map<Eigen::VectorXd> predictions, 
               Eigen::Map<Eigen::VectorXd> w0c,
               double eta, double cumulativeLoss, double loss_tau){
  
  const size_t T=experts.rows();
  const size_t N=experts.cols();
  
  using Row = Eigen::Array<double,1,Eigen::Dynamic>;
  
  auto w0=w0c.array().transpose();

  Row w = Row::Zero(N);
  Row lexp = Row::Zero(N);

  
  for (size_t t=0 ; t<T ; t++){
    const auto awaket = awake.row(t).array();
    w = ((eta * w0).exp() * awaket).unaryExpr(std::ptr_fun(truncate1));
    w /=w.sum();
    weights.row(t).array() = w;

    const double pred = experts.row(t) * w.matrix().transpose();//scalar product 
    predictions[t] = pred;
    
    const auto lf=LossPredFunctor<LT,false>(y[t],pred,loss_tau);
    cumulativeLoss += lf(pred);

    const auto lpf=LossPredFunctor<LT,G>(y[t],pred,loss_tau);
    const double lpred=lpf(pred);
    lexp = experts.row(t).array().unaryExpr(lpf);  
    
    w0 += awaket * (lpred-lexp);
  }
  
  return cumulativeLoss;
}
  

// [[Rcpp::export]]
double computeEWAEigen( Eigen::Map<Eigen::MatrixXd> awake, 
                      Eigen::Map<Eigen::MatrixXd> experts, 
                      Eigen::Map<Eigen::MatrixXd> weights,
                      Eigen::Map<Eigen::VectorXd> y, 
                      Eigen::Map<Eigen::VectorXd> predictions, 
                      Eigen::Map<Eigen::VectorXd> w0c,
                      double eta, double cumulativeLoss,
                      String loss_name,double loss_tau,bool loss_gradient ){
    
  
  
  const std::string cs=std::string(loss_name.get_cstring());
  
  if (loss_gradient){
    if (cs=="square" ) return EWAEigen<SquL,true>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau);
    if (cs=="absolute") return EWAEigen<AbsL,true>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau);
    if (cs=="percentage") return EWAEigen<PerL,true>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau);
    if (cs=="log") return EWAEigen<LogL,true>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau);
    if (cs=="pinball") return EWAEigen<PinL,true>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau);
  }
  else{
    if (cs=="square" ) return EWAEigen<SquL,false>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau);
    if (cs=="absolute") return EWAEigen<AbsL,false>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau);
    if (cs=="percentage") return EWAEigen<PerL,false>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau);
    if (cs=="log") return EWAEigen<LogL,false>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau);
    if (cs=="pinball") return EWAEigen<PinL,false>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau);
  }
  
  Rcout << "********** ERROR !!! " << cs << std::endl;
  return 0;
  
}

