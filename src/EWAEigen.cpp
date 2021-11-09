#include "closspred.h"
#include "progress.h"
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
               double eta, double cumulativeLoss, double loss_tau, bool quiet){
  
  const size_t T=experts.rows();
  const size_t N=experts.cols();
  
  using Row = Eigen::Array<double,1,Eigen::Dynamic>;
  
  auto w0=w0c.array().transpose();

  Row w = Row::Zero(N);
  Row lexp = Row::Zero(N);

  // init progress
  IntegerVector steps;
  if (! quiet) steps = init_progress_cpp(T);
  
  for (size_t t=0 ; t<T ; t++){
    // update progress
    if (! quiet) update_progress_cpp(t+1, steps);
    
    const auto awaket = awake.row(t).array();
    //w = ((eta * w0).exp() * awaket).unaryExpr(std::ptr_fun(truncate1));
    w = ((eta * w0).exp() * awaket).unaryExpr([](double c) {return truncate1(c);});
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
  // end progress
  if (! quiet) end_progress_cpp();
  
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
                      String loss_name,double loss_tau,bool loss_gradient, bool quiet){
    
  
  
  const std::string cs=std::string(loss_name.get_cstring());
  
  if (loss_gradient){
    if (cs=="square" ) return EWAEigen<SquL,true>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau,quiet);
    if (cs=="absolute") return EWAEigen<AbsL,true>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau,quiet);
    if (cs=="percentage") return EWAEigen<PerL,true>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau,quiet);
    if (cs=="log") return EWAEigen<LogL,true>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau,quiet);
    if (cs=="pinball") return EWAEigen<PinL,true>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau,quiet);
  }
  else{
    if (cs=="square" ) return EWAEigen<SquL,false>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau,quiet);
    if (cs=="absolute") return EWAEigen<AbsL,false>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau,quiet);
    if (cs=="percentage") return EWAEigen<PerL,false>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau,quiet);
    if (cs=="log") return EWAEigen<LogL,false>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau,quiet);
    if (cs=="pinball") return EWAEigen<PinL,false>(awake,experts,weights,y,predictions,w0c,eta,cumulativeLoss,loss_tau,quiet);
  }
  
  Rcout << "********** ERROR !!! " << cs << std::endl;
  return 0;
  
}

