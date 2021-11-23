#include "closspred.h"
#include "progress.h"
#include <RcppEigen.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]


template <class L, bool G>
double MLPolEigen( Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> eta, 
                 Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights,
                 Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> predictions, 
                 Eigen::Map<Eigen::VectorXd> Rc, Eigen::Map<Eigen::VectorXd> wc,
                 double B,double loss_tau, bool quiet){
  
  const size_t T=experts.rows();
  const size_t N=experts.cols();
  
  using Row = Eigen::Array<double,1,Eigen::Dynamic>;
  
  auto R=Rc.array().transpose();
  auto w=wc.array().transpose();
  
  Row p = Row::Zero(N);
  Row r = Row::Zero(N);
  Row lexp = Row::Zero(N);
  
  // init progress
  IntegerVector steps;
  if (! quiet) steps = init_progress_cpp(T);
  
  for (size_t t=0 ; t<T ; t++){
    // update progress
    if (! quiet) update_progress_cpp(t+1, steps);
    
    auto awaket = awake.row(t).array();
    
    // We check if there is at least one expert with positive weight
    if ((awaket*R).maxCoeff()>0){
      //w = eta.row(t).array() * R.unaryExpr(std::ptr_fun(ramp));
      w = eta.row(t).array() * R.unaryExpr([](double c) {return ramp(c);});
     // w /= w.sum(); // Necessary ?
    }
    else{
      w = 1.0;
    }
    // form the mixture and the prediction
    p = awaket * w.array();
    p = awake.row(t).array() * w.array();
    p /= p.sum();
    const double pred = (experts.row(t).array() * p).sum();
    //const double pred = experts.row(t) * p.matrix().transpose();//scalar product
    
    //save the mixture and the prediction
    weights.row(t).array() = p;
    predictions[t] = pred;

    // Observe losses
    const auto lpf=LossPredFunctor<L,G>(y[t],pred,loss_tau);
    const double lpred=lpf(pred);
    lexp = experts.row(t).array().unaryExpr(lpf);
    //Update the regret and the weight
    r = awaket * (lpred-lexp);
    R += r;
    //Update the learning rate
    const double newB = std::max(B,(r*r).maxCoeff());

    eta.row(t+1).array() = (eta.row(t).array().inverse() + r*r + newB -B).inverse();
    B = newB;
  }
  // end progress
  if (! quiet) end_progress_cpp();
  
  if (R.maxCoeff()>0){
    //w = eta.row(T).array() * R.unaryExpr(std::ptr_fun(ramp));
    w = eta.row(T).array() * R.unaryExpr([](double c) {return ramp(c);});
    w /= w.sum();
  }
  else{
    w = 1.0/double(N);
  }
  
  
  return B;
}


// [[Rcpp::export]]
double computeMLPolEigen( Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> eta, 
                          Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights,
                          Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> predictions, 
                          Eigen::Map<Eigen::VectorXd> R, Eigen::Map<Eigen::VectorXd> w,
                      double B,
                      String loss_name,double loss_tau,bool loss_gradient, bool quiet){
  
  
  const std::string cs=std::string(loss_name.get_cstring());
  
  
  if (loss_gradient){
    if (cs=="square" ) return MLPolEigen<SquL,true>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
    if (cs=="absolute") return MLPolEigen<AbsL,true>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
    if (cs=="percentage") return MLPolEigen<PerL,true>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
    if (cs=="log") return MLPolEigen<LogL,true>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
    if (cs=="pinball") return MLPolEigen<PinL,true>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
  }
  else{
    if (cs=="square" ) return MLPolEigen<SquL,false>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
    if (cs=="absolute") return MLPolEigen<AbsL,false>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
    if (cs=="percentage") return MLPolEigen<PerL,false>(awake,eta,experts,weights,y,predictions,w,R,B,loss_tau,quiet);
    if (cs=="log") return MLPolEigen<LogL,false>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
    if (cs=="pinball") return MLPolEigen<PinL,false>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
  }
  
  Rcout << "********** ERROR !!! " << cs << std::endl;
  return 0;
  
}

