#include "closspred.h"
#include <RcppEigen.h>
#include "progress.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

struct EtaBOAFunctor1{
  inline double operator() (double a, double b) const {return sqrt(std::log(1/a)/b);} 
};

struct EtaBOAFunctor2{
  inline double operator() (double a, double b) const {return std::min(1/a,b);}
};

template <class LT, bool G>
void BOAEigen( Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> eta, 
                 Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights,
                 Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> predictions, 
                 Eigen::Map<Eigen::VectorXd> wc, Eigen::Map<Eigen::VectorXd> w0c,
                 Eigen::Map<Eigen::VectorXd> Rc, Eigen::Map<Eigen::VectorXd> Regc,
                 Eigen::Map<Eigen::VectorXd> Bc, Eigen::Map<Eigen::VectorXd> Vc,
                 double loss_tau, bool quiet){
  
  const size_t T=experts.rows();
  const size_t N=experts.cols();
  
  using Row = Eigen::Array<double,1,Eigen::Dynamic>;
  
  auto w0=w0c.array().transpose();
  auto w=wc.array().transpose();
  auto R=Rc.array().transpose();
  auto Reg=Regc.array().transpose();
  auto B=Bc.array().transpose();
  auto V=Vc.array().transpose();
  
  Row p = Row::Zero(N);
  Row r = Row::Zero(N);
  Row lexp = Row::Zero(N);
  Row reg = Row::Zero(N);
  
  // init progress
  IntegerVector steps;
  if (! quiet) steps = init_progress_cpp(T);
  
  for (size_t t=0 ; t<T ; t++){
    // update progress
    if (! quiet) update_progress_cpp(t+1, steps);
    
    auto awaket = awake.row(t).array();
    p = awaket * w.array();
    p /= p.sum();
    const double pred = experts.row(t) * p.matrix().transpose();//scalar product
    weights.row(t).array() = p;
    predictions[t] = pred;
    
    // Observe losses
    const auto lpf=LossPredFunctor<LT,G>(y[t],pred,loss_tau);
    const double lpred=lpf(pred);
    lexp = experts.row(t).array().unaryExpr(lpf);
    
    // Instantaneous regret
    r = awaket * (lpred-lexp);
    // Update the learning rates
    B = B.binaryExpr(r,AbsMax());

    V += r*r;
    eta.row(t+1).array() = (B.unaryExpr([](double c) {return sq_log2(c);})).binaryExpr(w0.binaryExpr(V, EtaBOAFunctor1()), EtaBOAFunctor2());
    // Update the regret and the regularized regret used by BOA
    reg = 0.5 * (r - eta.row(t+1).array() * r*r + (B.unaryExpr([](double c) {return sq_log2(c);})) * (eta.row(t+1).array() * r).unaryExpr([](double c) {return sup_half(c);}));
    
    R+=r;
    Reg+=reg;
    //w = (eta.row(t+1).array().log()+w0.log()+eta.row(t+1).array()*Reg).exp().unaryExpr(std::ptr_fun(truncate1));
    w = (eta.row(t+1).array().log()+w0.log()+eta.row(t+1).array()*Reg).exp().unaryExpr([](double c) {return truncate1(c);});
    
  }
  // end progress
  if (! quiet) end_progress_cpp();
  
  return ;
}


// [[Rcpp::export]]
void computeBOAEigen( Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> eta, 
                      Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights,
                      Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> predictions, 
                      Eigen::Map<Eigen::VectorXd> wc, Eigen::Map<Eigen::VectorXd> w0c,
                      Eigen::Map<Eigen::VectorXd> Rc, Eigen::Map<Eigen::VectorXd> Regc,
                      Eigen::Map<Eigen::VectorXd> Bc, Eigen::Map<Eigen::VectorXd> Vc,
                      String loss_name,double loss_tau,bool loss_gradient, bool quiet){
  
  
  const std::string cs=std::string(loss_name.get_cstring());
  
  if (loss_gradient){
    if (cs=="square" ) return BOAEigen<SquL,true>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau,quiet);
    if (cs=="absolute") return BOAEigen<AbsL,true>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau,quiet);
    if (cs=="percentage") return BOAEigen<PerL,true>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau,quiet);
    if (cs=="log") return BOAEigen<LogL,true>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau,quiet);
    if (cs=="pinball") return BOAEigen<PinL,true>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau,quiet);
  }
  else{
    if (cs=="square" ) return BOAEigen<SquL,false>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau,quiet);
    if (cs=="absolute") return BOAEigen<AbsL,false>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau,quiet);
    if (cs=="percentage") return BOAEigen<PerL,false>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau,quiet);
    if (cs=="log") return BOAEigen<LogL,false>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau,quiet);
    if (cs=="pinball") return BOAEigen<PinL,false>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau,quiet);
  }
  
  Rcout << "********** ERROR !!! " << cs << std::endl;
  return ;
  
}

