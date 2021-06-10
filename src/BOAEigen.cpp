#include "closspred.h"
#include <RcppEigen.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]


struct EtaBOAFunctor{
  double logN_;
  EtaBOAFunctor(size_t N):logN_(std::log(N)){};
  inline double operator() (double a, double b) const {return std::min(std::min(1/(2*a),sqrt(logN_/b)),std::exp(350.));}
};


template <class LT, bool G>
void BOAEigen( Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> eta, 
                 Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights,
                 Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> predictions, 
                 Eigen::Map<Eigen::VectorXd> wc, Eigen::Map<Eigen::VectorXd> w0c,
                 Eigen::Map<Eigen::VectorXd> Rc, Eigen::Map<Eigen::VectorXd> Regc,
                 Eigen::Map<Eigen::VectorXd> Bc, Eigen::Map<Eigen::VectorXd> Vc,
                 double loss_tau){
  
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
  
  for (size_t t=0 ; t<T ; t++){
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
    eta.row(t+1).array() = B.binaryExpr(V,EtaBOAFunctor(N));
    if (eta.row(t+1).maxCoeff() > std::exp(300.)) {
      reg = r;
    }
    else{
      // if some losses still have not been observed
      reg = r - eta.row(t+1).array() * r*r;
    }
    // Update the regret and the regularized regret used by BOA
    R+=r;
    Reg+=reg;
    
    w = (w0.log()+eta.row(t+1).array()*Reg).exp().unaryExpr(std::ptr_fun(truncate1));
  }
  
  return ;
}


// [[Rcpp::export]]
void computeBOAEigen( Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> eta, 
                         Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights,
                         Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> predictions, 
                         Eigen::Map<Eigen::VectorXd> wc, Eigen::Map<Eigen::VectorXd> w0c,
                         Eigen::Map<Eigen::VectorXd> Rc, Eigen::Map<Eigen::VectorXd> Regc,
                         Eigen::Map<Eigen::VectorXd> Bc, Eigen::Map<Eigen::VectorXd> Vc,
                         String loss_name,double loss_tau,bool loss_gradient){
  
  
  const std::string cs=std::string(loss_name.get_cstring());
  
  if (loss_gradient){
    if (cs=="square" ) return BOAEigen<SquL,true>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau);
    if (cs=="absolute") return BOAEigen<AbsL,true>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau);
    if (cs=="percentage") return BOAEigen<PerL,true>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau);
    if (cs=="log") return BOAEigen<LogL,true>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau);
    if (cs=="pinball") return BOAEigen<PinL,true>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau);
  }
  else{
    if (cs=="square" ) return BOAEigen<SquL,false>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau);
    if (cs=="absolute") return BOAEigen<AbsL,false>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau);
    if (cs=="percentage") return BOAEigen<PerL,false>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau);
    if (cs=="log") return BOAEigen<LogL,false>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau);
    if (cs=="pinball") return BOAEigen<PinL,false>(awake,eta,experts,weights,y,predictions,wc,w0c,Rc,Regc,Bc,Vc,loss_tau);
  }
  
  Rcout << "********** ERROR !!! " << cs << std::endl;
  return ;
  
}
