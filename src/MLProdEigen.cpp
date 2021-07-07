#include "closspred.h"
#include <RcppEigen.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
struct EtaFunctor{
  double logN_;
  EtaFunctor(size_t N):logN_(std::log(N)){};
  inline double operator() (double a, double b) const {return std::min(1/(2*a),sqrt(logN_/b));}
};


template <class LT, bool G>
void MLProdEigen( Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> eta, 
                 Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights,
                 Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> predictions, 
                 Eigen::Map<Eigen::VectorXd> Rc, Eigen::Map<Eigen::VectorXd> Lc,
                 double dmaxloss,double loss_tau){
  
  const size_t T=experts.rows();
  const size_t N=experts.cols();
  
  auto R=Rc.array().transpose();
  auto L=Lc.array().transpose();
  
  using Row = Eigen::Array<double,1,Eigen::Dynamic>;
  
  Row w = Row::Zero(N);
  Row p = Row::Zero(N);
  Row r = Row::Zero(N);
  Row lexp = Row::Zero(N);
  Row neweta = Row::Zero(N);
  Row maxloss = Row::Ones(N)*dmaxloss;
  
  // init progress
  IntegerVector steps = init_progress_cpp(T); 
  
  for (size_t t=0 ; t<T ; t++){
    // update progress
    update_progress_cpp(t+1, steps);
    
    w = R.exp().unaryExpr(std::ptr_fun(truncate1));
    w *= eta.row(t).array() ;
    w /= w.sum();
    const auto awaket = awake.row(t).array();
    p = awaket * w.array();

    const double pred = experts.row(t) * p.matrix().transpose();//scalar product

    weights.row(t).array() = p;
    predictions[t] = pred;

    // Observe losses
    const auto lpf=LossPredFunctor<LT,G>(y[t],pred,loss_tau);
    const double lpred=lpf(pred);
    lexp = experts.row(t).array().unaryExpr(lpf);
    //Update the regret and the weight
    r = awaket * (lpred-lexp);
    L.array()+=(r*r);
    
    maxloss = maxloss.binaryExpr(r,AbsMax());
    neweta = maxloss.binaryExpr(L.array(),EtaFunctor(N));
    
    R = neweta*eta.row(t).array().inverse()*R+(1.0 + neweta*r).log();
    eta.row(t+1).array() = neweta; 
    
    if (std::isnan(R.sum())) Rcout << "Nan in R" << std::endl;
    
  }
  // end progress
  end_progress_cpp();
  
  return ;
}


// [[Rcpp::export]]
void computeMLProdEigen(Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> eta, 
                        Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights,
                        Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> predictions, 
                        Eigen::Map<Eigen::VectorXd> R, Eigen::Map<Eigen::VectorXd> L,
                        double maxloss,
                        String loss_name,double loss_tau,bool loss_gradient){
  
  
  const std::string cs=std::string(loss_name.get_cstring());
  
  if (loss_gradient){
    if (cs=="square" ) return MLProdEigen<SquL,true>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau);
    if (cs=="absolute") return MLProdEigen<AbsL,true>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau);
    if (cs=="percentage") return MLProdEigen<PerL,true>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau);
    if (cs=="log") return MLProdEigen<LogL,true>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau);
    if (cs=="pinball") return MLProdEigen<PinL,true>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau);
  }
  else{
    if (cs=="square" ) return MLProdEigen<SquL,false>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau);
    if (cs=="absolute") return MLProdEigen<AbsL,false>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau);
    if (cs=="percentage") return MLProdEigen<PerL,false>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau);
    if (cs=="log") return MLProdEigen<LogL,false>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau);
    if (cs=="pinball") return MLProdEigen<PinL,false>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau);
  }
  
  Rcout << "********** ERROR !!! " << cs << std::endl;
  return ;
  
}

