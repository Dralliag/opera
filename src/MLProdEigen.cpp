#include "closspred.h"
#include "progress.h"
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
                 Eigen::Map<Eigen::VectorXd> dmaxloss,double loss_tau, bool quiet){
  
  const size_t T=experts.rows();
  const size_t N=experts.cols();
  
  auto R=Rc.array().transpose();
  auto L=Lc.array().transpose();
  
  auto maxloss=dmaxloss.array().transpose();
  
  using Row = Eigen::Array<double,1,Eigen::Dynamic>;
  
  Row w = Row::Zero(N);
  Row p = Row::Zero(N);
  Row r = Row::Zero(N);
  Row lexp = Row::Zero(N);
  Row neweta = Row::Zero(N);
  //Row maxloss = Row::Ones(N)*dmaxloss;
  
  // init progress
  IntegerVector steps;
  if (! quiet) steps = init_progress_cpp(T); 
  
  for (size_t t=0 ; t<T ; t++){
    // update progress
    if (! quiet) update_progress_cpp(t+1, steps);
    
    //Rcout << "The value of R : " << R << "\n";
    //Rcout << "The value of eta : " << eta << "\n";
    
    const auto awaket = awake.row(t).array();
    
    //w = R.exp().unaryExpr(std::ptr_fun(truncate1)) * awaket;
    w = R.exp().unaryExpr([](double c) {return truncate1(c);}) * awaket;
    w *= eta.row(t).array() ;

    w /= w.sum();
    
    //Rcout << "The value of w : " << w << "\n";
  
    p = awaket * w.array();

    const double pred = experts.row(t) * p.matrix().transpose();//scalar product

    //Rcout << "The value of p : " << p << "\n";
    //Rcout << "The value of pred : " << pred << "\n";
    
    weights.row(t).array() = p;
    predictions[t] = pred;

    // Observe losses
    const auto lpf=LossPredFunctor<LT,G>(y[t],pred,loss_tau);
    const double lpred=lpf(pred);
    
    //Rcout << "The value of lpred : " << lpred << "\n";
    
    lexp = experts.row(t).array().unaryExpr(lpf);
    //Update the regret and the weight
    r = awaket * (lpred-lexp);
    L.array()+=(r*r);
    
    //Rcout << "The value of lexp : " << lexp << "\n";
    //Rcout << "The value of r : " << r << "\n";
    
    //Rcout << "The value of maxloss before : " << maxloss << "\n";
    
    maxloss = maxloss.binaryExpr(r,AbsMax());
    neweta = maxloss.binaryExpr(L.array(),EtaFunctor(N));
    
    //Rcout << "The value of maxloss after : " << maxloss << "\n";
    //Rcout << "The value of neweta : " << neweta << "\n";
    
    R = neweta*eta.row(t).array().inverse()*R+(1.0 + neweta*r).log();
    eta.row(t+1).array() = neweta; 
    
    if (std::isnan(R.sum())) Rcout << "Nan in R" << std::endl;
    
  }
  // end progress
  if (! quiet) end_progress_cpp();
  
  return ;
}


// [[Rcpp::export]]
void computeMLProdEigen(Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> eta, 
                        Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights,
                        Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> predictions, 
                        Eigen::Map<Eigen::VectorXd> R, Eigen::Map<Eigen::VectorXd> L,
                        Eigen::Map<Eigen::VectorXd> maxloss,
                        String loss_name,double loss_tau,bool loss_gradient, bool quiet){
  
  
  const std::string cs=std::string(loss_name.get_cstring());
  
  if (loss_gradient){
    if (cs=="square" ) return MLProdEigen<SquL,true>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau,quiet);
    if (cs=="absolute") return MLProdEigen<AbsL,true>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau,quiet);
    if (cs=="percentage") return MLProdEigen<PerL,true>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau,quiet);
    if (cs=="log") return MLProdEigen<LogL,true>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau,quiet);
    if (cs=="pinball") return MLProdEigen<PinL,true>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau,quiet);
  }
  else{
    if (cs=="square" ) return MLProdEigen<SquL,false>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau,quiet);
    if (cs=="absolute") return MLProdEigen<AbsL,false>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau,quiet);
    if (cs=="percentage") return MLProdEigen<PerL,false>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau,quiet);
    if (cs=="log") return MLProdEigen<LogL,false>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau,quiet);
    if (cs=="pinball") return MLProdEigen<PinL,false>(awake,eta,experts,weights,y,predictions,R,L,maxloss,loss_tau,quiet);
  }
  
  Rcout << "********** ERROR !!! " << cs << std::endl;
  return ;
  
}

