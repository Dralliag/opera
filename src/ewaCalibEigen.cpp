#include "closspred.h"
#include <RcppEigen.h>
using namespace Rcpp;
#include <limits>
// [[Rcpp::depends(RcppEigen)]]

template <class LT, bool G>
size_t ewaCalib( size_t tp1,double dbesteta,
                 Eigen::Map<Eigen::MatrixXd> awake,                 
                 Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights,
                 Eigen::Map<Eigen::MatrixXd> weta, Eigen::Map<Eigen::MatrixXd> w0,
                 Eigen::Map<Eigen::MatrixXd> grideta, 
                 Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> eta,
                 Eigen::Map<Eigen::VectorXd> cumulativeloss, Eigen::Map<Eigen::VectorXd> prediction,
                 double loss_tau, bool init_grid_eta){
  using Row = Eigen::Array<double,1,Eigen::Dynamic>;
  const size_t t=tp1-1;
  size_t besteta=std::floor(dbesteta)-1;
  
  const size_t T=experts.rows();
  const size_t N=experts.cols();
  const size_t neta=weta.cols();
  
  const auto awaket = awake.row(t).array();
  Row w = Row::Zero(N);
  Row lexp = Row::Zero(N);
  
  w = (weta.array().col(besteta).transpose() * awaket);
  w /= w.sum();
  weights.row(t).array() = w;
  prediction[t] = experts.row(t) * w.matrix().transpose();
  eta[t] = grideta(besteta,0);
  
  besteta=0;
  double mincl=std::numeric_limits<double>::max();
  
  double init_grideta=0.0;
  for (size_t l=0 ; l<neta ; l++){
    w = (weta.array().col(l).transpose() * awaket);
    w /=w.sum();
    const double pred = experts.row(t) * w.matrix().transpose();//scalar product 
    
    const auto lf=LossPredFunctor<LT,false>(y[t],pred,loss_tau);
    cumulativeloss[l] += lf(pred);
    if (mincl>cumulativeloss[l]){
      besteta=l;
      mincl=cumulativeloss[l];
    }
    const auto lpf=LossPredFunctor<LT,G>(y[t],pred,loss_tau);
    
    const double lpred=lpf(pred);
    lexp = experts.row(t).array().unaryExpr(lpf);  
    
    w0.col(l).array().transpose()+= awaket * (lpred-lexp);
    init_grideta+=w0.col(l).array().transpose().unaryExpr([](double c) {return abs(c);}).sum();
  }
  
  // init value of grid.eta if NULL at first step
  if (tp1 == 1 && init_grid_eta) {
    grideta(0,0) = 1 / (init_grideta / (w0.cols()*w0.rows()));
    eta[0]=R_NaN;
    for (size_t l2=1 ; l2<T ; l2++){
      eta[l2]=grideta(0,0);                
    }
  }
  
  
  for (size_t k=0 ; k < N ; k++){
    //weta.row(k)=(grideta.array().transpose() * w0.row(k).array()).array().exp().unaryExpr(std::ptr_fun(truncate1));
    weta.row(k)=(grideta.array().transpose() * w0.row(k).array()).array().exp().unaryExpr([](double c) {return truncate1(c);});
  }
  return besteta+1;
}

// [[Rcpp::export]]
size_t computeEWACalib( size_t tp1,double dbesteta,
                        Eigen::Map<Eigen::MatrixXd> awake,                 
                        Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights,
                        Eigen::Map<Eigen::MatrixXd> weta, Eigen::Map<Eigen::MatrixXd> w0,
                        Eigen::Map<Eigen::MatrixXd> grideta, 
                        Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> eta,
                        Eigen::Map<Eigen::VectorXd> cumulativeloss, Eigen::Map<Eigen::VectorXd> prediction,
                        String loss_name, double loss_tau, bool loss_gradient, bool init_grid_eta){
  
  
  const std::string cs=std::string(loss_name.get_cstring());
  
  if (loss_gradient){
    if (cs=="square" ) return ewaCalib<SquL,true>(tp1,dbesteta,awake,experts,weights,weta,w0,grideta,y,eta,cumulativeloss,prediction,loss_tau,init_grid_eta);
    if (cs=="absolute") return ewaCalib<AbsL,true>(tp1,dbesteta,awake,experts,weights,weta,w0,grideta,y,eta,cumulativeloss,prediction,loss_tau,init_grid_eta);
    if (cs=="percentage") return ewaCalib<PerL,true>(tp1,dbesteta,awake,experts,weights,weta,w0,grideta,y,eta,cumulativeloss,prediction,loss_tau,init_grid_eta);
    if (cs=="log") return ewaCalib<LogL,true>(tp1,dbesteta,awake,experts,weights,weta,w0,grideta,y,eta,cumulativeloss,prediction,loss_tau,init_grid_eta);
    if (cs=="pinball") return ewaCalib<PinL,true>(tp1,dbesteta,awake,experts,weights,weta,w0,grideta,y,eta,cumulativeloss,prediction,loss_tau,init_grid_eta);
  }
  else{
    if (cs=="square" ) return ewaCalib<SquL,false>(tp1,dbesteta,awake,experts,weights,weta,w0,grideta,y,eta,cumulativeloss,prediction,loss_tau,init_grid_eta);
    if (cs=="absolute") return ewaCalib<AbsL,false>(tp1,dbesteta,awake,experts,weights,weta,w0,grideta,y,eta,cumulativeloss,prediction,loss_tau,init_grid_eta);
    if (cs=="percentage") return ewaCalib<PerL,false>(tp1,dbesteta,awake,experts,weights,weta,w0,grideta,y,eta,cumulativeloss,prediction,loss_tau,init_grid_eta);
    if (cs=="log") return ewaCalib<LogL,false>(tp1,dbesteta,awake,experts,weights,weta,w0,grideta,y,eta,cumulativeloss,prediction,loss_tau,init_grid_eta);
    if (cs=="pinball") return ewaCalib<PinL,false>(tp1,dbesteta,awake,experts,weights,weta,w0,grideta,y,eta,cumulativeloss,prediction,loss_tau,init_grid_eta);
  }
  
  Rcout << "********** ERROR !!! " << cs << std::endl;
  return 0;
}





