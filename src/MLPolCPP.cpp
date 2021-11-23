#include "closspred.h"
#include "progress.h"


template <class L, bool G>
double MLPolCPP( NumericMatrix awake, NumericMatrix eta, NumericMatrix experts, NumericMatrix weights,
                 NumericVector y, NumericVector predictions, NumericVector R, NumericVector w,
                 double B,double loss_tau, bool quiet){
  
  const size_t T=experts.nrow();
  const size_t N=experts.ncol();
  
  NumericVector p(N);
  NumericVector lexp(N);
  NumericVector r(N);
  
  // init progress
  IntegerVector steps;
  if (! quiet) steps = init_progress_cpp(T);
  
  for (size_t t=0 ; t<T ; t++){
    // update progress
    if (! quiet) update_progress_cpp(t+1, steps);
    
    double mar=0;
    for (size_t k=0 ; k<N ; k++) mar = std::max(mar,awake(t,k)*R[k]);
    if (mar>0.0){
      double s=0.0;
      for (size_t k=0 ; k<N ; k++){
        w[k] = eta(t,k)*std::max(R[k],0.0);
        s+=w[k];
      }
      const double sm1=1.0/s;
      for (size_t k=0 ; k<N ; k++) w[k]*=sm1;
    }
    else{
      for (size_t k=0 ; k<N ; k++) w[k]=1.0;
    }
    
    //form the mixture and the prediction
    double s=0.0;
    for (size_t k=0 ; k<N ; k++){
      p[k] = awake(t,k)*w[k];
      s+=p[k];
    }
    const double sm1=1.0/s;
    for (size_t k=0 ; k<N ; k++) p[k]*=sm1;
    
    double pred=0.0;
    for (size_t k=0 ; k<N ; k++) pred+=experts(t,k)*p[k];
    
    //save the mixture and the prediction
    for (size_t k=0 ; k<N ; k++) weights(t,k) = p[k];
    predictions[t] = pred;
    
    
    // Observe losses
    double lpred=closspred<L,G>(pred,y[t],pred,loss_tau);
    for (size_t k=0 ; k<N ; k++) lexp[k]=closspred<L,G>(experts(t,k),y[t],pred,loss_tau);
    //Update the regret and the weight
    for (size_t k=0 ; k<N ; k++) r[k] = awake(t,k)*(lpred-lexp[k]);
    for (size_t k=0 ; k<N ; k++) R[k]+=r[k];
    //Update the learning rate
    double mar2=0;
    for (size_t k=0 ; k<N ; k++) mar2 = std::max(mar2,r[k]*r[k]);
    const double newB = std::max(B,mar2);
    for (size_t k=0 ; k<N ; k++) eta(t+1,k) = 1.0/((1.0/eta(t,k))+r[k]*r[k]+newB-B);
    B = newB;
  }
  // end progress
  if (! quiet) end_progress_cpp();
  
  double mar=0;
  for (size_t k=0 ; k<N ; k++) mar = std::max(mar,R[k]);
  if (mar>0){
    double s=0;
    for (size_t k=0 ; k<N ; k++){
      w[k]=eta(T,k)*std::abs(R[k]);
      s+=w[k];
    }
    const double sm1=1/s;
    for (size_t k=0 ; k<N ; k++) w[k]*=sm1;
  }
  else{
    const double nm1=1/double(N);
    for (size_t k=0 ; k<N ; k++) w[k]=nm1;
  }

  return B;
}


// [[Rcpp::export]]
double computeMLPolCPP( NumericMatrix awake, NumericMatrix eta, NumericMatrix experts, NumericMatrix weights,
                      NumericVector y, NumericVector predictions, NumericVector R, NumericVector w,
                      double B,
                      String loss_name,double loss_tau,bool loss_gradient, bool quiet){
  
  
  const std::string cs=std::string(loss_name.get_cstring());
  
  if (loss_gradient){
    if (cs=="square" ) return MLPolCPP<SquL,true>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
    if (cs=="absolute") return MLPolCPP<AbsL,true>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
    if (cs=="percentage") return MLPolCPP<PerL,true>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
    if (cs=="log") return MLPolCPP<LogL,true>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
    if (cs=="pinball") return MLPolCPP<PinL,true>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
  }
  else{
    if (cs=="square" ) return MLPolCPP<SquL,false>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
    if (cs=="absolute") return MLPolCPP<AbsL,false>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
    if (cs=="percentage") return MLPolCPP<PerL,false>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
    if (cs=="log") return MLPolCPP<LogL,false>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
    if (cs=="pinball") return MLPolCPP<PinL,false>(awake,eta,experts,weights,y,predictions,R,w,B,loss_tau,quiet);
  }
  
  Rcout << "********** ERROR !!! " << cs << std::endl;
  return 0;
  
}

