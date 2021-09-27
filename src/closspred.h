#pragma once
#include <Rcpp.h>
using namespace Rcpp;
#include <string>
#include <cmath>

struct SquL;
struct AbsL;
struct PerL;
struct LogL;
struct PinL;

inline double sign(double x) {
  if (x > 0.0) return 1.0;
  if (x < 0.0) return -1.0;
  return x;
}

inline double truncate1(double x){return std::min(std::max(x,std::exp(-700.0)),std::exp(700));}

inline double sq_log2(double x){return pow(2, ceil(log2(x)));}

inline double sup_half(double x){
  if (x > 0.5) return 1.0;
  if (x <= 0.5) return 0.0;
  return x;
}

inline double ramp(double x){return x>0.0 ? x : 0.0;}

template <class L, bool G>
inline double closspred(double x, double y, double p ,double t){
  Rcout << "********** ERROR !!! " << std::endl;
  return 0;
}
template<> 
inline double closspred<SquL,false>(double x, double y, double p ,double t){return (x-y)*(x-y);}
template<> 
inline double closspred<SquL,true>(double x, double y, double p ,double t){return 2*(p-y)*x;}
template<> 
inline double closspred<AbsL,false>(double x, double y, double p ,double t){return std::abs(x - y);}
template<> 
inline double closspred<AbsL,true>(double x, double y, double p ,double t){return sign(p - y) * x;}
template<> 
inline double closspred<PerL,false>(double x, double y, double p ,double t){return std::abs(x - y)/y;}
template<> 
inline double closspred<PerL,true>(double x, double y, double p ,double t){return x/y * sign(p - y);}
template<> 
inline double closspred<LogL,false>(double x, double y, double p ,double t){return -log(x);}
template<> 
inline double closspred<LogL,true>(double x, double y, double p ,double t){return -(x/p);}
template<> 
inline double closspred<PinL,false>(double x, double y, double p ,double t){return  ((y < x) - t) * (x - y);}
template<> 
inline double closspred<PinL,true>(double x, double y, double p ,double t){return ((y < p) - t) * x ;}

  

template <class L,bool G>
struct LossPredFunctor{
  double yt_,pred_,loss_tau_;
  LossPredFunctor(double yt, double pred, double loss_tau):yt_(yt),pred_(pred),loss_tau_(loss_tau){};
  inline double operator() (double a) const {return closspred<L,G>(a,yt_,pred_,loss_tau_);}
};

struct AbsMax{
  inline double operator() (double a, double b) const {return std::max(a,std::abs(b));} 
};

inline void show(NumericVector a,std::string as){
  const size_t s=a.size();
  Rcout << as<<"[ ";
  for (size_t k=0 ; k<s ; k++) {
    Rcout <<a[k]<<" ";
  }
  Rcout << "]"<<std::endl;
}

inline void show(double a , std::string as){
  Rcout << as<<"="<<a<<std::endl;
}

inline void show(NumericMatrix a, std::string as){
  const size_t T=a.nrow();
  const size_t N=a.ncol();
  for (size_t t=0 ; t<T ; t++) {
    Rcout << as<<" "<<std::to_string(t) <<" [" ;
    for (size_t k=0 ; k<N ; k++) {
      Rcout <<a(t,k)<<" ";
    }
    Rcout << "]"<<std::endl;
  }
}

inline void show(size_t a , std::string as){
  Rcout << as<<"="<<a<<std::endl;
}

