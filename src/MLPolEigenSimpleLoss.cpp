#include "closspred.h"
#include <RcppEigen.h>
#include "progress.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]


// inline double simplelosspredold(const char * loss_name, bool loss_gradient, double x, double y, double p ,double t){
//   if (!loss_gradient){
//     if (strcmp(loss_name,"square")==0){
//       return (x-y)*(x-y);
//     }
//     else if (strcmp(loss_name,"absolute")==0){
//       return std::abs(x - y);
//       }
//     else if (strcmp(loss_name,"percentage")==0){
//       return std::abs(x - y)/y;
//       }
//     else if (strcmp(loss_name,"log")==0){
//       return -log(x);
//       }
//     else if (strcmp(loss_name,"pinball")==0){
//       return  ((y < x) - t) * (x - y);
//     }
//     else{
//       Rcout << "********** ERROR !!! " << loss_name << std::endl;
//       }
//   }
//   else{
//     if (strcmp(loss_name,"square")==0){
//       return 2*(p-y)*x;
//     }
//     else if (strcmp(loss_name,"absolute")==0){
//       return sign(p - y) * x;
//     }
//     else if (strcmp(loss_name,"percentage")==0){
//       return x/y * sign(p - y);
//     }
//     else if (strcmp(loss_name,"log")==0){
//       return -(x/p);
//     }
//     else if (strcmp(loss_name,"pinball")==0){
//       return ((y < p) - t) * x ;
//     }
//     else{
//       Rcout << "********** ERROR !!! " << loss_name << std::endl;
//     }
//   }
//   return 0;
// }   
    
enum class LossType { Square, Absolute, Percentage, Log, Pinball};

LossType getLossType(const char * loss_name){
  if (strcmp(loss_name,"square")==0){
    return LossType::Square;
  }
  else if (strcmp(loss_name,"absolute")==0){
    return LossType::Absolute;
  }
  else if (strcmp(loss_name,"percentage")==0){
    return LossType::Percentage;
  }
  else if (strcmp(loss_name,"log")==0){
    return LossType::Log;
  }
  else if (strcmp(loss_name,"pinball")==0){
    return LossType::Pinball;
  }
  else{
    Rcout << "********** ERROR !!! " << loss_name << std::endl;
  }
  return LossType::Square;
}


inline double simplelosspred(LossType lt, bool loss_gradient, double x, double y, double p ,double t){
  if (!loss_gradient){
    switch(lt)
    {
      case LossType::Square:
        return (x-y)*(x-y);
        break;
      case LossType::Absolute:
        return  std::abs(x - y);
        break;
      case LossType::Percentage:
        return std::abs(x - y)/y;
        break;
      case LossType::Log:
        return -log(x);
        break;
      case LossType::Pinball:
        return  ((y < x) - t) * (x - y);
        break;
      default:
        Rcout << "********** ERROR !!! "  << std::endl;
        break;
    }
  }
  else{
    switch(lt)
    {
    case LossType::Square:
      return 2*(p-y)*x;
      break;
    case LossType::Absolute:
      return   sign(p - y) * x;
      break;
    case LossType::Percentage:
      return x/y * sign(p - y);
      break;
    case LossType::Log:
      return -(x/p);
      break;
    case LossType::Pinball:
      return  ((y < p) - t) * x ;
      break;
    default:
      Rcout << "********** ERROR !!! "  << std::endl;
    break;
    }
  }
  return 0;
}   

    
    
// [[Rcpp::export]]
double computeMLPolEigenSimpleLoss( Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> eta, 
                          Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights,
                          Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> predictions, 
                          Eigen::Map<Eigen::VectorXd> Rc, Eigen::Map<Eigen::VectorXd> wc,
                          double B,
                          String loss_name,double loss_tau,bool loss_gradient, bool quiet){
  
  
 
 const LossType lt = getLossType(loss_name.get_cstring());
  
  
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
    if ((awaket*R).maxCoeff()>0){
      //w = eta.row(t).array() * R.unaryExpr(std::ptr_fun(ramp));
      w = eta.row(t).array() * R.unaryExpr([](double c) {return ramp(c);});
      //w /= w.sum();
    }
    else{
      w = 1.0;
    }
    // form the mixture and the prediction
    p = awaket * w.array();
    p = awake.row(t).array() * w.array();
    p /= p.sum();
    const double pred = experts.row(t) * p.matrix().transpose();//scalar product
    
    //save the mixture and the prediction
    weights.row(t).array() = p;
    predictions[t] = pred;
    
    // Observe losses
    const double lpred=simplelosspred(lt,loss_gradient,pred,y[t],pred,loss_tau);
    for (size_t k=0 ; k<N ; k++) lexp[k]=simplelosspred(lt,loss_gradient,experts(t,k),y[t],pred,loss_tau);

    r = awaket * (lpred-lexp);
    R += r;

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
    
  
 
  


