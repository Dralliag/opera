#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int count_in(int x, IntegerVector y) {
  int count = 0;
  int ny = y.size();
  
  for (int i=0; i < ny; ++i) {
    if (y[i] == x)
      count += 1;
  }
  return count;
}

// [[Rcpp::export]]
IntegerVector init_progress_cpp(int T) {
  Rprintf("[---------|---------|----------|---------|---------|---------|---------|---------|---------|---------]\r");
  Rprintf("[");

  IntegerVector steps(100);

  for (int i=0; i<100; i++) {
    steps[i] = int(round((i+1)*T/100.0 + 0.49999));
  }
  
  return steps;
}

// [[Rcpp::export]]
void update_progress_cpp(int t, IntegerVector steps){
  
  int count = count_in(t, steps);

  if (count > 0) {
    for (int i=0; i<count; i++) {
      Rprintf("="); 
    } 
  }
}

// [[Rcpp::export]]
void end_progress_cpp(){
  Rprintf("]\n");
}