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
  REprintf("[---------|---------|----------|---------|---------|---------|---------|---------|---------|---------]\r");
  REprintf("[");
  IntegerVector steps = seq(1, 100);
  for (int i=1; i<=100; i++) {
    steps[i] = (steps[i]*T + 0.499)/100;
  }
  
  return steps;
}

// [[Rcpp::export]]
void update_progress_cpp(int t, IntegerVector steps){
  
  int count = count_in(t, steps);
  
  if (count > 0) {
    for (int i=0; i<count; i++) {
      REprintf("="); 
    } 
  }
}

// [[Rcpp::export]]
void end_progress_cpp(){
  REprintf("]\n");
}