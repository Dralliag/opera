#include <Rcpp.h>
using namespace Rcpp;

IntegerVector init_progress_cpp(int);
void update_progress_cpp(int,IntegerVector);
void end_progress_cpp();