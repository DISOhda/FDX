#include "old_helper.h"

// fast step function evaluation
// step function is represented by a single numeric vector under the conditions
// a) f(x) = x and b) 'x' and 'sfun' are sorted in ASCENDING order
// this is much faster than passing and evaluating R step function objects
NumericVector stepfun(const NumericVector &x, const NumericVector &sfun){
  // index variables and vector lengths
  int pos = 0, size = sfun.length(), len = x.length();
  // output vector of the same length as 'x'
  NumericVector out(len);
  
  // computing results
  for(int i = 0; i < len; i++){
    while(pos < size - 1 && sfun[pos] < x[i]) pos++;
    if(sfun[pos] == x[i]) out[i] = sfun[pos];
    else if(pos) out[i] = sfun[pos - 1]; else out[i] = 0;
  }
  
  return out;
}

// fast step function evaluation
// step function is represented by a single numeric vector under the conditions
// a) f(x) = x and b) 'x' and 'sfun' are sorted in DESCENDING order
// this is much faster than passing and evaluating R step function objects
NumericVector stepfun_leq(const NumericVector &x, const NumericVector &sfun){
  // index variables and vector lengths
  int size = sfun.length(), len = x.length(), pos = 0;
  // output vector of the same length as 'x'
  NumericVector out(len);
  
  // computing results
  for(int i = 0; i < len; i++){
    while(sfun[pos] > x[i] && pos < size) pos++;
    if(sfun[pos] == x[i]) out[i] = sfun[pos];
    else if(pos) out[i] = sfun[pos - 1]; else out[i] = 0;
  }
  
  return out;
}

// shortcut function that eliminates all values of a SORTED vector that
// are < limit, except the largest value <= limit
NumericVector short_eff(const NumericVector &x, const double limit){
  // length of the vector
  int len = x.length();
  // identify values <= limit
  NumericVector out = x[x <= limit];
  // eliminate values, but keep their maximum
  out = x[Range(which_max(out), len - 1)];
  
  return out;
}

// sort columns of a matrix in DESCENDING order
// using an intermediate numeric vector is necessary, because "in-column"
// sorting does not always work as expected, especially for large columns
void colsortdec(NumericMatrix &mat){
  // intermediate vector to store column values (necessary!)
  NumericVector vec;
  for(int i = 0; i < mat.ncol(); i++){
    // store values in a vector
    vec = NumericVector(mat(_, i));
    // sort values in DESCENDING order
    std::sort(vec.begin(), vec.end(), std::greater<double>());
    // write sorted values back to column
    mat(_, i) = vec;
  }
}

// sort columns of a matrix in ASCENDING order
// using an intermediate numeric vector is necessary, because "in-column"
// sorting does not always work as expected, especially for large columns
void colsortasc(NumericMatrix &mat){
  // intermediate vector to store column values (necessary!)
  NumericVector vec;
  for(int i = 0; i < mat.ncol(); i++){
    // store values in a vector
    vec = NumericVector(mat(_, i));
    // sort values in ASCENDING order
    std::sort(vec.begin(), vec.end());
    // write sorted values back to column
    mat(_, i) = vec;
  }
}
