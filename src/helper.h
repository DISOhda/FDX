#include <RcppArmadillo.h>
using namespace Rcpp;

inline void eval_pv(double &eval, double val, const NumericVector &vec, int len, int &pos){
  //if(val < 1){
  while(pos < len && vec[pos] <= 1 && vec[pos] <= val) pos++;
  if(pos) eval = vec[pos - 1];
  else eval = 0;
  //} else eval = 1;
}

inline void eval_pv_nolim(double &eval, double val, const NumericVector &vec, int len, int &pos){
  //if(val < 1){
  while(pos < len && vec[pos] <= val) pos++;
  if(pos) eval = vec[pos - 1];
  else eval = 0;
  //} else eval = 1;
}

inline void eval_pv_rev(double &eval, double val, const NumericVector &vec, int &pos){
  //if(val < 1){
  while(pos > 0 && (vec[pos] > val || vec[pos] > 1)) pos--;
  if(vec[pos] <= val) eval = vec[pos];
  else eval = 0;
  //} else eval = 1;
}

inline void eval_pv_rev_nolim(double &eval, double val, const NumericVector &vec, int &pos){
  //if(val < 1){
  while(pos > 0 && (vec[pos] > val)) pos--;
  if(vec[pos] <= val) eval = vec[pos];
  else eval = 0;
  //} else eval = 1;
}

inline int binary_search(const NumericVector &vec, const double value, const int len) {
  int pos_left = 0, pos_right = len - 1, pos_mid = len - 1;
  bool stop = false;
  while(!stop) {
    if(vec[pos_mid] > value) {
      if(pos_mid == 0) {
        stop = true;
      } else {
        pos_right = pos_mid;
        pos_mid = pos_left + (pos_right - pos_left) / 2;
      }
    } else if(vec[pos_mid] <= value) {
      if(vec[pos_mid] == value || pos_mid == len - 1 || pos_right - pos_mid == 1) {
        stop = true;
      } else{
        pos_left = pos_mid;
        pos_mid = pos_left + (pos_right - pos_left) / 2;
      }
    }
  }
  return pos_mid;
}

NumericVector sort_combine(const NumericVector &x, const NumericVector &y);



//////////// OLD, TO BE REMOVED!!!!! /////////////////

IntegerVector order(const NumericVector &x, bool descending);

NumericVector stepfun(const NumericVector &x, const NumericVector &sfun);
NumericVector stepfun_leq(const NumericVector &x, const NumericVector &sfun);

NumericVector short_eff(const NumericVector &x, const double limit);

void colsortdec(NumericMatrix &mat);
void colsortasc(NumericMatrix &mat);