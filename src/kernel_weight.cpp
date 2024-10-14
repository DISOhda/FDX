#include "kernel.h"

// geometric weighting function

// [[Rcpp::export]]
NumericVector geom_weight(const NumericVector &pvalues, const NumericVector &weights){
  // length of input vectors (equal)
  int len = pvalues.length();
  // results vector
  NumericVector res(len);
  // perform geometric weighting
  ////// does not work on older compilers (e.g. win_oldrelease)
  //std::transform(res.begin(), res.end(), weights.begin(), res.begin(), [](double x, double a){return std::pow(x, a);});
  //res = 1 - res;
  //////
  for(int i = 0; i < len; i++) {
    res[i] = 1 - std::pow(1 - pvalues[i], weights[i]);
    // was the i-th value rounded off to zero?
    if(res[i] <= 0) 
      // replace zero with first-order Taylor approximation
      res[i] = pvalues[i] * weights[i]; 
  }
  
  // return results
  return res;
}

NumericVector kernel_wLR_fast(const NumericVector &sorted_w_pv, const NumericVector &weights, double alpha, bool geom_weighting) {
  // length of input vectors (equal)
  int numTests = sorted_w_pv.length();
  // seqence 1 ... m
  IntegerVector seq_m = IntegerVector(seq_len(numTests));
  // sequence 1 * alpha ... m * alpha for adaptive procedure
  NumericVector seq_alpha = NumericVector(seq_m) * alpha;
  // sequence floor(1 * alpha) ... floor(m * alpha) for adaptive procedure
  IntegerVector a = IntegerVector(NumericVector(floor(seq_alpha)));
  // length of index sets for weighting
  IntegerVector numWeights = numTests - seq_m + a + 1;
  
  // compute transformations
  NumericVector pval_transf(numTests);
  for(int i = 0; i < numTests; i++) {
    // vector of probbilities
    NumericVector probs;
    if(geom_weighting)
      probs = geom_weight(NumericVector(numWeights[i], sorted_w_pv[i]), weights[Range(0, numWeights[i] - 1)]);
    else {
      probs = sorted_w_pv[i] * clone(weights)[Range(0, numWeights[i] - 1)];
    }
    // transformed values
    pval_transf[i] = sum(probs)/(a[i] + 1);
  }
  
  // output
  return pval_transf;
}
 
NumericVector kernel_wGR_fast(const NumericVector &sorted_w_pv, const NumericVector &weights, double alpha, bool geom_weighting) {
  // length of input vectors (equal)
  int numTests = sorted_w_pv.length();
  // seqence 1 ... m
  IntegerVector seq_m = IntegerVector(seq_len(numTests));
  // sequence 1 * alpha ... m * alpha for adaptive procedure
  NumericVector seq_alpha = NumericVector(seq_m) * alpha;
  // sequence floor(1 * alpha) ... floor(m * alpha) for adaptive procedure
  IntegerVector a = IntegerVector(NumericVector(floor(seq_alpha)));
  // length of index sets for weighting
  IntegerVector numWeights = numTests - seq_m + a + 1;
  
  // compute transformations
  NumericVector pval_transf(numTests);
  if(geom_weighting) {
    NumericVector running_weights = NumericVector(cumsum(weights)) / NumericVector(seq_m);
    running_weights = running_weights[numWeights - 1];
    NumericVector probs = geom_weight(sorted_w_pv, running_weights);
    for(int i = 0; i < numTests; i++) {
      pval_transf[i] = R::pbinom(a[i], numWeights[i], probs[i], 0, 0);
    }
  } else {
    for(int i = 0; i < numTests; i++) {
      NumericVector probs = sorted_w_pv[i] * clone(weights)[Range(0, numWeights[i] - 1)];
      pval_transf[i] = ppbinom_vec(a, probs, 2, false)[i];
    }
  }
  
  // output
  return pval_transf;
}
 
NumericVector kernel_wPB_fast(const NumericVector &sorted_w_pv, const NumericVector &weights, double alpha, bool geom_weighting, bool exact) {
  // length of input vectors (equal)
  int numTests = sorted_w_pv.length();
  // seqence 1 ... m
  IntegerVector seq_m = IntegerVector(seq_len(numTests));
  // sequence 1 * alpha ... m * alpha for adaptive procedure
  NumericVector seq_alpha = NumericVector(seq_m) * alpha;
  // sequence floor(1 * alpha) ... floor(m * alpha) for adaptive procedure
  IntegerVector a = IntegerVector(NumericVector(floor(seq_alpha)));
  // length of index sets for weighting
  IntegerVector numWeights = numTests - seq_m + a + 1;
  
  // compute transformations
  NumericVector pval_transf(numTests);
  for(int i = 0; i < numTests; i++) {
    // vector of probbilities
    NumericVector probs;
    if(geom_weighting)
      probs = geom_weight(NumericVector(numWeights[i], sorted_w_pv[i]), weights[Range(0, numWeights[i] - 1)]);
    else {
      probs = sorted_w_pv[i] * clone(weights)[Range(0, numWeights[i] - 1)];
      probs = pmin(probs, NumericVector(numTests, 1.0));
    }
    // transformed values
    pval_transf[i] = ppbinom_vec(a, probs, int(exact), false)[i];
  }
  
  // output
  return pval_transf;
}