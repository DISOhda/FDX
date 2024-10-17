#include "kernel.h"

List kernel_DGR_fast(const List &pCDFlist, const NumericVector &sorted_pv, const bool adaptive, const double alpha, const Nullable<IntegerVector> &pCDFcounts) {
  // number of tests
  int numTests = sorted_pv.length();
  // sequence 1 ... m
  IntegerVector seq_m = IntegerVector(seq_len(numTests));
  // sequence floor(1 * alpha) ... floor(m * alpha) for adaptive procedure
  IntegerVector a = IntegerVector(NumericVector(floor(NumericVector(seq_m) * alpha)));
  
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // get count of each unique p-value distribution
  IntegerVector CDFcounts;
  if(pCDFcounts.isNull()) CDFcounts = IntegerVector(numCDF, 1);
  else CDFcounts = pCDFcounts;
  
  // logarithms of pCDFlist
  NumericVector* log_sfuns = new NumericVector[numCDF];
  // lengths of CDFs
  int* lens = new int[numCDF];
  // gather CDFs and their lengths
  for(int i = 0; i < numCDF; i++) {
    log_sfuns[i] = NumericVector(-log(1 - as<NumericVector>(pCDFlist[i])));
    lens[i] = log_sfuns[i].length();
  }
  
  // vector to store transformed p-values
  NumericVector pval_transf(numTests);
  // vector of number of elements to be considered (number of summands)
  IntegerVector numSums;
  if(adaptive)
    numSums = numTests - seq_m + a + 1;
  else
    numSums = IntegerVector(numTests, numTests);
  
  // log p-value support
  NumericVector log_pv_list = -log(1 - sorted_pv);
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numCDF);
  // number of chunks
  int chunks = (numTests - 1) / size + 1;
  // number of p-values in chunk (differs in the last chunk)
  int numPV;
  // last evaluation positions of step functions
  int* pos = new int[numCDF]();
  // loop variables
  int i, j, k;
  // index of current p-value in pv_list(!)
  int idx_pval = 0;
  // sum
  double s;
  // number of remaining needed values (for adaptive sum)
  int rem;
  // sorting order (for adaptive sum)
  IntegerVector ord;
  // vector to store computed binomial probabilities
  NumericVector probs(numTests);
  
  for(i = 0; i < chunks; i++){
    checkUserInterrupt();
    // the min(..., numTests) is here for the last chunk
    NumericVector pv = log_pv_list[Range(i * size, std::min<int>((i + 1) * size, numTests) - 1)];
    // length of the vector
    numPV = pv.length();
    
    // rows:    indices from 1 to numCDF
    // columns: p-values
    NumericMatrix mat(numCDF, numPV);
    // compute columns \sum_{j=1}^numTests F_j(pv)
    for(k = 0; k < numCDF; k++)
      for(j = 0; j < numPV; j++)
        eval_pv_nolim(mat(k, j), pv[j], log_sfuns[k], lens[k], pos[k]);
    
    // compute transformed p-values
    for(j = 0; j < numPV; j++){
      checkUserInterrupt();
      // evaluated F_j
      //   (re-use variable "pv"; previous values are no longer needed)
      pv = NumericVector(mat(_, j));
      // compute sum
      s = 0;
      if(adaptive) {
        // get order
        if(pCDFcounts.isNull()){
          std::sort(pv.begin(), pv.end(), std::greater<double>());
          ord = IntegerVector(Range(0, numCDF - 1));
        } else ord = order(pv, true);
        // number of remaining needed values
        rem = numSums[idx_pval];
        // compute weighted sum
        for(k = 0; k < numCDF && CDFcounts[ord[k]] < rem; k++) {
          s += CDFcounts[ord[k]] * pv[ord[k]];
          rem -= CDFcounts[ord[k]];
        }
        s += rem * pv[ord[k]];
      } else {
        // compute weighted sum
        s = sum(NumericVector(CDFcounts) * pv);
      }
      // compute logarithms to avoid numerical problems
      // sum => log of product of probabilities
      probs[idx_pval] = 1 - std::exp(-s / numSums[idx_pval]);
      pval_transf[idx_pval] = R::pbinom(a[idx_pval], numSums[idx_pval], probs[idx_pval], false, false);
      
      // index of next p-value in pv_list (!)
      idx_pval++;
    }
  }
  
  // garbage collection
  delete[] log_sfuns;
  delete[] lens;
  delete[] pos;
  
  // output
  return List::create(Named("pval.transf") = pval_transf, Named("bin.probs") = probs);//pval_transf;
}

List kernel_DGR_crit(const List &pCDFlist, const NumericVector &support, const NumericVector &sorted_pv, const bool adaptive, const double alpha, const double zeta, const Nullable<IntegerVector> &pCDFcounts) {
  // number of tests
  int numTests = sorted_pv.length();
  // sequence 1 ... m
  IntegerVector seq_m = IntegerVector(seq_len(numTests));
  // sequence 1 * alpha ... m * alpha for adaptive procedure
  NumericVector seq_alpha = NumericVector(seq_m) * alpha;
  // sequence floor(1 * alpha) ... floor(m * alpha)
  IntegerVector a = IntegerVector(NumericVector(floor(seq_alpha)));
  
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // get count of each unique p-value distribution
  IntegerVector CDFcounts;
  if(pCDFcounts.isNull()) CDFcounts = IntegerVector(numCDF, 1);
  else CDFcounts = pCDFcounts;
  
  // logarithms of pCDFlist
  NumericVector* log_sfuns = new NumericVector[numCDF];
  // last evaluation positions of step functions
  int* pos = new int[numCDF]();
  // gather CDFs and their lengths and starting positions
  for(int i = 0; i < numCDF; i++) {
    //log_sfuns[i] = NumericVector(-log(1 - as<NumericVector>(pCDFlist[i])));
    //pos[i] = log_sfuns[i].length() - 1;
    log_sfuns[i] = clone(as<NumericVector>(pCDFlist[i]));
    pos[i] = log_sfuns[i].length() - 1;
    while(pos[i] > 0 && log_sfuns[i][pos[i]] > 1) pos[i]--;
    log_sfuns[i] = log_sfuns[i][Range(0, pos[i])];
    for(int j = pos[i]; j >= 0; j--)
      log_sfuns[i][j] = -std::log(1 - log_sfuns[i][j]);
  }
  
  // vector to store critical values indices
  IntegerVector crit(numTests);
  // vector to store transformed p-values
  NumericVector pval_transf(numTests);
  // vector of number of elements to be considered (number of summands)
  IntegerVector numSums;
  
  // Continuous Guo-Romano critical values
  NumericVector log_tau;
  NumericVector tau_GR(numTests);
  if(adaptive) {
    for(int i = 0; i < numTests; i++) 
      tau_GR[i] = R::qbeta(zeta, a[i] + 1, numTests - i, true, false);
    log_tau = -log(1 - tau_GR);
    numSums = numTests - seq_m + a + 1;
    log_tau = log_tau * NumericVector(numSums);
  } else {
    for(int i = 0; i < numTests; i++) 
      tau_GR[i] = R::qbeta(zeta, a[i] + 1, numTests - a[i], true, false);
    log_tau = -log(1 - tau_GR) * numTests;
    numSums = IntegerVector(numTests, numTests);
  }
  
  // reduce and revert p-value support
  NumericVector pv_list = rev(support);
  //pv_list = rev(sort_combine(pv_list, sorted_pv));
  // log p-value support
  NumericVector log_pv_list = -log(1 - pv_list);
  // number of p-values in support
  int numValues = pv_list.length();
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numCDF);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  // actual number of p-values in current chunk
  int numPV = size;
  // loop variables
  int i, j, k;
  // index of current p-value
  int idx_pval = 0;
  // index of current critical value
  int idx_crit = numTests - 1;
  // index of current raw p-value to be transformed
  int idx_transf = numTests - 1;
  // current p-value in chunk
  double pv_current = 0;
  // sum
  double s;
  // number of remaining needed values (for adaptive sum)
  int rem = 0;
  // sorting order (for adaptive sum)
  IntegerVector ord;
  // vector to store computed binomial probabilities
  NumericVector probs(numTests);
  // evaluated F_j
  NumericVector temp;
  
  // compute critical values and transformed raw p-values
  for(i = 0; i < chunks && (idx_transf >= 0 || idx_crit >= 0); i++){
    checkUserInterrupt();
    // the min(..., numValues) is here for the last chunk
    NumericVector pv = log_pv_list[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    numPV = pv.length();
    
    // rows:    indices from 1 to numCDF
    // columns: p-values
    NumericMatrix mat(numCDF, numPV);
    // compute columns \sum_{j=1}^numTests F_j(pv)
    for(k = 0; k < numCDF; k++)
      for(j = 0; j < numPV; j++)
        eval_pv_rev_nolim(mat(k, j), pv[j], log_sfuns[k], pos[k]);
    
    // compute transformed p-value support (as in pv_list)
    // start with first p-value of the current chunk
    j = 0;
    while(j < numPV && (idx_transf >= 0 || idx_crit >= 0)){
      checkUserInterrupt();
      // evaluated F_j in ascending order
      NumericVector temp = NumericVector(mat(_, j));
      // current p-value
      pv_current = pv_list[idx_pval];
      
      if(adaptive) {
        // get order
        if(pCDFcounts.isNull()){
          std::sort(temp.begin(), temp.end(), std::greater<double>());
          ord = IntegerVector(Range(0, numCDF - 1));
        } else ord = order(temp, true);
      }
      
      // compute sum
      s = 0;
      if(!adaptive || idx_crit >= 0){
        // sum of logarithms to avoid numerical problems
        //s = sum(temp[Range(0, numSums[idx_crit] - 1)]);
        if(adaptive) {
          // number of remaining needed values
          rem = numSums[idx_crit];
          // compute weighted sum
          for(k = 0; k < numCDF && CDFcounts[ord[k]] < rem; k++) {
            s += CDFcounts[ord[k]] * temp[ord[k]];
            rem -= CDFcounts[ord[k]];
          }
          s += rem * temp[ord[k]];
        } else {
          // compute weighted sum
          s = sum(NumericVector(CDFcounts) * temp);
        }
      }
      // check satisfaction of condition
      if(idx_crit >= 0 && s <= log_tau[idx_crit]){
        // current p-value satisfies condition
        // => save index of current p-value as critical value
        crit[idx_crit] = idx_pval;
        // go to next critical value index to search for
        idx_crit--;
      } else {
        // current p-value does not satisfy condition
        // compute transformed raw p-value for step-down, if there is at least
        // one equal to current support p-value
        while(idx_transf >= 0 && pv_current < sorted_pv[idx_transf]) idx_transf--;
        while(idx_transf >= 0 && pv_current == sorted_pv[idx_transf]){
          //s = sum(temp[Range(0, numSums[idx_transf] - 1)]);
          if(adaptive) {
            // number of remaining needed values
            rem = numSums[idx_transf];
            // compute weighted sum
            s = 0;
            for(k = 0; k < numCDF && CDFcounts[ord[k]] < rem; k++) {
              s += CDFcounts[ord[k]] * temp[ord[k]];
              rem -= CDFcounts[ord[k]];
            }
            s += rem * temp[ord[k]];
          }
          probs[idx_transf] = 1 - std::exp(-s / numSums[idx_transf]);
          pval_transf[idx_transf] = R::pbinom(a[idx_transf], numSums[idx_transf], probs[idx_transf], false, false);
          idx_transf--;
        }
        // go to next p-value in this chunk
        j++;
        idx_pval++;
      }
    }
    if(idx_transf < 0 && idx_crit < 0) break;
  }
  
  // garbage collection
  delete[] log_sfuns;
  //delete[] lens;
  delete[] pos;
  
  // output
  return List::create(Named("crit.consts") = pv_list[crit], Named("pval.transf") = pval_transf, Named("bin.probs") = probs);
}