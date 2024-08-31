#include "dist_fun.h"
#include "helper.h"

//' @title Kernel functions
//' @name kernel
//' 
//' @description
//' Kernel functions transform observed p-values or their support according to
//' \[HLR\], \[PB\] and \[HGR\]. The output is used by [`discrete.LR()`],
//' [`discrete.PB()`] and [`discrete.GR()`], respectively.
//' For each procedure, there is a kernel for fast computation and one for
//' calculation of critical values. Kernels followed by ".crit", e.g.
//' `kernel.DGR.crit`, compute and return these critical values, while
//' kernels ending in ".fast" only transform p-values and are therefore faster.
//' The end user should not use these functions directly.
//' 
//' @seealso
//' [`FDX`][`FDX-package`], [`discrete.LR()`]
//' [`discrete.GR()`], [`discrete.PB()`],
//' [`weighted.LR()`], [`weighted.GR()`],
//' [`discrete.PB()`]
//'
//' @templateVar pCDFlist TRUE
//' @templateVar pvalues TRUE
//' @templateVar sorted_pv TRUE
//' @templateVar adaptive TRUE
//' @templateVar alpha TRUE
//' @templateVar stepUp TRUE
//' @templateVar zeta TRUE
//' @templateVar support TRUE
//' @templateVar exact TRUE
//' @templateVar raw.pvalues FALSE
//' @templateVar direction FALSE
//' @templateVar ret.crit.consts FALSE
//' @templateVar weights TRUE
//' @template param 
//' @param geom_weighting    a boolean specifying whether to conduct geometric
//'                          (`TRUE`) or arithmetic (`FALSE`)
//'                          weighting.
//' @param qvalues           a numeric vector. Contains weighted raw p-values.
//' 
//' @template example
//' @examples
//' 
//' alpha <- 0.05
//' 
//' # If not searching for critical constants, we use only the observed p-values
//' sorted.pvals <- sort(raw.pvalues)
//' y.DLR.fast <- kernel_DLR_fast(pCDFlist, sorted.pvals, TRUE)
//' y.NDGR.fast <- kernel_DGR_fast(pCDFlist, sorted.pvals, FALSE)$pval.transf
//' # transformed values
//' y.DLR.fast
//' y.NDGR.fast
//' 
//' # compute support
//' pv.list <- sort(unique(unlist(pCDFlist)))
//' y.DGR.crit <- kernel_DGR_crit(pCDFlist, pv.list, sorted.pvals, TRUE)
//' y.NDPB.crit <- kernel_DPB_crit(pCDFlist, pv.list, sorted.pvals, FALSE)
//' # critical constants
//' y.DGR.crit$crit.consts
//' y.NDPB.crit$crit.consts
//' # transformed values
//' y.DGR.crit$pval.transf
//' y.NDPB.crit$pval.transf
//' 
//' @return
//' For ".fast" kernels, a vector of transformed p-values is returned; ".crit"
//' kernels return a list object with critical constants (`$crit.consts`)
//' and transformed p-values (`$pval.transf`).
//' 

double DLR_su_tau_m2(const List &pCDFlist, const int numTests, const NumericVector &pvalues, const bool adaptive = true, const double alpha = 0.05, const double zeta = 0.5) {
  // "m(l)" for tau_m (see 16)
  int a = (int)std::floor(alpha * numTests) + 1;
  
  // number of p-values to be transformed
  int numValues = pvalues.length();
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numTests);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  
  // set of p-values to transform
  NumericVector pv_list = rev(pvalues);
  // index of tau.m
  int idx_tau = 0;
  // bool variable to store/check end of search
  bool stop = false;
  // start computations
  for(int i = 0; i < chunks; i++){
    checkUserInterrupt();
    // the min( , numValues) is here for the last chunk
    NumericVector pv = pv_list[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int len = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numTests, len);
    // compute columns \sum_{j=1}^numTests F_j(pv) / (1 - F_j(pv))
    for(int j = 0; j < numTests; j++){
      NumericVector f_eval = rev(stepfun(rev(pv), as<NumericVector>(pCDFlist[j])));
      mat(j, _) = f_eval / (1 - f_eval);
    }
    // sort columns in descending order
    if(adaptive) colsortdec(mat);
    // compute transformed p-values
    // start with first p-value of the current chunk
    int j = 0;
    // stop loop, when either the last p-value of the chunk is reached or the
    // critical value is found
    double s;
    while(j < len && !stop){
      checkUserInterrupt();
      // evaluated F_j in descending order
      // (re-use variable "pv"; previous values are no longer needed)
      pv = NumericVector(mat(_, j));
      // compute sum
      if(adaptive)
        s = sum(pv[Range(0, a - 1)]);
      else
        s = sum(pv);
      
      // check satisfaction of condition
      if(s <= zeta * a){
        // current p-value satisfies condition => mark index
        idx_tau = i * size + j;
        // tau.m is found, stop searching
        stop = true;
        break;
      }else{
        // current p-value does not satisfy condition
        // go to next p-value of this chunk
        j++;
      }
    }
    // if search is stopped, there is need for further chunk computations
    if(stop) break;
  }
  //Rcout << pv_list[idx_tau] << "\n";
  return pv_list[idx_tau];
}

// [[Rcpp::export]]
NumericVector kernel_DLR_fast2(const List &pCDFlist, const NumericVector &pvalues, const bool adaptive = true, const double alpha = 0.05, const bool stepUp = false, const double zeta = 0.5, const NumericVector &support = 0) {
  // number of tests
  int numTests = pCDFlist.length();
  // seqence 1 ... m
  NumericVector seq_m = NumericVector(IntegerVector(seq_len(numTests)));
  // sequence floor(1 * alpha) ... floor(m * alpha) for adaptive procedure
  IntegerVector a = IntegerVector(NumericVector(floor(seq_m * alpha)));
  // vector to store F_i(tau_m) for SU case only
  NumericVector f_denom(numTests);
  
  // set of p-values to transform
  NumericVector pv_list;
  // reduce p-value set (if possible)
  if(!stepUp){
    // SD case, see (18)
    // do not reduce p-value set
    pv_list = pvalues;
  }else{
    // SU case, see (16, 17)
    // compute tau_m
    double tau_m = DLR_su_tau_m2(pCDFlist, numTests, support, adaptive, alpha, zeta);
    // search the values of the vector <= numTests * alpha
    // restrict attention to these values, because tau.k needs to be <= tau.m
    pv_list = pvalues[pvalues <= tau_m];
    // pre-compute F_i(tau_m) to avoid unnecessary re-computing
    for(int i = 0; i < numTests; i++) 
      f_denom[i] = 1 - stepfun(NumericVector(1, tau_m), as<NumericVector>(pCDFlist[i]))[0];
  }
  
  // number of p-values to be transformed
  int numValues = pv_list.length();
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numTests);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  // vector to store transformed p-values
  NumericVector pval_transf(numValues);
  
  for(int i = 0; i < chunks; i++){
    checkUserInterrupt();
    // the min( , numValues) is here for the last chunk
    NumericVector pv = pv_list[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int len = pv.length();
    
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numTests, len);
    // compute columns \sum_{j=1}^numTests F_j(pv)
    if(!stepUp)
      // SD case, see (18)
      for(int j = 0; j < numTests; j++)
        mat(j, _) = stepfun(pv, as<NumericVector>(pCDFlist[j]));
    else
      // SU case, see (17)
      for(int j = 0; j < numTests; j++)
        mat(j, _) = stepfun(pv, as<NumericVector>(pCDFlist[j])) / f_denom[j];
    // sort columns in descending order
    if(adaptive) colsortdec(mat);
    // compute transformed p-values
    for(int j = 0; j < len; j++){
      checkUserInterrupt();
      // index of current value in pv_list (!)
      int idx_pval = i * size + j;
      // evaluated F_j in descending order
      // (re-use variable "pv"; previous values are no longer needed)
      pv = NumericVector(mat(_, j));
      // compute sum
      if(adaptive)
        pval_transf[idx_pval] = sum(pv[Range(0, numTests - idx_pval + a[idx_pval] - 1)]);
      else
        pval_transf[idx_pval] = sum(pv);
    }
  }
  
  // output
  return pval_transf;
}

// [[Rcpp::export]]
List kernel_DLR_crit2(const List &pCDFlist, const NumericVector &pvalues, const NumericVector &sorted_pv, const bool adaptive = true, const double alpha = 0.05, const double zeta = 0.5, const bool stepUp = false){
  // number of tests
  int numTests = pCDFlist.length();
  // seqence 1 ... m
  NumericVector seq_m = NumericVector(IntegerVector(seq_len(numTests)));
  // sequence floor(1 * alpha) ... floor(m * alpha)
  IntegerVector a = IntegerVector(NumericVector(floor(seq_m * alpha)));
  // vector to store critical values indices
  IntegerVector crit(numTests);
  // vector to store transformed p-values
  NumericVector pval_transf(numTests);
  // tau.1 for reducing support for SD case only (????????????)
  double tau_1 = zeta * (std::floor(alpha) + 1) / (numTests + std::floor(alpha));
  // vector to store F_i(tau_m) for SU case only
  NumericVector f_denom(numTests);
  
  // set of p-values to transform
  NumericVector pv_list;
  // reduce p-value support
  if(!stepUp){
    // SD case
    // apply the shortcut drawn from ???, that is
    // c.1 >= the effective critical value associated to zeta * (floor(alpha) + 1) / (numTests + floor(alpha))
    pv_list = short_eff(pvalues, tau_1);
    // then re-add the observed p-values (needed to compute the adjusted p-values),
    // because we may have removed some of them by the shortcut
    pv_list = NumericVector(sort_combine(sorted_pv, pv_list));
    // set minimum critical values indices to the one of the largest value <= tau_1
  }
  else{
    // SU case, see (16, 17)
    // compute tau.m
    double tau_m = DLR_su_tau_m2(pCDFlist, numTests, pvalues, adaptive, alpha, zeta);
    // search the values of the vector <= numTests * alpha
    // restrict attention to these values, because tau.k needs to be <= tau.m
    pv_list = pvalues[pvalues <= tau_m];
    // apply the shortcut drawn from Lemma ???, that is
    // c.1 >= the effective critical value associated to zeta * (floor(alpha) + 1) / (numTests + floor(alpha))
    //pv_list = short_eff(pv_list, tau_1);
    // pre-compute F_i(tau_m) to avoid unnecessary re-computing
    for(int i = 0; i < numTests; i++) 
      f_denom[i] = 1 - stepfun(NumericVector(1, tau_m), as<NumericVector>(pCDFlist[i]))[0];
  }
  
  // number of p-values to be transformed
  int numValues = pv_list.length();
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numTests);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  
  // index of current critical value
  int idx_crit = 0;
  // index of current raw p-value to be transformed
  int idx_transf = 0;
  // compute critical values (and transformed raw p-values for step-down)
  for(int i = 0; i < chunks; i++){
    checkUserInterrupt();
    // the min( , numValues) is here for the last chunk
    NumericVector pv = pv_list[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int len = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numTests, len);
    // compute columns
    if(stepUp)
      // SU case, see (17)
      // F_j(pv)/F_j(tau.m)
      for(int j = 0; j < numTests; j++)
        mat(j, _) = NumericVector(stepfun(pv, as<NumericVector>(pCDFlist[j])) / f_denom[j]);
    else
      // SD case, see (18)
      // F_j(pv)
      for(int j = 0; j < numTests; j++)
        mat(j, _) = stepfun(pv, as<NumericVector>(pCDFlist[j]));
    // sort columns in descending order
    if(adaptive) colsortdec(mat);
    // compute transformed p-value support (as in pv_list)
    // start with first p-value of the current chunk
    int j = 0;
    // stop loop, when either the last p-value of the chunk is reached or the
    // last critical value is found
    while(j < len && ((!stepUp && (idx_transf < numTests || idx_crit < numTests)) || (stepUp && idx_crit < numTests))){
      checkUserInterrupt();
      // evaluated F_j in descending order
      NumericVector temp = NumericVector(mat(_, j));
      // compute sum
      double s = 0, pv_current = pv[j];
      if(adaptive){
        if(idx_crit < numTests) s = sum(temp[Range(0, numTests - idx_crit + a[idx_crit] - 1)]);
      }else s = sum(temp);
      if(idx_crit < numTests){
        // check satisfaction of condition
        if(s <= zeta * (a[idx_crit] + 1)){
          // current p-value does not satisfy condition
          // current p-value satisfies condition
          // => save index of current p-value as critical value
          crit[idx_crit] = i * size + j;
          // go to next p-value in this chunk
          j++;
        }else{
          // current p-value does not satisfy condition
          if(idx_crit > 0 && i * size + j > 0 && crit[idx_crit - 1] == i * size + j - 1)
            crit[idx_crit] = i * size + j - 1;
          // go to next critical value index to search for
          idx_crit++;
        }
      }
      // compute transformed raw p-value for step-down, if there is at least
      // one equal to current support p-value
      if(!stepUp){
        while(idx_transf < numTests && pv_current > sorted_pv[idx_transf]) idx_transf++;
        while(idx_transf < numTests && pv_current == sorted_pv[idx_transf]){
          if(adaptive)
            pval_transf[idx_transf] = sum(temp[Range(0, numTests - idx_transf + a[idx_transf] - 1)]);
          else
            pval_transf[idx_transf] = s;
          idx_transf++;
        }
        if(idx_crit == numTests) j++;
      }
    }
    if((!stepUp && idx_transf == numTests && idx_crit == numTests) || (stepUp && idx_crit  == numTests)) break;
  }
  if(idx_crit < numTests){
    for(int i = idx_crit + 1; i < numTests; i++) crit[i] = crit[idx_crit];
  }
  
  // output
  if(stepUp)
    return List::create(Named("crit.consts") = pv_list[crit]);
  else
    return List::create(Named("crit.consts") = pv_list[crit], Named("pval.transf") = pval_transf);
}

// [[Rcpp::export]]
List kernel_DGR_fast2(const List &pCDFlist, const NumericVector &pvalues, const bool adaptive = true, const double alpha = 0.05) {
  // number of tests
  int numTests = pCDFlist.length();
  // form logarithms of pCDFlist
  List logCDFs(numTests);
  for(int i = 0; i < numTests; i++) logCDFs[i] = NumericVector(log(1 - as<NumericVector>(pCDFlist[i])));
  // seqence 1 ... m
  NumericVector seq_m = NumericVector(IntegerVector(seq_len(numTests)));
  // sequence floor(1 * alpha) ... floor(m * alpha) for adaptive procedure
  IntegerVector a = IntegerVector(NumericVector(floor(seq_m * alpha)));
  // vector to store transformed p-values
  NumericVector pval_transf(numTests);
  // vector to store computed binomial probabilities
  NumericVector probs(numTests);
  
  // log p-value support
  NumericVector log_pv_list = log(1 - pvalues);
  // number of p-values to be transformed
  int numValues = pvalues.length();
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numTests);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  
  for(int i = 0; i < chunks; i++){
    checkUserInterrupt();
    // the min( , numValues) is here for the last chunk
    NumericVector pv = log_pv_list[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int len = pv.length();
    
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numTests, len);
    // compute columns \sum_{j=1}^numTests F_j(pv)
    for(int j = 0; j < numTests; j++)
      mat(j, _) = stepfun_leq(pv, as<NumericVector>(logCDFs[j]));
    // sort columns in descending order
    if(adaptive) colsortasc(mat);
    // compute transformed p-values
    for(int j = 0; j < len; j++){
      checkUserInterrupt();
      // index of current value in pv_list (!)
      int idx_pval = i * size + j;
      // number of summands
      int numSum;
      // compute sum
      if(adaptive){
        // number of summands
        numSum = numTests - idx_pval + a[idx_pval];
        // evaluated F_j in descending order
        // (re-use variable "pv"; previous values are no longer needed)
        pv = NumericVector(mat(_, j))[Range(0, numSum - 1)];
      }else{
        // number of summands
        numSum = numTests;
        // evaluated F_j in descending order
        // (re-use variable "pv"; previous values are no longer needed)
        pv = NumericVector(mat(_, j));
      }
      // compute logarithms to avoid numerical problems
      // sum => log of product of probabilities
      double s = sum(pv) / numSum;
      probs[idx_pval] = 1 - std::exp(s);
      pval_transf[idx_pval] = R::pbinom(a[idx_pval], numSum, probs[idx_pval], false, false);
    }
  }
  
  // output
  return List::create(Named("pval.transf") = pval_transf, Named("bin.probs") = probs);//pval_transf;
}

// [[Rcpp::export]]
List kernel_DGR_crit2(const List &pCDFlist, const NumericVector &pvalues, const NumericVector &sorted_pv, const bool adaptive = true, const double alpha = 0.05, const double zeta = 0.5) {
  // number of tests
  int numTests = pCDFlist.length();
  // form logarithms of pCDFlist
  List log_sfuns(numTests);
  for(int i = 0; i < numTests; i++) log_sfuns[i] = NumericVector(log(1 - as<NumericVector>(pCDFlist[i])));
  // seqence 1 ... m
  IntegerVector seq_m = IntegerVector(seq_len(numTests));
  // sequence 1 * alpha ... m * alpha for adaptive procedure
  NumericVector seq_alpha = NumericVector(seq_m) * alpha;
  // sequence floor(1 * alpha) ... floor(m * alpha) for adaptive procedure
  IntegerVector a = IntegerVector(NumericVector(floor(seq_alpha)));
  // vector to store critical values indices
  IntegerVector crit(numTests);
  // vector to store transformed p-values
  NumericVector pval_transf(numTests);
  // vector of number of elements to be considered (number of summands)
  IntegerVector numSums;
  
  // Continuous Guo-Romano critical values
  NumericVector log_tau;
  NumericVector tau_GR(numTests);
  if(adaptive){
    for(int i = 0; i < numTests; i++) tau_GR[i] = R::qbeta(zeta, a[i] + 1, numTests - i, true, false);
    log_tau = log(1 - tau_GR);
    numSums = numTests - seq_m + a + 1;
    log_tau = log_tau * NumericVector(numSums);
  }else{
    for(int i = 0; i < numTests; i++) tau_GR[i] = R::qbeta(zeta, a[i] + 1, numTests - a[i], true, false);
    log_tau = log(1 - tau_GR) * numTests;
    numSums = IntegerVector(numTests, numTests);
  }
  
  // reduce and revert p-value support
  NumericVector pv_list = rev(pvalues);// short_eff(pvalues, tau_GR[0]);
  //pv_list = rev(sort_combine(pv_list, sorted_pv));
  // log p-value support
  NumericVector log_pv_list = log(1 - pv_list);
  // number of p-values to be transformed
  int numValues = pv_list.length();
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numTests);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  
  // index of current p-value
  int idx_pval = 0;
  // index of current critical value
  int idx_crit = numTests - 1;
  // index of current raw p-value to be transformed
  int idx_transf = numTests - 1;
  // compute critical values (and transformed raw p-values for step-down)
  for(int i = 0; i < chunks && (idx_transf >= 0 || idx_crit >= 0); i++){
    checkUserInterrupt();
    // the min( , numValues) is here for the last chunk
    NumericVector pv = log_pv_list[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int len = pv.length();
    
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numTests, len);
    // compute columns \sum_{j=1}^numTests F_j(pv)
    for(int j = 0; j < numTests; j++)
      mat(j, _) = rev(stepfun_leq(rev(pv), as<NumericVector>(log_sfuns[j])));
    // sort columns in ascending order
    if(adaptive) colsortasc(mat);//colsortdec(mat);
    
    // compute transformed p-value support (as in pv_list)
    // start with first p-value of the current chunk
    int j = 0;
    // sum
    double s = 0;
    // stop loop, when either the last p-value of the chunk is reached or the
    // last critical value is found
    while(j < len && (idx_transf >= 0 || idx_crit >= 0)){
      checkUserInterrupt();
      // evaluated F_j in ascending order
      NumericVector temp = NumericVector(mat(_, j));
      if(idx_crit >= 0){
        // sum of logarithms to avoid numerical problems
        s = sum(temp[Range(0, numSums[idx_crit] - 1)]);
      }
      // check satisfaction of condition
      if(idx_crit >= 0 && s >= log_tau[idx_crit]){
        // current p-value satisfies condition
        // => save index of current p-value as critical value
        crit[idx_crit] = idx_pval;
        // go to next critical value index to search for
        idx_crit--;
      }else{
        // current p-value does not satisfy condition
        // compute transformed raw p-value for step-down, if there is at least
        // one equal to current support p-value
        while(idx_transf >= 0 && pv_list[idx_pval] < sorted_pv[idx_transf]) idx_transf--;
        while(idx_transf >= 0 && pv_list[idx_pval] == sorted_pv[idx_transf]){
          s = sum(temp[Range(0, numSums[idx_transf] - 1)]);
          pval_transf[idx_transf] = R::pbinom(a[idx_transf], numSums[idx_transf], 1 - std::exp(s / numSums[idx_transf]), false, false);
          idx_transf--;
        }
        // go to next p-value in this chunk
        j++;
        idx_pval++;
      }
    }
    if(idx_transf < 0 && idx_crit < 0) break;
  }
  
  // output
  return List::create(Named("crit.consts") = pv_list[crit], Named("pval.transf") = pval_transf);
}

// [[Rcpp::export]]
NumericVector kernel_DPB_fast2(const List &pCDFlist, const NumericVector &pvalues, const bool adaptive, const double alpha, const bool exact) {
  // number of tests
  int numTests = pCDFlist.length();
  // seqence 1 ... m
  NumericVector seq_m = NumericVector(IntegerVector(seq_len(numTests)));
  // sequence floor(1 * alpha) ... floor(m * alpha) for adaptive procedure
  IntegerVector a = IntegerVector(NumericVector(floor(seq_m * alpha)));
  // method of Poisson-Binomial Distribution
  int method = (int)exact;
  
  // number of p-values to be transformed
  int numValues = pvalues.length();
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numTests);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  // vector to store transformed p-values
  NumericVector pval_transf(numValues);
  // index of current value in pvalues
  int idx_pval = 0;
  
  for(int i = 0; i < chunks; i++){
    checkUserInterrupt();
    // the min( , numValues) is here for the last chunk
    NumericVector pv = pvalues[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int len = pv.length();
    
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numTests, len);
    // compute columns \sum_{j=1}^numTests F_j(pv)
    for(int j = 0; j < numTests; j++)
      mat(j, _) = stepfun(pv, as<NumericVector>(pCDFlist[j]));
    
    // sort columns in descending order
    if(adaptive) colsortdec(mat);
    // compute transformed p-values
    for(int j = 0; j < len; j++){
      checkUserInterrupt();
      // evaluated F_j in descending order
      // (re-use "pv"; previous values are no longer needed)
      pv = NumericVector(mat(_, j));
      // compute transformation
      if(adaptive){
        //NumericVector temp = ppbinom(a, pv[Range(0, numTests - idx_pval + a[idx_pval] - 1)], method, false);
        pval_transf[idx_pval] = ppbinom_vec(a, pv[Range(0, numTests - idx_pval + a[idx_pval] - 1)], method, false)[idx_pval];
      }else{
        //NumericVector temp = ppbinom(a, pv, method, false);
        pval_transf[idx_pval] = ppbinom_vec(a, pv, method, false)[idx_pval];
      }
      // go to next p-value
      idx_pval++;
    }
  }
  
  // output
  return pval_transf;
}

// [[Rcpp::export]]
List kernel_DPB_crit2(const List &pCDFlist, const NumericVector &pvalues, const NumericVector &sorted_pv, const bool adaptive = true, const double alpha = 0.05, const double zeta = 0.5, const bool exact = true) {
  // number of tests
  int numTests = pCDFlist.length();
  // seqence 1 ... m
  NumericVector seq_m = NumericVector(IntegerVector(seq_len(numTests)));
  // sequence floor(1 * alpha) ... floor(m * alpha)
  IntegerVector a = IntegerVector(NumericVector(floor(seq_m * alpha)));
  // vector to store transformed p-values
  NumericVector pval_transf(numTests);
  // set of p-values to transform
  NumericVector pv_list;
  // apply the shortcut drawn from ???, that is
  // c.1 >= the effective critical value of the continuous GR approach
  pv_list = short_eff(pvalues, R::qbeta(zeta, a[0] + 1, numTests, true, false));
  // then re-add the observed p-values (needed to compute the adjusted p-values),
  // because we may have removed some of them by the shortcut
  pv_list = sort_combine(sorted_pv, pv_list);
  // number of p-values to be transformed
  int numValues = pv_list.length();
  
  // method of Poisson-Binomial Distribution
  int method = (int)exact;
  // vector to store critical values indices
  IntegerVector crit(numTests);// = match(crit_GRstar, pv_list) - 1;
  // GR* critical values
  if(exact){
    NumericVector crit_GRstar = as<NumericVector>(kernel_DGR_crit2(pCDFlist, pvalues, sorted_pv, adaptive, alpha, zeta)["crit.consts"]);
    crit = match(crit_GRstar, pv_list) - 1;
  }
  
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numTests);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  
  // index of current critical value
  int idx_crit = 0;
  // index of current raw p-value to be transformed
  int idx_transf = 0;
  // index of current p-value of the support
  int idx_pval = 0;
  // compute critical values (and transformed raw p-values for step-down)
  for(int i = 0; i < chunks; i++){
    checkUserInterrupt();
    // the min( , numValues) is here for the last chunk
    NumericVector pv = pv_list[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int len = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numTests, len);
    // compute columns
    for(int j = 0; j < numTests; j++)
      mat(j, _) = stepfun(pv, as<NumericVector>(pCDFlist[j]));
    
    // sort columns in descending order
    if(adaptive) colsortdec(mat);
    // compute transformed p-value support (as in pv_list)
    // start with first p-value of the current chunk
    int j = 0;
    // boolean to indicate if transformation of a raw p-value is to be computed
    bool raw_transf;
    // stop loop, when either the last p-value of the chunk is reached or the
    // last critical value is found
    while(j < len && (idx_transf < numTests || idx_crit < numTests)){
      checkUserInterrupt();
      // is there a raw p-value equal to the current support p-value?
      if(idx_transf < numTests){
        while(pv[j] > sorted_pv[idx_transf] && idx_transf < numTests - 1) idx_transf++;
        if(pv[j] == sorted_pv[idx_transf])
          raw_transf = true;
        else
          raw_transf = false;
      }else raw_transf = false;
      // vector for transformation(s)
      NumericVector s(numTests);
      // evaluated F_j in descending order
      NumericVector temp = NumericVector(mat(_, j));
      // performance shortcut for non-adaptive procedure
      if(!adaptive && ((idx_crit < numTests && idx_pval > crit[idx_crit]) || raw_transf))
        s = ppbinom_vec(a, temp, method, false);
      // is it still necessary to search for critical values?
      if(idx_crit < numTests && idx_pval > crit[idx_crit]){
        if(adaptive)
          s = ppbinom_vec(a, temp[Range(0, numTests - idx_crit + a[idx_crit] - 1)], method, false);
        
        if(s[idx_crit] <= zeta){
          // current p-value satisfies condition
          // => save index of current p-value as critical value
          crit[idx_crit] = idx_pval;
        }else{
          // current p-value does not satisfy condition
          // go to next critical value index to search for
          idx_crit++;
          continue;
        }
      }
      
      // compute transformed raw p-value for step-down, if there is at least
      // one equal to current support p-value
      while(raw_transf && (pv_list[idx_pval] == sorted_pv[idx_transf])){
        if(adaptive)
          s = ppbinom_vec(a, temp[Range(0, numTests - idx_transf + a[idx_transf] - 1)], method, false);
        
        pval_transf[idx_transf] = s[idx_transf];
        if(++idx_transf == numTests) raw_transf = false;
      }
      // go to next p-value in this chunk
      j++;
      idx_pval++;
    }
    if(idx_transf == numTests && idx_crit == numTests) break;
  }
  // critical values are non-decreasing, but optimized computation left some untouched (=0)
  crit = IntegerVector(cummax(crit));
  
  // output
  return List::create(Named("crit.consts") = pv_list[crit], Named("pval.transf") = pval_transf);
}