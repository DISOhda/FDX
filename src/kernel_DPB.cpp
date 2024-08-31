#include "kernel.h"

NumericVector kernel_DPB_fast(const List &pCDFlist, const NumericVector &sorted_pv, const bool adaptive, const double alpha, const bool exact, const Nullable<IntegerVector> &pCDFcounts) {//, const Nullable<List> &pCDFindices) {
  // number of tests
  int numTests = sorted_pv.length();
  // sequence 1 ... m
  NumericVector seq_m = NumericVector(IntegerVector(seq_len(numTests)));
  // sequence floor(1 * alpha) ... floor(m * alpha) for adaptive procedure
  IntegerVector a = IntegerVector(NumericVector(floor(seq_m * alpha)));
  // method of Poisson-Binomial Distribution
  int method = (int)exact;
  
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // array of unique p-value CDF indices
  //IntegerVector* sfuns_idx = new IntegerVector[numCDF];
  // get count of each unique p-value distribution
  IntegerVector CDFcounts;
  //if(pCDFcounts.isNull()) CDFcounts = IntegerVector(numCDF, 1);
  //else CDFcounts = pCDFcounts;
  if(pCDFcounts.isNotNull()) {
    CDFcounts = pCDFcounts;
    //if(!adaptive) {
    //  for(int i = 0; i < numCDF; i++)
    //    sfuns_idx[i] = as<IntegerVector>(List(pCDFindices)[i]);
    //} else delete[] sfuns_idx;
  }
  
  // array of p-value CDF vectors
  NumericVector* sfuns = new NumericVector[numCDF];
  // lengths of CDFs
  int* lens = new int[numCDF];
  // gather CDFs and their lengths
  for(int i = 0; i < numCDF; i++) {
    sfuns[i] = as<NumericVector>(pCDFlist[i]);
    lens[i] = sfuns[i].length();
  }
  
  // number of p-values to be transformed
  int numValues = sorted_pv.length();
  // vector to store transformed p-values
  NumericVector pval_transf(numValues);
  
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numCDF);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  // number of p-values in chunk (differs in the last chunk)
  int numPV;
  // last evaluation positions of step functions
  int* pos = new int[numCDF]();
  // loop variables
  int i, j, k, l, m;
  // index of current p-value
  int idx_pval = 0;
  // number of remaining needed values (for adaptive sum)
  int rem;
  // sorting order (for adaptive sum)
  IntegerVector ord;
  // vector to store computed binomial probabilities
  NumericVector probs;
  
  for(i = 0; i < chunks; i++) {
    checkUserInterrupt();
    // the min( , numValues) is here for the last chunk
    NumericVector pv = sorted_pv[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    numPV = pv.length();
    
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numCDF, numPV);
    // compute columns \sum_{j=1}^numTests F_j(pv)
    for(k = 0; k < numCDF; k++)
      for(j = 0; j < numPV; j++)
        eval_pv(mat(k, j), pv[j], sfuns[k], lens[k], pos[k]);
    
    // compute transformed p-values
    for(int j = 0; j < numPV; j++) {
      checkUserInterrupt();
      // evaluated F_j in descending order
      // (re-use "pv"; previous values are no longer needed)
      pv = NumericVector(mat(_, j));
      
      // extract probabilities for Poisson-Binomial-Distribution
      if(adaptive) {
        // number of remaining needed values
        rem = numTests - idx_pval + a[idx_pval];
        // get order and extract probabilities
        if(pCDFcounts.isNull()) {
          // sort in descending order;
          std::sort(pv.begin(), pv.end(), std::greater<double>());
          // extract probabilities
          probs = pv[Range(0, rem - 1)];
        } else {
          // get order
          ord = order(pv, true);
          // vector of probabilities
          probs = NumericVector(rem);
          // extract probabilities
          m = 0;
          for(k = 0; k < numCDF && m < rem; k++) {
            for(l = 0; l < CDFcounts[ord[k]] && m < rem; l++) {
              probs[m++] = pv[ord[k]];
            }
          }
        }
      } else {
        // vector of probabilities
        probs = NumericVector(numTests);
        // extract probabilities
        if(pCDFcounts.isNull()) {
          probs = pv;
        } else {
          m = 0;
          for(k = 0; k < numCDF; k++) {
            for(l = 0; l < CDFcounts[k]; l++) {
              //probs[sfuns_idx[k][l] - 1] = pv[k];
              probs[m++] = pv[k];
            }
          }
        }
      }
      
      // compute FDX
      pval_transf[idx_pval] = ppbinom(a[idx_pval], probs, method, false);
      
      // go to next p-value
      idx_pval++;
    }
  }
  
  //garbage collection
  delete[] sfuns;
  delete[] lens;
  delete[] pos;
  //if(!adaptive && pCDFcounts.isNotNull()) 
  //  delete[] sfuns_idx;
  
  // output
  return pval_transf;
}

List kernel_DPB_crit(const List &pCDFlist, const NumericVector &support, const NumericVector &sorted_pv, const bool adaptive, const double alpha, const double zeta, const bool exact, const Nullable<IntegerVector> &pCDFcounts) {//, const Nullable<List> &pCDFindices) {
  // number of tests
  int numTests = sorted_pv.length();
  // sequence 1 ... m
  NumericVector seq_m = NumericVector(IntegerVector(seq_len(numTests)));
  // sequence floor(1 * alpha) ... floor(m * alpha)
  IntegerVector a = IntegerVector(NumericVector(floor(seq_m * alpha)));
  // method of Poisson-Binomial Distribution
  int method = (int)exact;
  
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // array of unique p-value CDF indices
  //IntegerVector* sfuns_idx = nullptr;
  // get count of each unique p-value distribution
  IntegerVector CDFcounts;
  //if(pCDFcounts.isNull()) CDFcounts = IntegerVector(numCDF, 1);
  //else CDFcounts = pCDFcounts;
  if(pCDFcounts.isNotNull()) {
    CDFcounts = pCDFcounts;
    //if(!adaptive) {
    //  sfuns_idx = new IntegerVector[numCDF];
    //  for(int i = 0; i < numCDF; i++)
    //    sfuns_idx[i] = as<IntegerVector>(List(pCDFindices)[i]);
    //}
  }
  
  // array of p-value CDF vectors
  NumericVector* sfuns = new NumericVector[numCDF];
  // lengths of CDFs
  int* lens = new int[numCDF];
  // gather CDFs and their lengths
  for(int i = 0; i < numCDF; i++) {
    sfuns[i] = as<NumericVector>(pCDFlist[i]);
    lens[i] = sfuns[i].length();
  }
  
  // number of p-values in the support
  int numValues = support.length();
  // vector to store transformed p-values
  NumericVector pval_transf(numTests);
  // set of p-values to transform
  NumericVector pv_list;
  
  // apply the shortcut drawn from ???, that is
  // c.1 >= the effective critical value of the continuous GR approach
  int idx = binary_search(support, R::qbeta(zeta, a[0] + 1, numTests, true, false), numValues);
  pv_list = support[Range(idx, numValues - 1)];
  // then re-add the observed p-values (needed to compute the adjusted p-values),
  // because we may have removed some of them by the shortcut
  pv_list = sort_combine(sorted_pv, pv_list);
  // number of p-values to be transformed
  numValues = pv_list.length();
  
  // vector to store critical values indices
  IntegerVector crit(numTests);// = match(crit_GRstar, pv_list) - 1;
  // GR* critical values
  if(exact){
    NumericVector crit_GRstar = as<NumericVector>(kernel_DGR_crit(pCDFlist, support, sorted_pv, adaptive, alpha, zeta, pCDFcounts)["crit.consts"]);
    crit = match(crit_GRstar, pv_list) - 1;
  }
  
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numCDF);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  // actual number of p-values in current chunk
  int numPV;
  // last evaluation positions of step functions
  int* pos = new int[numCDF]();
  // loop variables
  int i, j, k, l, m;
  // index of current p-value of the support
  int idx_pval = 0;
  // index of current critical value
  int idx_crit = 0;
  // index of current raw p-value to be transformed
  int idx_transf = 0;
  // Poisson-binomial CDF values 
  NumericVector s;
  // number of remaining needed values (for adaptive sum)
  int rem = 0, rem_c = 0, rem_t = 0;
  // sorting order (for adaptive sum)
  IntegerVector ord;
  // vector to store computed binomial probabilities
  NumericVector probs;
  // evaluated F_j
  NumericVector temp;
  // current probability
  double prob = 0;
  
  // compute critical values (and transformed raw p-values for step-down)
  for(i = 0; i < chunks; i++){
    checkUserInterrupt();
    // the min( , numValues) is here for the last chunk
    NumericVector pv = pv_list[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    numPV = pv.length();
    // rows:    indices from 1 to numCDF
    // columns: p-values
    NumericMatrix mat(numCDF, numPV);
    // compute columns
    for(k = 0; k < numCDF; k++)
      for(j = 0; j < numPV; j++)
        eval_pv(mat(k, j), pv[j], sfuns[k], lens[k], pos[k]);
    
    // compute transformed p-value support (as in pv_list)
    // start with first p-value of the current chunk
    j = 0;
    // boolean to indicate if transformation of a raw p-value is to be computed
    bool raw_transf;
    // stop loop, when either the last p-value of the chunk is reached or the
    // last critical value is found
    while(j < numPV && (idx_transf < numTests || idx_crit < numTests)){
      checkUserInterrupt();
      
      // is there a raw p-value equal to the current support p-value?
      if(idx_transf < numTests) {
        while(pv[j] > sorted_pv[idx_transf] && idx_transf < numTests - 1) idx_transf++;
        if(pv[j] == sorted_pv[idx_transf])
          raw_transf = true;
        else
          raw_transf = false;
      } //else raw_transf = false;
      
      // evaluated F_j in descending order
      temp = NumericVector(mat(_, j));
      
      // extract probabilities for Poisson-Binomial-Distribution
      if(adaptive) {
        // number of remaining needed values
        if(idx_crit < numTests)
          rem_c = numTests - idx_crit   + a[idx_crit];
        if(idx_transf < numTests)
          rem_t = numTests - idx_transf + a[idx_transf];
        rem = std::max<int>(rem_c, rem_t);
        
        // get order and extract probabilities
        if(pCDFcounts.isNull()) {
          // sort in descending order;
          std::sort(temp.begin(), temp.end(), std::greater<double>());
          // extract probabilities
          probs = temp[Range(0, rem - 1)];
        } else {
          // get order
          ord = order(temp, true);
          // vector of probabilities
          probs = NumericVector(rem);
          // extract probabilities
          m = 0;
          for(k = 0; k < numCDF && m < rem; k++) {
            for(l = 0; l < CDFcounts[ord[k]] && m < rem; l++) {
              probs[m++] = temp[ord[k]];
            }
          }
        }
      } else {
        // vector of probabilities
        probs = NumericVector(numTests);
        // extract probabilities
        if(pCDFcounts.isNull()) {
          probs = temp;
        } else {
          m = 0;
          for(k = 0; k < numCDF; k++) {
            for(int l = 0; l < CDFcounts[k]; l++) {
              //probs[sfuns_idx[k][l] - 1] = temp[k];
              probs[m++] = temp[k];
            }
          }
        }
      }
      
      // performance shortcut for non-adaptive procedure
      if(!adaptive && ((idx_crit < numTests && idx_pval > crit[idx_crit]) || raw_transf)) {
        s = ppbinom_vec(a, probs, method, false);
      }
      // is it still necessary to search for critical values?
      if(idx_crit < numTests && idx_pval > crit[idx_crit]) {
        if(adaptive)
          prob = ppbinom(a[idx_crit], probs[Range(0, rem_c - 1)], method, false);
        else
          prob = s[idx_crit];
        
        if(prob <= zeta) {
          // current p-value satisfies condition
          // => save index of current p-value as critical value
          crit[idx_crit] = idx_pval;
        } else {
          // current p-value does not satisfy condition
          // go to next critical value index to search for
          idx_crit++;
          continue;
        }
      }
      
      // compute transformed raw p-value for step-down, if there is at least
      // one equal to current support p-value
      while(raw_transf && idx_transf < numTests && pv[j] == sorted_pv[idx_transf]) {
        if(adaptive) {
          rem_t = numTests - idx_transf + a[idx_transf];
          pval_transf[idx_transf] = ppbinom(a[idx_transf], probs[Range(0, rem_t - 1)], method, false);
        } else 
          pval_transf[idx_transf] = s[idx_transf];
        
        idx_transf++;
      }
      if(idx_transf == numTests) raw_transf = false;
      // go to next p-value in this chunk
      j++;
      idx_pval++;
    }
    if(idx_transf == numTests && idx_crit == numTests) break;
  }
  // critical values are non-decreasing, but optimized computation left some untouched (=0)
  crit = IntegerVector(cummax(crit));
  
  //garbage collection
  delete[] sfuns;
  delete[] lens;
  delete[] pos;
  //if(!adaptive && pCDFcounts.isNotNull()) 
  //  delete[] sfuns_idx;
  
  // output
  return List::create(Named("crit.consts") = pv_list[crit], Named("pval.transf") = pval_transf);
}