#include "kernel.h"

tau_m_results DLR_su_tau_m(const NumericVector* sfuns, const IntegerVector &CDFcounts, const int numCDF, const int* lens, const int numTests, const NumericVector &support, const bool adaptive = true, const double alpha = 0.05, const double zeta = 0.5) {
  // "m(l)" for tau_m (see 16)
  int a = (int)std::floor(alpha * numTests) + 1;
  // number of p-values to be transformed
  int numValues = support.length();
  
  // index of tau_m
  int idx_tau = 0;
  // tau_m
  double tau_m = 1;
  // step function evaluations
  std::vector<double> tau_m_eval(numCDF);
  
  // last evaluation positions of step functions
  int* pos = new int[numCDF];
  for(int i = 0; i < numCDF; i++) pos[i] = lens[i] - 1;
  
  // positions for binary search
  int idx_left = 0, idx_right = numValues - 1, idx_mid = numValues - 1;
  // bool variable to indicate end of search
  bool stop = false;
  // traverse step functions increasingly
  bool sf_incr = false;
  // loop index
  int i = 0;
  // sum
  double sum = 0;
  // number of remaining needed values for adaptive sum
  int rem = 0;
  // sorting order (for adaptive sum)
  IntegerVector ord;
  
  // start binary search
  while(!stop) {
    // check if user wishes to abort computations
    checkUserInterrupt();
    // evaluate pCDFs
    if(sf_incr)
      for(i = 0; i < numCDF; i++)
        eval_pv(tau_m_eval[i], support[idx_mid], sfuns[i], lens[i], pos[i]);
    else
      for(i = 0; i < numCDF; i++)
        eval_pv_rev(tau_m_eval[i], support[idx_mid], sfuns[i], pos[i]);
    
    // compute sum in (21)
    sum = 0;
    if(adaptive) {
      NumericVector fracs(numCDF);
      for(i = 0; i < numCDF; i++)
        fracs[i] = tau_m_eval[i] / (1 - tau_m_eval[i]);
      
      if(numCDF == numTests) {
        std::sort(fracs.begin(), fracs.end(), std::greater<double>());
        ord = IntegerVector(Range(0, numCDF - 1));
      } else ord = order(fracs, true);
      
      i = 0;
      rem = a;
      while(i < numCDF && CDFcounts[ord[i]] <= rem){
        sum += CDFcounts[ord[i]] * fracs[ord[i]];
        rem -= CDFcounts[ord[i]];
        i++;
      }
      if(rem > 0) sum += rem * fracs[ord[i]];
    } else {
      for(i = 0; i < numCDF; i++)
        sum += CDFcounts[i] * tau_m_eval[i] / (1 - tau_m_eval[i]);
    }
    
    if(sum > zeta * a) {
      if(idx_mid == 0) {
        // no p-value can satisfy condition
        idx_tau = idx_mid;
        tau_m = support[idx_tau];
        for(i = 0; i < numCDF; i++) 
          eval_pv_rev(tau_m_eval[i], support[idx_tau], sfuns[i], pos[i]);
        stop = true;
      } else if(idx_mid - idx_left == 1) {
        stop = true;
        // left p-value is the last one to satisfy condition
        idx_tau = idx_left;
        tau_m = support[idx_left];
        for(i = 0; i < numCDF; i++) 
          eval_pv_rev(tau_m_eval[i], support[idx_left], sfuns[i], pos[i]);
      } else {
        // tau_m MUST be smaller than the current p-value
        idx_right = idx_mid;
        idx_mid = idx_left + (idx_mid - idx_left) / 2;
        sf_incr = false;
      }
    } else {
      // tau_m COULD be larger than the current p-value
      if(idx_mid == numValues - 1 ||  sum == zeta * a || idx_right - idx_mid == 1) {
        // if difference between mid and right position equals 1 or the largest
        //   p-value satisfies the condition or sum equals threshold, then we
        //   found tau_m
        stop = true;
        idx_tau = idx_mid;
        tau_m = support[idx_mid];
      } else {
        idx_left = idx_mid;
        idx_mid = idx_left + (idx_right - idx_mid + 1) / 2;
        sf_incr = true;
      }
    }
  }
  
  // garbage collection
  delete[] pos;
  //Rcout << tau_m << "\n";
  return {tau_m, idx_tau, tau_m_eval};
}

NumericVector kernel_DLR_fast(const List &pCDFlist, const NumericVector &sorted_pv, const bool adaptive, const double alpha, const bool stepUp, const double zeta, const NumericVector &support, const Nullable<IntegerVector> &pCDFcounts) {
  // number of tests
  int numTests = sorted_pv.length();
  // sequence 1 ... m
  NumericVector seq_m = NumericVector(IntegerVector(seq_len(numTests)));
  // sequence floor(1 * alpha) ... floor(m * alpha) for adaptive procedure
  IntegerVector a = IntegerVector(NumericVector(floor(seq_m * alpha)));
  // vector to store F_i(tau_m) for SU case only
  std::vector<double> f_denom;
  
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // get count of each unique p-value distribution
  IntegerVector CDFcounts;
  if(pCDFcounts.isNull()) CDFcounts = IntegerVector(numCDF, 1);
  else CDFcounts = pCDFcounts;
  
  // array of p-value CDF vectors
  NumericVector* sfuns = new NumericVector[numCDF];
  // lengths of CDFs
  int* lens = new int[numCDF];
  // gather CDFs and their lengths
  for(int i = 0; i < numCDF; i++) {
    sfuns[i] = as<NumericVector>(pCDFlist[i]);
    lens[i] = sfuns[i].length();
  }
  
  // set of p-values to transform
  NumericVector pv_list;
  // reduce p-value set (if possible)
  if(!stepUp) {
    // SD case, see (18)
    // do not reduce p-value set
    pv_list = sorted_pv;
  } else {
    // SU case, see (16, 17)
    // compute tau_m
    tau_m_results tau_m = DLR_su_tau_m(sfuns, CDFcounts, numCDF, lens, numTests, support, adaptive, alpha, zeta);
    // restrict attention to values <= tau_m
    pv_list = sorted_pv[Range(0, binary_search(sorted_pv, tau_m.value, numTests))];
    // vector to store F_i(tau_m) for SU case only
    f_denom = tau_m.eval;
  }
  
  // number of p-values to be transformed
  int numValues = pv_list.length();
  // vector to store transformed p-values
  NumericVector pval_transf(numValues);
  // nothing to compute if there are no p-values
  if(numValues){
    // possibly large data size requires to use chunks
    // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
    int numPV = std::max<int>(1, std::pow(2.0, 26.0) / numCDF);
    // number of chunks
    int chunks = (numValues - 1) / numPV + 1;
    // last evaluation positions of step functions
    int* pos = new int[numCDF]();
    // loop variables
    int i, j, k;
    // index of current p-value in pv_list(!)
    int idx_pval = 0;
    // number of remaining needed values (for adaptive sum)
    int rem = 0;
    // sorting order (for adaptive sum)
    IntegerVector ord;
    
    for(i = 0; i < chunks; i++) {
      // the min(..., numValues) is here for the last chunk
      NumericVector pv = pv_list[Range(i * numPV, std::min<int>((i + 1) * numPV, numValues) - 1)];
      // length of the vector (differs in the last chunk)
      numPV = pv.length();
      // rows:    indices from 1 to numCDF
      // columns: p-values
      NumericMatrix mat(numCDF, numPV);
      
      // compute columns \sum_{j=1}^numTests F_j(pv)
      for(j = 0; j < numCDF; j++) {
        for(k = 0; k < numPV; k++)
          eval_pv(mat(j, k), pv[k], sfuns[j], lens[j], pos[j]);
        
        if(stepUp) // SU case, see (17)
          mat(j, _) = mat(j, _) / (1 - f_denom[j]);
      }
      
      // compute transformed p-values
      for(j = 0; j < numPV; j++) {
        checkUserInterrupt();
        // evaluated F_j
        //   (re-use variable "pv"; previous values are no longer needed)
        pv = NumericVector(mat(_, j));
        
        // compute sum
        if(adaptive) {
          // number of remaining needed values
          rem = numTests - idx_pval + a[idx_pval];
          
          if(pCDFcounts.isNull()) {
            // sort F_i evaluations for current p-value
            std::sort(pv.begin(), pv.end(), std::greater<double>());
            // compute weighted sum
            for(k = 0; k < rem; k++)
              pval_transf[idx_pval] += pv[k];
          } else {
            // get order
            ord = order(pv, true);
            // compute weighted sum
            for(k = 0; k < numCDF && CDFcounts[ord[k]] < rem; k++) {
              pval_transf[idx_pval] += CDFcounts[ord[k]] * pv[ord[k]];
              rem -= CDFcounts[ord[k]];
            }
            pval_transf[idx_pval] += rem * pv[ord[k]];
          }
        } else {
          if(pCDFcounts.isNull()) {
            for(k = 0; k < numCDF; k++)
              pval_transf[idx_pval] += pv[k];
          } else
            // compute weighted sum
            pval_transf[idx_pval] = sum(NumericVector(CDFcounts) * pv);
        }
        // increase index of current p-value in pv_list(!)
        idx_pval++;
      }
    }
    // garbage collection
    delete[] pos;
  }
  
  // garbage collection
  delete[] sfuns;
  delete[] lens;
  
  // output
  return pval_transf;
}

List kernel_DLR_crit(const List &pCDFlist, const NumericVector &support, const NumericVector &sorted_pv, const bool adaptive, const double alpha, const double zeta, const bool stepUp, const Nullable<IntegerVector> &pCDFcounts) {
  // number of tests
  int numTests = sorted_pv.length();
  // sequence 1 ... m
  NumericVector seq_m = NumericVector(IntegerVector(seq_len(numTests)));
  // sequence floor(1 * alpha) ... floor(m * alpha)
  IntegerVector a = IntegerVector(NumericVector(floor(seq_m * alpha)));
  
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // get count of each unique p-value distribution
  IntegerVector CDFcounts;
  if(pCDFcounts.isNull()) CDFcounts = IntegerVector(numCDF, 1);
  else CDFcounts = pCDFcounts;
  
  // array of p-value CDF vectors
  NumericVector* sfuns = new NumericVector[numCDF];
  // lengths of CDFs
  int* lens = new int[numCDF];
  // gather CDFs and their lengths
  for(int i = 0; i < numCDF; i++) {
    sfuns[i] = as<NumericVector>(pCDFlist[i]);
    lens[i] = sfuns[i].length();
  }
  
  // vector to store critical values indices
  IntegerVector crit(numTests);
  // support size
  int numValues = support.length();
  // set of p-values to transform
  NumericVector pv_list;
  // vector to store F_i(tau_m) for SU case only
  std::vector<double>  f_denom;
  // reduce p-value support
  if(!stepUp) {
    // SD case
    // tau.1 for reducing support
    double tau_1 = zeta * (std::floor(alpha) + 1) / (numTests + std::floor(alpha));
    // apply shortcut, that is tau_1 >= the effective critical value associated 
    //   to zeta * (floor(alpha) + 1) / (numTests + floor(alpha))
    //pv_list = short_eff(support, tau_1);
    int idx = binary_search(support, tau_1, numValues);
    pv_list = support[Range(idx, numValues - 1)];
    // then re-add the observed p-values (needed to compute the adjusted p-values),
    // because we may have removed some of them by the shortcut
    pv_list = NumericVector(sort_combine(sorted_pv, pv_list));
    // set minimum critical values indices to the one of the largest value <= tau_1
    crit.fill(idx);
  } else {
    // SU case, see (16, 17)
    // compute tau_m
    tau_m_results tau_m = DLR_su_tau_m(sfuns, CDFcounts, numCDF, lens, numTests, support, adaptive, alpha, zeta);
    // restrict attention to values <= tau_m
    pv_list = support[Range(0, tau_m.index)];
    // vector to store F_i(tau_m) for SU case only
    f_denom = tau_m.eval;
    // last critical value found
    crit[numTests - 1] = tau_m.index;
  }
  
  // number of p-values in support
  numValues = pv_list.length();
  // number of p-values to be transformed
  //int numTransf = sorted_pv.length();
  // vector to store transformed p-values
  NumericVector pval_transf(numTests);
  
  // nothing to compute if there are no p-values
  if(numValues) {
    // possibly large data size requires to use chunks
    // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
    int size = std::max<int>(1, std::pow(2.0, 26.0) / numCDF);
    // number of chunks
    int chunks = (numValues - 1) / size + 1;
    // actual number of p-values in current chunk
    int numPV = size;
    // last evaluation positions of step functions
    int* pos = new int[numCDF]();
    // loop variables
    int i, j, k;
    // index of current critical value
    int idx_crit = 0;
    // index of current raw p-value to be transformed
    int idx_transf = 0;
    // current p-value in chunk
    double pv_current = 0;
    // sum
    double s;
    // number of remaining needed values (for adaptive sum)
    int rem = 0;
    // sorting order (for adaptive sum)
    IntegerVector ord;
    // evaluated F_j
    NumericVector temp;
    
    // compute critical values (and transformed raw p-values for step-down)
    for(i = 0; i < chunks; i++){
      checkUserInterrupt();
      // the min(..., numValues) is here for the last chunk
      NumericVector pv = pv_list[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
      // actual length of the vector (differs in the last chunk)
      numPV = pv.length();
      // rows:    indices from 1 to numCDF
      // columns: p-values
      NumericMatrix mat(numCDF, numPV);
      
      // compute columns \sum_{j=1}^numTests F_j(pv)
      for(k = 0; k < numCDF; k++) {
        for(j = 0; j < numPV; j++)
          eval_pv(mat(k, j), pv[j], sfuns[k], lens[k], pos[k]);
        
        if(stepUp) // SU case, see (17)
          mat(k, _) = mat(k, _) / (1 - f_denom[k]);
      }
      
      // compute transformed p-value support (as in pv_list)
      //   stop loop, when either the last p-value of the chunk is reached
      //   or the last critical value is found
      j = 0;
      //while(j < numPV && (idx_transf < numTests || ((!stepUp && idx_crit < numTests) || (stepUp && idx_crit < numTests - 1)))){
      while(j < numPV && ((!stepUp && (idx_transf < numTests || idx_crit < numTests)) || (stepUp && idx_crit < numTests - 1))){
        checkUserInterrupt();
        // evaluated F_j
        temp = NumericVector(mat(_, j));
        // current p-value
        pv_current = pv[j];
        // compute sum
        s = 0;
        if(adaptive) {
          // get order
          if(pCDFcounts.isNull()) {
            std::sort(temp.begin(), temp.end(), std::greater<double>());
            //ord = IntegerVector(Range(0, numCDF - 1));
          } else ord = order(temp, true);
          
          if((!stepUp && idx_crit < numTests) || (stepUp && idx_crit < numTests - 1)) {
            // number of remaining needed values
            rem = numTests - idx_crit + a[idx_crit];
            // compute weighted sum
            if(pCDFcounts.isNull()) {
              for(k = 0; k < rem && s <= zeta * (a[idx_crit] + 1); k++)
                s+= temp[k];
            } else {
              for(k = 0; k < numCDF && CDFcounts[ord[k]] < rem && s <= zeta * (a[idx_crit] + 1); k++) {
                s += CDFcounts[ord[k]] * temp[ord[k]];
                rem -= CDFcounts[ord[k]];
              }
              s += rem * temp[ord[k]];
            }
          }
        } else {
          // compute weighted sum
          if(pCDFcounts.isNull()) {
            s = sum(temp);
          } else {
            s = sum(NumericVector(CDFcounts) * temp);
          }
        }
        
        if((!stepUp && idx_crit < numTests) || (stepUp && idx_crit < numTests - 1)) {
          // check satisfaction of condition
          if(s <= zeta * (a[idx_crit] + 1)){
            // current p-value does not satisfy condition
            // current p-value satisfies condition
            // => save index of current p-value as critical value
            crit[idx_crit] = i * size + j;
            // go to next p-value in this chunk
            j++;
          } else {
            // current p-value does not satisfy condition
            //   previous p-value is critical value
            if(idx_crit > 0 && i * size + j > 0 && crit[idx_crit - 1] == i * size + j - 1)
              crit[idx_crit] = i * size + j - 1;
            // go to next critical value index to search for
            idx_crit++;
          }
        }
        // compute transformed raw p-value for step-down, if there is at least
        // one equal to current support p-value
        if(!stepUp) {
          while(idx_transf < numTests && pv_current > sorted_pv[idx_transf]) idx_transf++;
          while(idx_transf < numTests && pv_current == sorted_pv[idx_transf]){
            if(adaptive) {
              // number of remaining needed values
              rem = numTests - idx_transf + a[idx_transf];
              // compute weighted sum
              if(pCDFcounts.isNull()) {
                for(k = 0; k < rem; k++)
                  pval_transf[idx_transf] += temp[k];
              } else {
                for(k = 0; k < numCDF && CDFcounts[ord[k]] < rem; k++) {
                  pval_transf[idx_transf] += CDFcounts[ord[k]] * temp[ord[k]];
                  rem -= CDFcounts[ord[k]];
                }
                pval_transf[idx_transf] += rem * temp[ord[k]];
              }
            } else {
              pval_transf[idx_transf] = s;
            }
            idx_transf++;
          }
          if(idx_crit == numTests) j++;
        }
      }
      if((!stepUp && idx_transf == numTests && idx_crit == numTests) || (stepUp && idx_crit  == numTests)) break;
    }
    
    if((!stepUp && idx_crit < numTests - 1) || (stepUp && idx_crit < numTests - 2)){
      for(i = idx_crit + 1; i < numTests; i++) crit[i] = crit[idx_crit];
    }
    // garbage collection
    delete[] pos;
  }
  
  // garbage collection
  delete[] sfuns;
  delete[] lens;
  
  // output
  if(stepUp)
    return List::create(Named("crit.consts") = pv_list[crit]);
  else
    return List::create(Named("crit.consts") = pv_list[crit], Named("pval.transf") = pval_transf);
}