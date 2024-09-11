#' @name match.pvals
#' 
#' @title
#' Matching Raw P-Values with Supports
#' 
#' @keywords internal
#' 
#' @description
#' `r lifecycle::badge('deprecated')`
#' 
#' Constructs the observed p-values from the raw observed p-values, by rounding
#' them to their nearest neighbor matching with the supports of their
#' respective CDFs (as in function `p.discrete.adjust()` of package
#' `discreteMTP`, which is no longer available on CRAN).
#' 
#' **Note**: This is an internal function and has to be called directly via
#' `:::`, i.e. `FDX:::match.pvals()`.
#' 
#' @details
#' Well computed raw p-values should already belong to their respective CDF
#' support. So this function is called at the beginning of [`discrete.GR()`],
#' [`discrete.LR()`], [`discrete.PB()`] and their respective wrappers, just in
#' case raw $p$-values may be biased.
#'
#' For each raw p-value that needs to be rounded, a warning is issued.
#'
#' @seealso
#' [`discrete.GR()`], [`discrete.LR()`], [`discrete.PB()`]
#'
#' @templateVar pCDFlist TRUE
#' @templateVar raw.pvalues TRUE
#' @templateVar pCDFlist.indices TRUE
#' @templateVar weights FALSE
#' @template param
#' 
#' @examples \dontrun{
#' toyList <- list(c(0.3,0.7,1),c(0.1,0.65,1))
#' toyRaw1 <- c(0.3,0.65)
#' match.pvals(toyList,toyRaw1)
#' toyRaw2 <- c(0.31,0.6)
#' match.pvals(toyList,toyRaw2)
#' }
#'
#' @return
#' A vector where each raw p-value has been replaced by its nearest neighbor, if
#' necessary.
#'
match.pvals <- function(pCDFlist, raw.pvalues, pCDFlist.indices = NULL){
  m <- length(raw.pvalues)
  if(!is.null(pCDFlist.indices)){
    idx <- unlist(pCDFlist.indices)
    counts <- sapply(pCDFlist.indices, length)
    pCDFlist <- rep(pCDFlist, counts)[order(idx)]
  }
  n <- length(pCDFlist)
  if(m > 0 && m == n){
    pvec <- raw.pvalues
    in.CDF <- numeric(m)
    for (k in seq_len(m)) {
      in.CDF[k] <- match(pvec[k], pCDFlist[[k]])
      if (is.na(in.CDF[k])){
        in.CDF[k] <- which.min(abs(pCDFlist[[k]] - pvec[k]))
        pvec[k] <- pCDFlist[[k]][in.CDF[k]]
        ordinal <- "th"
        if (k%%10==1) ordinal <- "st"
        if (k%%10==2) ordinal <- "nd"
        if (k%%10==3) ordinal <- "rd"
        if (k%%100-k%%10==10) ordinal <- "th"
        warning("Since ", raw.pvalues[k], 
                " is not a value of the CDF of the ",k,ordinal ," p-value,\n  the p-value is rounded to be ", 
                pCDFlist[[k]][in.CDF[k]], call. = F)
      }
    }
    return(pvec)
  }else{
    stop("'pCDFlist' and 'raw.pvalues' must have the same non-zero length")
  }
}

discrete.fdx.int <- function(
    pvec,
    pCDFlist,
    pCDFlist.indices,
    method = "GR",
    alpha = 0.05,
    zeta = 0.5,
    adaptive = TRUE,
    method.parameter = NULL,
    crit.consts = FALSE,
    threshold = 1,
    data.name = NULL
) {
  # original number of hypotheses
  n <- length(pvec)
  
  #--------------------------------------------
  #       prepare output object
  #--------------------------------------------
  input.data <- list()
  switch(
    EXPR = method, 
    GR = {
      alg <- "Discrete Guo-Romano procedure"
    },
    LR = {
      alg <- paste(
        "Discrete Lehmann-Romano procedure",
        ifelse(method.parameter, "(step-up)", "(step-down)")
      )
    },
    PB = {
      alg <- paste(
        "Discrete Poisson-Binomial procedure",
        ifelse(method.parameter, "(exact)", "(normal approximation)")
      )
    }
  )
  input.data$Method <- ifelse(!adaptive, paste("Non-adaptive", alg), alg)
  input.data$Raw.pvalues <- pvec
  if(length(pCDFlist) == n) {
    input.data$pCDFlist <- pCDFlist
  } else {
    idx <- unlist(pCDFlist.indices)
    pCDFlist.counts <- sapply(pCDFlist.indices, length)
    input.data$pCDFlist <- rep(pCDFlist, pCDFlist.counts)[order(idx)]
  }
  input.data$FDP.threshold <- alpha
  input.data$Exceedance.probability <- zeta
  input.data$Adaptive <- adaptive
  input.data$Data.name <- ifelse(
    !is.null(data.name),
    data.name,
    paste(deparse(substitute(pvec)), "and", deparse(substitute(pCDFlist)))
  )
  
  #--------------------------------------------
  #       apply p-value selection
  #--------------------------------------------
  if(threshold < 1) {
    # which p-values are not above threshold?
    select <- which(pvec <= threshold)
    # number of selected p-values
    m <- length(select)
    # filter pCDFs, indices and counts of selected p-values
    if(is.null(pCDFlist.indices)) {
      pCDFlist.indices <- as.list(select)
      pCDFlist.counts  <- rep(1, m)
      pCDFlist         <- pCDFlist[select]
    } else{
      pCDFlist.indices <- lapply(pCDFlist.indices, setdiff, seq_len(n)[-select])
      pCDFlist.counts  <- sapply(pCDFlist.indices, length)
      idx.nonempty     <- which(pCDFlist.counts > 0)
      pCDFlist         <- pCDFlist[idx.nonempty]
      pCDFlist.indices <- pCDFlist.indices[idx.nonempty]
      pCDFlist.counts  <- pCDFlist.counts[idx.nonempty]
    }
    pCDFlist.idx <- order(unlist(pCDFlist.indices))
    # rescale pCDFs
    F.thresh <- sapply(
      X = pCDFlist, 
      FUN = function(X){
        t <- which(X <= threshold)
        if(length(t)) X[max(t)] else 0
      }
    )
    pCDFlist <- sapply(
      X = seq_along(pCDFlist),
      FUN = function(k) unique(pmin(pCDFlist[[k]] / F.thresh[k], 1))
    )
    # rescale selected p-values
    pvec <- pvec[select] / rep(F.thresh, pCDFlist.counts)[pCDFlist.idx]
  } else {
    # all p-values were selected
    select <- seq_len(n)
    # number of selected p-values
    m <- n
    # use original counts (or 1 for all, if all pCDFs are unique)
    pCDFlist.counts <- if(is.null(pCDFlist.indices)) {
      rep(1, m)
    } else
      sapply(pCDFlist.indices, length)
    # F_i(1) = 1 for all i = 1, ..., n
    F.thresh <- rep(1.0, n)
  }
  
  #--------------------------------------------
  #       determine sort order and do sorting
  #--------------------------------------------
  o <- order(pvec)
  sorted.pvals <- pvec[o]
  
  #--------------------------------------------
  #       construct the vector of the overall support of the p-values
  #--------------------------------------------
  support <- unique(sort(pmin(as.numeric(unlist(pCDFlist)), 1.0)))
  
  #--------------------------------------------
  #        compute significant p-values, their
  #        indices and the number of rejections
  #--------------------------------------------
  if(crit.consts) {
    switch(
      EXPR = method,
      GR = {
        # compute critical constants
        res <- kernel_DGR_crit(pCDFlist, support, sorted.pvals, adaptive, alpha, zeta, pCDFlist.counts)
        # extract critical constants
        crit.constants <- res$crit.consts 
        # determine significant p-values
        idx.rej <- which(sorted.pvals > crit.constants)
      },
      LR = {
        # compute critical constants
        res <- kernel_DLR_crit(pCDFlist, support, sorted.pvals, adaptive, alpha, zeta, method.parameter, pCDFlist.counts)
        # extract critical constants
        crit.constants <- res$crit.consts
        # determine significant p-values
        idx.rej <- if(method.parameter) 
          which(sorted.pvals <= crit.constants) else 
            which(sorted.pvals > crit.constants)
      },
      PB = {
        # compute critical constants
        res <- kernel_DPB_crit(pCDFlist, support, sorted.pvals, adaptive, alpha, zeta, method.parameter, pCDFlist.counts)
        # extract critical constants
        crit.constants <- res$crit.consts 
        # determine significant p-values
        idx.rej <- which(sorted.pvals > crit.constants)
      }
    )
  } else {
    switch(
      EXPR = method,
      GR = {
        # compute transformed sorted p-values
        res <- kernel_DGR_fast(pCDFlist, sorted.pvals, adaptive, alpha, pCDFlist.counts)$pval.transf
        # determine significant (transformed) p-values
        idx.rej <- which(res > zeta)
      },
      LR = {
        # compute transformed sorted p-values
        res <- kernel_DLR_fast(pCDFlist, sorted.pvals, adaptive, alpha, method.parameter, zeta, support, pCDFlist.counts)
        # determine significant (transformed) p-values
        idx.rej <- if(method.parameter) 
          which(res <= zeta * (floor(seq_along(res) * alpha) + 1)) else 
            which(res > zeta * (floor(seq_len(m) * alpha) + 1)) 
      },
      PB = {
        # compute transformed sorted p-values
        res <- kernel_DPB_fast(pCDFlist, sorted.pvals, adaptive, alpha, method.parameter, pCDFlist.counts)
        # determine significant (transformed) p-values
        idx.rej <- which(res > zeta)
      }
    )
  }
  
  k <- length(idx.rej)
  
  if(method == "LR" && method.parameter) {
    if(k > 0) {
      m.rej <- max(idx.rej)
      # determine significant (observed) p-values in sorted.pvals
      idx.rej  <- which(pvec <= sorted.pvals[m.rej]) 
      pvec.rej <- input.data$Raw.pvalues[select][idx.rej]
    } else {
      m.rej    <- 0
      idx.rej  <- integer(0)
      pvec.rej <- numeric(0)
    }
  } else {
    if(k > 0) {
      m.rej <- min(idx.rej) - 1
      if(m.rej) {
        # determine significant (observed) p-values in sorted.pvals
        idx.rej  <- which(pvec <= sorted.pvals[m.rej])
        pvec.rej <- input.data$Raw.pvalues[select][idx.rej]
      } else {
        idx.rej  <- numeric(0)
        pvec.rej <- numeric(0)
      }
    } else {
      m.rej    <- m
      idx.rej  <- seq_len(m)
      pvec.rej <- input.data$Raw.pvalues[select]
    }
  }
  
  #--------------------------------------------
  #       create output object
  #--------------------------------------------
  # rejections
  output <- list(
    Rejected     = pvec.rej,
    Indices      = select[idx.rej],
    Num.rejected = m.rej
  )
  
  # adjusted p-values
  if(method != "LR" || (method == "LR" && !method.parameter)){
    if(crit.consts){
      res <- res$pval.transf
    }
    # compute adjusted p-values
    pv.adj <- switch(
      EXPR = method,
      GR = cummax(pmin(res, 1)),
      LR = cummax(pmin(res / (floor(alpha * seq_len(m)) + 1), 1)),
      PB = cummax(pmin(res, 1))
    ) 
    # add adjusted p-values to output list
    ro <- order(o)
    output$Adjusted          <- numeric(n)
    output$Adjusted[select]  <- pv.adj[ro]
    output$Adjusted[-select] <- NA
  }
  # add critical values to output list
  if(crit.consts) 
    output$Critical.values <- c(crit.constants, rep(NA, n - m))
  
  # include selection data, if selection was applied
  if(threshold < 1) {
    output$Select <- list()
    output$Select$Threshold <- threshold
    output$Select$Effective.Thresholds <- F.thresh
    output$Select$Pvalues <- input.data$Raw.pvalues[select]
    output$Select$Indices <- select
    output$Select$Scaled <- pvec
    output$Select$Number <- m
  }
  
  # original test data (often included, e.g. when using 'binom.test()')
  output$Data <- input.data
  
  class(output) <- c("FDX", class(output))
  return(output)
}

continuous.fdx.int <- function(
    pvec,
    method      = "GR",
    alpha       = 0.05,
    zeta        = 0.5,
    adaptive    = TRUE,
    crit.consts = FALSE,
    exact       = TRUE,
    threshold   = 1,
    data.name   = NULL
) {
  # original number of tests
  n <- length(pvec)
  
  #--------------------------------------------
  #       prepare output object
  #--------------------------------------------
  input.data <- list()
  alg <- paste0(
    "Continuous ",
    ifelse(method == "GR", "Guo", "Lehmann"),
    "-Romano procedure"
  )
  input.data$Method <- ifelse(!adaptive, paste("Non-adaptive", alg), alg)
  input.data$Raw.pvalues <- pvec
  input.data$FDP.threshold <- alpha
  input.data$Exceedance.probability <- zeta
  input.data$Adaptive <- adaptive
  input.data$Data.name <- ifelse(
    !is.null(data.name),
    data.name,
    deparse(substitute(pvec))
  )
  
  #--------------------------------------------
  #       apply p-value selection
  #--------------------------------------------
  if(threshold < 1) {
    # which p-values are not above threshold?
    select <- which(pvec <= threshold)
    # number of selected p-values
    m <- length(select)
    # rescale selected p-values
    F.thresh <- rep(threshold, m)
    pvec <- pvec[select] / F.thresh
  } else {
    # all p-values were selected
    select <- seq_len(n)
    # number of selected p-values
    m <- n
    # F_i(1) = 1 for all i = 1, ..., n
    F.thresh <- rep(1.0, n)
  }
  # [alpha * k], k = 1, ..., m
  a <- floor(alpha * seq_len(m)) + 1
  
  #--------------------------------------------
  #       determine sort order and do sorting
  #--------------------------------------------
  o <- order(pvec)
  sorted.pvals <- pvec[o]
  
  #--------------------------------------------
  #        compute significant p-values, their
  #        indices and the number of rejections
  #--------------------------------------------
  # compute transformed sorted p-values
  switch(
    EXP = method,
    GR =  {
      res <- if(adaptive) {
        cummax(pbinom(a - 1, m - seq_len(m) + a, pmin(1, sorted.pvals), lower.tail = FALSE))
        # equivalent: cummax(pbeta(pmin(1, sorted.pvals), a, m - seq_len(m) + 1))
      } else {
        cummax(pbinom(a - 1, m, pmin(1, sorted.pvals), lower.tail = FALSE))
        # equivalent: cummax(pbeta(pmin(1, sorted.pvals), a, m - a + 1))
      }
    },
    LR = {
      res <- if(adaptive) {
        pmin(1, cummax(sorted.pvals * (m - seq_len(m) + a) / a))
      } else {
        pmin(1, cummax(sorted.pvals * m / a))
      }
    }
  )
  # find critical constants
  if(crit.consts) {
    switch(
      EXPR = method,
      GR = {
        crit.constants <- if(adaptive) 
          qbeta(zeta, a, m - seq_len(m) + 1) else 
            qbeta(zeta, a, m - a + 1)
      },
      LR = {
        crit.constants <- if(adaptive)
          zeta * a / (m - 1:m + a) else
            zeta * a / m
      }
    )
  }
  
  # determine significant (transformed) p-values
  idx.rej <- which(res > zeta)
  
  if(length(idx.rej)) {
    m.rej <- min(idx.rej) - 1
    if(m.rej) {
      # determine significant (observed) p-values in sorted.pvals
      idx.rej  <- which(pvec <= sorted.pvals[m.rej])
      pvec.rej <- input.data$Raw.pvalues[select][idx.rej]
    } else {
      idx.rej  <- numeric(0)
      pvec.rej <- numeric(0)
    }
  } else {
    m.rej    <- m
    idx.rej  <- seq_len(m)
    pvec.rej <- input.data$Raw.pvalues[select]
  }
  
  #--------------------------------------------
  #       create output object
  #--------------------------------------------
  # rejections
  output <- list(
    Rejected     = pvec.rej,
    Indices      = select[idx.rej],
    Num.rejected = m.rej
  )
  
  # adjusted p-values
  # add adjusted p-values to output list
  ro <- order(o)
  output$Adjusted          <- numeric(n)
  output$Adjusted[select]  <- res[ro]
  output$Adjusted[-select] <- NA
  
  # add critical values to output list
  if(crit.consts) 
    output$Critical.values <- c(crit.constants, rep(NA, n - m))
  
  # include selection data, if selection was applied
  if(threshold < 1) {
    output$Select <- list()
    output$Select$Threshold <- threshold
    output$Select$Effective.Thresholds <- F.thresh
    output$Select$Pvalues <- input.data$Raw.pvalues[select]
    output$Select$Indices <- select
    output$Select$Scaled  <- pvec
    output$Select$Number  <- m
  }
  
  # original test data (often included, e.g. when using 'binom.test()')
  output$Data <- input.data
  
  class(output) <- c("FDX", class(output))
  return(output)
}

#' @importFrom PoissonBinomial ppbinom
#' @importFrom pracma fzero
#' @importFrom stats qbeta
weighted.fdx.int <- function(
    pvec,
    weights,
    method      = "GR",
    weight.meth = "AM",
    alpha       = 0.05,
    zeta        = 0.5,
    exact       = TRUE,
    crit.consts = FALSE,
    threshold   = 1,
    data.name   = NULL
) {
  #--------------------------------------------
  #       remove NAs from raw p-values
  #--------------------------------------------
  non.missing <- !is.na(pvec) & !is.na(weights)
  pvec.ok     <- pvec[non.missing]
  weights.ok  <- weights[non.missing]
  
  #--------------------------------------------
  #       prepare output object
  #--------------------------------------------
  input.data <- list()
  input.data$Method <- switch(
    EXPR = method, 
    GR = "Weighted Guo-Romano procedure",
    LR = "Weighted Lehmann-Romano procedure",
    PB = paste(
      "Weighted Poisson-Binomial procedure",
      ifelse(exact, "(exact)", "(normal approximation)")
    )
  )
  input.data$Raw.pvalues <- pvec.ok
  input.data$Weights <- weights.ok
  input.data$FDP.threshold <- alpha
  input.data$Exceedance.probability <- zeta
  input.data$Weighting <- ifelse(weight.meth == "AM", "Arithmetic Mean", "Geometric Mean")
  input.data$Data.name <- ifelse(
    !is.null(data.name),
    data.name,
    paste(deparse(substitute(pvec)), "and", deparse(substitute(weights)))
  )
  
  # original number of tests
  n   <- length(pvec.ok)
  
  #--------------------------------------------
  #       apply p-value selection
  #--------------------------------------------
  if(threshold < 1) {
    # which p-values are not above threshold?
    select <- which(pvec.ok <= threshold)
    # number of selected p-values
    m <- length(select)
    # rescale selected p-values
    F.thresh <- rep(threshold, m)
    pvec <- pvec.ok[select] / F.thresh
    # select weights
    weights <- weights.ok[select]
  } else {
    # all p-values were selected
    select <- seq_len(n)
    # number of selected p-values
    m <- n
    # use all (OK) p-values
    pvec <- pvec.ok
    # F_i(1) = 1 for all i = 1, ..., n
    F.thresh <- rep(1.0, n)
    # select weights
    weights <- weights.ok
  }
  # [alpha * k], k = 1, ..., m
  a <- floor(alpha * seq_len(m)) + 1
  # adaptive size
  m.l <- m - seq_len(m) + a
  
  #--------------------------------------------
  #       Rescale weights
  #--------------------------------------------
  weights.rescaled   <- weights / mean(weights)
  weights.decreasing <- sort(weights.rescaled, decreasing = TRUE)
  
  #--------------------------------------------
  #       Compute weighted p-values
  #--------------------------------------------
  qvalues <- if(weight.meth == "AM")
    pmin(pvec / weights.rescaled, 1) else
      geom_weight(pvec, 1 / weights.rescaled)
  
  #--------------------------------------------
  #       Determine sort order and do sorting
  #--------------------------------------------
  o <- order(qvalues)
  ro <- order(o)
  qvalues.sorted <- qvalues[o]
  
  #--------------------------------------------
  #        Compute [wLR-AM] or [wLR-GM] significant p-values,
  #        their indices and the number of rejections
  #--------------------------------------------
  # Compute transformed weighted p-values
  res <- switch(
    EXPR = method,
    GR = pmin(1, cummax(kernel_wGR_fast(qvalues.sorted, weights.decreasing, alpha, weight.meth == "GM"))),
    LR = pmin(1, cummax(kernel_wLR_fast(qvalues.sorted, weights.decreasing, alpha, weight.meth == "GM"))),
    PB = pmin(1, cummax(kernel_wPB_fast(qvalues.sorted, weights.decreasing, alpha, weight.meth == "GM", exact)))
  )
  
  # determine significant (transformed) p-values
  idx.rej <- which(res > zeta)
  if(length(idx.rej)) {
    m.rej <- min(idx.rej) - 1
    if (m.rej) {
      # determine significant (observed) p-values in sorted.pvals
      idx.rej  <- which(qvalues <= qvalues.sorted[m.rej])
      pvec.rej <- pvec.ok[select][idx.rej]
    } else {
      idx.rej  <- numeric(0)
      pvec.rej <- numeric(0)
    }
  } else {
    m.rej    <- m
    idx.rej  <- seq_len(m)
    pvec.rej <- pvec.ok[select]
  }
  
  # find critical values
  if(crit.consts) {
    switch(
      EXPR = method,
      GR = {
        if(weight.meth == "AM") {
          crit.constants <- numeric(m)
          for(l in seq_len(m)) 
            crit.constants[l] <- fzero(
              fun = function(x, w, b, z) {
                ppbinom(b, pmin(w * x, 1), method = "GeoMeanCounter", lower.tail = FALSE) - z
              }, 
              x = c(crit.constants[max(1, l - 1)], 1),#c(max(crit.constants), 1),
              w = weights.decreasing[seq_len(m.l[l])],
              b = a[l] - 1,
              z = zeta,
              tol = .Machine$double.neg.eps
            )$x
        } else {
          crit.constants <- geom_weight(
            qbeta(zeta, a, m.l - a + 1),
            (seq_len(m) / cumsum(weights.decreasing))[m.l]
          )
        }
      },
      LR = {
        if(weight.meth == "AM") {
          crit.constants <- zeta * a / cumsum(weights.decreasing)[m.l]
        } else {
          crit.constants <- numeric(m)
          for(l in seq_len(m)) 
            crit.constants[l] <- fzero(
              fun = function(x, w, b, z) {
                sum(1 - (1 - x)^w) / b - z
              },
              x = c(crit.constants[max(1, l - 1)], 1),#c(max(crit.constants), 1),
              w = weights.decreasing[seq_len(m.l[l])],
              b = a[l],
              z = zeta,
              tol = .Machine$double.neg.eps
            )$x
        }
      },
      PB = {
        meth.pb <- ifelse(exact, "DivideFFT", "RefinedNormal")
        crit.constants <- numeric(m)
        if(weight.meth == "AM") {
          for(l in seq_len(m)) 
            crit.constants[l] <- fzero(
              fun = function(x, w, b, z){
                ppbinom(b, pmin(w * x, 1), method = meth.pb, lower.tail = FALSE) - z
              },
              x = c(crit.constants[max(1, l - 1)], 1),#c(max(crit.constants) * 0.9, 1),
              w = weights.decreasing[seq_len(m.l[l])],
              b = a[l] - 1,
              z = zeta,
              tol = .Machine$double.neg.eps
            )$x
        } else {
          for(l in seq_len(m))
            crit.constants[l] <- fzero(
              fun = function(x, w, b, z){
                ppbinom(b, geom_weight(rep(x, length(w)), w), method = meth.pb, lower.tail = FALSE) - z
              },
              x = c(crit.constants[max(1, l - 1)], 1),#c(max(crit.constants) * 0.9, 1),
              w = weights.decreasing[seq_len(m.l[l])],
              b = a[l] - 1,
              z = zeta,
              tol = .Machine$double.neg.eps
            )$x
        }
      }
    )
  }
  
  #--------------------------------------------
  #       create output object
  #--------------------------------------------
  # rejections
  output <- list(
    Rejected     = pvec.rej,
    Indices      = select[idx.rej],
    Num.rejected = m.rej
  )
  
  # add adjusted p-values and weighted q-values to output list
  ro <- order(o)
  output$Adjusted          <- numeric(n)
  output$Adjusted[select]  <- res[ro]
  output$Adjusted[-select] <- NA
  output$Weighted          <- numeric(n)
  output$Weighted[select]  <- qvalues
  output$Weighted[-select] <- NA
  
  # add critical values to output list
  if(crit.consts) 
    output$Critical.values <- c(crit.constants, rep(NA, n - m))
  
  # include selection data, if selection was applied
  if(threshold < 1) {
    output$Select <- list()
    output$Select$Threshold <- threshold
    output$Select$Effective.Thresholds <- F.thresh
    output$Select$Pvalues <- input.data$Raw.pvalues[select]
    output$Select$Indices <- select
    output$Select$Scaled  <- pvec
    output$Select$Number  <- m
  }
  
  # original test data (often included, e.g. when using 'binom.test()')
  output$Data <- input.data
  
  class(output) <- c("FDX", class(output))
  return(output)
}