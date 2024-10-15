discrete.GR <- function(raw.pvalues, pCDFlist, alpha = 0.05, zeta = 0.5, adaptive = TRUE, critical.values = FALSE) {
  #--------------------------------------------
  #       check arguments
  #--------------------------------------------
  if(is.null(alpha) || is.na(alpha) || !is.numeric(alpha) || alpha < 0 || alpha > 1)
    stop("'alpha' must be a probability between 0 and 1!")
  
  if(is.null(zeta) || is.na(zeta) || !is.numeric(zeta) || zeta < 0 || zeta > 1)
    stop("'zeta' must be a probability between 0 and 1!")
  
  m <- as.integer(length(raw.pvalues))
  if(m != length(pCDFlist)) stop("The lengths of 'raw.pvalues' and 'pCDFlist' must be equal!") 
  #--------------------------------------------
  #       prepare p-values for processing
  #--------------------------------------------
  pvec <- match.pvals(raw.pvalues, pCDFlist)
  #--------------------------------------------
  #       Determine sort order and do sorting
  #--------------------------------------------
  o <- order(pvec)
  sorted.pvals <- pvec[o]
  #--------------------------------------------
  #       construct the vector of all values of all supports of the p-values
  #--------------------------------------------
  pv.list.all <- sort(unique(as.numeric(unlist(pCDFlist))))
  #--------------------------------------------
  #        Compute [HSU] or [HSD] significant p-values,
  #        their indices and the number of rejections
  #--------------------------------------------
  if(critical.values){
    # compute transformed support
    y <- kernel_DGR_crit(pCDFlist, pv.list.all, sorted.pvals, adaptive, alpha, zeta)
    # find critical constants
    crit.constants <- y$crit.consts
    idx <- which(sorted.pvals > crit.constants)
  }else{
    # compute transformed sorted p-values
    y <- kernel_DGR_fast(pCDFlist, sorted.pvals, adaptive, alpha)$pval.transf
    # determine significant (transformed) p-values
    idx <- which(y > zeta)
  }
  if(length(idx)){
    m.rej <- min(idx) - 1
    if(m.rej){
      # determine significant (observed) p-values in sorted.pvals
      idx <- which(pvec <= sorted.pvals[m.rej])
      pvec.rej <- raw.pvalues[idx]
    }
    else{
      idx <- numeric(0)
      pvec.rej <- numeric(0)
    }
  }
  else{
    m.rej <- m
    idx <- 1:m
    pvec.rej <- raw.pvalues
  }
  #--------------------------------------------
  #       Create output S3 object
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = idx, Num.rejected = m.rej)
  
  if(critical.values){
    y <- y$pval.transf
  }
  # compute adjusted p-values
  pv.adj = cummax(pmin(y, 1))
  # add adjusted p-values to output list
  ro <- order(o)
  output$Adjusted = pv.adj[ro]
  
  # add critical values to output list
  if(critical.values) output$Critical.values <- crit.constants
  
  # include details of the used algorithm as strings
  alg <- "Discrete Guo-Romano procedure"
  output$Method <- if(!adaptive) paste("Non-Adaptive", alg) else alg
  output$FDP.threshold <- alpha
  output$Exceedance.probability <- zeta
  output$Adaptive <- adaptive
  
  # original test data (often included, e.g. when using 'binom.test()')
  output$Data <- list()
  output$Data$raw.pvalues <- raw.pvalues
  output$Data$pCDFlist <- pCDFlist
  # object names of the data as strings
  output$Data$data.name <- paste(deparse(substitute(raw.pvalues)), "and", deparse(substitute(pCDFlist)))
  
  class(output) <- "FDX"
  return(output)
}

discrete.LR <- function(raw.pvalues, pCDFlist, alpha = 0.05, zeta = 0.5, direction = "sd", adaptive = TRUE, critical.values = FALSE){
  #--------------------------------------------
  #       check arguments
  #--------------------------------------------
  if(is.null(alpha) || is.na(alpha) || !is.numeric(alpha) || alpha < 0 || alpha > 1)
    stop("'alpha' must be a probability between 0 and 1!")
  
  if(is.null(zeta) || is.na(zeta) || !is.numeric(zeta) || zeta < 0 || zeta > 1)
    stop("'zeta' must be a probability between 0 and 1!")
  
  m <- as.integer(length(raw.pvalues))
  if(m != length(pCDFlist)) stop("The lengths of 'raw.pvalues' and 'pCDFlist' must be equal!") 
  #--------------------------------------------
  #       prepare p-values for processing
  #--------------------------------------------
  pvec <- match.pvals(raw.pvalues, pCDFlist)
  #--------------------------------------------
  #       Determine sort order and do sorting
  #--------------------------------------------
  o <- order(pvec)
  sorted.pvals <- pvec[o]
  #--------------------------------------------
  #       construct the vector of all values of all supports of the p-values
  #--------------------------------------------
  pv.list.all <- sort(unique(as.numeric(unlist(pCDFlist))))
  #--------------------------------------------
  #        Compute [HSU] or [HSD] significant p-values,
  #        their indices and the number of rejections
  #--------------------------------------------
  direction <- match.arg(direction, c("su", "sd"))
  if(direction == "su"){
    # SU case
    if(critical.values){
      # compute transformed support
      y <- kernel_DLR_crit(pCDFlist, pv.list.all, sorted.pvals, adaptive, alpha, zeta, TRUE)
      # find critical constants
      crit.constants <- y$crit.consts
      idx <- which(sorted.pvals <= crit.constants)
    }
    else{
      # compute transformed observed p-values
      y <- kernel_DLR_fast(pCDFlist, sorted.pvals, adaptive, alpha, TRUE, zeta, pv.list.all)
      # determine significant (transformed) p-values
      if(length(y)){
        idx <- which(y <= zeta * (floor(1:length(y) * alpha) + 1))
      }else{
        idx <- integer(0)
      }
    }
    if(length(idx)){
      m.rej <- max(idx)
      # determine significant (observed) p-values in sorted.pvals
      idx <- which(pvec <= sorted.pvals[m.rej]) 
      pvec.rej <- raw.pvalues[idx]
    }
    else{
      m.rej <- 0
      idx <- integer(0)
      pvec.rej <- numeric(0)
    }
  }
  else{
    # SD case
    if(critical.values){
      # compute transformed support
      y <- kernel_DLR_crit(pCDFlist, pv.list.all, sorted.pvals, adaptive, alpha, zeta, FALSE)
      # find critical constants
      crit.constants <- y$crit.consts
      idx <- which(sorted.pvals > crit.constants)
    }
    else{
      # compute transformed sorted p-values
      y <- kernel_DLR_fast(pCDFlist, sorted.pvals, adaptive, alpha, FALSE)
      # determine significant (transformed) p-values
      idx <- which(y > zeta * (floor(1:m * alpha) + 1)) 
    }
    if(length(idx)){
      m.rej <- min(idx) - 1
      if(m.rej){
        # determine significant (observed) p-values in sorted.pvals
        idx <- which(pvec <= sorted.pvals[m.rej])
        pvec.rej <- raw.pvalues[idx]
      }
      else{
        idx <- numeric(0)
        pvec.rej <- numeric(0)
      }
    }
    else{
      m.rej <- m
      idx <- 1:m
      pvec.rej <- raw.pvalues
    }
  }
  #--------------------------------------------
  #       Create output S3 object
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = idx, Num.rejected = m.rej)
  if(direction == "sd"){
    if(critical.values){
      y <- y$pval.transf
    }
    # compute adjusted p-values
    pv.adj = cummax(pmin(y / (floor(alpha * 1:m) + 1), 1))
    # add adjusted p-values to output list
    ro <- order(o)
    output$Adjusted = pv.adj[ro]
  }
  # add critical values to output list
  if(critical.values) output$Critical.values = crit.constants
  
  # include details of the used algorithm as strings
  alg <- "Discrete Lehmann-Romano procedure"
  alg <- if(!adaptive) paste("Non-Adaptive", alg) else alg
  output$Method <- paste(alg, switch(direction, su = "(step-up)", sd = "(step-down)"))
  output$FDP.threshold <- alpha
  output$Exceedance.probability <- zeta
  output$Adaptive <- adaptive
  
  # original test data (often included, e.g. when using 'binom.test()')
  output$Data <- list()
  output$Data$raw.pvalues <- raw.pvalues
  output$Data$pCDFlist <- pCDFlist
  # object names of the data as strings
  output$Data$data.name <- paste(deparse(substitute(raw.pvalues)), "and", deparse(substitute(pCDFlist)))
  
  class(output) <- "FDX"
  return(output)
}

discrete.PB <- function(raw.pvalues, pCDFlist, alpha = 0.05, zeta = 0.5, adaptive = TRUE, critical.values = FALSE, exact = TRUE){
  #--------------------------------------------
  #       check arguments
  #--------------------------------------------
  if(is.null(alpha) || is.na(alpha) || !is.numeric(alpha) || alpha < 0 || alpha > 1)
    stop("'alpha' must be a probability between 0 and 1!")
  
  if(is.null(zeta) || is.na(zeta) || !is.numeric(zeta) || zeta < 0 || zeta > 1)
    stop("'zeta' must be a probability between 0 and 1!")
  
  m <- as.integer(length(raw.pvalues))
  if(m != length(pCDFlist)) stop("The lengths of 'raw.pvalues' and 'pCDFlist' must be equal!") 
  #--------------------------------------------
  #       prepare p-values for processing
  #--------------------------------------------
  pvec <- match.pvals(raw.pvalues, pCDFlist)
  #--------------------------------------------
  #       Determine sort order and do sorting
  #--------------------------------------------
  o <- order(pvec)
  sorted.pvals <- pvec[o]
  #--------------------------------------------
  #       construct the vector of all values of all supports of the p-values
  #--------------------------------------------
  pv.list.all <- sort(unique(as.numeric(unlist(pCDFlist))))
  #--------------------------------------------
  #        Compute [HSU] or [HSD] significant p-values,
  #        their indices and the number of rejections
  #--------------------------------------------
  if(critical.values){
    # compute transformed support
    y <- kernel_DPB_crit(pCDFlist, pv.list.all, sorted.pvals, adaptive, alpha, zeta, exact)
    # find critical constants
    crit.constants <- y$crit.consts
    idx <- which(sorted.pvals > crit.constants)
  }else{
    # compute transformed sorted p-values
    y <- kernel_DPB_fast(pCDFlist, sorted.pvals, adaptive, alpha, exact)
    # determine significant (transformed) p-values
    idx <- which(y > zeta)
  }
  if(length(idx)){
    m.rej <- min(idx) - 1
    if(m.rej){
      # determine significant (observed) p-values in sorted.pvals
      idx <- which(pvec <= sorted.pvals[m.rej])
      pvec.rej <- raw.pvalues[idx]
    }
    else{
      idx <- numeric(0)
      pvec.rej <- numeric(0)
    }
  }
  else{
    m.rej <- m
    idx <- 1:m
    pvec.rej <- raw.pvalues
  }
  #--------------------------------------------
  #       Create output S3 object
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = idx, Num.rejected = m.rej)
  
  if(critical.values){
    y <- y$pval.transf
  }
  # compute adjusted p-values
  pv.adj = cummax(pmin(y, 1))
  # add adjusted p-values to output list
  ro <- order(o)
  output$Adjusted = pv.adj[ro]
  
  # add critical values to output list
  if(critical.values) output$Critical.values = crit.constants
  
  # include details of the used algorithm as strings
  alg <- "Discrete Poisson-Binomial procedure"
  alg <- if(exact) paste(alg, "(exact)") else paste(alg, "(normal approximation)")
  output$Method <- if(!adaptive) paste("Non-Adaptive", alg) else alg
  output$FDP.threshold <- alpha
  output$Exceedance.probability <- zeta
  output$Adaptive <- adaptive
  
  # original test data (often included, e.g. when using 'binom.test()')
  output$Data <- list()
  output$Data$raw.pvalues <- raw.pvalues
  output$Data$pCDFlist <- pCDFlist
  # object names of the data as strings
  output$Data$data.name <- paste(deparse(substitute(raw.pvalues)), "and", deparse(substitute(pCDFlist)))
  
  class(output) <- "FDX"
  return(output)
}
