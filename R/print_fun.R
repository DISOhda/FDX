#' @name print.FDX
#' 
#' @title
#' Printing FDX results
#' 
#' @description
#' Prints the results of discrete FDX analysis, stored in a `FDX`
#' S3 class object.
#' 
#' @param x          object of class `FDX`.
#' @param ...        further arguments to be passed to or from other methods.
#'                   They are ignored in this function.
#' 
#' @return
#' The respective input object is invisibly returned via `invisible(x)`. 
#' 
#' @template example
#' @examples
#' 
#' DPB.crit <- DPB(raw.pvalues, pCDFlist, critical.values = TRUE)
#' print(DPB.crit)
#' 
#' @method print FDX
#' @importFrom stats p.adjust
#' @export
## S3 method for class 'FDX'
print.FDX <- function(x, ...){
  if(!any(c("FDX", "summary.FDX") %in% class(x)))
    return(print(x))
  
  # determine if selection was performed
  select <- exists('Select', x)
  if(select) m <- x$Select$Number
  
  # number of tests
  n <- length(x$Data$Raw.pvalues)
  # number of rejected null hypotheses
  k <- x$Num.rejected
  
  if(grepl("Lehmann", x$Data$Method)) {
    k.o <- continuous.LR(
      x$Data$Raw.pvalues,
      x$Data$FDP.threshold,
      x$Data$Exceedance.probability,
      TRUE,
      FALSE
    )$Num.rejected
    orig <- "Lehmann-Romano"
  } else {
    k.o <- continuous.GR(
      x$Data$Raw.pvalues,
      x$Data$FDP.threshold,
      x$Data$Exceedance.probability,
      TRUE,
      FALSE
    )$Num.rejected
    orig <- "Guo-Romano"
  }
  
  # print title (i.e. algorithm)
  cat("\n")
  cat("\t", x$Data$Method, "\n")
  
  # print dataset name(s)
  cat("\n")
  cat("Data: ", x$Data$Data.name, "\n")
  
  # print short results overview
  if(!select) {
    cat("Number of tests =", n, "\n")
  } else {
    cat("Number of selected tests =", m, "out of", n, "\n")
    cat("Selection threshold =", x$Select$Threshold, "\n")
  }
    
  cat("Number of rejections = ", k, " when controlling FDP at level ", x$Data$FDP.threshold, " with probability ",
      x$Data$Exceedance.probability, ",\n", paste(rep(" ", 24 + nchar(as.character(k))), collapse = ""),
      "i.e. P(FDP > ", x$Data$FDP.threshold, ") <= ", x$Data$Exceedance.probability, "\n", sep = "")
  
  if(!grepl("Continuous", x$Data$Method))
    cat("Original", orig, "rejections =", k.o, "\n")
  
  cat("Original Benjamini-Hochberg rejections =", sum(p.adjust(x$Data$Raw.pvalues, "BH") <= x$Data$FDP.threshold),
      "at global FDR level", x$Data$FDP.threshold, "\n")
  
  if(k && !select) {
    if(!grepl("Weighted", x$Data$Method))
      cat("Largest rejected p value: ", max(x$Rejected), "\n")
    else
      cat("Largest rejected weighted p value: ", max(x$Weighted[x$Indices]), "\n")
  }
  
  cat("\n")
  invisible(x)
}
