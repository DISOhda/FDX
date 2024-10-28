#' @name summary.FDX
#' 
#' @title
#' Summarizing Discrete FDX Results
#' 
#' @description
#' `summary` method for class `FDX`
#'
#' @param object       object of class "`FDX`".
#' @param x            object of class "`summary.FDX`".
#' @param max          numeric or `NULL`, specifying the maximum number of
#'                     *rows* of the p-value table to be printed; if `NULL`
#'                     (the default), `getOption("max.print")` is used.
#' @param ...          further arguments passed to or from other methods.
#' 
#' @details
#' `summary.FDX` objects include all data of an `FDX` class
#' object, but also include an additional table which includes the raw p-values,
#' their indices, the respective critical values (if present), the adjusted
#' p-values (if present) and a logical column to indicate rejection. The table
#' is sorted in ascending order by the raw p-values.
#' 
#' `print.summary.FDX` simply prints the same output as
#' `print.FDX`, but also prints the p-value table.
#'
#' @return
#' `summary.FDX` computes and returns a list that includes all the
#' data of an input `FDX`, plus
#' \item{Table}{a `data.frame`, sorted by the raw p-values, that contains
#'              the indices, that raw p-values themselves, their respective
#'              critical values (if present), their adjusted p-values (if
#'              present) and a logical column to indicate rejection.}
#'              
#' `print.summary.FDX` returns that object invisibly.
#' 
#' @template example
#' @examples
#' 
#' # DGR with critical values; using test results object
#' DGR.crit <- DGR(test.results, critical.values = TRUE)
#' # create summary
#' DGR.crit.summary <- summary(DGR.crit)
#' # print summary
#' print(DGR.crit.summary)
#' 
#' @method summary FDX
#' @export
## S3 method for class 'FDX'
summary.FDX <- function(object, ...) {
  if(!("FDX" %in% class(object)))
    return(summary(object))
  
  # determine if selection was performed
  select <- exists('Select', object)
  if(select) m <- object$Select$Number
  
  # determine if weighting was performed
  weight <- grepl("Weighted", object$Data$Method)
  
  # number of tests
  n <- length(object$Data$Raw.pvalues)
  # ordered indices
  i <- seq_len(n)
  # determine for each p-value if its corresponding null hypothesis is rejected
  r <- i %in% object$Indices
  
  # create summary table
  tab <- data.frame(
    Index = i,
    P.value = object$Data$Raw.pvalues
  )
  
  # add selection (T/F) and scaled p-values
  if(select) {
    tab$Selected <- i %in% object$Select$Indices
    tab$Scaled <- NA
    tab$Scaled[tab$Selected] <- object$Select$Scaled
  }
  
  # add weights and weighted p-values; determine order
  if(weight) {
    # add weights and weighted p-values
    tab$Weights <- out$Data$Weights
    tab$Weighted <- out$Weighted
    # determine order of weighted p-values
    o <- order(tab$Weighted, object$Data$Raw.pvalues)
  } else {
    if(select) {
      # determine order of scaled selected p-values
      o <- order(tab$Scaled, object$Data$Raw.pvalues)
    } else
      # determine order of raw p-values
      o <- order(object$Data$Raw.pvalues)
  }
  
  # sort rows in ascending order
  tab <- tab[o, ]
  
  # add critical constants (if present)
  if(exists('Critical.values', where = object)) {
    tab$Critical.value <- object$Critical.values
  }
  
  # add adjusted p-values (if present)
  if(exists('Adjusted', where = object)) {
    tab$Adjusted <- object$Adjusted[o]
  }
  
  # add rejection decisions
  tab <- data.frame(tab, 'Rejected' = r[o])
  # if row names are numbers, rearrange them to represent sort order
  if(all(rownames(tab) == tab$Index)) rownames(tab) <- i
  
  # return output object
  out <- c(object, list(Table = tab))
  class(out) <- c("summary.FDX", class(object)) # basically a 'FDX' object, but with a summary table (just like 'lm' and 'summary.lm' classes)
  return(out)
}

#' @rdname summary.FDX
#' @method print summary.FDX
#' @export
## S3 method for class 'summary.FDX'
print.summary.FDX <- function(x, max = NULL, ...) {
  if(!("summary.FDX" %in% class(x)))
    return(print(x))
  
  # determine number of tests
  m <- ncol(x$Table)
  
  # print 'FDX' part of the object
  print.FDX(x)
  
  # rows to print: number of rejections + 5 (if not requested otherwise)
  max <- if(!is.null(max)) m * max else getOption("max.print")
  
  # print additional summary table
  print(x$Table, max = max, ...)
  
  cat("\n")
  invisible(x)
}
