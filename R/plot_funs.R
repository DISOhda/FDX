#' @name hist.FDX
#' 
#' @title
#' Histogram of Raw P-Values
#' 
#' @description
#' Computes a histogram of the raw p-values of a `FDX` object.
#' 
#' @param x          object of class `FDX`.
#' @param breaks     as in [`graphics::hist()`]; here, the Friedman-Diaconis
#'                   algorithm (`"FD"`) is used as default.
#' @param mode       single character string specifying for which $p$-values the
#'                   histogram is to be generated; must be one of `"raw"`,
#'                   `"selected"` or `"weighted"`.
#' @param ...        further arguments to [`graphics::hist()`] or 
#'                   [`graphics::plot.histogram()`], respectively.
#' 
#' @details
#' If `x` does not contain results of a weighting or selection approach, a
#' warning is issued and a histogram of the raw $p$-values is drawn. 
#' 
#' @return
#' An object of class `histogram`.
#' 
#' @template example
#' @examples
#' 
#' # DGR
#' DGR <- DGR(raw.pvalues, pCDFlist)
#' # histogram of raw p-values
#' hist(DGR)
#' 
#' # arithmetic-weighted GR (using 1 - raw.pvalues as weights)
#' wGR <- wGR.AM(raw.pvalues, 1 - raw.pvalues)
#' # histogram of raw p-values
#' hist(wGR)
#' # histogram of weighted p-values
#' hist(wGR, mode = "weighted")
#' 
#' @importFrom graphics hist
#' @export
hist.FDX <- function(
    x,
    breaks = "FD",
    mode = c("raw", "selected", "weighted"),
    ...
) {
  if(!("FDX" %in% class(x)))
    stop("'x' must be an object of class FDX")
  
  mode <- match.arg(tolower(mode), c("raw", "selected", "weighted"))
  
  # determine if selection was performed
  select <- exists('Select', x)
  select.mode <- mode == "selected"
  # warn if histogram is for selected p-values that were never computed
  if(select.mode && !select)
    warning("No selected p-values present in object. Using raw p-values.")
  
  # determine if weighting was performed
  weight <- exists('Weighted', x)
  weight.mode <- mode == "weighted"
  # warn if histogram is for weighted p-values that were never computed
  if(weight.mode && !weight)
    warning("No weighted p-values present in object. Using raw p-values.")
  
  # get values of ...-arguments
  lst <- list(...)
  
  # p-value type
  pv.type <- ifelse(
    test = select.mode, 
    yes = ifelse(select, "Selected and Scaled", "Raw"),
    no = ifelse(
      test = weight.mode,
      yes = ifelse(
        test = weight, 
        yes = paste0(
          ifelse(select, "Selected, Scaled and ", ""),
          "Weighted"
        ), 
        no = "Raw"),
      no = "Raw"
    )
  )
  
  # labels
  if(!exists('main', where = lst)) 
    lst$main <- paste("Histogram of", pv.type, "P-Values")
  if(!exists('xlab', where = lst)) 
    lst$xlab <- paste(pv.type, "P-Values")
  
  # add 'x' and 'breaks' to 'lst' (for 'do.call')
  lst$breaks <- breaks
  lst$x <- switch(
    EXPR = mode,
    raw = x$Data$Raw.pvalues,
    selected = if(select) x$Select$Scaled else x$Data$Raw.pvalues,
    weighted = if(weight) x$Weighted else x$Data$Raw.pvalues
  )
  
  # call 'hist'
  r <- do.call(hist, lst)
  
  r$xname <- deparse(substitute(x))
  
  plot <- if(exists('plot', lst)) lst$plot else TRUE
  
  if(plot) return(invisible(r)) else return(r)
}


#' @name plot.FDX
#' 
#' @title
#' Plot Method for `FDX` objects
#' 
#' @description
#' Plots raw $p$-values of a `FDX` object and highlights rejected and
#' non-rejected $p$-values. If present, the critical values are plotted, too.
#' 
#' @param x          an object of class "`FDX`".
#' @param col        numeric or character vector of length 3 indicating the
#'                   colors of the \enumerate{
#'                     \item rejected $p$-values
#'                     \item non-rejected $p$-values
#'                     \item critical values (if present).
#'                   }
#' @param pch        numeric or character vector of length 3 indicating the
#'                   point characters of the \enumerate{
#'                     \item rejected $p$-values
#'                     \item non-rejected $p$-values
#'                     \item critical values (if present and `type.crit`
#'                           is a plot type like `'p'`, `'b'` etc.).
#'                   }
#' @param lwd        numeric vector of length 3 indicating the thickness of the
#'                   points and lines; defaults to current `par()$lwd` setting
#'                   for all components.
#' @param type.crit  single character giving the type of plot desired for the
#'                   critical values (e.g.: `'p'`, `'l'` etc; see
#'                   [`graphics::plot.default()`]).
#' @param legend     if `NULL`, no legend is plotted; otherwise expecting a
#'                   character string like `"topleft"` etc. or a numeric vector
#'                   of two elements indicating (x, y) coordinates.
#' @param cex        numeric vector of length 3 indicating the size of point
#'                   characters or lines of the \enumerate{
#'                     \item rejected p-values
#'                     \item accepted p-values
#'                     \item critical values (if present).
#'                   }
#'                   defaults to current `par()$cex` setting for all components.
#' @param ...        further arguments to [`graphics::plot.default()`].
#' 
#' @details
#' If `x` contains results of a weighted approach, the Y-axis of the plot
#' is derived from the weighted p-values. Otherwise, it is constituted by the
#' raw ones. 
#' 
#' @template example
#' @examples
#' 
#' # DLR without critical values; using extracted p-values and supports
#' DLR.sd.fast <- DLR(raw.pvalues, pCDFlist)
#' # plot with default settings
#' plot(DLR.sd.fast)
#'
#' # DLR (step-up) with critical values; using test results object
#' DLR.su.crit <- DLR(test.results, direction = "su", critical.values = TRUE)
#' # limited plot range
#' plot(DLR.su.crit, xlim = c(1, 5), ylim = c(0, 0.4))
#'
#' # DPB without critical values; using test results object
#' DPB.fast <- DPB(test.results)
#' # limited plot range, custom colors, line widths and point symbols, top-left legend 
#' plot(DPB.fast, col = c(2, 4), pch = c(2, 3), lwd = c(2, 2), 
#'      legend = "topleft", xlim = c(1, 5), ylim = c(0, 0.4))
#' 
#' # DGR with critical values; using extracted p-values and supports
#' DGR.crit <- DGR(raw.pvalues, pCDFlist, critical.values = TRUE)
#' # additional customized plot parameters
#' plot(DGR.crit, col = c(2, 4, 1), pch = c(1, 1, 4), lwd = c(1, 1, 2), 
#'      type.crit = 'o', legend = c(1, 0.4), lty = 1, xlim = c(1, 5), 
#'      ylim = c(0, 0.4), cex = c(3, 3, 2))
#' 
#' @importFrom graphics legend lines par plot points
#' @importFrom checkmate assert assert_string check_character check_choice check_numeric
#' @importFrom stats na.omit
#' @export
plot.FDX <- function(
    x,
    col = c(2, 4, 1),
    pch = c(20, 20, 17),
    lwd = rep(par()$lwd, 3),
    cex = rep(par()$cex, 3),
    type.crit = 'b',
    legend = NULL,
    ...
){
  if(!("FDX" %in% class(x)))
    stop("'x' must be an object of class FDX")
  
  # make sure 'col' includes integers or color strings
  assert(
    check_character(col, min.len = 1, max.len = 3, null.ok = TRUE),
    check_numeric(col, min.len = 1, max.len = 3, null.ok = TRUE)
  )
  
  # make sure 'pch' includes integers or color strings
  assert(
    check_character(col, min.len = 1, max.len = 3, null.ok = TRUE),
    check_numeric(col, upper = 25, min.len = 1, max.len = 3, null.ok = TRUE)
  )
  
  # make sure 'type.crit' is a single character string
  assert_string(type.crit, null.ok = TRUE)
  
  # make sure 'legend' is a single string or a numerical vector of two values
  assert(
    check_choice(
      x = legend, 
      choices = c("bottomright", "bottom", "bottomleft", "left",
                  "topleft", "top", "topright", "right", "center"),
      null.ok = TRUE
    ),
    check_numeric(legend, len = 2, any.missing = FALSE)
  )
  
  # determine if selection was performed
  select <- exists('Select', x)
  if(select) m <- x$Select$Number
  
  # determine if weighting was performed
  weight <- exists('Weighted', x)
  
  # determine if critical constants were computed
  critical <- exists('Critical.values', x)
  
  # number of tests
  n <- length(x$Data$Raw.pvalues)
  # number of rejected null hypotheses
  k <- x$Num.rejected
  
  # replace NAs in plot parameters with current par() settings
  col[is.na(col)] <- par()$col
  pch[is.na(pch)] <- par()$pch
  lwd[is.na(lwd)] <- par()$lwd
  cex[is.na(cex)] <- par()$cex
  
  # replicate shorter plot parameter vectors to avoid errors
  len <- 2L + as.integer(critical)
  col <- if(!is.null(col)) rep_len(col, len) else c( 2,  4,  1)[seq_len(len)]
  pch <- if(!is.null(pch)) rep_len(pch, len) else c(20, 20, 17)[seq_len(len)]
  lwd <- if(!is.null(lwd)) rep_len(lwd, len) else rep(par()$lwd, len)
  cex <- if(!is.null(cex)) rep_len(cex, len) else rep(par()$cex, len)
  
  # get values of ...-arguments
  lst <- list(...)
  
  # labels
  if(!exists('main', where = lst)) lst$main <- x$Data$Method
  if(!exists('xlab', where = lst)) {
    lst$xlab <- "Index"
    if(select) lst$xlab <- paste(lst$xlab, "(selected)")
  }
  if(!exists('ylab', where = lst)) {
    pv.type <- ifelse(
      test = select, 
      yes = paste0("selected (scaled", ifelse(weight, " and weighted)", ")")),
      no = ifelse(weight, "weighted", "raw")
    )
    lst$ylab <- paste0(
      paste("Sorted", pv.type, "p-values"),
      ifelse(critical, " / Critical values", "")
    )
  }
  
  # p-values that are to be plotted
  lst$x <- if(weight) {
    na.omit(x$Weighted)
  } else if(select) {
    x$Select$Scaled 
  } else {
    x$Data$Raw.pvalues
  }
  
  # start plotting with empty area
  lst$col <- "NA"
  do.call(plot, lst)
  
  # plot accepted p-values
  if(select) {
    idx <- which(!(x$Select$Indices %in% x$Indices))
    if(length(idx)) {
      y_acc <- if(weight) sort(na.omit(x$Weighted)[idx]) else sort(x$Select$Scaled[idx])
      x_acc <- (m - length(y_acc) + 1):m
    }
  } else {
    idx <- setdiff(seq_len(n), x$Indices)
    if(length(idx)) {
      y_acc <- if(weight) sort(na.omit(x$Weighted)[idx]) else sort(x$Data$Raw.pvalues[idx])
      x_acc <- (k + 1):n #setdiff(seq_len(n), seq_len(k))
    }
  }
  if(length(idx))
    points(x_acc, y_acc, col = col[2], pch = pch[2],
           lwd = lwd[2], cex = cex[2], ...)
  
  # plot rejected p-values
  if(x$Num.rejected) {
    if(!length(idx)) idx <- x$Num.rejected + 1
    y_rej <- if(select) {
      if(weight) sort(na.omit(x$Weighted)[-idx]) else sort(x$Select$Scaled[-idx])
    } else if(weight) {
      sort(na.omit(x$Weighted)[-idx])
    } else sort(x$Rejected)
    points(1:x$Num.rejected, y_rej, col = col[1], pch = pch[1],
           lwd = lwd[1], cex = cex[1], ...)
  }
  
  # plot critical values (if present and plotting is requested by the user)
  if(exists('Critical.values', where = x) && type.crit != 'n'){
    lines(x$Critical.values, col = col[3], lwd = lwd[3], pch = pch[3],
          type = type.crit, cex = cex[3], ...)
  }
  
  # plot legend
  if(!is.null(legend)) {
    idx <- 1:2
    lt <- rep(0, 3)
    if(critical && type.crit != "n") {
      idx <- c(idx, 3)
      if(!(type.crit %in% c("p", "o", "b"))) pch[3] <- NA
      lt[3] <- if(type.crit %in% c('b', 'l', 'o')) {
        if(exists('lty', where = lst)) lst$lty else 1
      } else 0
    }
    len <- length(legend)
    if(len <= 2 & len >= 1) {
      if(len == 1){
        x <- legend
        y <- NULL
      }else{
        x <- legend[1]
        y <- legend[2]
      }
      legend(x, y, c("Rejected", "Not rejected", "Critical values")[idx], col = c(col[idx], "darkgrey"), pch = pch[idx], lty = lt[idx], lwd = lwd[idx])
    } else warning("Expecting character string or numeric vector of one or two elements for creating a legend")
  }
}


#' @name rejection.path
#' @title Rejection Path Plot (for `FDX` objects)
#' 
#' @description
#' Displays the number of rejections of the raw p-values in a `FDX`
#' object in dependence of the exceedance probability `zeta`.
#' 
#' @param x                 object of class "`FDX`".
#' @param xlim              x axis limits of the plot. If `NULL` (default),
#'                          the \[0, 1\] range is used.
#' @param ylim              the y limits of the plot. If `NULL` (default),
#'                          the double of the median of the number of possible
#'                          rejections is used as upper limit.
#' @param main              main title. If `NULL` (default), a description
#'                          string is used.
#' @param xlab,ylab         labels for x and y axis.
#' @param verticals         logical; if `TRUE`, draw vertical lines at
#'                          steps.
#' @param pch               jump point character.
#' @param ref.show          logical; if `TRUE` a vertical reference line
#'                          is plotted, whose height is the number of
#'                          rejections of the original Benjamini-Hochberg (BH)
#'                          procedure.
#' @param ref.col           color of the reference line.
#' @param ref.lty,ref.lwd   line type and thickness for the reference line.
#' @param ...               further arguments to [`stats::plot.stepfun()`].
#' 
#' @return
#' Invisibly returns a `stepfun` object that computes the number of
#' rejectionsin dependence on the exceedance probability `zeta`.
#' 
#' @template example
#' @examples
#' 
#' # DLR without critical values; using extracted p-values and supports
#' DLR <- DLR(raw.pvalues, pCDFlist)
#' # plot number of rejections dependent on the exceedance probability zeta
#' rejection.path(DLR, xlim = c(0, 1), ref.show = TRUE, ref.col = "green", ref.lty = 4)
#' 
#' # None-adaptive DLR without critical values; using test results object
#' NDLR <- NDLR(test.results)#' 
#' # add plot for non-adaptive procedure (in red)
#' rejection.path(NDLR, col = "red", add = TRUE)
#' 
#' @importFrom graphics plot abline mtext
#' @importFrom stats stepfun ecdf p.adjust plot.stepfun
#' @importFrom methods is
#' @export
rejection.path <- function(x, xlim = NULL, ylim = NULL, main = NULL, xlab = expression(zeta), ylab = "Number of Rejections", verticals = FALSE, pch = 19, ref.show = FALSE, ref.col = "gray", ref.lty = 2, ref.lwd = 2, ...){
  if(!is(x, "FDX")) stop("'x' must be an object of class 'FDX'!")
  # number of hypotheses
  m <- length(x$Data$Raw.pvalues)
  # number of BH rejections
  num.rejections.BH <- sum(p.adjust(x$Data$Raw.pvalues, "BH") <= x$Data$FDP.threshold)
  
  # create step function
  ecdf.env <- environment(ecdf(x$Adjusted))
  stepfun.x <- stepfun(ecdf.env$x, c(0, m * ecdf.env$y))
  
  lst <- list(...)
  subt <- NULL
  
  # plot
  if(is.null(xlim)) xlim <- c(0, 1)
  if(is.null(ylim)) ylim <- c(0, min(m, 2 * stepfun.x(0.5)))
  if(is.null(main)){
    main <- bquote(bold("Rejection path for"~alpha==.(as.character(x$Data$FDP.threshold))))
    if(!exists('add', where = lst) || (exists('add', where = lst) && !lst$add)) subt <- x$Data$Method
  }
  
  plot.stepfun(stepfun.x, xlim = xlim, ylim = ylim, main = main, xlab = xlab, ylab = ylab, verticals = verticals, pch = pch, ...)
  if(!is.null(subt)) mtext(subt, line = 0.3)
  
  if(ref.show) abline(h = num.rejections.BH, col = ref.col, lty = ref.lty, lwd = ref.lwd)
  
  return(invisible(stepfun.x))
}
