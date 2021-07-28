#' <%=ifelse(exists("pCDFlist") && pCDFlist, "@param pCDFlist a list of the supports of the CDFs of the p-values. Each support is represented by a vector that must be in increasing order.","") %>
#' <%=ifelse(exists("pvalues") && pvalues, "@param pvalues a numeric vector. Contains all values of the p-values supports if we search for the critical constants. If not, contains only the observed p-values. Must be sorted in increasing order!","") %>
#' <%=ifelse(exists("stepUp") && stepUp, "@param stepUp a numeric vector. Identical to \\code{pvalues} for a step-down procedure. Equals \\code{c.m} for a step-up procedure.","") %>
#' <%=ifelse(exists("support") && support, "@param support a numeric vector. Contains all values of the p-values supports. Ignored, if \\code{stepUp = FALSE}. Must be sorted in increasing order!","") %>
#' <%=ifelse(exists("alpha") && alpha, "@param alpha the target FDP, a number strictly between 0 and 1. For \\code{*.fast} kernels, it is only necessary, if \\code{stepUp = TRUE}.","") %>
#' <%=ifelse(exists("zeta") && zeta, "@param zeta the target probability of not exceeding the desired FDP, a number strictly between 0 and 1. If \\code{zeta=NULL} (the default), then \\code{zeta} is chosen equal to \\code{alpha}.","") %>
#' <%=ifelse(exists("sorted_pv") && sorted_pv, "@param sorted_pv a vector of observed p-values, in increasing order.","") %>
#' <%=ifelse(exists("raw.pvalues") && raw.pvalues, "@param raw.pvalues vector of the raw observed p-values, as provided by the end user and before matching with their nearest neighbor in the CDFs supports.","") %>
#' <%=ifelse(exists("direction") && direction, "@param direction a character string specifying whether to conduct a step-up (\\code{direction=\"su\"}, the default) or step-down procedure (\\code{direction=\"sd\"}).","") %>
#' <%=ifelse(exists("critical.values") && critical.values, "@param critical.values a boolean. If \\code{TRUE}, critical constants are computed and returned (this is computationally intensive).","") %>
#' <%=ifelse(exists("adaptive") && adaptive, "@param adaptive a boolean specifying whether to conduct an adaptive procedure or not.","") %>
#' <%=ifelse(exists("exact") && exact, "@param exact a boolean specifying whether to compute the Poisson-Binomial distribution exactly or by a normal approximation.","") %>
#' <%=ifelse(exists("weights") && weights, "@param weights a numeric vector. Contains the weights of the p-values.","") %>
#' <%=ifelse(exists("weighting.method") && weighting.method, "@param weighting.method a character string specifying whether to conduct arithmetic (\\code{direction=\"AM\"}, the default) or geometric weighting (\\code{direction=\"GM\"}) of p-values.","") %>
