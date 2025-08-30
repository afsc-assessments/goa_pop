# note that the 2nd time period in the selectivity blocks is not scaled to one
# so have this add on to account for that 
to_one <- function(x) {
  x / max(x)
}


#' standard MESA logistic selectivity
#'
#' @param age age to examine
#' @param a50 inflection point
#' @param delta steepness of slope
#' @param adj based upon index that starts at 1 - unless adjusted to start at recruitment age
#' doesn't really do anything other than scale
sel_logistic <- function(age, a50, delta, adj) {
  sel = 1 / (1 + exp(-log(19) * ((age + adj) - a50) / delta))
  sel #= sel / max(sel)
}

#' standard MESA double logistic selectivity
#'
#' @param age age to examine
#' @param a50 inflection point
#' @param delta steepness of slope
#' @param adj based upon index that starts at 1 - unless adjusted to start at recruitment age
#' doesn't really do anything other than scale
sel_double_logistic <- function(age, a50, delta, adj) {
  expa50 = exp(a50)  # log scale for dome
  denom = 0.5 * (sqrt(expa50^2 + 4 * delta^2) - expa50)
  sel = (((age + adj) / expa50)^(expa50 / denom)) * exp((expa50 - (age + adj)) / denom)
  sel #= sel / max(sel)
}

sel_double_normal <- function(age, a50, delta, adj = 0) {
  sel <- exp(-((age + adj - a50)^2) / (2 * delta^2))
  sel
}

# Helper function for Dirichlet-multinomial likelihood
dirmultinom <- function(x, prob, theta, log = FALSE) {
  # x: observed proportions
  # prob: predicted proportions  
  # theta: effective sample size
  n = sum(x)
  alpha = prob * theta
  ll = lgamma(n + 1) + lgamma(theta) - lgamma(n + theta) +
    sum(lgamma(x + alpha) - lgamma(alpha) - lgamma(x + 1))
  if(log) return(ll) else return(exp(ll))
}

filt <- function(item, data) {
  newd = data
  newd[[item]] = 0
  newd
}
