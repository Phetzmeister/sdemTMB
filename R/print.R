
#' Summary of sdem class
#' 
#' @description 'summary.sdem' prints parameter estimates, standard deviations,
#' the slope of the objective function in parameter directions etc.
#' @param object an object of class 'sdem'.
#' @param extended logical. if 'TRUE' the parameter correlation matrix is also printed.
#' @param ... further arguments passed to/from other methods.
#' 
#' @return a lot of information from the optimisation. parameter estimates,
#' state estimates etc. 
#' 

summary.sdem = function(object,extended=FALSE,...) {
  object$summary(extended,...)
}
