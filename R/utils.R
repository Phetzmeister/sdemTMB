#######################################################
# SIMPLIFY FORMULA
#######################################################

simplify_formula = function(form) {

  form = as.formula(paste(
    form[[2]],
    # paste(deparse(Deriv::Simplify(form[[2]])),collapse=""),
    paste(deparse(Deriv::Simplify(form[[3]])),collapse=""),
    sep="~"
  ))

}

#######################################################
# EXTRA FUNCTIONS FOR DERIV THAT ALSO WORKS IN TMB
#######################################################

# this only works because logit and invlogit is already defined in TMB
# if not they must be defined in the c++ file

logit = function(x) log(x/(1-x))
invlogit = function(x) 1/(1+exp(-x))

#######################################################
# CHANGE FROM R POWER NOTATION TO C++
#######################################################

# Changes powers in expressions from R-style (e.g. x^2) to pow(x,2) which is interpretable by C++
# from https://stackoverflow.com/questions/40606723/substitute-the-power-symbol-with-cs-pow-syntax-in-mathematical-expression

hat2pow <- function(e) {
  #check if you are at the end of the tree's branch
  if (is.name(e) || is.atomic(e)) {
    #replace ^
    if (e == quote(`^`)) return(quote(pow))
    return(e)
  }
  #follow the tree with recursion
  for (i in seq_along(e)) e[[i]] <- hat2pow(e[[i]])
  return(e)
}

#######################################################
# Print - S3 Method
#######################################################


#' Basic print of objects of class 'sdem'
#' @returns A huge amount of information
#' @export
print.sdem = function(object,...) {

  print("I'm using the s3 method")

  # For calling summary on model object
  if (inherits(object,"sdem") & inherits(object,"R6")) {


    obj = object$print()
    return(invisible(obj))

    # for calling summary on object returned by model$estimate(data) called 'fit'
  }
  if (inherits(object,"sdem")) {

    fit = object
    mat = cbind(fit$par.fixed,fit$sd.fixed,fit$tvalue,fit$Pr.tvalue)
    colnames(mat) = c("Estimate","Std. Error","t value","Pr(>|t|)")
    cat("Coefficent Matrix \n")
    stats::printCoefmat(mat)

    returnlist = list(parameters=mat)
    class(returnlist) = "sdem"

    return(invisible(returnlist))

  }

  return(invisible(NULL))
}

#######################################################
# Summary - S3 Method
#######################################################

#' Basic summary of objects of class 'sdem'
#' @returns A huge amount of information
#' @export
summary.sdem = function(object,...) {

  print("I'm using the s3 method")

  # For calling summary on model object
  if (inherits(object,"sdem") & inherits(object,"R6")) {

    obj = object$summary()
    return(invisible(obj))

    # for calling summary on object returned by model$estimate(data) called 'fit'
  }
  if (inherits(object,"sdem")) {

    fit = object
    mat = cbind(fit$par.fixed,fit$sd.fixed,fit$tvalue,fit$Pr.tvalue)
    colnames(mat) = c("Estimate","Std. Error","t value","Pr(>|t|)")
    cat("Coefficent Matrix \n")
    stats::printCoefmat(mat)

    returnlist = list(parameters=mat)
    class(returnlist) = "sdem"

    return(invisible(returnlist))

  }

}

#######################################################
# Plot - S3 Method
#######################################################


#' Basic summary of objects of class 'sdem'
#' @param extended logical. if TRUE additional information is printed
#' @returns A huge amount of information
#' @export
plot.sdem = function(object,extended=FALSE,...) {
  object$plot(extended,...)
}


