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

  obj = object$print()
  #
  return(invisible(obj))

}

#' Basic print of objects of class 'sdem'
#' @returns A huge amount of information
#' @export
print.sdem.fit = function(fit) {

  mat = cbind(fit$par.fixed,fit$sd.fixed,fit$tvalue,fit$Pr.tvalue)
  colnames(mat) = c("Estimate","Std. Error","t value","Pr(>|t|)")
  cat("Coefficent Matrix \n")
  stats::printCoefmat(mat)

  return(invisible(mat))

}


#######################################################
# Summary - S3 Method
#######################################################

#' Basic summary of objects of class 'sdem'
#' @returns A huge amount of information
#' @export
summary.sdem = function(object,correlation=FALSE) {

  obj = object$summary(correlation)
  #
  return(invisible(obj))

}

#' @returns summary of fit object from \code{obj$estimate()}
#' @export
summary.sdem.fit = function(fit,correlation=FALSE) {

  if (!is.logical(correlation)) {
    stop("correlation must be logical")
  }

  mat = cbind(fit$par.fixed,fit$sd.fixed,fit$tvalue,fit$Pr.tvalue)
  colnames(mat) = c("Estimate","Std. Error","t value","Pr(>|t|)")
  cat("Coefficent Matrix \n")
  stats::printCoefmat(mat)
  if (correlation){
    cat("\nCorrelation Matrix\n")
    S = diag(1/private$fit$sd.fixed)
    cor = S %*% private$fit$cov.fixed %*% S
    rownames(cor) = names(fit$par.fixed)
    colnames(cor) = names(fit$par.fixed)
    cor <- format(round(cor, 2), nsmall = 2, digits = 6)
    cor[!lower.tri(cor,diag=T)] <- ""
    print(cor, drop = FALSE, quote = FALSE)
    # cor[!lower.tri(cor)] <- ""
    # print(cor[-1, -(dim(cor)[1]), drop = FALSE], quote = FALSE)
  }

  return(invisible(list(parameters=mat)))
}

#######################################################
# GGPLOT2 FUNCTIONS FOR USE IN PLOTS
#######################################################

getggplot2theme = function() {
  mytheme =
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = element_text("Avenir Next Condensed",size=12),
      legend.text = element_text(size=12),
      axis.text = element_text(size=12),
      strip.text = element_text(face="bold",size=12),
      # panel.grid.major = element_blank(),
      # panel.grid.minor = element_blank(),
      legend.box = "vertical",
      legend.position = "top",
      plot.title = element_text(hjust=0.5)
    )
  return(mytheme)
}

getggplot2colors = function(n) {
  hues = seq(15, 375, length = n + 1)
  ggcolors = hcl(h = hues, l = 65, c = 100)[1:n]
  return(ggcolors)
}

#######################################################
# Plot - S3 Method
#######################################################

#' Basic summary of objects of class 'sdem'
#' @param extended logical. if TRUE additional information is printed
#' @returns A huge amount of information
#' @export
plot.sdem = function(object,plot.obs=1,extended=FALSE,use.ggplot=FALSE,ggtheme=getggplot2theme()) {

  object$plot(plot.obs=plot.obs,use.ggplot=use.ggplot,extended=extended,ggtheme=ggtheme)
  return(invisible(NULL))
}

#' Basic summary of objects of class 'sdem'
#' @param extended logical. if TRUE additional information is printed
#' @returns A huge amount of information
#' @export
plot.sdem.fit = function(fit, plot.obs=1,use.ggplot=FALSE,extended=FALSE,ggtheme=getggplot2theme()){

  if (!is.logical(use.ggplot)) {
    stop("use.ggplot must be logical")
  }
  if (!is.logical(extended)) {
    stop("extended must be logical")
  }
  if (!(inherits(ggtheme,"theme") & inherits(ggtheme,"gg"))) {
    stop("ggtheme must be a ggplot2 theme")
  }
  if (!is.numeric(plot.obs)) {
    stop("plot.obs must be an integer")
  }

  l = length(private$fit$residuals$normalized)-1
  nams = names(private$fit$residuals$normalized)[-1]
  mycolor = getggplot2colors(2)[2]
  plot.obs = min(l,plot.obs)

  if (use.ggplot) {
    plots = list()
    for (i in 1:l) {
      t = private$fit$residuals$normalized[["t"]]
      e = private$fit$residuals$normalized[[i+1]]
      e = na.omit(e)
      id = attr(e,"na.action")
      t = t[-id]
      nam = nams[i]

      # time vs residuals
      plot.res =
        ggplot2::ggplot(data=data.frame(t,e)) +
        ggplot2::geom_point(aes(x=t,y=e),color=mycolor) +
        ggtheme +
        ggplot2::labs(
          title = paste("Time Series of Residuals: "),
          y = "",
          x = "Time"
        )
      # quantile plots
      plot.qq =
        ggplot2::ggplot(data=data.frame(e)) +
        ggplot2::stat_qq(aes(sample=e),color=mycolor) +
        ggplot2::stat_qq_line(aes(sample=e),lty="dashed") +
        ggtheme +
        ggplot2::labs(
          title = paste("Normal Q-Q Plot: "),
          y = "Sample Quantiles",
          x = "Theoretical Quantiles"
        )
      # histogram
      plot.hist =
        ggplot2::ggplot(data=data.frame(e)) +
        ggplot2::geom_histogram(aes(x=e,y=..density..),bins=100,color="black",fill=mycolor) +
        ggtheme +
        ggplot2::labs(
          title = paste("Histogram: "),
          y = "",
          x = ""
        )
      # acf
      myacf = stats::acf(e,na.action=na.pass,plot=FALSE)
      plot.acf =
        ggplot2::ggplot(data=data.frame(lag=myacf$lag[-1],acf=myacf$acf[-1])) +
        ggplot2::geom_errorbar(aes(x=lag,ymax=acf,ymin=0),width=0,color=mycolor) +
        ggplot2::geom_hline(yintercept=0) +
        ggplot2::geom_hline(yintercept=c(-2/sqrt(myacf$n.used),2/sqrt(myacf$n.used)),color="blue",lty="dashed") +
        ggtheme +
        ggplot2::labs(
          title = paste("Auto-Correlation: "),
          y = "",
          x = "Lag"
        )
      # cpgram
      plot.cpgram =
        ggfortify::ggcpgram(e,colour=mycolor) +
        ggtheme +
        ggplot2::labs(
          title = paste("Cumulative Periodogram: "),
          y = "",
          x = "Lag"
        )
      #
      templist = list(
        plot.res,
        plot.hist,
        plot.qq,
        plot.acf,
        plot.cpgram
      )
      # save plot
      plots[[i]] = (plot.qq + plot.hist) / (plot.acf + plot.cpgram) +
        plot_annotation(title=paste("Residuals for ", nam),theme=ggtheme + ggplot2::theme(text=element_text("Avenir Next Condensed",size=18,face="bold")),)
    }
    # return plot list, and print one of the plots
    print(plots[[plot.obs]])
    return(invisible(plots))
  }

  # if we aren't using ggplot we can't return the plots, we just print one plot
  graphics::par(mfrow=c(2,2))
  e = private$fit$residuals$normalized[[plot.obs+1]]
  nam = nams[plot.obs]
  # qqnorm and line
  stats::qqnorm(e, main=paste("Normal Q-Q Plot:"),axes=FALSE)
  graphics::grid()
  graphics::axis(1,-3:3)
  graphics::axis(2)
  graphics::grid()
  stats::qqline(e)
  # histogram
  graphics::hist(e,main=paste("Histogram:"),breaks=50)
  graphics::grid()
  graphics::axis(1)
  graphics::axis(2)
  # acf
  max.lag = round(10*log10(length(stats::na.omit(e))))
  plot(acf(e,plot=F,na.action=na.pass)[1:max.lag],main=paste("Auto-Correlation:"),axes=FALSE)
  graphics::grid()
  graphics::axis(1,seq(0,max.lag,round(max.lag/10)))
  graphics::axis(2)
  # cpgram
  x <- e[!is.na(e)]
  x <- spec.taper(scale(x, TRUE, FALSE), p = 0.1)
  y <- Mod(fft(x))^2/length(x)
  y[1L] <- 0
  n <- length(x)
  x <- (0:(n/2)) * frequency(e)/n
  if (length(x)%%2 == 0) {
    n <- length(x) - 1
    y <- y[1L:n]
    x <- x[1L:n]
  } else {
    y <- y[seq_along(x)]
  }
  xm <- frequency(e)/2
  mp <- length(x) - 1
  crit <- 1.358/(sqrt(mp) + 0.12 + 0.11/sqrt(mp))
  plot(x, cumsum(y)/sum(y), type = "s", xlim = c(0, xm), ylim = c(0,1),
       xaxs = "i", yaxs = "i", xlab = "frequency", ylab = "",
       main=paste("Cumulative Periodogram:"), axes=FALSE)
  lines(c(0, xm * (1 - crit)), c(crit, 1), col = "blue", lty = 2)
  lines(c(xm * crit, xm), c(0, 1 - crit), col = "blue", lty = 2)
  # stats::cpgram(e,main=paste("Cumulative Periodogram:"),axes=FALSE)
  graphics::grid()
  graphics::axis(1)
  graphics::axis(2)
  # overall title
  graphics::title(main=paste("Residuals for ", nam),outer=TRUE,line=-1.25,cex.main=1.5)

  # return
  return(invisible(NULL))
}
