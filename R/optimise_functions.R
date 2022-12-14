
#######################################################
# CHECK AND SET DATA BEFORE OPTIMIZATION
#######################################################

check_and_set_data = function(data, self, private) {

  # is data a list of data.frame
  if (!is.list(data) & !is.data.frame(data)) {
    stop("The data should be a list or a data.frame")
  }
  # if data.frame convert to list
  if (is.data.frame(data)){
    nams = names(data)
    datalist = list()
    for (i in 1:length(data)){
      datalist[[nams[i]]] = data[[nams[i]]]
    }
    data = datalist
  }

  # check that all required entries are provided
  required.names = c(private$input.names,private$obs.names)
  given.names = names(data)
  bool = required.names %in% given.names
  if (any(!bool)){
    stop("The following required entries were not supplied in the data: \n\t",
         paste(required.names[!bool],collapse=", "))
  }
  bool = given.names %in% required.names
  if (any(!bool)){
    # stop("The following provided data entries are not used by the model: \n\t",
    # paste(given.names[!bool],collapse=", "))
    warning("The following provided data entries are redundant since they are not used in the model: \n\t",
            paste(given.names[!bool],collapse=", "))
  }

  # check time entry
  # must be increasing
  if (any(diff(data$t)<=0)) {
    ids = which(diff(data$t)<=0)
    stop(sprintf("The time-vector is non-increasing at the following indice(s) %s",paste(ids,collapse=", ")))
  }
  # must not contain non-numeric values
  if (any(is.na(data$t))) {
    ids = which(is.na(data$t))
    stop(sprintf("The time-vector is NA at the following indice(s) %s",paste(ids,collapse=", ")))
  }

  # add data to private
  private$data = data

  # time-step must be smaller than smallest time diff

  # set ODE time-steps
  if (is.null(private$ode.N)) {
    private$ode.N = rep(1,length(data$t)-1)
  } else if (length(private$ode.N)==1) {
    private$ode.N = rep(private$ode.N,length(data$t)-1)
  } else {
    if (!(length(private$ode.N) == length(private$data$t)-1)) {
      stop("You must provide time steps factor of length 1 or length(data$t)-1")
    }
  }
  private$ode.dt = diff(data$t)/private$ode.N
  private$ode.Ncumsum = c(0,cumsum(private$ode.N))

  #
  return(invisible(self))
}

#######################################################
# WHEN WE USE TMB THE STATES SHOULD BE ADDED
#######################################################

if_method_is_tmb_add_states_to_parameters = function(self, private) {

  if (!any(private$method == c("tmb","tmb_exact"))) {
    return(invisible(self))
  }

  if (!is.null(private$tmb.initial.state)) {
    # repeat each state guess N times (except the last state value)
    for (i in 1:length(private$tmb.initial.state)) {
      private$tmb.initial.state[[i]] = rep(private$tmb.initial.state[[i]],times=c(private$ode.N,1))
    }
    return(invisible(self))
  }

  message("No initial state values provided for TMB. I will use the provided
          initial state values for all time points.")
  # The number of states is sum(private$ode.N) + 1
  for (i in 1:private$n) {
    private$tmb.initial.state[[i]] = rep(private$initial.state$mean[i],sum(private$ode.N)+1)
  }
  names(private$tmb.initial.state) = private$state.names

  return(invisible(self))
}

#######################################################
# CONSTRUCT AD-FUN AND RUN OPTIMIZATION
#######################################################

construct_and_optimise = function(self ,private) {

  # Find free parameters
  bool = private$parameter.names %in% names(private$fixed.pars)
  private$free.pars = private$parameters[!bool]

  ################################################
  # Data
  ################################################

  # add mandatory entries to data
  .data1 = list(
    X0__ = private$initial.state$mean,
    P0__ = private$initial.state$cov,
    dt__ = private$ode.dt,
    N__ = private$ode.N,
    Ncumsum__ = private$ode.Ncumsum,
    which_loss__ = private$loss$loss,
    loss_c_value__ = private$loss$c,
    tukey_pars__ = private$tukey.pars
  )
  # add MAP entries to data if MAP is active
  if (!is.null(private$map)) {
    .data2 = list(
      map_bool__ = 1,
      map_mean__ = private$map$mean,
      map_cov__ = private$map$cov,
      map_ints__ = as.numeric(!bool),
      sum_map_ints__ = sum(as.numeric(!bool))
    )
    private$map$mean = private$map$mean[!bool]
    private$map$cov = private$map$cov[!bool,!bool]
  } else {
    .data2 = list(
      map_bool__ = 0
    )
  }
  # add constants to data
  .data3 = lapply(private$constants, function(x) x$rhs)

  # construct final data list
  data = c( private$data , .data1, .data2, .data3)

  ################################################
  # Parameters
  ################################################

  parameters = c(private$tmb.initial.state, lapply(private$parameters, function(x) x$init))

  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  .random = private$state.names
  if (any(private$method==c("ekf","ukf"))) {
    .random = NULL
  }

  message("Constructing function template...")

  nll = tryCatch(
    TMB::MakeADFun(data = data,
                   parameters = parameters,
                   map = private$fixed.pars,
                   DLL = private$modelname,
                   silent = private$silent,
                   random = .random)
    ,
    error = function(e) {
      invisible(structure(e,class="try-error"))
    },
    warning = function(w) {
      invisible(structure(w,class="try-warning"))
    }
  )
  if (inherits(nll,c("try-error","try-warning"))) {
    stop("Pre-compiled model parameters and the supplied do not match. You need to recompile your model. Use $set_compile(TRUE)")
  }
  private$nll = nll

  ################################################
  # Optimise
  ################################################

  lb = unlist(lapply(private$free.pars, function(x) x$lb))
  ub = unlist(lapply(private$free.pars, function(x) x$ub))

  # Need to do something with tryCatch errors here
  # maybe add private$didOptimConverge for flag to use in estimate
  # function instead of !is.null(private$opt)

  message("Estimation initiated - minimizing the objective function...")
  if (any(private$method==c("ekf","ukf"))) {
    if (private$use.hessian) {
      comptime = system.time(
        opt <- stats::nlminb(start = private$nll$par,
                             objective = private$nll$fn,
                             gradient = private$nll$gr,
                             he = private$nll$he,
                             lower = lb,
                             upper = ub)
      )
    } else {
      comptime = system.time(
        opt <- stats::nlminb(start = private$nll$par,
                             objective = private$nll$fn,
                             gradient = private$nll$gr,
                             lower = lb,
                             upper = ub)
      )
    }
  }
  if (private$method =="tmb") {
    comptime = system.time(
      opt <- stats::nlminb(start = private$nll$par,
                           objective = private$nll$fn,
                           gradient = private$nll$gr,
                           lower = lb,
                           upper = ub)
    )
  }

  if (inherits(opt,"error")) {
    message("The optimisation failed due to an error: \n\n\t",opt)
  }
  if (inherits(opt,"warning")) {
    message("The optimisation failed due to a warning: \n\n\t",opt)
  }
  if (is.null(opt$convergence)) {
    private$opt = NULL
    return(invisible(self))
  } else {
    private$opt = opt
  }

  # Success / Failure message
  # comptime = round(as.numeric(comptime["elapsed"])*1e2)/1e2
  # if(private$opt$convergence == 1){
  #   message("Failure! The estimation encountered numerical issues: elapsed time ", comptime, " seconds.")
  # } else {
  #   message("Success! The estimation was completed: elapsed time: ", comptime, " seconds.")
  # }

  if (private$method=="tmb") {
    private$sdr = TMB::sdreport(private$nll)
  }

  return(invisible(self))
}

#######################################################
# MAKE RETURN DATA NICE AFTER OPTIMIZATION
#######################################################

create_return_fit = function(self, private) {

  if (is.null(private$opt)) {
    return(NULL)
  }

  private$fit$convergence = private$opt$convergence

  ################################################
  # FOR KALMAN FILTERS
  ################################################

  if (private$method == "ekf" | private$method == "ukf") {


    # Objective and Gradient
    private$fit$nll.value = private$opt$objective
    private$fit$nll.gradient = tryCatch( as.vector(private$nll$gr(private$opt$par)),
                                         error=function(e) NA,
                                         warning=function(w) NA
    )
    names(private$fit$nll.gradient) = names(private$free.pars)

    # Hessian
    private$fit$nll.hessian = tryCatch(private$nll$he(private$opt$par),
                                       error=function(e) NA,
                                       warning=function(w) NA
    )
    rownames(private$fit$nll.hessian) = names(private$free.pars)
    colnames(private$fit$nll.hessian) = names(private$free.pars)


    # Parameter Estimate
    private$fit$par.fixed = private$opt$par
    private$fit$sd.fixed = tryCatch(sqrt(diag(solve(private$nll$he(private$opt$par)))),
                                    error=function(e) NA,
                                    warning=function(w) NA
    )

    # Parameter Covariance
    private$fit$cov.fixed = tryCatch(solve(private$nll$he(private$opt$par)),
                                     error=function(e) NA,
                                     warning=function(w) NA
    )

    # Prior States
    rep = private$nll$report()
    temp.states = c()
    temp.sd = c()
    for (i in 1:length(private$data$t)) {
      temp.states = rbind(temp.states, rep$xPrior[[i]])
      temp.sd = rbind(temp.sd, diag(rep$pPrior[[i]]))
    }
    temp.states = cbind(private$data$t, temp.states)
    temp.sd = cbind(private$data$t, temp.sd)
    colnames(temp.states) = c("t",private$state.names)
    colnames(temp.sd) = c("t",private$state.names)
    private$fit$states$mean$prior = as.data.frame(temp.states)
    private$fit$states$sd$prior = as.data.frame(temp.sd)

    # Posterior States
    temp.states = c()
    temp.sd = c()
    for (i in 1:length(private$data$t)) {
      temp.states = rbind(temp.states, rep$xPost[[i]])
      temp.sd = rbind(temp.sd, diag(rep$pPost[[i]]))
    }
    temp.states = cbind(private$data$t, temp.states)
    temp.sd = cbind(private$data$t, temp.sd)
    colnames(temp.states) = c("t",private$state.names)
    colnames(temp.sd) = c("t",private$state.names)
    private$fit$states$mean$posterior = as.data.frame(temp.states)
    private$fit$states$sd$posterior = as.data.frame(temp.sd)

    # Full Covariance Matrix - Prior and Posterior States
    private$fit$states$cov$prior = rep$pPost
    private$fit$states$cov$posterior = rep$pPrior
    names(private$fit$states$cov$prior) = paste("t=",private$data$t,sep="")
    names(private$fit$states$cov$posterior) = paste("t=",private$data$t,sep="")

    # Residuals
    rowNAs = as.matrix(!is.na(do.call(cbind,private$data[private$obs.names]))[-1,])
    sumrowNAs = rowSums(rowNAs)

    innovation = rep$Innovation
    innovation[[1]] = NULL
    innovation.cov = rep$InnovationCovariance
    innovation.cov[[1]] = NULL

    temp.res = matrix(nrow=length(private$data$t)-1,ncol=private$m)
    temp.sd =  matrix(nrow=length(private$data$t)-1,ncol=private$m)
    for (i in 1:(length(private$data$t)-1)) {
      if (sumrowNAs[i] > 0) {
        temp.res[i,rowNAs[i,]] = innovation[[i]]
        temp.sd[i,rowNAs[i,]] = diag(innovation.cov[[i]])
      }
    }
    temp.res = cbind(private$data$t[-1],temp.res)
    temp.sd = cbind(private$data$t[-1],temp.sd)

    names(innovation.cov) = paste("t=",private$data$t[-1],sep="")

    # should we remove the empty matrices?
    # innovation.cov = innovation.cov[sumrowNAs!=0]

    colnames(temp.res) = c("t",private$obs.names)
    colnames(temp.sd) = c("t",private$obs.names)
    private$fit$residuals$mean = as.data.frame(temp.res)
    private$fit$residuals$sd = as.data.frame(temp.sd)
    private$fit$residuals$normalized = as.data.frame(temp.res)
    private$fit$residuals$normalized[,-1] = private$fit$residuals$normalized[,-1]/temp.sd[,-1]
    private$fit$residuals$cov = innovation.cov


    # Observations?
    # private$fit$observations$mean =

    # t-values and Pr( t > t_test )
    private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
    private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)

  }

  ################################################
  # FOR TMB
  ################################################

  if (private$method == "tmb") {

    # Objective and Gradient
    private$fit$nll.value = private$opt$objective
    private$fit$nll.gradient = private$sdr$gradient.fixed
    names(private$fit$nll.gradient) = names(private$free.pars)

    # Parameter Estimate
    private$fit$par.fixed = private$opt$par
    private$fit$sd.fixed = diag(private$sdr$cov.fixed)

    # Parameter Covariance
    private$fit$cov.fixed = private$sdr$cov.fixed

    # Posterior States (Smoothed)
    # Only pick out correct number of states
    states.at.timepoints = cumsum(c(1,private$ode.N))
    temp = cbind(private$sdr$par.random,sqrt(private$sdr$diag.cov.random))[states.at.timepoints,]
    temp.states = cbind(private$data$t, matrix(temp[,1],nrow=length(private$data$t)))
    temp.sd = cbind(private$data$t, matrix(temp[,2],nrow=length(private$data$t)))
    #
    private$fit$states$mean$posterior = as.data.frame(temp.states)
    private$fit$states$sd$posterior = as.data.frame(temp.sd)
    colnames(private$fit$states$sd$posterior) = c("t",private$state.names)
    colnames(private$fit$states$mean$posterior) = c("t",private$state.names)

    # Residuals
    rowNAs = as.matrix(!is.na(do.call(cbind,private$data[private$obs.names]))[-1,])
    sumrowNAs = rowSums(rowNAs)

    # still need to figure out how to get residuals
    # TMB::oneStepPredict(private$nll,
    #                     observation.name=private$obs.names,
    #                     method="fullGaussian",
    #                     discrete=FALSE)

    # t-values and Pr( t > t_test )
    private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
    private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)

  }

  class(private$fit) = "sdem"

  return(invisible(self))
}

