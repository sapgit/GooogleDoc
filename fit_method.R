#' Title
#'
#' @description Fits the ZINB model to the dataset with grouped covariates using gooogle function or zipath for "EMLasso"
#' @param dataset Simulation dataset
#' @param yvar Name of the outcome variable
#' @param xvars Name of the predictor variables of the count model
#' @param zvars Name of the predictor variables of the zero model
#' @param method Different penalties
#' @param group Vector containing grouping structure of the covariates
#'
#' @return  fitted coefficients for count model and zero model along with the AIC, BIC and loglikelihood of the ZINB model

fit.method<-function(dataset,yvar,xvars,zvars,method,group)
{
  fit.formula<-as.formula(paste(yvar,"~",paste(paste(xvars,collapse="+"),"|",paste(zvars,collapse="+"),sep="")))
  if (method=="EMLasso") {
    ptm<-proc.time()
    fit<-try(eval(call(paste("func_",method,sep=""),formula=fit.formula,data=dataset,dist="negbin")),silent=T)
    time.taken<-round((proc.time()-ptm)[3],3)
  } else {
    ptm<-proc.time()
    fit<-try(eval(call(paste("func_",method,sep=""),data=dataset,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist="negbin")),silent=T)
    time.taken<-round((proc.time()-ptm)[3],3)
  }

  # if there is an error in calling func_methods
  if (class(fit)=="try-error")
  {
    fit<-NA
  }
  return(list(fit=fit,time=time.taken))
}
