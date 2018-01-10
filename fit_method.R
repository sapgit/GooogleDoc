#' Title
#'
#' @description Fits the ZINB model to the dataset with grouped covariates using gooogle function or zipath for "EMLasso"
#' @param dataset Simulation dataset
#' @param yvar Name of the outcome variable
#' @param xvars Name of the predictor variables of the count model
#' @param zvars Name of the predictor variables of the zero model
#' @param method Different penalties
#' @param dist Distribution of the count model which is Negative Binomial in this case
#' @param group Vector containing grouping structure of the covariates
#'
#' @return  fitted coefficients for count model and zero model along with the AIC, BIC and loglikelihood of the ZINB model

fit.method<-function(dataset,yvar,xvars,zvars,method,dist,group)
{
  fit.formula<-as.formula(paste(yvar,"~",paste(paste(xvars,collapse="+"),"|",paste(zvars,collapse="+"),sep="")))
  if (method=="EMLasso") {
    ptm<-proc.time()
    fit<-try(eval(call(paste("func_",method,sep=""),formula=fit.formula,data=dataset,dist=dist)),silent=T)
  } else {
    ptm<-proc.time()
    fit<-try(eval(call(paste("func_",method,sep=""),data=dataset,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist=dist)),silent=T)
  }

  # if there is an error in calling func_methods
  if (class(fit)=="try-error")
  {
    fit<-NA
  }
  return(fit)
}
