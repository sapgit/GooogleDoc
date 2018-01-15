b.mean <- function(vec, num, na.rm=T) {
  if(na.rm) vec<-vec[!is.na(vec)]
  vec<-rm.outlier(vec)
  resamples <- lapply(1:num, function(i) sample(vec, replace=T))
  r.mean <- sapply(resamples, mean)
  std.err <- sqrt(var(r.mean))
  return(round(std.err,2))
}


#' Title
#'
#' @description For a given training and test dataset and the fitted coefficients this function calculates MAE and MASE
#' @param train Training dataset obtained from the original dataset
#' @param test Complement of the training dataset used for prediction
#' @param fit The output of the function fit.method
#' @param yvar Name of the outcome variable
#' @param xvars Name of the predictor variables for the count model
#' @param zvars Name of the predictor variables for the zero model
#'
#' @return The predictive measures MAE and MASE calculated from the function accuracy of the forecast package

measures.func<-function(train,test,fit,yvar,xvars,zvars)
{
  if (is.na(fit))
  {
    measures<-c(NA,NA)
  } else {
    betahat<-fit$coefficients$count
    gammahat<-fit$coefficients$zero
    z.test<-as.matrix(cbind(1,test[,zvars]))
    phi.hat<-1/(1+exp(-z.test%*%gammahat)) #
    x.test<-as.matrix(cbind(1,test[,xvars]))
    lam.hat<-exp(x.test%*%betahat)
    y.pred<-(1-phi.hat)*lam.hat
    y.test<-test[,yvar]
    y.train<-train[,yvar]
    forecast <- structure(list(mean=y.pred, fitted=y.test, x=y.train), class='forecast')
    measures<-c(round(accuracy(forecast,y.test)[2,c(3,6)],4))
  }
  return(measures)
}

#' Title
#'
#' @description This function fits the ZINB model with different penalties to the training part of a dataset and calculates MAE and MASE from the test set. It outputs the median of MAE and MASE over all the simulated datasets
#' @param n.train Sample size in the training dataset
#' @param n.test Sample size in the test dataset
#' @param data.list output of datagen.sim.all
#' @param method Different penalties
#' @param ITER Number of simulations
#' @param group Vector containing grouping structure of the covariates
#'
#' @return The median of MAE and MASE, calculated over all the simulated datasets

measures.summary<-function(n.train,n.test,data.list,method,ITER,group)
{
  options(warn=-1)
  measures.mat<-NULL
  for (i in 1:ITER)
  {
    dataset<-data.list[[i]]
    train<-dataset$data[1:n.train,]
    test<-dataset$data[-(1:n.train),]

    yvar<-dataset$yvar
    xvars<-dataset$xvars
    zvars<-dataset$zvars

    fit.summary<-fit.method(dataset=train,yvar=yvar,xvars=xvars,zvars=zvars,method=method,group=group)
    fit<-fit.summary$fit
    time.taken<-fit.summary$time

    predict.measures<- measures.func(train=train,test=test,fit=fit,yvar=yvar,xvars=xvars,zvars=zvars)
    measures.mat<-rbind(measures.mat,c(predict.measures,time.taken))
  }
  
  measures.se<-t(apply(apply(measures.mat[,-3],2,function(x) return(as.numeric(x))),2,b.mean,num=1000,na.rm=T))
  
  measures.median<-apply(measures.mat,2, function(x){median(x, na.rm=T)})
  
  mae.summary<-paste(round(measures.median[1],3),"(",round(measures.se[1],3),")",sep="")

  mase.summary<-paste(round(measures.median[2],3),"(",round(measures.se[2],3),")",sep="")
  
  output<-c(MAE=mae.summary,MASE=mase.summary,time.taken=measures.median[3])
  
  options(warn=0)
  return(output)
}
