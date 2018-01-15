# Required Libraries
require(forecast)
require(cvTools)
require(Gooogle)
require(mpath)

#' Title
#'
#' @param train Training data
#' @param test Test data
#' @param yvar The name of the outcome variable
#' @param xvars The name of the predictor variables for the count model
#' @param zvars The name of the predictor variables for the zero model
#' @param group The vector specifying grouping structure of the covariates
#' @param penalty Penalty
#'
#' @return Predicted value of the outcome variable and outcome variable from the test and the training set

fit.out<-function(train,test,yvar,xvars,zvars,group,penalty)
{
  if(penalty=="EMLasso")
  {
    fit.formula<-as.formula(paste(yvar,"~",paste(paste(xvars,collapse="+"),"|",paste(zvars,collapse="+")),sep=""))

    fit.em<-zipath(formula=fit.formula,data=train,family="negbin")

    bic.idx<-which.min(fit.em$bic)
    gammahat<-fit.em$coefficients$zero[,bic.idx]
    betahat<-fit.em$coefficients$count[,bic.idx]
    z.test<-as.matrix(cbind(1,test[,zvars]))
    phi.hat<-1/(1+exp(-z.test%*%gammahat))
    x.test<-as.matrix(cbind(1,test[,xvars]))
    lam.hat<-exp(x.test%*%betahat)
    y.pred<-(1-phi.hat)*lam.hat
    y.test<-test[,yvar]
    y.train<-train[,yvar]
  } else {
    fit<-gooogle(data=train,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist="negbin",penalty=penalty)

  z.test<-as.matrix(cbind(1,test[,zvars]))
  phi.hat<-1/(1+exp(-z.test%*%fit$coefficients$zero))
  x.test<-as.matrix(cbind(1,test[,xvars]))
  lam.hat<-exp(x.test%*%fit$coefficients$count)
  y.pred<-(1-phi.hat)*lam.hat
  y.test<-test[,yvar]
  y.train<-train[,yvar]
  }

  return(list(y.train=y.train,y.test=y.test,y.pred=y.pred))
}

data<-read.table(".\\Realdata_cluster\\Data\\docvisits_spline.txt",header=T,sep=",")
group=c(rep(1:5,each=3),(6:14))

yvar<-names(data)[1]
xvars<-names(data)[-1]
zvars<-xvars

#' Title
#'
#' @param data The real German healthcare dataset 
#' @param yvar @param yvar The name of the outcome variable
#' @param xvars The name of the predictor variables for the count model
#' @param zvars The name of the predictor variables for the zero model
#' @param cv.iter Number of iteration for calculating MASE using 5-fold CV mechanism
#' @param k.fold Number of CV folds (here it is 5)
#' @param seedval 
#' @param penalty Penalty
#'
#' @return Calculates MAE and MASE using the function accuracy from the forecast package

realdata.func<-function(data,yvar,xvars,zvars,cv.iter,k.fold,seedval,penalty)
{
  ptm<-proc.time()
  set.seed(seedval)
  folds <- cvFolds(nrow(data), K = k.fold, R = cv.iter)
  cv<-cbind.data.frame(folds$which,folds$subsets)
  names(cv)<-c("fold",paste("iter",1:cv.iter,sep=""))

  measures<-NULL

  for (i in 1:cv.iter)
  {
    cat("\n")
    print(paste("CV iter",i))
    cat("\n")

    y.testv<-NULL
    y.predv<-NULL
    y.trainv<-NULL

    for (k in 1:k.fold)
    {
      print(paste("fold",k))

      train.idx<-folds$subsets[folds$which!=k,i]
      train<-data[train.idx,]
      test<-data[-train.idx,]

      temp<-fit.out(train=train,test=test,yvar=yvar,xvars=xvars,zvars=zvars,group=group,penalty=penalty)

      y.test<-temp$y.test
      y.pred<-temp$y.pred
      y.train<-temp$y.train

      y.testv<-c(y.testv,y.test)
      y.predv<-c(y.predv,y.pred)
      y.trainv<-c(y.trainv,y.train)
    }

    forecast <- structure(list(mean=y.predv, fitted=y.testv, x=y.trainv), class='forecast')

    measures<-rbind(measures,c(iter=i,accuracy(forecast,y.test)[2,c(3,6)]))

  }
  measures<-apply(measures[,-1],2,function(x) return(round(median(x,na.rm=T),3)))
  time.taken<-round((proc.time()-ptm)[3],3)
  return(list(measures=measures,time=time.taken))
}

# result<-realdata.func(data=data,yvar=yvar,xvars=xvars,zvars=zvars,cv.iter=5,k.fold=5,seedval=123,penalty="EMLasso")
# result

