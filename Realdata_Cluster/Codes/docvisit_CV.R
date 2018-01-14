func_main<-function(NUM,ITER)
{
  penalties<-c("grLasso", "grMCP", "grSCAD", "gBridge")

  group=c(rep(1:5,each=3),(6:14))
  
  yvar<-names(data)[1]
  xvars<-names(data)[-1]
  zvars<-xvars


  ################## compute MAE, MASE for different penalties #############
  measures.list<-list()
  coeff.list<-list()
  
  for(penalty in penalties)
  {
      measures<-NULL

      for (k in (NUM*ITER+1):((NUM+1)*ITER))
      {
        print(c(penalty,k))
        
        y.testv<-NULL
        y.predv<-NULL
        y.trainv<-NULL

        time.taken<-0
        for (i in 1:5)
        {
          train.idx<-cv.result[cv.result$fold!=i,(k+1)]
          train<-data[train.idx,]
          test<-data[-train.idx,]
          
          ptm<-proc.time()
          fit<-gooogle(data=train,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist="negbin",penalty=penalty)
          time.taken<-time.taken+round((proc.time()-ptm)[3],3)
          
          z.test<-as.matrix(cbind(1,test[,zvars]))
          phi.hat<-1/(1+exp(-z.test%*%fit$coefficients$zero))
          x.test<-as.matrix(cbind(1,test[,xvars]))
          lam.hat<-exp(x.test%*%fit$coefficients$count)
          y.pred<-(1-phi.hat)*lam.hat
          y.test<-test[,yvar]
          y.train<-train[,yvar]

          y.testv<-c(y.testv,y.test)
          y.predv<-c(y.predv,y.pred)
          y.trainv<-c(y.trainv,y.train)
        }
        
        forecast <- structure(list(mean=y.predv, fitted=y.testv, x=y.trainv), class='forecast')

        measures<-rbind(measures,c(iter=k,accuracy(forecast,y.test)[2,c(3,6)],time.taken=time.taken))
      }

      measures.list[[penalty]]<-measures
  }

  
  ################## ------ Compute same measures for EM ------ ############
  
  fit.formula<-as.formula(paste(yvar,"~",paste(paste(xvars,collapse="+"),"|",paste(zvars,collapse="+")),sep=""))
  
  measures_EM<-NULL
  
  for (k in (NUM*ITER+1):((NUM+1)*ITER))
  {
    print(c("EM_Lasso",k))
    
    y.testv<-NULL
    y.predv<-NULL
    y.trainv<-NULL

    time.taken<-0
    for (i in 1:5)
    {
      train.idx<-cv.result[cv.result$fold!=i,(k+1)]
      train<-data[train.idx,]
      test<-data[-train.idx,]
      
      ptm<-proc.time()
      fit.em<-zipath(formula=fit.formula,data=train,family="negbin")
      time.taken<-time.taken+round((proc.time()-ptm)[3],3)
      
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

      y.testv<-c(y.testv,y.test)
      y.predv<-c(y.predv,y.pred)
      y.trainv<-c(y.trainv,y.train)
    }
    
    forecast <- structure(list(mean=y.predv, fitted=y.testv, x=y.trainv), class='forecast')

    measures_EM<-rbind(measures_EM,c(iter=k,accuracy(forecast,y.test)[2,c(3,6)],time=time.taken))
  }

  measures.list[["EM"]]<-measures_EM
  return(measures.list)
}
