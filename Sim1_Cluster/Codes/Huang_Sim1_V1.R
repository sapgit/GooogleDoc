func_main<-function(rho,phi,n,ITER)
{
# data generating function
data.func.sim1<-function(n.train,n.test,p,ngrp,beta,gamma,rho,family)
{
  n<-n.train+n.test
  
  R<-matrix(rnorm(n*p),n,p)
  V<-matrix(0,ngrp,ngrp)
  for(i in 1:ngrp)
  {
    for(j in 1:ngrp)
    {
      V[i,j]=rho^(abs(i-j))
    }
  }
  Z<-mvrnorm(n,mu=rep(0,ngrp),Sigma=V)
  
  X<-matrix(0,n,p)
  size=rep(p/ngrp,ngrp)
  for(g in 1:ngrp)
  {
    for (j in 1:size[g])
    {
      X[,(g-1)*size[g]+j]<-(Z[,g]+R[,(g-1)*size[g]+j])/sqrt(2)
    }
  }
  X<-scale(X)
  colnames(X)<-paste("X",c(1:ncol(X)),sep="")
  xvars=colnames(X)
  zvars=xvars

  ## capture zero inflation ##
  rn<-round(runif(1)*10^5)
  set.seed(rn)
  y<-rzi(n,x=X,z=X,a=beta,b=gamma,family=family)
  
  set.seed(rn)
  ## capture from console output ##
  yout<-capture.output(rzi(n,x=X,z=X,a=beta,b=gamma,family=family))
  zeroinfl<-as.numeric(substring(yout[1],15))
  
  data<-cbind.data.frame(y,X)
  return(list(data=data,yvar="y",xvars=xvars,zvars=zvars,zeroinfl=zeroinfl))
}

# list of all methods compared
methods<-c("EMLasso","grLasso", "grMCP", "grSCAD", "gBridge")

# true parameter values used for count part 
beta<-c(-1, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 0.75, rep(0.2,8), rep(0,24))
beta<-c(5,beta) # beta0 changed from 1

# true parameter values used for zero part (excluding intercept)
gamma<-c(-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, rep(0.2,8), rep(0,24))

# determine gamma0 according to zero inflation (phi)
if (phi==0.3) gamma<-c(-1,gamma)
if (phi==0.4) gamma<-c(-0.5,gamma)
if (phi==0.5) gamma<-c(0,gamma)

# size for each group
size=rep(8,5)

# number of groups
ngrp<-length(size)

# generate a list of 200 datasets (containing both training and testing)
  data.list<-NULL
  for (i in 1:ITER) {
    dataset<-data.func.sim1(n.train=n,n.test=n,p=sum(size),ngrp=length(size),beta=beta,gamma=gamma,rho=rho,family="negbin")
    data.list<-c(data.list,list(dataset))
  }
  
  result.list<-list()
  
  for(method in methods)
  {
    measures.mat<-NULL
    for (i in 1:ITER)
    {
      dataset<-data.list[[i]]
      train<-dataset$data[1:n,]
      test<-dataset$data[-(1:n),]
      yvar<-dataset$yvar
      xvars<-dataset$xvars
      zvars<-dataset$zvars
      
      fit.formula<-as.formula(paste(yvar,"~",paste(paste(xvars,collapse="+"),"|",paste(zvars,collapse="+"),sep="")))
      
      # percentage of zero outcome
      pzero<-dataset$zeroinfl

      if (method=="EMLasso") {
        ptm<-proc.time()
        fit<-try(eval(call(paste("func_",method,sep=""),formula=fit.formula,data=train,dist="negbin")),silent=T) 
        time.taken<-round((proc.time()-ptm)[3],3)
      } else {
        ptm<-proc.time()
        fit<-try(eval(call(paste("func_",method,sep=""),data=train,yvar=yvar,xvars=xvars,zvars=zvars,group=rep(1:5,each=8),dist="negbin")),silent=T) 
        time.taken<-round((proc.time()-ptm)[3],3)
      }
      
      # if there is an error in calling func_methods
      if (class(fit)=="try-error")
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
        measures<-c(round(accuracy(forecast,y.test)[2,c(3,6)],4),time.taken=time.taken)
        measures.mat<-rbind(measures.mat,measures)
      }
    }
    print(paste(method,"complete"))
    result.list[[method]]<-measures.mat
  }
  return(result.list)
}

