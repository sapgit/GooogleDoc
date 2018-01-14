func_main<-function(rho,phi,n,ITER)
{
  # data generating function
  data.func.sim2<-function(n.train,n.test,size,beta,gamma,rho,family) 
  {
    n<-n.train+n.test
    p1=6;p2=5;ngrp1=6;ngrp2=5;
    R<-matrix(rnorm(n*p1),n,p1)
    w<- rnorm(n,0,1)
    x1<-matrix(0,n,p1)
    
    for(g in 1:ngrp1)
    {
      x1[,g]<- (R[,g]+w)/(sqrt(2))
    }
    xc<-cbind(x1[,1],x1[,2],x1[,3],(x1[,3])^2,(x1[,3])^3,x1[,4],x1[,5],x1[,6],(x1[,6])^2,(x1[,6])^3)
    
    V<-matrix(0,ngrp2,ngrp2)
    
    for(i in 1:ngrp2)
    {
      for(j in 1:ngrp2)
      {
        V[i,j]=rho^(abs(i-j))
      }
    }
    
    x2<-matrix(mvrnorm(n*p2,mu=rep(0,ngrp2),Sigma=V),n,p2)
    
    for(i in 1:n)
    {
      for (j in 1:p2)
      {
        if (x2[i,j] < qnorm(0.2,0,1)){
          x2[i,j]=0
        } else if (qnorm(0.2,0,1) < x2[i,j] && x2[i,j] < qnorm(0.4,0,1)){
          x2[i,j]=1
        } else if (qnorm(0.4,0,1) < x2[i,j] && x2[i,j] < qnorm(0.6,0,1)){
          x2[i,j]=2
        } else if (qnorm(0.6,0,1) < x2[i,j] && x2[i,j] < qnorm(0.8,0,1)){
          x2[i,j]=3
        } else {
          x2[i,j]=4
        }
      }
    }
    
    xd<-NULL
    
    for(j in 1:p2)
    {
      xd<-cbind(xd,dummy(x2[,j]))
    }
    xd<-xd[,-seq(p2,p2*ngrp2,ngrp2)]
    
    X<-cbind(scale(xc),xd)
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

  betag1<-c(0)
  betag2<-c(0)
  betag3<-c(-0.1,0.2,0.1)
  betag4<-c(0)
  betag5<-c(0)
  betag6<-c(2/3,-1,1/3)
  betag7<-c(-2,-1,1,2)
  betag8<-c(0,0,0,0)
  betag9<-c(0,0,0,0)
  betag10<-rep(0,4)
  betag11<-c(0,0,0,0)
  
    
  beta<-c(5,betag1,betag2,betag3,betag4,betag5,betag6,betag7,betag8,betag9,betag10,betag11)
  
  # determine gamma0 according to zero inflation (phi)
  if (phi==0.3) gamma<-c(-1.4,beta[-1])
  if (phi==0.4) gamma<-c(-0.7,beta[-1])
  if (phi==0.5) gamma<-c(-0.15,beta[-1])
  
  # size for each group
  size=c(1,1,3,1,1,3,4,4,4,4,4)
  ngrp<-length(size)
  
  # group
  group<-NULL
  for (k in 1:ngrp)
  {
    group<-c(group,rep(k,size[k]))
  }
  
  # generate a list of 200 datasets (containing both training and testing)
  data.list<-NULL
  for (i in 1:ITER) {
    dataset<-data.func.sim2(n.train=n.train,n.test=n.test,beta=beta,gamma=gamma,rho=rho,family="negbin")
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
        fit<-try(eval(call(paste("func_",method,sep=""),data=train,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist="negbin")),silent=T) 
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