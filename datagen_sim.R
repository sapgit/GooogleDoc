#' Title
#' @description This function generates one dataset according to the specified parameters for simulation 1
#'
#' @param n.train Sample size in the training dataset
#' @param n.test Sample size in the test dataset
#' @param beta True regression coeffficients in the count model
#' @param gamma True regression coeffficients in the zero model
#' @param rho Correlation parameter of AR(1) covariance matrix for multivariate normal distribution
#' @param phi Zero-inflation parameter
#' @param family Distribution of the count model which is Negative Binomial in this case
#'
#' @return Dataset along with the variable names and the proportion of zero inflation
datagen.sim1.func<-function(n.train,n.test,beta,gamma,rho,phi,family)
{
  n=n.train+n.test
  p=40
  ngrp=5

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
  y<-rzi(n=n,x=X,z=X,a=beta,b=gamma,family=family)

  set.seed(rn)
  ## capture from console output ##
  yout<-capture.output(rzi(n,x=X,z=X,a=beta,b=gamma,family=family))
  zeroinfl<-as.numeric(substring(yout[1],15))

  data<-cbind.data.frame(y,X)
  return(list(data=data,yvar="y",xvars=xvars,zvars=zvars,zeroinfl=zeroinfl))
}

#' Title
#' @description This function generates one dataset according to the specified parameters for simulation 2
#'
#' @param n.train Sample size in the training dataset
#' @param n.test Sample size in the test dataset
#' @param beta True regression coeffficients in the count model
#' @param gamma True regression coeffficients in the zero model
#' @param rho Correlation parameter of AR(1) covariance matrix for multivariate normal distribution
#' @param phi Zero-inflation parameter
#' @param family Distribution of the count model which is Negative Binomial in this case
#'
#' @return Dataset along with the variable names and the proportion of zero inflation
#'
datagen.sim2.func<-function(n.train,n.test,beta,gamma,rho,phi,family)
{
  n=n.train+n.test
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


#' Title
#' @description This function generates one dataset according to the specified parameters for simulation 1
#'
#' @param n.train Sample size in the training dataset
#' @param n.test Sample size in the test dataset
#' @param beta True regression coeffficients in the count model
#' @param gamma True regression coeffficients in the zero model
#' @param rho Correlation parameter of AR(1) covariance matrix for multivariate normal distribution
#' @param phi Zero-inflation parameter
#' @param family Distribution of the count model which is Negative Binomial in this case
#' @param sim Simulation number
#' @param ITER Iteration number
#' @return ITER number of datasets in the list form `data.list`
datagen.sim.all<-function(n.train,n.test,beta,gamma,rho,phi,family,sim,ITER)
{
  data.list<-NULL

  if(sim==1)
  {
    for (i in 1:ITER) {
      dataset<-datagen.sim1.func(n.train=n.train,n.test=n.test,beta=beta,gamma=gamma,rho=rho,phi=phi,family=family)
      data.list<-c(data.list,list(dataset))
    }
  }else{for (i in 1:ITER) {
    dataset<-datagen.sim2.func(n.train=n.train,n.test=n.test,beta=beta,gamma=gamma,rho=rho,phi=phi,family=family)
    data.list<-c(data.list,list(dataset))
  }
  }

  return(data.list)
}
