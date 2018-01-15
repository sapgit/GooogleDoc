#' Title
#' @description This function generates one dataset according to the specified parameters for simulation 1
#'
#' @param n.train Sample size in the training dataset
#' @param n.test Sample size in the test dataset
#' @param grpsize Group size
#' @param rho Correlation parameter of AR(1) covariance matrix for multivariate normal distribution
#' @param phi Zero-inflation parameter
#'
#' @return Dataset along with the variable names and the proportion of zero inflation
datagen.sim1.func<-function(n.train,n.test,grpsize,rho,phi,seedval)
{
  # true parameter values used for count part 
  beta<-c(-1, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 0.75, rep(0.2,8), rep(0,24))
  beta<-c(5,beta) 
  
  # true parameter values used for zero part (excluding intercept)
  gamma<-c(-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, rep(0.2,8), rep(0,24))
  
  # determine gamma0 according to zero inflation (phi)
  if (phi==0.3) gamma<-c(-1,gamma)
  if (phi==0.4) gamma<-c(-0.5,gamma)
  if (phi==0.5) gamma<-c(0,gamma)
  
  if (!is.null(seedval)) set.seed(seedval)

  n=n.train+n.test

  p=sum(grpsize)
  ngrp=length(grpsize)

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
  for(g in 1:ngrp)
  {
    for (j in 1:grpsize[g])
    {
      X[,(g-1)*grpsize[g]+j]<-(Z[,g]+R[,(g-1)*grpsize[g]+j])/sqrt(2)
    }
  }
  X<-scale(X)
  colnames(X)<-paste("X",c(1:ncol(X)),sep="")
  xvars=colnames(X)
  zvars=xvars

  ## capture zero inflation ##
  if(is.null(seedval))
  {
    rn<-round(runif(1)*10^5)
    set.seed(rn)
    y<-rzi(n=n,x=X,z=X,a=beta,b=gamma,family="negbin")

    set.seed(rn)
    ## capture from console output ##
    yout<-capture.output(rzi(n,x=X,z=X,a=beta,b=gamma,family="negbin"))
    zeroinfl<-as.numeric(substring(yout[1],15))
  } else {
    y<-rzi(n=n,x=X,z=X,a=beta,b=gamma,family="negbin")

    ## capture from console output ##
    yout<-capture.output(rzi(n,x=X,z=X,a=beta,b=gamma,family="negbin"))
    zeroinfl<-as.numeric(substring(yout[1],15))
  }

  data<-cbind.data.frame(y,X)
  return(list(data=data,yvar="y",xvars=xvars,zvars=zvars,zeroinfl=zeroinfl))
}

#' Title
#' @description This function generates one dataset according to the specified parameters for simulation 2
#'
#' @param n.train Sample size in the training dataset
#' @param n.test Sample size in the test dataset
#' @param rho Correlation parameter of AR(1) covariance matrix for multivariate normal distribution
#' @param grpsize Group size
#' @param phi Zero-inflation parameter
#'
#' @return Dataset along with the variable names and the proportion of zero inflation
#'
datagen.sim2.func<-function(n.train,n.test,rho,grpsize,phi,seedval)
{
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
  
  if (!is.null(seedval)) set.seed(seedval)
  
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
  if(is.null(seedval))
  {
    rn<-round(runif(1)*10^5)
    set.seed(rn)
    y<-rzi(n=n,x=X,z=X,a=beta,b=gamma,family="negbin")

    set.seed(rn)
    ## capture from console output ##
    yout<-capture.output(rzi(n,x=X,z=X,a=beta,b=gamma,family="negbin"))
    zeroinfl<-as.numeric(substring(yout[1],15))
  } else {
    y<-rzi(n=n,x=X,z=X,a=beta,b=gamma,family="negbin")

    ## capture from console output ##
    yout<-capture.output(rzi(n,x=X,z=X,a=beta,b=gamma,family="negbin"))
    zeroinfl<-as.numeric(substring(yout[1],15))
  }

  data<-cbind.data.frame(y,X)
  return(list(data=data,yvar="y",xvars=xvars,zvars=zvars,zeroinfl=zeroinfl))
}


#' Title
#' @description This function generates one dataset according to the specified parameters for simulation 1
#'
#' @param n.train Sample grpsize in the training dataset
#' @param n.test Sample grpsize in the test dataset
#' @param rho Correlation parameter of AR(1) covariance matrix for multivariate normal distribution
#' @param phi Zero-inflation parameter
#' @param grpsize Group size
#' @param sim Simulation number
#' @param ITER Iteration number
#' @return ITER number of datasets in the list form `data.list`
datagen.sim.all<-function(n.train,n.test,rho,phi,grpsize,sim,ITER,seedval)
{
  data.list<-NULL

  if (!is.null(seedval))
  {
    set.seed(seedval)
    seedvals<-round(runif(ITER)*10^5)
  } else {
    seedvals<-rep(NULL,ITER)
  }

  if(sim==1)
  {
    for (i in 1:ITER) {
      dataset<-datagen.sim1.func(n.train=n.train,n.test=n.test,grpsize=grpsize,rho=rho,phi=phi,seedval=seedvals[i])
      data.list<-c(data.list,list(dataset))
    }
  }else{for (i in 1:ITER) {
    dataset<-datagen.sim2.func(n.train=n.train,n.test=n.test,grpsize=grpsize,rho=rho,phi=phi,seedval=seedvals[i])
    data.list<-c(data.list,list(dataset))
  }
  }

  return(data.list)
}
