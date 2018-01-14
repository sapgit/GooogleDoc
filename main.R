# setwd("C:\\Users\\schowdhury\\Dropbox\\Gooogle\\Codes\\Github\\GooogleDoc")
source("penalty.R")
source("datagen_sim.R")
source("fit_method.R")
source("predict_measures.R")
source("master.R")

require(MASS)
require(forecast)
require(Gooogle)
require(mpath)
require(dummies)

############################ Sim 1 ########################################
beta<-c(-1, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 0.75, rep(0.2,8), rep(0,24))
beta<-c(5,beta) # beta0 changed from 1

# true parameter values used for zero part (excluding intercept)
gamma<-c(-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, rep(0.2,8), rep(0,24))
if (phi==0.3) gamma<-c(-1,gamma)
if (phi==0.4) gamma<-c(-0.5,gamma)
if (phi==0.5) gamma<-c(0,gamma)

# specify group size
grpsize=rep(8,5)

# Example
master.func(n.train=200,n.test=200,beta=beta,gamma=gamma,grpsize=grpsize,rho=0,phi=0.3,family="negbin",method="grLasso",sim=1,ITER=10,seedval=123)

############################ Sim 2 ########################################

# true parameter values used for count part
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

# specify group size
grpsize=c(1,1,3,1,1,3,4,4,4,4,4)


# Example
master.func(n.train=200,n.test=200,beta=beta,gamma=gamma,grpsize=grpsize,rho=0,phi=0.3,family="negbin",method="grLasso",sim=2,ITER=10,seedval=NULL)
