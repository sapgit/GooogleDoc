# Supplementary Codes for the Manuscript "Group Regularization for Zero-inflated Regression Models with Application to Healthcare Demand in Germany"


### Example 1
In this example, we generate 40 continuous predictors which are grouped into 5 groups consisting of 8 predictors each. The covariates are generated from multivariate normal distribution with `rho` as the correlation parameter of the AR(1) covariance matrix. The wrapper function `master.func`  generates a pre-specified number of datasets for a given training and test data size `n.train` and `n.test`, the correlation parameter `rho` and zero abundance parameter `phi` from ZINB model and fits the `gooogle` function on the training data and calculates median MAE and MASE for both the count and zero model based on the test datasets. The values of the regression coefficients for the count model (beta) and those for the zero model (gamma) are provided inside the function. 

```
## Install Gooogle package
install.packages("devtools")
devtools::install_github("himelmallick/Gooogle")
library(Gooogle)
```


```
## Load other packages

require(MASS)
require(forecast)
require(mpath)
require(dummies)
```

```
 ## Specify the true parameter values
  
 count model: beta<-c(5, -1, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 0.75, rep(0.2,8), rep(0,24))    
    
 zero model: 
 For phi=0.3, gamma<-c(-1,-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, rep(0.2,8), rep(0,24))
 For phi=0.4, gamma<-c(-0.5,-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, rep(0.2,8), rep(0,24))
 For phi=0.5, gamma<-c(0,-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, rep(0.2,8), rep(0,24))
```
    
```
 ## The main function for running the simulation and calculating the predictive measures
 
 master.func<-function(n.train,n.test,rho,phi,family,method,sim,ITER)
{
  if(sim==1)
  {
    # size for each group
    size=rep(8,5)
    p=sum(size)
    # number of groups
    ngrp<-length(size) 
  }else{
    size=c(1,1,3,1,1,3,4,4,4,4,4)
    p=sum(size)
    ngrp<-length(size)
  }
  group<-NULL
  for (k in 1:ngrp)
  {
    group<-c(group,rep(k,size[k]))
  }
  ### Generate the data ###
  data.list<-datagen.sim.all(n.train=n.train,n.test=n.test,beta=beta,gamma=gamma,rho=rho,phi=phi,family=family,sim=sim,ITER)
  
  ### Calculate the predictive measures ###
  measures<-measures.summary(n.train=n.train,n.test=n.test,data.list=data.list,method=method,ITER=ITER,group=group,family=family)
  
  return(measures)
}
```

Below are some examples comparing MAE and MASE of group Lasso and EM Lasso based on 10 replications for different zero inflation proportions (`phi`).

```
## phi=0.3
## group lasso 
master.func(n.train=200,n.test=200,beta=beta,gamma=gamma,rho=0,phi=0.3,family="negbin",method="grLasso",sim=1,ITER=10)
## EM lasso
master.func(n.train=200,n.test=200,beta=beta,gamma=gamma,rho=0,phi=0.3,family="negbin",method="EMLasso",sim=1,ITER=10)
 
 ```
 |Method | MAE     | MASE  |
 |-------|:-------:|------:|
 |grLasso|148.7548 |0.7341 | 
 |EMLasso|232.3708 |1.1582 |
 
 ```
 ## phi=0.4
## group Lasso
master.func(n.train=200,n.test=200,beta=beta,gamma=gamma,rho=0,phi=0.4,family="negbin",method="grLasso",sim=1,ITER=10)
## EM lasso
master.func(n.train=200,n.test=200,beta=beta,gamma=gamma,rho=0,phi=0.4,family="negbin",method="EMLasso",sim=1,ITER=10)

 ```
 |Method | MAE     | MASE  |
 |-------|:-------:|------:|
 |grLasso|137.4921 |0.6958 | 
 |EMLasso|192.1370 |0.8552 |
 
 ```
 ## phi=0.5
## group Lasso
master.func(n.train=200,n.test=200,beta=beta,gamma=gamma,rho=0,phi=0.5,family="negbin",method="grLasso",sim=1,ITER=10)
## EM lasso
master.func(n.train=200,n.test=200,beta=beta,gamma=gamma,rho=0,phi=0.5,family="negbin",method="EMLasso",sim=1,ITER=10)
 ```
 |Method | MAE     | MASE  |
 |-------|:-------:|------:|
 |grLasso|100.7777 |0.7067 | 
 |EMLasso|160.3597 |0.8542 |
           

### Example 2

In example 2 we generate 6 continuous covariates each of which forms a singleton group. Four more variables (two each from X3 and X6) are further polynomially constructed, giving rise to a total of 10 continuous predictors. For constructing the categorical variables, we generate 5 continuous variables from a multivariate normal distribution and quantile-discretize each of them into 5 new variables based on their quantiles. This leads to a combination of 20 categorical predictors with 5 non-overlapping groups of equal size.

```
## The true parameter values

For count model: 
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

For zero model:
phi=0.3: gamma<-c(-1.4,beta[-1])
phi=0.4: gamma<-c(-.7,beta[-1])
phi=0.5: gamma<-c(-.15,beta[-1])
```

Below are some examples illustrating MAE and MASE comparisons between Group Lasso and EM Lasso across a range of zero abundance values. 

```
 ## phi=0.3
 master.func(n.train=200,n.test=200,beta=beta,gamma=gamma,rho=0,phi=0.3,family="negbin",method="grLasso",sim=2,ITER=10)
 master.func(n.train=200,n.test=200,beta=beta,gamma=gamma,rho=0,phi=0.3,family="negbin",method="EMLasso",sim=2,ITER=10)
 ```
 |Method | MAE     | MASE  |
 |-------|:-------:|------:|
 |grLasso|145.8202 |0.5172 | 
 |EMLasso|225.4927 |0.9596|
 
 ```
 ## phi=0.4
 master.func(n.train=200,n.test=200,beta=beta,gamma=gamma,rho=0,phi=0.4,family="negbin",method="grLasso",sim=2,ITER=10)
 master.func(n.train=200,n.test=200,beta=beta,gamma=gamma,rho=0,phi=0.4,family="negbin",method="EMLasso",sim=2,ITER=10)
 ```
 |Method | MAE     | MASE  |
 |-------|:-------:|------:|
 |grLasso|127.3018 |0.8545 | 
 |EMLasso|170.8427 |0.9443 |
 
 ```
 ## phi=0.5
 master.func(n.train=200,n.test=200,beta=beta,gamma=gamma,rho=0,phi=0.5,family="negbin",method="grLasso",sim=2,ITER=10)
 master.func(n.train=200,n.test=200,beta=beta,gamma=gamma,rho=0,phi=0.5,family="negbin",method="EMLasso",sim=2,ITER=10)
 ```
 |Method | MAE     | MASE  |
 |-------|:-------:|------:|
 |grLasso|85.2794  |0.7185 | 
 |EMLasso|150.7711 |1.1151 |
 
   
## Real Data Example

We compare the performance of Gooogle method with EM Lasso on German Health Care demand dataset (link). For the sake of meaningful comparison we implement 5 fold cross validation where we randomly partition the dataset into five equal folds - four folds to train the model and one fold for test. Prediction errors are calculated based on the comprehensive set of held out samples. To employ this method we use the `realdata.func` function which takes the dataset, the outcome and predictor variable names, number of CV iterations, number of folds, seed value and the concerning method as input and returns the MASE and MAE. 

Below is an example code for calculating the median of MAE and MASE over 5 iterations using 5-fold CV for the Gooogle method using "gBridge" penalty and for "EMLasso". 

```
result.gooogle<-realdata.func(data=data,yvar=yvar,xvars=xvars,zvars=zvars,cv.iter=5,k.fold=5,seedval=123,method="Gooogle")

result.em<-realdata.func(data=data,yvar=yvar,xvars=xvars,zvars=zvars,cv.iter=5,k.fold=5,seedval=123,method="EMLasso")
```

 |Method | MAE     | MASE  |
 |-------|:-------:|------:|
 |gBridge|2.9844   |0.9616 | 
 |EMLasso|3.1874   |1.0270 |
 
 
