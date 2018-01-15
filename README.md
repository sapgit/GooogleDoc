# Supplementary Codes for the Manuscript "Group Regularization for Zero-inflated Regression Models with Application to Healthcare Demand in Germany"


## Simulation 1
In this simulation, we generate 40 continuous predictors which are grouped into 5 groups consisting of 8 predictors each. The covariates are generated from multivariate normal distribution with `rho` as the correlation parameter of the AR(1) covariance matrix. The wrapper function `master.func` in [master.R](master.R) generates a pre-specified number of datasets for a given training and test data size `n.train` and `n.test`, the correlation parameter `rho` and zero abundance parameter `phi` from ZINB model using the data generation functions from [datagen.R](datagen.R). This function then fits the one of the five methods (apart from EMLasso, Gooogle method with four penalties - grLasso, grMCP, grSCAD and gBridge) on the training data using [fit_method.R](fit_method.R). Lastly, the prediction functions from [predict_measures.R](predict_measures.R) are used to calculate median of MAE and MASE along with bootstrap standard error based on the test datasets. Apart from the prediction measures, the master function also outputs the median run time per iteration. 


### Install Gooogle package
```
install.packages("devtools")
devtools::install_github("himelmallick/Gooogle")
library(Gooogle)
```

### Load other packages
```
require(MASS)
require(forecast)
require(mpath)
require(cvTools)
require(dummies)
```

### Specify the true parameter values
```
## Train and test sample size
n.train<-n.test<-c(200,500,1000)
## AR(1) correlation
rho<-c(0,0.4,0.8)
## methods considered
method<-c("grLasso","grMCP","grSCAD","gBridge","EMLasso")
## zero inflation
phi<-c(0.3,0.4,0.5)

## True model coefficients
 count model: beta<-c(5, -1, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 0.75, rep(0.2,8), rep(0,24))    
    
 zero model: 
 For phi=0.3, gamma<-c(-1,-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, rep(0.2,8), rep(0,24))
 For phi=0.4, gamma<-c(-0.5,-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, rep(0.2,8), rep(0,24))
 For phi=0.5, gamma<-c(0,-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, rep(0.2,8), rep(0,24))
```
    
### Examples

Below are some examples of the implementation on the simulation functions to compare MAE and MASE between EM Lasso and Gooogle method with a range of penalties based on 10 replications across different zero inflation proportions.

```
## phi=0.3
## group lasso 
master.func(n.train=500,n.test=500,rho=0,phi=0.3,method="grLasso",sim=1,ITER=10,seedval=123) 
## EM lasso
master.func(n.train=500,n.test=500,rho=0,phi=0.3,method="EMLasso",sim=1,ITER=10,seedval=123)  
 ```
 |Method | MAE  | MASE  | Time/iteration (sec)|
 |:-------|:-------:|:------:|---------:|
 |Gooogle (grLasso)|167.227(8.02) |0.784(0.06) |1.08 |
 |EMLasso|236.041(10.2) |1.095(0.07) |12.93|
 
 ```
 ## phi=0.4
## group SCAD 
master.func(n.train=500,n.test=500,rho=0,phi=0.4,method="grSCAD",sim=1,ITER=10,seedval=123) 
## EM lasso
master.func(n.train=500,n.test=500,rho=0,phi=0.4,method="EMLasso",sim=1,ITER=10,seedval=123)  

 ```
 |Method | MAE  | MASE  | Time/iteration (sec)|
 |:-------|:-------:|:------:|---------:|
 |Gooogle (grSCAD)|134.357(7.91) |0.752(0.07) |0.77 |
 |EMLasso|210.695(8.66) |1.191(0.08) |14.045|
 ```
 
 ## phi=0.5
## group Bridge
master.func(n.train=500,n.test=500,rho=0,phi=0.5,method="gBridge",sim=1,ITER=10)
## EM lasso
master.func(n.train=500,n.test=500,rho=0,phi=0.5,method="EMLasso",sim=1,ITER=10)
 ```
 |Method | MAE  | MASE  | Time/iteration (sec)|
 |:-------|:-------:|:------:|---------:|
 |Gooogle (gBridge)|86.318(4.8) |0.703(0.03) |0.725|
 |EMLasso|173.184(10.67)|1.466(0.09)|12.31|

 
 ### Original Codes
The original simulation was performed in a high performance cluster because of the high computation cost. The codes are given in the [Sim1_Cluster](Sim1_Cluster) folder. HOWEVER, PLEASE NOTE THAT THE COMPUTATION TIME MIGHT BE HUGE.

### Simulation 2

In simulation 2 we generate 6 continuous covariates each of which forms a singleton group. Four more variables (two each from X3 and X6) are further polynomially constructed, giving rise to a total of 10 continuous predictors. For constructing the categorical variables, we generate 5 continuous variables from a multivariate normal distribution and quantile-discretize each of them into 5 new variables based on their quantiles. This leads to a combination of 20 categorical predictors with 5 non-overlapping groups of equal size. 

The following zero and count model coefficients have been used to generate the data under this simulation setting.

```
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

Similar to Simulation 1, below are some examples illustrating MAE and MASE comparisons between Gooogle and EM Lasso across a range of zero abundance values. 

```
## phi=0.3
## group lasso 
master.func(n.train=500,n.test=500,rho=0,phi=0.3,method="grLasso",sim=2,ITER=10,seedval=123) 
## EM lasso
master.func(n.train=500,n.test=500,rho=0,phi=0.3,method="EMLasso",sim=2,ITER=10,seedval=123)  
 ```
 |Method | MAE  | MASE  | Time/iteration (sec)|
 |:-------|:-------:|:------:|---------:|
 |Gooogle (grLasso)|179.172(6.49) |0.737(0.04) |0.48 |
 |EMLasso|236.041(10.2) |1.02(0.04) |17.885|
 
 ```
 ## phi=0.4
## group SCAD 
master.func(n.train=500,n.test=500,rho=0,phi=0.4,method="grSCAD",sim=2,ITER=10,seedval=123) 
## EM lasso
master.func(n.train=500,n.test=500,rho=0,phi=0.4,method="EMLasso",sim=2,ITER=10,seedval=123)  

 ```
 |Method | MAE  | MASE  | Time/iteration (sec)|
 |:-------|:-------:|:------:|---------:|
 |Gooogle (grSCAD)|118.423(2.99) |0.678(0.03) |0.385 |
 |EMLasso|194.301(5.56) |1.08(0.05) |16.04|
 
 ```
 ## phi=0.5
## group Bridge
master.func(n.train=500,n.test=500,rho=0,phi=0.5,method="gBridge",sim=2,ITER=10)
## EM lasso
master.func(n.train=500,n.test=500,rho=0,phi=0.5,method="EMLasso",sim=2,ITER=10)
 ```
 |Method | MAE  | MASE  | Time/iteration (sec)|
 |:-------|:-------:|:------:|---------:|
 |Gooogle (gBridge)|80.228(3.48)) |0.632(0.04) |0.48|
 |EMLasso|158.726(6.67) |1.18(0.08) |23.555|
 
  ### Original Codes
The original simulation was performed in a high performance cluster because of the high computation cost. The codes are given in the [Sim2_Cluster](Sim2_Cluster) folder. HOWEVER, PLEASE NOTE THAT THE COMPUTATION TIME MIGHT BE HUGE.

   
## Real Data 

We compare the performance of Gooogle method with EM Lasso on German Health Care demand [dataset](Realdata_Cluster/Data/). For the sake of meaningful comparison we implement 5 fold cross validation where we randomly partition the dataset into five equal folds - four folds to train the model and one fold for test. Prediction errors are calculated based on the comprehensive set of held out samples. To employ this method we use the function `realdata.func` in  [realdata.R](realdata.R) which takes the dataset, the outcome and predictor variable names, number of CV iterations, number of folds, seed value and the selected method as input and returns the MASE and MAE along with total computation time. 

Below is an example code for calculating the median of MAE and MASE over 5 iterations using 5-fold CV for the Gooogle method using "gBridge" penalty and for "EMLasso". 

```
source("realdata.R")

## Gooogle (grLasso)
result.gooogle<-realdata.func(data=data,yvar=yvar,xvars=xvars,zvars=zvars,cv.iter=5,k.fold=5,seedval=123,method="grLasso")
result.gooogle

## EM Lasso
result.em<-realdata.func(data=data,yvar=yvar,xvars=xvars,zvars=zvars,cv.iter=5,k.fold=5,seedval=123,method="EMLasso")
result.em
```

 |Method | MAE     | MASE  | Time (sec) |
 |-------|:-------:|:------:|----------:|
 |Gooogle (grLasso)|2.9844   |0.9616 | 20.83|
 |EMLasso|3.1874   |1.0270 | 480.18     |
 
 
 ### Original Codes
The original codes were run in a high performance cluster because of the high computation cost. The codes are given in the [Realdata_Cluster](Realdata_Cluster) folder. HOWEVER, PLEASE NOTE THAT THE COMPUTATION TIME MIGHT BE HUGE.
