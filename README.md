# Group Regularization for Zero-inﬂated Count Regression Models with An Application to Healthcare Demand in Germany

## Introduction

In many biomedical applications, covariates are naturally grouped, with variables in the same group being systematically related or statistically correlated. Under such setting [Gooogle](https://github.com/himelmallick/Gooogle) provides a unified algorithm for group regularization in presence of zero abundant count outcome data using a least squares approximation of the mixture likelihood and a variety of group-wise penalties on the coefﬁcients. We investigate the finite sample performance of these methods through extensive simulation experiments. In this repository, we provide the R functions of the simulation and some basic examples on those functions. In addition to the simulation settings, we applied our method on German Healthcare demand dataset and compared the performance with EM Lasso. The associated R codes are also provided in this repository.   


## Simulation 
We conduct two extensive simulation studies to investigate the statistical properties of our proposed group regularization approach across a wide range of scenarios. To effectively evaluate the group selection performance of various methods, we prefix the number of groups to a moderately large number, and vary the group sizes from equal to unequal, and within each group the coefficients are set all zero, all nonzero or a mixture of both. We partition each generated synthetic data into a training set and a test set; models are fitted on the training set and MAE and MASE's are calculated on the test set. We also note the percentage of times correct groups are selected by each method. We consider the sample size in training and test group as {200,200}, {500,500} and {1000,1000}, the zero inflation parameter 0.3, 0.4 and 0.5, the correlation parameter of the covariance matrix of the covariates 0, 0.4 and 0.8. Below we explain each function that we have used to generate the data and calculate the predictive measures that are presented in our simulation study in the manuscript.

### The R Codes

#### Function for data generation: `datagen_sim.R`
This R file contains three functions for generating simulated data:

1. `datagen.sim1.func` & `datagen.sim2.func`: These are the data generating functions for simulation 1 and 2 respectively which generate one dataset according to the specified parameters. 

* **n:** sample size in each of training and testing group
* __p:__ Number of covariates
* __ngrp:__ Number of covariates in each group (constant for this example)
* __beta:__ True regression coefficients in the count model
* __gamma:__ True regression coefficients in the zero model
* __rho:__ Correlation parameter of AR(1) of the covariance matrix for multivariate normal distribution
* __sim:__ "1" for simulation 1, else takes "2"
* __family:__ The distribution of the count model ("negbin"). 

These functions generate the dataset, zero-inflated outcome variable, covariates for the count model (X) which are assumed to be equal to that of the zero model (Z), and the proportion of zero inflation.

2. `datagen.sim.all`: The predictive measures require the simulation datasets to be generated multiple times. The function `datagen.sim.all` uses the function `datagen.sim1.func` to generate a prespecified number ("ITER") of the simulated data set for a given set of parameters. It outputs all the datasets in the list form (`data.list`).

#### Function to fit the gooogle method on the training dataset: `fit_method`
This R file contains the function `fit.method` which takes the the following as arguments:

* __dataset:__ generated dataset using datagen.sim1.func
* __yvar,xvars,zvars:__ corresponding variables
* __method:__ "grLasso" (group Lasso), "gBridge" (group bridge), "grSCAD" (group SCAD), "grMCP" (group MCP) or "EMLasso" (EM Lasso) 
* __group:__ The vector specifying grouping structure of the covariates
* __dist:__ "negbin"

Depending on the method specified in the argument, this function calls the `gooogle` function from the `Gooogle` package to fit the zero inflated dataset with grouped covariates or calls the function `func_EMLasso` (also present in this R file) which fits the data with Lasso penalty using the zipath function. The functon `fit.method` outputs the fitted coefficients for count model and zero model along with the AIC, BIC and loglikelihood of the zero-inflated model.

### Calculate the predictive measure: `predict_measures.R`
This R file contains two functions described below:

1. `measures.func`: This function takes the following arguments:

* __train:__ Training dataset obtained from the original dataset 
* __test:__ Complement of the training dataset used for prediction
* __fit:__ The output of the function fit.method
* __yvar,xvars,zvars:__ Corresponding variables

For a given training and test dataset and the fitted coefficients this function calculates and returns the predictive measure MAE and MASE by using the `accuracy` function from the `forecast` package.

2. `measures.summary`: This function takes the following arguments:

* __n:__  Sample size in the training dataset
* __data.list:__ output of `datagen.sim1.all`
* __method:__ "grLasso" (group Lasso), "gBridge" (group bridge), "grSCAD" (group SCAD), "grMCP" (group MCP) or "EMLasso" (EM Lasso).
* __ITER:__ number of simulations
* __group:__ grouping of the covariates
* __family:__ "negbin"

This function fits the gooogle method to the training part of a dataset in data.list and calculate the predicted value using the fitted coefficients from the test data. It calls the function `measures.func` to calculate MASE. This is repeated for "ITER" number of datasets and outputs the median MAE and MASE.

## Example 1
In this example we generate 40 continuous predictors which are grouped into 5 groups consisting of 8 predictors each. The covariates are generated from multivariate normal distribution with `rho` as the corrrelation parameter of the AR(1) covariance matrix. The following function master.func generates 10 datasets for a given n, rho and phi from ZIP model and fits the google function and calculates median MASE as well percentage of correct group selection for both the count and zero model. The values of the regression coefficients for the count model (beta) and those for the zero model (gamma) are given inside the function. 

```
## Load the packages

require(MASS)
require(forecast)
require(Gooogle)
require(mpath)
require(dummies)
```

```
 ## Specify the true parameter values
  
 count model: beta<-c(5,-1, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 0.75, rep(0.2,8), rep(0,24))    
    
 zero model: gamma<-c(-1,-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, rep(0.2,8), rep(0,24))
 For phi=0.3, gamma_0 (intercept) of the zero model is -1
 
 group: (8,8,8,8,8)
   ```
    
  ```
 ## Generate the list of datasets
  
  data.list<-datagen.sim.all(n=200,beta=beta,gamma=gamma,rho=0.4,phi=0.3,family="negbin",sim=1,ITER=10)
  ```
  
  ```
 ## Calculate MASE
  
  measures<-measures.summary(n=200,data.list=data.list,method="gBridge",ITER=10,group=group,family="negbin")
  ```
  
  ```
  ## Calculate the Percentage of correct group selection for both the count and zero model
  
  betahat<-measures$betahat
  gammahat<-measures$gammahat
  
  grp.count<- grpresult.func(betahat=betahat,betatrue=beta,sim=1)
  
  grp.zero<- grpresult.func(betahat=gammahat,betatrue=gamma,sim=1)
  
  pgrp.corr.count<- grp.count$pgrp.correct
  pgrp.corr.zero<- grp.zero$pgrp.correct
  
  result<-(data.frame(measures$output,pgrp.corr.count,pgrp.corr.zero))
  names(result)<- c("MASE","corr_group_count","corr_group_zero")
  return(result)
  ```

 ```
 ## Output 
      MASE corr_group_count corr_group_zero
  0.99965                1             0.8
```
## Real Data Example

We illustrate our proposal by re-analyzing the auto insurance claim dataset from SAS Enterprise Miner database. The response variable of interest (y) is the aggregate claim loss of an auto insurance policy. Considering only policy records corresponding to the new customers the reduced dataset has 2,812 observations with 56 predictors being grouped into 21 groups of different group sizes. For the comparison of Gooogle methods with EMLasso we employ a repeated 5-fold cross validation (CV) procedure in which the dataset is randomly partitioned into 5 equal folds, iteratively taking each fold as the test set and the remaining set as the trainng set. We fit the models on the training sets and predictions are based on the test sets. We calculate average and median of the Mean Absolute Scaled Error (MASE) as the metric of evaluation, calculated over 100 iterations. 

## Function

The R file realdata.R contains two functions which are described below:

1) The function fit.out fits gBridge/EMLasso on the training dataset and obtains the regression coefficients for the count as well as the zero model corresponding to the minimum BIC. It takes the following arguments:

* train: Training data
* test: Test data
* yvar: The name of the outcome variable
* xvars: The name of the predictor variables for the count model
* zvars: The name of the predictor variables for the zero model
* group: The vector specifying grouping structure of the covariates
* penalty: "gBridge" or "EMLasso"

The function outputs the predicted value of the outcome variable and outcome variable from the test and the training set.

2) The function realdata.R takes the following arguments:

* data: The real auto insurance dataset 
* yvar, xvars, zvars: The name of the corresponding variables
* cv.iter: Number of iteration for calculating MASE using 5-fold CV mechanism
* k.fold: Number of CV folds (here it is 5)
* penalty: "gBridge" or "EMLasso"

This function forms the training and test data set iteratively in the 5-fold CV for every dataset and calls the function fit.out and using the output of the fit.out it calculates MASE using the function accuracy{forecast}.

The following codes can be used to obtain MASE of Gooogle and EM methods respecively.

```
result.gooogle<-realdata.func(data=data,yvar=yvar,xvars=xvars,zvars=zvars,cv.iter=5,k.fold=5,seedval=123,method="Gooogle")
result.em<-realdata.func(data=data,yvar=yvar,xvars=xvars,zvars=zvars,cv.iter=5,k.fold=5,seedval=123,method="EMLasso")
```




