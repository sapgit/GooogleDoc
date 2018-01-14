#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
rho<-as.numeric(args[1])
phi<-as.numeric(args[3])
n<-as.numeric(args[4])
ITER<-as.numeric(args[5])
NUM<-as.numeric(args[6])

# Required Libraries
require(forecast); require(MASS); require(pscl); require(grpreg); 
require(mpath); require(zic); require(SGL); require(stringr); 
require(gamlss.dist); require(dummies); require(Gooogle)

source('./Codes/penalty_func.R')
source('./Codes/Park_Sim2_V1.R')

result<-func_main(rho=rho, phi=phi, n=n, ITER=ITER)  
save(result, file=paste('./Results/Sim2_', rho, '_', phi, '_', n, '_', NUM, '.RData', sep=''))
