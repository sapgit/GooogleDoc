#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
ITER<-as.numeric(args[2])
NUM<-as.numeric(args[3])

# Required Libraries
require(forecast); require(MASS); require(pscl); require(grpreg); require(mpath); require(zic); require(SGL); require(stringr); require(gamlss.dist); require(Gooogle)

source('./Codes/docvisit_CV.R')
data<-read.table("./Data/docvisits_spline.txt",header=T,sep=",")

## load cv index data
load("./Codes/CV_result_doc.RData")
final<-func_main(NUM=NUM,ITER=ITER)  
save(final, file=paste('./Results/docvisit_CV_', NUM, '.RData', sep=''))


