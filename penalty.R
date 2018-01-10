# function for each method used in simulation
func_EMLasso<-function(formula,data,dist){fit<-zipath(formula=formula,data=data,family=dist); idx<-which.min(fit$bic); coeff.final=list(count=fit$coefficients$count[,idx],zero=fit$coefficients$zero[,idx]); return(list(coefficients=coeff.final,aic=fit$aic[idx], bic=fit$bic[idx], loglik=fit$loglik[idx]))}

func_grLasso<-function(data,yvar,xvars,zvars,group,dist){return(gooogle(data=data,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist=dist,penalty="grLasso"))}

func_grMCP<-function(data,yvar,xvars,zvars,group,dist){return(gooogle(data=data,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist=dist,penalty="grMCP"))}

func_grSCAD<-function(data,yvar,xvars,zvars,group,dist){return(gooogle(data=data,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist=dist,penalty="grSCAD"))}

func_gBridge<-function(data,yvar,xvars,zvars,group,dist){return(gooogle(data=data,yvar=yvar,xvars=xvars,zvars=zvars,group=group,dist=dist,penalty="gBridge"))}
