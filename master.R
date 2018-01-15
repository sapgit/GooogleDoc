########## Generate data and calculate the measures ##############

master.func<-function(n.train,n.test,rho,phi,method,sim,ITER=10,seedval=NULL)
{
  methods<-c("EMLasso","grLasso", "grMCP", "grSCAD","gBridge")
  if (sum(method == methods)==0) stop("The method argument can only take 'EMLasso', 'grLasso', 'grMCP', 'grSCAD' or 'gBridge'. ")
  
  # group size
  
  if (sim==1)
  {
    grpsize<-rep(8,5)
  } else {
    grpsize<-c(1,1,3,1,1,3,4,4,4,4,4)
  }
  
  # number of groups
  ngrp<-length(grpsize)

  # group formation
  group<-NULL
  for (k in 1:ngrp)
  {
    group<-c(group,rep(k,grpsize[k]))
  }

  data.list<-datagen.sim.all(n.train=n.train,n.test=n.test,rho=rho,phi=phi,grpsize=grpsize,sim=sim,ITER=ITER,seedval=seedval)

  measures<-measures.summary(n.train=n.train,n.test=n.test,data.list=data.list,method=method,ITER=ITER,group=group)
  names(measures)[3]<-"Time(sec)"

  return(measures)
}
