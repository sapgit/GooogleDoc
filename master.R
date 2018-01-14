########## Generate data and calculate the measures ##############

master.func<-function(n.train,n.test,beta,gamma,rho,grpsize,phi,family,method,sim,ITER=10,seedval=NULL)
{
  # number of groups
  ngrp<-length(grpsize)

  # group
  group<-NULL
  for (k in 1:ngrp)
  {
    group<-c(group,rep(k,grpsize[k]))
  }

  data.list<-datagen.sim.all(n.train=n.train,n.test=n.test,rho=rho,phi=phi,beta=beta,gamma=gamma,grpsize=grpsize,family=family,sim=sim,ITER=ITER,seedval=seedval)

  measures<-measures.summary(n.train=n.train,n.test=n.test,data.list=data.list,method=method,ITER=ITER,group=group,family=family)
  return(measures)
}
