setwd("./Results/")

rhos<-c(0,0.4,0.8)
phis<-c(0.3,0.4,0.5)
ns<-c(200,500,1000)
nums<-1:10

# list of all methods compared
methods<-c("EMLasso","grLasso", "grMCP", "grSCAD", "gBridge")

mase.summary.all<-NULL
mae.summary.all<-NULL

for (rho in rhos)
{
  for (phi in phis)
  {
      for (n in ns)
      {
        measures.method.all<-NULL
        for (method in methods)
        {
          m<-which(methods==method)
          measures.mat<-NULL
          
          files<-list.files(pattern=paste("Sim2_",rho,"_",phi,"_",n,sep=""))
          for (file in files)
          {
            load(file)
            result.method<-result[[m]]
            row.names(result.method)<-NULL

            colnames(result.method)[3]<-"time.taken"
            measures.mat<-rbind.data.frame(measures.mat,result.method)
          }
          measures.method<-t(apply(apply(measures.mat,2,function(x) return(as.numeric(x))),2,median,na.rm=T))
          measures.method<-cbind.data.frame(rho=rho,phi=phi,n=n,method=method,measures.method)#*#
          measures.method.all<-rbind.data.frame(measures.method.all,measures.method)

          print(paste("rho=",rho,"; phi=",phi,"; n=",n,"; method=",method," is complete",sep=""))
        }
        
        mase.summary<-c(rho=rho,phi=phi,n=n, round(measures.method.all[1:5,"MASE"],3))
        mase.summary.all<-rbind(mase.summary.all,mase.summary)
        
        mae.summary<-c(rho=rho, phi=phi,n=n,round(measures.method.all[1:5,"MAE"],3))
        mae.summary.all<-rbind(mae.summary.all,mae.summary)

        write.table(measures.method.all,file=paste("Sim_1_",rho,"_",dist,"_",phi,"_",n,"_measures.summary.csv",sep=""),sep=",",row.names = F)
        
      }
    }
  }

colnames(mase.summary.all)[-(1:3)]<-methods[1:5]
write.table(mase.summary.all,file="mase.summary.csv",sep=",",row.names=F)

colnames(mae.summary.all)[-(1:3)]<-methods[1:5]
write.table(mae.summary.all,file="mae.summary.csv",sep=",",row.names=F)