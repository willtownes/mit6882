#utility function used by other scripts
library(car)

mu<-c(1,1)
Sigma<-matrix(c(1,.5,.5,1),nrow=2)
ellipse(mu,Sigma,radius=1,add=FALSE,grid=FALSE,col="blue",xlim=c(-1,2),ylim=c(-1,2))

mu<-c(0,0)
Sigma<-matrix(c(1,-.5,-.5,1),nrow=2)
ellipse(mu,Sigma,radius=1,add=TRUE,grid=FALSE,col="red")
