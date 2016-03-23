#utility function used by other scripts
library(car)

mu<-c(1,1)
Sigma<-matrix(c(1,.5,.5,1),nrow=2)
ellipse(mu,Sigma,radius=1,add=FALSE)

mu<-c(1,1)
Sigma<-matrix(c(1,-.5,-.5,1),nrow=2)
ellipse(mu,Sigma,radius=1,add=FALSE)
