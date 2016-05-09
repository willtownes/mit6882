### General-purpose useful functions
mse<-function(x,y){
  #compute mean squared error betwen object x and object y
  mean((x-y)^2)
}
trace<-function(x){
  #return the trace of matrix x
  sum(diag(x))
}
list_mean<-function(x,f=NULL){
  # input: x, a list where each element in the list is some numeric array
  # optional: f, a function, to be applied to each element of x
  # output: arithmetic mean computed across elements of x
  # ie, [f(x[[1]])+f(x[[2]])+....f(x[[n]])]/n
  if(is.null(f)){
    return(Reduce("+",x)/length(x))
  } else {
    return(Reduce("+",lapply(x,f))/length(x))
  }
}

### graphing ellipses
# library(car)
# mu<-c(1,1)
# Sigma<-matrix(c(1,.5,.5,1),nrow=2)
# ellipse(mu,Sigma,radius=1,add=FALSE,grid=FALSE,col="blue",xlim=c(-1,2),ylim=c(-1,2))
# mu<-c(0,0)
# Sigma<-matrix(c(1,-.5,-.5,1),nrow=2)
# ellipse(mu,Sigma,radius=1,add=TRUE,grid=FALSE,col="red")
