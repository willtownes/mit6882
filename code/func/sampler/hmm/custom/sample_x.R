sample_x <-function(y,log_lik_func=NULL,z,m=NULL,x,p0=NULL,pk=NULL,theta,lambda,hyper=NULL){
  #convenience wrapper to match argument signature expected by HMM_HDP main function
  #function rLDS is under func/sampler/lds/lds_util.R
  #the provided x is only needed in order to pre-pend x0 to the top of the matrix
  return(rbind(x[1,],rLDS(1,y,z,theta,lambda)[[1]]))
}