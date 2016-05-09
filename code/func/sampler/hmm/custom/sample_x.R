sample_x <-function(y,log_lik_func=NULL,z,m=NULL,x=NULL,p0=NULL,pk=NULL,theta,lambda,hyper=NULL){
  #convenience wrapper to match argument signature expected by HMM_HDP main function
  #function rLDS is under func/sampler/lds/lds_util.R
  rLDS(1,y,z,theta,lambda)[[1]]
}