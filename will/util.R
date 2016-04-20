#utility function used by other scripts
library(matrixStats) #logSumExp functions
library(RSpectra) #faster eigenvalue decomposition

### Hidden Markov Model Functions
get_stationary_dist<-function(transition,logScale=TRUE){
  # given a transition matrix, compute the implied stationary distribution
  # stopifnot(all(rowSums(transition)==1))
  eigen(t(transition))$vectors[,1]
}

comp_bk_msg<-function(transition,elps,res=NULL){
  # transition is the L by L transition matrix for the hidden state markov chain
  # assume transition rows sum to one (standard from for stochastic matrix)
  # elps ("emission log probabilities") is a Tmax by L matrix...
  # ...where element (t,i) is log-probability...
  # ...of emitting the observed value x[t] conditional on the hidden state being "i"
  # optionally provide pre-initialized res matrix to save computation
  # RETURNS: a matrix of dimension Tmax by L of unnormalized backward messages on log-scale
  # add these to the output from comp_fwd_msg to get the (unnormalized) log probabilities 
  # of hidden states at each time point
  Tmax<-nrow(elps)
  L<-ncol(transition)
  # stopifnot(all(rowSums(transition)==1))
  tlps<-log(transition)
  # original transition[i,j] is P(z_{t+1}=j|z_t=i) (standard form)
  # in Murphy book this is called Psi
  ttlps<-t(tlps)
  #now ttlps[i,j] is log P(z_{t+1}=i|z_t=j)
  #taking transpose facilitates R's vector*matrix recycling
  if(is.null(res)) res<-matrix(0.0,nrow=Tmax,ncol=L)
  res[Tmax,]<-0 #messages on log-scale
  for(t in Tmax:2){
    #below formula modified from Murphy p. 611, underflowed/failed
    #res[(t-1),]<-transition %*% (exp(elps[t,])*res[t,])
    # more stable? log-scale version, also underflowed
    # res[(t-1),]<-transition %*% exp(elps[t,]+log(res[t,]))
    # even more stable approach, use logsumexp function
    res[(t-1),]<-colLogSumExps(res[t,]+elps[t,]+ttlps)
    #note the use of vector recycling!!
    # res[t,] is 2x1, elps[t,] is 2x1, tlps is 2x2!!
  }
  return(res)
}

comp_fwd_msg<-function(transition,elps,res=NULL){
  # transition is the L by L transition matrix for the hidden state markov chain
  # assume transition rows sum to one (standard from for stochastic matrix)
  # elps ("emission log probabilities") is a Tmax by L matrix...
  # ...where element (t,i) is log-probability...
  # ...of emitting the observed value x[t] conditional on the hidden state being "i"
  # optionally provide pre-initialized res matrix to save computation
  # RETURNS: a matrix of dimension Tmax by L of unnormalized forward messages on log-scale
  # add these to the output from comp_fwd_msg to get the (unnormalized) log probabilities 
  # of hidden states at each time point
  Tmax<-nrow(elps)
  L<-ncol(transition)
  # stopifnot(all(rowSums(transition)==1))
  tlps<-log(transition)
  if(is.null(res)) res<-matrix(0.0,nrow=Tmax,ncol=L)
  #compute stationary distribution as marginal log probs to initialize
  mlpz<-log(get_stationary_dist(transition)) 
  res[1,]<-elps[1,]+mlpz
  for(t in 2:Tmax){
    res[t,]<-elps[t,]+colLogSumExps(res[(t-1),]+tlps) #vector recycling!
  }
  return(res)
}

forward_backward<-function(transition,elps,logscale=TRUE,normalize=FALSE,res=NULL){
  # transition is the L by L transition matrix for the hidden state markov chain
  # assume transition rows sum to one (standard from for stochastic matrix)
  # elps ("emission log probabilities") is a Tmax by L matrix...
  # ...where element (t,i) is log-probability...
  # ...of emitting the observed value x[t] conditional on the hidden state being "i"
  # optionally provide pre-initialized res matrix to save computation
  # RETURNS: a matrix of dimension Tmax by L of "gammas", i.e. marginal probs of hidden states
  # res[t,i] = Pr(z_t=i | x_{1:T}) where z is hidden state and x is observed state
  # if logscale=TRUE, values are returned on the log-scale
  # if normalize=TRUE, the values are probabilities, otherwise they are just "messages"
  # note that normalization is not necessary to determine the most likely state at each time
  Tmax<-nrow(elps)
  L<-ncol(transition)
  # stopifnot(all(rowSums(transition)==1))
  msg_b<-comp_bk_msg(transition,elps,res=res)
  msg_a<-comp_fwd_msg(transition,elps,res=res)
  msg_gam<-msg_a+msg_b #unnormalized log-probs
  if(normalize){
    msg_gam_Z<-rowLogSumExps(msg_gam) #calculate normalizing constant
    msg_gam<- msg_gam - msg_gam_Z #normalize via log-scale
  }
  if(logscale){
    return(msg_gam)
  } else {
    return(exp(msg_gam))
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
