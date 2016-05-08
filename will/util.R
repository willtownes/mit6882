### Utility function used by other scripts
library(matrixStats) #logSumExp functions
library(mvtnorm) #multivariate normal sampling
library(MCMCpack) #inverse wishart sampling
if(require(RSpectra)){
  FASTEIG=TRUE
} else {
  FASTEIG=FALSE
}

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

### Multivariate and Matrix Normal Sampling
# see matnorm.Rmd for examples/ visual tests
# see also speed_tests.Rmd
rmvnorm_info<-function(n,theta,Lambda,D=length(theta)){
  # sample n variates from multivariate normal in "information form"
  # theta="offset" and Lambda="information matrix"
  # Lambda must be positive definite
  # standard params, mu=Lambda^{-1}*theta, Sigma=Lambda^{-1}
  Q<-chol(Lambda)
  Z<-matrix(rnorm(n*D),nrow=D,ncol=n)
  #using vector recycling tricks below
  #Q^{-1}*Z+Lambda^{-1}*theta == Q^{-1}*(Z+(Q')^{-1}*theta)
  backsolve(Q,Z+drop(backsolve(Q,theta,transpose=TRUE)))
}

rmatnorm<-function(n,M,V,K,foxpar=TRUE){
  # return n samples from matrix normal(M,V,K) distribution
  # M = mean matrix
  # V = among-row covariance matrix
  # if foxpar==TRUE
  # K = inverse of among-column covariance matrix
  # if foxpar==FALSE
  # K= among-column covariance matrix
  m<-nrow(M)
  p<-ncol(M)
  A<-t(chol(V))
  R<-chol(K)
  fn1<-function(){A%*%matrix(rnorm(m*p),nrow=m,ncol=p)}
  AZs<-replicate(n,fn1(),simplify=FALSE)
  if(foxpar){
    # case of K is inverse covariance matrix
    fn2<-function(AZ){t(backsolve(R,t(AZ)))+M}
  } else {
    # case of K is covariance matrix
    B<-R
    fn2<-function(AZ){AZ%*%B+M}
  }
  lapply(AZs,fn2)
}

bayes_mlinreg_post<-function(Y,X,hyper=list(),nBurn=100,nSample=100){
  #Bayesian Multivariate Linear Regression
  #Given hyperparameters and data, run Gibbs Sampler to obtain samples from the posterior of 
  #input: Y= (d by n) outcome matrix
  #X= (d by n) input matrix
  #model: Y=AX+B+E with E~MN(0,Sigma,I)
  #hyper is a list of hyperparameters. Elements:
  #nu0,Delta0=hyperparameters of Sigma~Inv.Wishart()
  #M_A,K_A=hyperparameters of transition matrix A~MatNorm()
  #M_B,kappa0=hyperparameters of additive matrix B~Norm()
  #Any hyperparameters not provided will be set to following defaults "empirical bayes":
  #nu0=nrow(Y)+2 (most diffuse value while still proper)
  #Delta0=0.75*empirical covariance of Y rows
  #kappa0=0 uninformative
  #M_A= (d by d) zero matrix 
  #M_B= zero vector of length d
  #K_A= (d by d) zero matrix (uninformative prior)
  #output: list of length nSample
  #each element has three components, "A","B","Sigma"
  #corresponding to samples from posterior of these parameters.
  d<-nrow(Y)
  n<-ncol(Y)
  M_A<-hyper[["M_A"]]
  M_B<-hyper[["M_B"]]
  K_A<-hyper[["K_A"]]
  kappa0<-hyper[["kappa0"]]
  nu0<-hyper[["nu0"]]
  Delta0<-hyper[["Delta0"]]
  #set defaults if not specified
  if(is.null(M_A)) M_A<-matrix(0,nrow=d,ncol=d)
  if(is.null(M_B)) M_B<-rep(0,d)
  if(is.null(K_A)) K_A<-matrix(0,nrow=d,ncol=d)
  if(is.null(kappa0)) kappa0<-.01
  if(is.null(nu0)) nu0<-d+2
  if(is.null(Delta0)) Delta0<-.75*cov(t(Y))
  #sufficient statistics
  Ymn<-rowMeans(Y)
  Xmn<-rowMeans(X)
  Saxx<-X%*%t(X)+K_A
  YXt<-Y%*%t(X)
  MaKa<-M_A%*%K_A
  kappan<-n+kappa0
  nun<-nu0+n
  wt_b<-n/kappan
  J<-rep(1,n)
  #initialize data containers
  #resS<-resA<-matrix(0,nrow=d,ncol=d)
  #resB<-rep(0,d)
  res<-replicate(nSample,list(),simplify=FALSE)
  A<-M_A
  B<-rowMeans(Y)
  Sigma<-Delta0 #/(nu0-d-1)
  for(t in (1-nBurn):nSample){
    mu_b<-wt_b*(Ymn-A%*%Xmn)+(1-wt_b)*M_B
    B<-mvrnorm(1,mu_b,Sigma/kappan)
    Sayx<-YXt-B%*%t(n*Xmn)+MaKa
    A<-rmatnorm(1,t(solve(Saxx,t(Sayx))),Sigma,Saxx,foxpar=TRUE)[[1]]
    resids<-Y-A%*%X-B%o%J
    Sigma<-riwish(nun,Delta0+resids%*%t(resids))
    if(t>0){
      res[[t]][["A"]]<-A
      res[[t]][["B"]]<-B
      res[[t]][["Sigma"]]<-Sigma
    }
  }
  return(res)
}

### Hidden Markov Model Functions
# refer to HMM.Rmd for examples/ visual tests
get_stationary_dist<-function(transition){
  # given a transition matrix, compute the implied stationary distribution
  # stopifnot(all(rowSums(transition)==1))
  if(FASTEIG && nrow(transition)>2){
    evec<-eigs(t(transition),1) #fast method using RSpectra
  } else {
    evec<-eigen(t(transition))$vectors[,1] #slow method using base
  }
  evec/sum(evec) #normalize
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

hidden_int2char<-function(hidden_ints,categs){
  # given a sequence of integers representing a hidden state path,
  # and a character vector "categs" with hidden state labels
  # returns the same hidden state path in character form for easier interpretation
  vapply(hidden_ints,function(q){categs[q]},FUN.VALUE="S")
}

sample_hidden_states<-function(msg_b,elps,transition,pi0=NULL,categs=NULL){
  # inputs: msg_b is a Tmax by L matrix of backward messages on log scale
  # Tmax is number of time steps, L is number of hidden states
  # elps is a Tmax by L matrix of emission log-probabilities
  # transition is the transition matrix for the hidden state markov chain
  # pi0 is a L-vector of marginal hidden state log-probabilities (optional)...
  # ...if not provided it is calculated from transition matrix stationary distribution
  # if categs provided, converts sample path from integer to category label
  Tmax<-nrow(elps)
  L<-ncol(elps)
  tlps<-log(transition)
  if(is.null(pi0)) pi0<-log(get_stationary_dist(transition))
  zs<-rep(0,Tmax) #will use numeric index
  #compute sampled hidden state sequence
  msg_a<-elps[1,]+pi0+msg_b[1,]
  msg_a<-msg_a - logSumExp(msg_a) #normalized, still on log scale
  zs[1]<-sample.int(L,1,prob=exp(msg_a))
  for(t in 2:Tmax){
    msg_a<- tlps[zs[t-1],] + elps[t,] + msg_b[t,]
    msg_a<-msg_a - logSumExp(msg_a)
    zs[t]<-sample.int(L,1,prob=exp(msg_a))
  }
  if(is.null(categs)){
    return(zs)
  } else {
    return(hidden_int2char(zs,categs))
  }
}

### Linear Dynamical System Functions
# test code/ examples are in LDS.Rmd
gen_lds<-function(z,x0,theta,lambda){
  # create simulated data for linear 
  # z is sequence of hidden mode indicators (integers)
  # x1 is initial position (a vector)
  # Given dynamical parameters eg theta[[1]]=list(A,B,Sigma)
  # theta[z[t]] tells dynamical regime at time t
  # emissions noise is in lambda. lambda[["R"]] is the matrix of interest
  # observation matrix is in lambda. lambda[["C"]].
  # returns a data matrix x with latent trajectory
  # also returns data matrix y with observed trajectory
  nIter<-length(z)
  xdim<-length(x0)
  R<-lambda[["R"]]
  C<-lambda[["C"]]
  ydim<-nrow(C)
  x<-matrix(NA,nrow=nIter,ncol=xdim) #eg, position and velocity
  y<-matrix(NA,nrow=nIter,ncol=ydim) #eg, observed position
  for(t in 1:nIter){
    theta_t<-theta[[z[t]]]
    A<-theta_t[["A"]]
    B<-theta_t[["B"]]
    Sigma<-theta_t[["Sigma"]]
    if(t==1){
      x_prev<-x0
    } else {
      x_prev<-x[(t-1),]
    }
    x[t,]<-rmvnorm(1,drop(A%*%x_prev+B),Sigma)
    y[t,]<-rmvnorm(1,drop(C%*%x[t,]),R)
  }
  return(list("x"=x,"y"=y))
}

backward_kalman_msgs<-function(y,z,pars,C,R){
  # inputs:
  # y is a matrix where each row is an observation
  # z is a vector of mode indicators. For time-invarying dynamics, set z=rep(1,ncol(y))
  # pars is a list, one element for each unique mode (==length(unique(z)))
  # pars[[i]] is a list containing time-varying parameters
  # "A" is transition matrix for linear dynamical system (LDS)
  # "B" is additive term for transition in LDS
  # "Sigma" is noise for transition in LDS
  # C and R are (time-invariant) parameters for observation model
  Tmax<-nrow(y)
  #D<-ncol(y)
  eye<-diag(ncol(A))
  #note if u=chol(R) then R=u'u, NOT uu'.
  #R^{-1} = u^{-1}u^{-T}
  #implementing Algorithm 19 from Emily Fox MIT dissertation
  #U<-chol(R)
  #UtiC<-solve(t(U),C)
  #CtRiC<-t(UtiC)%*%UtiC
  CtRiC<-t(C)%*%solve(R,C)
  #initialize messages
  Lambda_tt<-lapply(1:Tmax,function(x){CtRiC})
  #theta_tt<-t(apply(y,1,function(yt){t(UtiC)%*%solve(t(U),yt)}))
  theta_tt<-t(apply(y,1,function(yt){t(C)%*%solve(R,yt)}))
  for(t in Tmax:1){
    Lambda<-Lambda_tt[[t]]
    par_t<-pars[[z[t]]]
    A<-par_t[["A"]]
    B<-par_t[["B"]]
    Sigma<-par_t[["Sigma"]]
    Si<-solve(Sigma)
    Jt<-t(solve(Lambda+Si,Lambda))
    Lt<-eye - Jt
    Lam<-t(A)%*%(Lt%*%Lambda%*%t(Lt)+Jt%*%solve(Sigma,t(Jt)))%*%A
    theta<-t(A)%*%Lt%*%(theta_tt[t,] - Lambda%*%B)
    if(t>1){
      Lambda_tt[[t-1]]<-Lam+CtRiC #or Lam+Lambda_tt[[t-1]]
      theta_tt[t-1,]<-theta+theta_tt[t-1,]
    } else {
      Lambda_00<-Lam
      theta_00<-theta
    }
  }
  res<-list()
  res[["info_matrix_msgs"]]<-Lambda_tt #latent state dim x latent state dim
  res[["offset_msgs"]]<-theta_tt #Tmax by dimension of latent state
  res[["info_matrix_init"]]<-Lambda_00
  res[["offset_msg_init"]]<-theta_00
  return(res)
}

fwd_kalman_sample<-function(z,pars,Lambda_00,Lambda_tt,theta_00,theta_tt,xs=NULL){
  # z,pars are as described in backward_kalman_msgs
  # remaining parameters are form output of backward_kalman_msgs:
  # Lambda_tt is list of all the information matrix messages from backward Kalman Filter
  # Lambda_00 is initial information matrix
  # theta_00 is initial offset parameter vector
  # theta_tt is matrix where row t is offset parameter for step t in the time series
  # can save time by pre-allocating result matrix xs
  Tmax<-length(z)
  if(is.null(xs)){
    D<-ncol(pars[[1]][["Sigma"]])
    xs<-matrix(NA,nrow=Tmax,ncol=D)
  }
  xs0<-rmvnorm_info(1,theta_00,Lambda_00) #initialize x0, not sure if right
  for(t in 1:Tmax){
    par_t<-pars[[z[t]]]
    A<-par_t[["A"]]
    B<-par_t[["B"]]
    Sigma<-par_t[["Sigma"]]
    info_mat<-solve(Sigma)+Lambda_tt[[t]]
    if(t==1){
      offset<-solve(Sigma,A%*%xs0+B)+theta_tt[t,]
    } else {
      offset<-solve(Sigma,A%*%xs[t-1,]+B)+theta_tt[t,]
    }
    xs[t,]<-rmvnorm_info(1,offset,info_mat)
  }
  return(xs)
}
rLDS<-function(n,y,z,theta,lambda){
  # convenience wrapper for combining backward_kalman_msgs with fwd_kalman_sample
  # n is number of desired samples
  # each sample is a matrix
  # returns a list of length n
  # special case n=1, returns a single matrix
  # dim(matrix) is nrow=number of steps in time series, ncol=dim(latent state)
  C<-lambda[["C"]]
  R<-lambda[["R"]]
  msg_b<-backward_kalman_msgs(y,z,theta,C,R)
  L_tt<-msg_b[["info_matrix_msgs"]] #latent state dim x latent state dim
  th_tt<-msg_b[["offset_msgs"]] #Tmax by dimension of latent state
  L_00<-msg_b[["info_matrix_init"]]
  th_00<-msg_b[["offset_msg_init"]]
  res<-replicate(n,fwd_kalman_sample(z,theta,L_00,L_tt,th_00,th_tt),simplify=FALSE)
  return(res)
}
rLDS_melt<-function(LDS_samples,cnames=NULL){
  # takes output of rLDS and converts into a giant data frame
  # facilitates making plots with ggplot()
  #for(n in 1:length(LDS_samples)){
  #  LDS_samples[[n]]<-cbind(LDS_samples[[n]],"id"=n)
  #}
  Tmax<-nrow(LDS_samples[[1]])
  D<-ncol(LDS_samples[[1]])
  n<-length(LDS_samples)
  res<-do.call("rbind",LDS_samples)
  res<-cbind(res,"id"=rep(1:n,each=Tmax))
  res<-as.data.frame(res)
  if(is.null(cnames)){
    return(res)
  } else {
    colnames(res)<-c(cnames,"id")
    return(res)
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
