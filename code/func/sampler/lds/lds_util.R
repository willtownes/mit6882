### Functions used by Linear Dynamical System models
library(MASS) #multivariate normal sampling
library(MCMCpack) #inverse wishart sampling

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

bayes_mlinreg_post<-function(Y,X,hyper,nBurn=100,nSample=100){
  #Bayesian Multivariate Linear Regression
  #Given hyperparameters and data, run Gibbs Sampler to obtain samples from the posterior of 
  #input: Y= (d by n) outcome matrix
  #X= (d by n) input matrix
  #model: Y=AX+B+E with E~MN(0,Sigma,I)
  #hyper is a list of hyperparameters. Elements:
  #S0_df,S0=hyperparameters of Sigma~Inv.Wishart()
  #M_A,K_A=hyperparameters of transition matrix A~MatNorm()
  #M_B,kappa0=hyperparameters of additive matrix B~Norm()
  #Any hyperparameters not provided will be set to following defaults "empirical bayes":
  #S0_df=nrow(Y)+2 (most diffuse value while still proper)
  #S0=0.75*empirical covariance of Y rows
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
  S0_df<-hyper[["S0_df"]]
  S0<-hyper[["S0"]]
  #set defaults if not specified
  # if(is.null(M_A)) M_A<-matrix(0,nrow=d,ncol=d)
  # if(is.null(M_B)) M_B<-rep(0,d)
  # if(is.null(K_A)) K_A<-matrix(0,nrow=d,ncol=d)
  # if(is.null(kappa0)) kappa0<-.01
  # if(is.null(S0_df)) S0_df<-d+2
  # if(is.null(S0)) S0<-.75*cov(t(Y))
  #sufficient statistics
  Ymn<-rowMeans(Y)
  Xmn<-rowMeans(X)
  Saxx<-X%*%t(X)+K_A
  YXt<-Y%*%t(X)
  MaKa<-M_A%*%K_A
  kappan<-n+kappa0
  nun<-S0_df+n
  wt_b<-n/kappan
  J<-rep(1,n)
  #initialize data containers
  #resS<-resA<-matrix(0,nrow=d,ncol=d)
  #resB<-rep(0,d)
  res<-replicate(nSample,list(),simplify=FALSE)
  A<-M_A
  B<-rowMeans(Y)
  Sigma<-S0 #/(S0_df-d-1)
  for(t in (1-nBurn):nSample){
    mu_b<-wt_b*(Ymn-A%*%Xmn)+(1-wt_b)*M_B
    B<-mvrnorm(1,mu_b,Sigma/kappan)
    Sayx<-YXt-B%*%t(n*Xmn)+MaKa
    A<-rmatnorm(1,t(solve(Saxx,t(Sayx))),Sigma,Saxx,foxpar=TRUE)[[1]]
    resids<-Y-A%*%X-B%o%J
    Sigma<-riwish(nun,S0+resids%*%t(resids))
    if(t>0){
      res[[t]][["A"]]<-A
      res[[t]][["B"]]<-B
      res[[t]][["Sigma"]]<-Sigma
    }
  }
  return(res)
}

### Linear Dynamical System Functions
# test code/ examples are in LDS.Rmd
gen_lds<-function(z,x0, theta, lambda){
  # create simulated data for linear dynamical system
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
    x[t,]<-mvrnorm(1,drop(A%*%x_prev+B),Sigma)
    y[t,]<-mvrnorm(1,drop(C%*%x[t,]),R)
  }
  return(list("x"=x,"y"=y))
}

backward_kalman_msgs<-function(y, z, pars, C, R){
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
  eye<-diag(ncol(C))
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

fwd_kalman_sample<-function(z,pars,Lambda_00,Lambda_tt,theta_00,theta_tt,xs=NULL,x0_return=FALSE){
  # z,pars are as described in backward_kalman_msgs
  # remaining parameters are from output of backward_kalman_msgs:
  # Lambda_tt is list of all the information matrix messages from backward Kalman Filter
  # Lambda_00 is initial information matrix
  # theta_00 is initial offset parameter vector
  # theta_tt is matrix where row t is offset parameter for step t in the time series
  # can save time by pre-allocating result matrix xs
  # if x0_return==TRUE, returns the sample of the "0^th" state in row 1, otherwise only rows 1:T
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
  if(x0_return==TRUE) xs<-rbind(t(xs0),xs)
  return(xs)
}
rLDS<-function(n,y,z,theta,lambda,x0_return=FALSE){
  # convenience wrapper for combining backward_kalman_msgs with fwd_kalman_sample
  # n is number of desired samples
  # each sample is a matrix
  # returns a list of length n
  # dim(matrix) is nrow=number of steps in time series, ncol=dim(latent state)
  # if x0_return==TRUE, returns the sample of the "0^th" state in row 1, otherwise only rows 1:T
  C<-lambda[["C"]]
  R<-lambda[["R"]]
  msg_b<-backward_kalman_msgs(y,z,theta,C,R)
  L_tt<-msg_b[["info_matrix_msgs"]] #latent state dim x latent state dim
  th_tt<-msg_b[["offset_msgs"]] #Tmax by dimension of latent state
  L_00<-msg_b[["info_matrix_init"]]
  th_00<-msg_b[["offset_msg_init"]]
  res<-replicate(n,fwd_kalman_sample(z,theta,L_00,L_tt,th_00,th_tt,x0_return=x0_return),simplify=FALSE)
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

### Functions for initializing Dynamical Parameters Priors from Data

hyper_init<-function(y,D=ncol(y),hyper=list()){
  #given a matrix of observations "y" (each row is an observation vector) nrow=number of time points
  #specify the desired dimension of the latent state x
  #returns a list of hyperparameters for the matrix normal inverse wishart prior over dynamical parameters
  # the hyperparameters are:
  # S0_df = prior degrees of freedom on Sigma, the transition model covariance/noise
  # S0 = prior scale matrix for Sigma
  # R0_df = prior degrees of freedom on R, the emission model covariance/noise (observation noise)
  # R0 = prior scale matrix for R
  # M_B = (vector) prior mean for additive term B in transition model
  # kappa0 = prior sample size for B
  # M_A = (matrix) prior mean for the transition matrix A of the latent states
  # K_A = inverse of "column covariance" matrix in prior for A
  # Any matching params in the input list "hyper" are returned as-is
  # otherwise, the params are set based on statistics of the observed data y
  # lowess_f only needed for S0 or R0 estimation. Small value=wigglier estimated latent x values
  # large value= more linear latent x values path
  ydim<-ncol(y)
  stopifnot(D>=ydim)
  if(is.null(hyper[["S0_df"]])) hyper[["S0_df"]]<-D+2
  if(is.null(hyper[["R0_df"]])) hyper[["R0_df"]]<-ydim+2
  if(is.null(hyper[["M_B"]])) hyper[["M_B"]]<-rep(0,D)
  if(is.null(hyper[["kappa0"]])) hyper[["kappa0"]]<-.01
  if(is.null(hyper[["M_A"]])) hyper[["M_A"]]<-diag(D)
  if(is.null(hyper[["K_A"]])) hyper[["K_A"]]<-.01*diag(D)
  if(!is.null(hyper[["S0"]]) && !is.null(hyper[["R0"]])) return(hyper)
  S0_df<-hyper[["S0_df"]]
  R0_df<-hyper[["R0_df"]]
  #case where either R0 or S0 not specified, estimate from covariance of data:
  Tmax<-nrow(y)
  time<-1:Tmax
  lowess_f<-10/Tmax #fraction of data points used to estimate x from y, essentially 20 nearest neighbors
  lws_func<-function(u,f){
    # infer latent x column for one of y's columns
    #print(str(u))
    lowess(time,u,f=f)$y
  }
  xhat<-apply(y,2,function(u){lws_func(u,lowess_f)})
  if(is.null(hyper[["R0"]])) hyper[["R0"]]<-cov(y-xhat) #*(R0_df-ydim-1)
  #R0<-hyper[["R0"]]
  xtrend<-apply(xhat,2,function(u){lws_func(u,20/Tmax)})
  if(is.null(hyper[["S0"]])){
    S0<-matrix(0,nrow=D,ncol=D)
    S0_topleft<-cov(xhat-xtrend)
    #S0_topleft<-cov(xhat[2:Tmax,]-xhat[1:(Tmax-1),]) 
    #S0_topleft<-cov(xhat)
    S0[1:ydim,1:ydim]<-S0_topleft
    if(D>ydim){
      S0[(ydim+1):D,(ydim+1):D]<-diag(D-ydim)*det(S0_topleft)/(D-ydim) #formula from E. Fox thesis p. 160
    }
    hyper[["S0"]]<- S0 #*(S0_df-D-1)
  }
  return(hyper)
}
# sample a set of parameters theta from the prior defined by the passed-in hyperparameters
theta_init<-function(hp){
  # hp is a list of hyperparameters (see init_mniw_hyper)
  # this is a sample from the prior for initialization, use bayes_mlinreg_post to sample from full conditional
  Sigma<-riwish(hp$S0_df,hp$S0)
  B<-mvrnorm(1,hp$M_B,Sigma/hp$kappa0)
  A<-rmatnorm(1,hp$M_A,Sigma,hp$K_A,foxpar=TRUE)[[1]]
  return(list(Sigma=Sigma,A=A,B=B))
}
lambda_init<-function(hp){
  # C matrix regarded as fixed and known
  # R matrix random but constant across all dynamical modes
  R<-riwish(hp$R0_df,hp$R0)
  ydim<-ncol(R)
  D<-ncol(hp$S0)
  C<-matrix(0,nrow=ydim,ncol=D)
  C[1:ydim,1:ydim]<-diag(ydim)
  return(list(R=R,C=C))
}
x0_init<-function(x,time=1:10){
  # given other xvalues, initialize x0 by linear extrapolation
  # x is a matrix with each row = an observation vector
  # time is the index vector of observations used to extrapolate to time zero
  x<-x[time,]
  coef(lm(x~time))[1,]
}