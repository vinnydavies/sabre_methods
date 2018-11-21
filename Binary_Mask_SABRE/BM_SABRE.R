library("mvtnorm")
library("Rlab")
library("coda")
library("msm")
library("pscl")
library("lme4")
library("dummies")
library("doBy")

source("./Binary_Mask_SABRE/BM_SABRE_gamma_density.R")

is.invert <- function(m,XTX,gamma) class(try(solve(solve(m) + XTX*(gamma%*%t(gamma))),silent=T))=="matrix"

BM_SABRE<-function(y,X,Z,iters,it_count,
                    gamma0=c(1,rbern((ncol(X)-1),pi_a/(pi_a + pi_b))),
                    w0=rnorm(ncol(X),c(max(y),rep(min(0,((mean(y)-max(y))/(sum(gamma0)-1))),ncol(X)-1)),sd(y)/sqrt(sum(gamma0))),
                    sigsq_e0=runif(1,0,var(y)),
                    sigsq_w0=c(100,runif((length(sigsq_w_change)),0,var(y))),
                    b0=rnorm(ncol(Z),(((y%*%Z)/colSums(Z))-mean((y%*%Z)/colSums(Z))),sd((((y%*%Z)/colSums(Z))-mean((y%*%Z)/colSums(Z))))),
                    sigsq_b0=runif(c,0,var(b0)),
                    sigsq_w_change=c(ncol(X)),
                    mu_w0=c(max(y),rtnorm(length(sigsq_w_change),rep((mean(y)-max(y))/length(sigsq_w_change),length(sigsq_w_change)),sd(y)/sqrt(length(sigsq_w_change)))),
                    var_mu_w=rep(100,length(sigsq_w_change)),
                    alpha_e=0.001,
                    beta_e=0.001,
                    alpha_w=rep(0.001,length(sigsq_w_change)+1),
                    beta_w=rep(0.001,length(sigsq_w_change)+1),
                    alpha_b=rep(0.001,numcolZ),
                    beta_b=rep(0.001,numcolZ),
                    pi0=rbeta(1,pi_a+sum(gamma0)-1,pi_b + length(gamma0)-sum(gamma0)),
                    pi_a=1,
                    pi_b=4,
                    pi_prop=0.3,
                    f=round((ncol(X)-1)/10),
                    mean_sig_prop=0.05,
                    prop_var_sig=0.00001){
  time.check="TRUE"
  method="log-MH"
  pi.vary=FALSE
  if(time.check=="TRUE"){
    t0<-Sys.time()
    t<-0
  }
  if(sum(sigsq_w_change)!=ncol(X)){
    stop("Not all variables are assigned a sigsq_w")
  } 
  if(sum(is.na(pi.vary)==c("TRUE","FALSE"))!=1){
    stop("pi.vary is neither True or False (requires double speak marks)")
  }
  
  PP<-c()
  
  X<-cbind(rep(1,length(y)),X)
  XTX<-t(X)%*%X
  y<-c(y)
  NY<-length(y)
  
  gamma<-matrix(rep(0,ncol(X)*iters),ncol=ncol(X))
  gamma[1,]<-gamma0
  gamma_curr<-gamma0
  
  NE_prop<-NA
  NE_curr<-sum(gamma0)!=1
  
  wf<-matrix(rep(0,ncol(X)*iters),ncol=ncol(X))
  wf[1,]<-w0
  
  sigsq_e<-c()
  sigsq_e[1]<-sigsq_e0
  
  sigsq_w<-matrix(rep(0,iters*(length(sigsq_w_change)+1)),ncol=(length(sigsq_w_change)+1))
  sigsq_w[1,]<-sigsq_w0
  change_start<-c(1)
  change_end<-c(1)
  if(length(sigsq_w_change)>1){
    change_start<-cbind(change_start,2)
    for(q in 2:length(sigsq_w_change)){
      change_start<-cbind(change_start,sum(sigsq_w_change[1:(q-1)])+2)
      change_end<-cbind(change_end,sum(sigsq_w_change[1:(q-1)])+1)
    }
    change_end<-cbind(change_end,ncol(X))
  }
  else{
    change_start<-c(change_start,2)
    change_end<-cbind(change_end,max(sigsq_w_change)+1)
  }
  sigsq_W<-rep(sigsq_w0,(change_end-change_start+1))
  
  mu_w<-matrix(rep(0,iters*(length(sigsq_w_change)+1)),ncol=(length(sigsq_w_change)+1))
  mu_w[1,]<-mu_w0
  Mu_w<-rep(mu_w0,(change_end-change_start+1))
  Var_mu_w<-rep(c(sigsq_w[1,1],var_mu_w),(change_end-change_start+1))  
  numcolZ<-ncol(data.matrix(Z))
  
  levelcolZ<-c()
  sumprevcolZ<-c()
  Z<-data.matrix(Z)
  c=ncol(Z)
  ZZ<-matrix(rep(0,length(y)),ncol=1)
  for(d in 1:c){
    sumprevcolZ[d]<-ncol(ZZ)-1
    ZZ<-cbind(ZZ,dummy(Z[,d]))
    levelcolZ[d]<-ncol(dummy(Z[,d]))
  }
  Z<-ZZ[,-1]+0
  
  b<-matrix(rep(0,ncol(data.matrix(Z))*iters),ncol=ncol(data.matrix(Z)))
  b[1,]<-b0
  sigsq_b<-matrix(rep(0,c*iters),ncol=c)
  sigsq_b[1,]<-sigsq_b0
  sigsq_B<-rep(sigsq_b0,levelcolZ)
  
  
  XG<-data.matrix(sweep(X,2,gamma0,"*"))
  XTXg<-XTX*(gamma0%*%t(gamma0))
  
  pi<-c()
  pi[1]<-pi0
  
  # Acceptance Rate Counters
  Success<-0
  Success2<-0
  Total<-0
  
  # Parameters with fixed values
  a_b<-0.5*levelcolZ + alpha_b  
  ZTZ<-t(Z)%*%Z
  NX<-ncol(X)
  
  M0<-rep(c(mu_w[1,1],rep(0,length(sigsq_w_change))),(change_end-change_start+1))
  M<-dummy(rep(1:(length(sigsq_w_change)+1),c(change_end-change_start+1)))
  V0<- M%*%diag(c(sigsq_w[1,1],var_mu_w))%*%t(M)
  V0[1,1]<-0
  for(i in 2:iters){
    #wf
    var_w<-data.matrix(solve(XTXg + diag((1/sigsq_W)) ))
    mean_w<-var_w%*%((t(XG)%*%(y - Z%*%b[i-1,])) + Mu_w/sigsq_W)
    wf[i,]<-rmvnorm(1,mean_w,sigsq_e[i-1]*var_w)
    
    
    sigsq_w[i,1]<-sigsq_w[i-1,1]
    for(p in 2:(length(sigsq_w_change)+1)){
      ALPHA_sigsq_w<-(length(change_start[p]:change_end[p]) + 2*alpha_w[p])/2
      BETA_sigsq_w<-(beta_w[p] + 0.5*sum((wf[i,change_start[p]:change_end[p]]-mu_w[i-1,p])^2)/sigsq_e[i-1])
      sigsq_w[i,p]<-rigamma(1,ALPHA_sigsq_w,BETA_sigsq_w)
    }  
    sigsq_W<-rep(sigsq_w[i,],(change_end-change_start+1))
    
    #b
    var_b<-data.matrix(solve(ZTZ/sigsq_e[i-1] + diag(1/sigsq_B)))
    mu_b<-var_b%*%(t(Z)%*%(y - XG%*%wf[i,]))/sigsq_e[i-1]
    b[i,]<-rmvnorm(1,mu_b,var_b)
    
    # sigsq_b
    for(k in 1:numcolZ){
      b_b<-beta_b[k] + 0.5*sum((b[i,((sumprevcolZ[k]+1):(sumprevcolZ[k]+levelcolZ[k]))])^2)
      sigsq_b[i,k]<-rigamma(1,a_b[k],b_b)
    }
    sigsq_B<-rep(sigsq_b[i,],levelcolZ)
    
    # pi
    if(pi.vary=="TRUE"){
      pi[i]<-rbeta(1,sum(gamma[i-1,-1])+pi_a,NX-1-sum(gamma[i-1,-1])+pi_b)
    }
    else{
      pi[i]<-pi[i-1]
    }
    
    #gamma
    YZB<-y - Z%*%b[i,]
    if(method=="fixed"){
      gamma[i,]<-gamma[i-1,]
      m<-SSdensity.BM(NY,X,XTX,gamma[i,],alpha_e,beta_e,YZB,sigsq_W,pi_a,pi_b,M0,V0)
      XTXg<-m$XTXg
      XG<-m$XG
      ySy<-m$YSY
      INV.SIGMA<-m$inv.sigma
    }
    else if(method=="log-MH"){
      if(pi_prop=="pi"){
        PIP1<-pi[i]
        PIP2<-pi[i-1]
      }
      else{
        PIP1<-pi_prop
        PIP2<-pi_prop
      }
      gamma[i,1]<-gamma[i-1,1]
      
      sample_gam<-sample((2:NX),(NX-1))
      if((ncol(X)-1)%%f==0){
        fact_samp<-c(rep(f,(NX-1)%/%f))
      }
      else{
        fact_samp<-c(rep(f,(NX-1)%/%f),(NX-1)%%f)      
      }
      # j = 1
      gamma_prop <- gamma_curr
      gamma_prop[sample_gam[1:fact_samp[1]]]<-rbern(fact_samp[1],PIP1)
      prop<-gamma_prop[sample_gam[1:fact_samp[1]]]
      curr<-gamma_curr[sample_gam[1:fact_samp[1]]]
      
      NE_prop<-sum(gamma_prop)!=1
      NE_curr<-sum(gamma_curr)!=1
      
      q_curr<-sum(dbern(curr,PIP2,log=TRUE))
      post_prop<-SSdensity.BM(NY,X,XTX,gamma_prop,alpha_e,beta_e,YZB,sigsq_W,pi_a,pi_b,M0,V0)$logmarg
      q_prop<-sum(dbern(prop,PIP1,log=TRUE))
      post_curr<-SSdensity.BM(NY,X,XTX,gamma_curr,alpha_e,beta_e,YZB,sigsq_W,pi_a,pi_b,M0,V0)$logmarg
      
      logMH<-q_curr + post_prop - q_prop - post_curr
      if(is.na(logMH>-Inf)=="FALSE"){
        logu<-log(runif(1))
        
        if(logMH>logu){
          if(sum(gamma_curr==gamma_prop)==length(gamma_curr)){
            Success2<-Success2 + 1
          }
          gamma_curr<-gamma_prop
          post_curr<-post_prop
          Success<-Success + 1
        }            
      }
      Total<-Total + 1
      
      
      for(j in 2:length(fact_samp)){
        gamma_prop <- gamma_curr
        gamma_prop[sample_gam[(((j-1)*f + 1):((j-1)*f + fact_samp[j]))]]<-rbern(fact_samp[j],PIP1)
        prop<-gamma_prop[sample_gam[(((j-1)*f + 1):((j-1)*f + fact_samp[j]))]]
        curr<-gamma_curr[sample_gam[(((j-1)*f + 1):((j-1)*f + fact_samp[j]))]]
        
        
        q_curr<-sum(dbern(curr,PIP2,log=TRUE))
        post_prop<-SSdensity.BM(NY,X,XTX,gamma_prop,alpha_e,beta_e,YZB,sigsq_W,pi_a,pi_b,M0,V0)$logmarg
        q_prop<-sum(dbern(prop,PIP1,log=TRUE))
        
        logMH<-q_curr + post_prop - q_prop - post_curr
        if(is.na(logMH>-Inf)=="FALSE"){
          logu<-log(runif(1))
          
          if(logMH>logu){
            if(mean(gamma_curr!=gamma_prop)){
              Success2<-Success2 + 1
            }
            gamma_curr<-gamma_prop
            post_curr<-post_prop
            Success<-Success + 1
          }            
        }
        Total<-Total + 1
      }
      NE_curr<-sum(gamma_curr)!=1
      gamma[i,]<-gamma_curr
      gind<-(1:NX)*gamma[i,]
      g_ind<-gind[gind!=0]
      m<-SSdensity.BM(NY,X,XTX,gamma[i,],alpha_e,beta_e,YZB,sigsq_W,pi_a,pi_b,M0,V0)
      XTXg<-m$XTXg
      XG<-m$XG
      INV.SIGMA<-m$inv.sigma
      ySy<-m$YSY
    }
    else{
      stop("Invalid method selected")
    }
    
    #sigsq_e
    a_e<-0.5*length(y) + alpha_e
    b_e<-beta_e + 0.5*ySy 
    sigsq_e[i]<-rigamma(1,a_e,b_e)
    
    
    # sigsq_w
    mu_w[i,1]<-mu_w[i-1,1]
    XM<-(XG%*%M[,-1])
    Var_mu_w1<-solve(diag(1,NY) + XG%*%((diag(sigsq_W)))%*%t(XG))
    Var_mu_w<-solve(diag((var_mu_w^(-1)),length(var_mu_w)) + t(XM)%*%Var_mu_w1%*%XM)
    MU_mu_w<-Var_mu_w%*%(t(XM)%*%Var_mu_w1%*%(y - mu_w[i,1] - Z%*%b[i,]))
    mu_w[i,-1]<-rmvnorm(1,MU_mu_w,sigsq_e[i]*Var_mu_w)
    Mu_w<-rep(mu_w[i,],(change_end-change_start+1))
    
    # Iteration Counter
    if(i/it_count==floor(i/it_count)){
      print(i)
    }
    if(time.check=="TRUE"){
      t<-c(t,difftime(Sys.time(),t0,u="secs")) 
    }
  }
  if(method=="Bern"){
    accept<-c(1,"Gibbs Sampling")
    accept2<-c(1,"Gibbs Sampling")
  }
  else{
    accept<-Success/Total
    accept2<-Success2/Total
  }
  return(list(t=t,wf=wf,sigsq_e=sigsq_e,sigsq_w=sigsq_w,b=b,sigsq_b=sigsq_b,gamma=gamma,pi=pi,iters=iters,mu_w=mu_w,accept=accept,accept2=accept2))
}
