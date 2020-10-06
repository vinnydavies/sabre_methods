library("mvtnorm")
library("Rlab")
library("coda")
library("msm")
library("pscl")
library("lme4")
library("dummies")
library("doBy")
library("matrixStats")

source("./eSABRE/eSABRE_gamma_density.R")

eSABRE<-function(y,X,Z,Challenge,Protective,iters,it_count=floor(iters/10),pi_a=1,pi_b=4,
                 gamma0=c(1,rbern((ncol(X)-1),pi_a/(pi_a + pi_b))),
                 w0=rnorm(ncol(X),c(max(y),rep(min(0,((mean(y)-max(y))/(sum(gamma0)-1))),ncol(X)-1)),sd(y)/sqrt(sum(gamma0))),
                 sigsq_e0=runif(1,0,var(y)),
                 sigsq_w0=c(100,runif(1,0,var(y))),
                 b0=rnorm(ncol(Z),(((y%*%Z)/colSums(Z))-mean((y%*%Z)/colSums(Z))),sd((((y%*%Z)/colSums(Z))-mean((y%*%Z)/colSums(Z))))),
                 sigsq_b0=runif(c,0,var(b0)),
                 mu_intercept0=max(y),
                 mu_w0=0,
                 var_mu_w0=100,
                 alpha_e=0.001,
                 beta_e=0.001,
                 alpha_w=0.001,
                 beta_w=0.001,
                 alpha_b=rep(0.001,numcolZ),
                 beta_b=rep(0.001,numcolZ),
                 pi0=rbeta(1,pi_a+sum(gamma0)-1,pi_b + length(gamma0)-sum(gamma0)),
                 pi_prop=0.3,
                 f=round((ncol(X)-1)/10),
                 alpha_y=0.001,
                 beta_y=0.001){
  
  pi.vary=FALSE
  time.check=TRUE
  method="log-MH"
  if(is.na(as.vector(Z)[1])){
    
    if(time.check=="TRUE"){
      t0<-Sys.time()
      t<-0
    }
    
    y<-c(y)
    NY<-length(y)
    
    sigsq_y0<-var(y)
    sigsq_y<-c()
    sigsq_y[1]<-sigsq_y0
    
    XX<-c()
    M<-c()
    Mvec<-c()
    Virus1<-c()
    Virus2<-c()
    LEVELS<-levels(factor(c(as.vector(Challenge),as.vector(Protective))))
    level_number<-1
    for(i in 1:length(LEVELS)){
      for(j in i:length(LEVELS)){
        if(sum((as.vector(Protective)==as.vector(LEVELS[i]))*(as.vector(Challenge)==as.vector(LEVELS[j])))>0 | sum((as.vector(Protective)==as.vector(LEVELS[j]))*(as.vector(Challenge)==as.vector(LEVELS[i])))>0 ){
          Mcol<-rep(0,NY)
          if(as.vector(LEVELS[i])!=as.vector(LEVELS[j])){
            Mcol<-Mcol + (as.vector(Protective)==as.vector(LEVELS[i]))*(as.vector(Challenge)==as.vector(LEVELS[j])) + (as.vector(Protective)==as.vector(LEVELS[j]))*(as.vector(Challenge)==as.vector(LEVELS[i]))
          }
          else{
            Mcol<-Mcol + (as.vector(Protective)==as.vector(LEVELS[i]))*(as.vector(Challenge)==as.vector(LEVELS[j]))
          }
          M<-cbind(M,Mcol)
          XX<-rbind(XX,X[which(Mcol==1)[1],])
          Virus1<-c(Virus1,as.vector(Challenge[which(Mcol==1)[1]]))
          Virus2<-c(Virus2,as.vector(Protective[which(Mcol==1)[1]]))
          Mvec[which(Mcol==1)]<-level_number
          level_number<-level_number + 1
        }
      }
    }
    X<-cbind(1,XX)
    MTMc<-colSums(M)
    NM<-ncol(M)
    
    a_e<-0.5*NM + alpha_e
    XTX<-t(X)%*%X
    mu_y<-matrix(ncol=nrow(X),nrow=iters)
    mu_y[1,]<-rep(mean(y),nrow(X))
    
    gamma<-matrix(rep(0,ncol(X)*iters),ncol=ncol(X))  
    gamma[1,]<-gamma0
    gamma_curr<-gamma0
    if(sum(gamma0)==1){
      wf<-matrix(rep(0,ncol(X)*iters),ncol=ncol(X))
      wf[1,1]<-w0[1]
    }
    else{
      wf<-matrix(rep(0,ncol(X)*iters),ncol=ncol(X))
      wf[1,]<-w0*gamma0
    }
    sigsq_e<-c()
    sigsq_e[1]<-sigsq_e0
    
    sigsq_w<-c()
    sigsq_w[1]<-sigsq_w0[2]
    
    mu_w<-c()
    mu_w[1]<-mu_w0
    
    Z<-matrix(1,nrow=length(y),ncol=1)
    
    c=ncol(Z)
    
    b<-matrix(rep(0,iters),ncol=1)
    b[1,]<-0
    sigsq_b<-matrix(rep(0,c*iters),ncol=c)
    sigsq_b[1,]<-0
    sigsq_B<-0
    
    gind<-(1:ncol(X))*gamma[1,]
    g_ind<-gind[gind!=0]
    Xg<-data.matrix(X[,g_ind])
    XTXg<-XTX[g_ind,g_ind]
    
    pi<-c()
    pi[1]<-pi0
    
    # Acceptance Rate Counters
    Success<-0
    Success2<-0
    Total<-0
    
    # Parameters with fixed values
    NX<-ncol(X)
    J<-length(gamma0[-1])
    
    Mu_W<-c(mu_intercept0,rep(mu_w[1],J))
    sigsq_W<-c(sigsq_w0[1],rep(sigsq_w[1],J))
    yZb<-y
    Xw<-Xg%*%wf[1,g_ind]
    M<-Matrix(M)
    DB<-ncol(Z)
    alpha_sigsq_y<-NY/2 + alpha_y
    
    m0<-c(mu_intercept0,rep(mu_w0,(NX-1)))
    V0<-rbind(0,cbind(0,matrix(var_mu_w0,nrow=NX-1,ncol=NX-1)))
    ZTZ<-data.matrix(t(Z)%*%Z)

    for(i in 2:iters){    
      var_w<-data.matrix(solve(XTXg + diag((1/sigsq_W[g_ind]),length(g_ind)) ))
      mean_w<-var_w%*%((t(Xg)%*%mu_y[i-1,]) + Mu_W[g_ind]/sigsq_W[g_ind])
      wf[i,g_ind]<-rmvnorm(1,mean_w,sigsq_e[i-1]*var_w)
      Xw<-Xg%*%wf[i,g_ind]
      
      # mu_y
      var_mu_y<-(1/(rep((1/sigsq_e[i-1]),NM) + (1/sigsq_y[i-1])*MTMc))
      mean_mu_y<-var_mu_y * (Xw/sigsq_e[i-1] + (1/sigsq_y[i-1])*c(matrix(t(M)%*%(y-Z%*%b[i-1,]))))
      mu_y[i,]<-rnorm(NM,mean_mu_y,sqrt(var_mu_y))
      
      # sigsq_y
      beta_sigsq_y<-0.5*sum((y - c(mu_y[i,Mvec]) - c(matrix(Z%*%b[i-1,])))^2) + beta_y
      sigsq_y[i]<-rigamma(1,alpha_sigsq_y,beta_sigsq_y)
      
      # b
      b[i,]<-b[i-1,]
      
      # mu_w
      MW<-solve(diag(1,NM) + Xg%*%diag(sigsq_W[g_ind],length(sigsq_W[g_ind]))%*%t(Xg))
      var_mu_w<-(1/((1/var_mu_w0) + t(rowSums(data.matrix(Xg[,-1])))%*%MW%*%rowSums(data.matrix(Xg[,-1]))))
      mean_mu_w<-var_mu_w*(mu_w0/var_mu_w0 + t(rowSums(data.matrix(Xg[,-1])))%*%MW%*%(mu_y[i,] - mu_intercept0))
      mu_w[i]<-rnorm(1,mean_mu_w,sqrt(sigsq_e[i-1]*var_mu_w)) # mean((REGRES_PAIRS[,which(REGRES_PAIRS[,,12]!=0),12])[-1]) # 
      Mu_W<-c(mu_intercept0,rep(mu_w[i],J))
      
      #sigsq_w
      if(sum(gamma[i-1,-1])>1){
        alpha_sigsq_w<-sum(gamma[i-1,-1])/2 + alpha_w
        beta_sigsq_w<-(1/(2*sigsq_e[i-1]))*sum((wf[i,g_ind[-1]]-mu_w[i])^2) + beta_w
        sigsq_w[i]<-rigamma(1,alpha_sigsq_w,beta_sigsq_w)
      }
      else{
        sigsq_w[i]<-0.001
      }
      sigsq_W<-c(sigsq_w0[1],rep(sigsq_w[i],J))
      
      # pi
      alpha_pi<-sum(gamma[i-1,-1])+pi_a
      beta_pi<-sum(gamma[i-1,-1]==0)+pi_b
      pi[i]<-rbeta(1,alpha_pi,beta_pi)
      
      #gamma
      if(method=="fixed"){
        gamma[i,]<-gamma[i-1,]
        m<-gamma_density(gamma[i,],pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)
        XTXg<-m$XTXg
        Xg<-m$Xg
        ySy<-m$YSY
        INV.SIGMA<-m$inv.sigma
      }
      else if(method=="Bern"){
        gamma[i,1]<-gamma[i-1,1]
        samp<-sample((2:NX),(NX-1))
        gam1 <- gamma_curr
        gam1[samp[1]]<-1
        gam0 <- gamma_curr
        gam0[samp[1]] <-0
        
        P1<-gamma_density(gam1,pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)$logmarg
        P0<-gamma_density(gam0,pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)$logmarg
        prob<-(1/(1 + exp(P0-P1)))
        
        if(is.na(prob>-1)=="FALSE"){
          gamma_curr[samp[1]]<-rbern(1,prob)
          gamma[i,samp[1]]<-gamma_curr[samp[1]]
          if(gamma_curr[samp[1]]==1){
            PP<-P1
          }
          else{
            PP<-P0
          }
        }
        else{
          gamma[i,samp[1]]<-gamma[i-1,samp[1]]
        }
        
        for(j in samp[2:(NX-1)]){
          gam1 <- gamma_curr
          gam1[j]<-1
          gam0 <- gamma_curr
          gam0[j] <-0
          if(all(gamma_curr==gam1)){
            P1<-PP
            P0<-gamma_density(gam0,pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)$logmarg
          }
          else{
            P1<-gamma_density(gam1,pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)$logmarg
            P0<-PP
          }
          prob<-(1/(1 + exp(P0-P1)))
          if(is.na(prob>-1)=="FALSE"){
            gamma_curr[j]<-rbern(1,prob)
            gamma[i,j]<-gamma_curr[j]
            if(gamma_curr[j]==1){
              PP<-P1
            }
            else{
              PP<-P0
            }
          }
          else{
            gamma[i,j]<-gamma[i-1,j]
          }
        }
        gind<-(1:NX)*gamma[i,]
        g_ind<-gind[gind!=0]
        m<-gamma_density(gamma[i,],pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)
        XTXg<-m$XTXg
        Xg<-m$Xg
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
        
        q_curr<-sum(dbern(curr,PIP2,log=TRUE))    
        post_prop<-gamma_density(gamma_prop,pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)$logmarg
        q_prop<-sum(dbern(prop,PIP1,log=TRUE))
        post_curr<-gamma_density(gamma_curr,pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)$logmarg
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
          post_prop<-gamma_density(gamma_prop,pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)$logmarg
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
        gamma[i,]<-gamma_curr
        gind<-(1:NX)*gamma[i,]
        g_ind<-gind[gind!=0]
        m<-gamma_density(gamma[i,],pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)
        XTXg<-m$XTXg
        Xg<-m$Xg
        ySy<-m$YSY
      }
      else{
        stop("Invalid method selected")
      }
      
      # sigsq_e
      b_e<-beta_e + 0.5*ySy 
      sigsq_e[i]<-rigamma(1,a_e,b_e)
      
      # Iteration Counter
      if(i/it_count==floor(i/it_count)){
        print(i)
      }
      if(time.check=="TRUE"){
        t<-c(t,difftime(Sys.time(),t0,u="secs")) 
      }
    }
    
    return(list(t=t,X=X,y=y,M=M,Virus1=Virus1,Virus2=Virus2,Z=Z,Mvec=Mvec,gamma=gamma,wf=wf,mu_w=mu_w,mu_y=mu_y,b=b,sigsq_e=sigsq_e,sigsq_b=sigsq_b,sigsq_w=sigsq_w,sigsq_y=sigsq_y,pi=pi))
    
    
  }
  else{
    
    if(time.check=="TRUE"){
      t0<-Sys.time()
      t<-0
    }
    
    y<-c(y)
    NY<-length(y)
    
    sigsq_y0<-var(y)
    sigsq_y<-c()
    sigsq_y[1]<-sigsq_y0
    
    XX<-c()
    M<-c()
    Mvec<-c()
    Virus1<-c()
    Virus2<-c()
    LEVELS<-levels(factor(c(as.vector(Challenge),as.vector(Protective))))
    level_number<-1
    for(i in 1:length(LEVELS)){
      for(j in i:length(LEVELS)){
        if(sum((as.vector(Protective)==as.vector(LEVELS[i]))*(as.vector(Challenge)==as.vector(LEVELS[j])))>0 | sum((as.vector(Protective)==as.vector(LEVELS[j]))*(as.vector(Challenge)==as.vector(LEVELS[i])))>0 ){
          Mcol<-rep(0,NY)
          if(as.vector(LEVELS[i])!=as.vector(LEVELS[j])){
            Mcol<-Mcol + (as.vector(Protective)==as.vector(LEVELS[i]))*(as.vector(Challenge)==as.vector(LEVELS[j])) + (as.vector(Protective)==as.vector(LEVELS[j]))*(as.vector(Challenge)==as.vector(LEVELS[i]))
          }
          else{
            Mcol<-Mcol + (as.vector(Protective)==as.vector(LEVELS[i]))*(as.vector(Challenge)==as.vector(LEVELS[j]))
          }
          M<-cbind(M,Mcol)
          XX<-rbind(XX,X[which(Mcol==1)[1],])
          Virus1<-c(Virus1,as.vector(Challenge[which(Mcol==1)[1]]))
          Virus2<-c(Virus2,as.vector(Protective[which(Mcol==1)[1]]))
          Mvec[which(Mcol==1)]<-level_number
          level_number<-level_number + 1
        }
      }
    }
    X<-cbind(1,XX)
    MTMc<-colSums(M)
    NM<-ncol(M)
    
    a_e<-0.5*NM + alpha_e
    XTX<-t(X)%*%X
    mu_y<-matrix(ncol=nrow(X),nrow=iters)
    mu_y[1,]<-rep(mean(y),nrow(X))
    
    gamma<-matrix(rep(0,ncol(X)*iters),ncol=ncol(X))  
    gamma[1,]<-gamma0
    gamma_curr<-gamma0
    if(sum(gamma0)==1){
      wf<-matrix(rep(0,ncol(X)*iters),ncol=ncol(X))
      wf[1,1]<-w0[1]
    }
    else{
      wf<-matrix(rep(0,ncol(X)*iters),ncol=ncol(X))
      wf[1,]<-w0*gamma0
    }
    sigsq_e<-c()
    sigsq_e[1]<-sigsq_e0
    
    sigsq_w<-c()
    sigsq_w[1]<-sigsq_w0[2]
    
    mu_w<-c()
    mu_w[1]<-mu_w0
    
    if(is.list(Z)==FALSE){
      print("WARNING: Z must be a list")
    }
    
    Zname<-list()
    ZX<-c()
    levelcolZ<-c()
    numcolZ<-dim(array(Z))
    for(i in 1:numcolZ){
      Zname[[i]]<-levels(factor(Z[[i]]))
      ZX<-cbind(ZX,Z[[i]])
      levelcolZ<-c(levelcolZ,length(levels(factor(Z[[i]]))))
    }  
    Z<-data.matrix(ZX)
    
    if(numcolZ==1){
      sumprevcolZ<-0
    }
    else{
      sumprevcolZ<-c()
      for(l in 1:numcolZ){
        sumprevcolZ[l]<-sum(levelcolZ[1:l])-levelcolZ[l]
      }
    }
    
    Z<-data.matrix(ZX)
    c=ncol(Z)
    ZZ<-c()
    for(d in 1:c){
      ZZ<-cbind(ZZ,dummy(Z[,d]))
    }
    Z<-ZZ
    colnames(Z)<-NULL
    
    b<-matrix(rep(0,ncol(data.matrix(Z))*iters),ncol=ncol(data.matrix(Z)))
    b[1,]<-b0
    sigsq_b<-matrix(rep(0,c*iters),ncol=c)
    sigsq_b[1,]<-sigsq_b0
    sigsq_B<-rep(sigsq_b0,levelcolZ)
    
    gind<-(1:ncol(X))*gamma[1,]
    g_ind<-gind[gind!=0]
    Xg<-data.matrix(X[,g_ind])
    XTXg<-XTX[g_ind,g_ind]
    
    pi<-c()
    pi[1]<-pi0
    
    # Acceptance Rate Counters
    Success<-0
    Success2<-0
    Total<-0
    
    # Parameters with fixed values
    a_b<-0.5*levelcolZ + alpha_b  
    ZTZc<-colSums(Z)
    NX<-ncol(X)
    J<-length(gamma0[-1])
    
    Mu_W<-c(mu_intercept0,rep(mu_w[1],J))
    sigsq_W<-c(sigsq_w0[1],rep(sigsq_w[1],J))
    yZb<-y-Z%*%b[1,]
    Xw<-Xg%*%wf[1,g_ind]
    Z<-Matrix(Z)
    M<-Matrix(M)
    DB<-ncol(Z)
    alpha_sigsq_y<-NY/2 + alpha_y
    
    m0<-c(mu_intercept0,rep(mu_w0,(NX-1)))
    V0<-rbind(0,cbind(0,matrix(var_mu_w0,nrow=NX-1,ncol=NX-1)))
    ZTZ<-data.matrix(t(Z)%*%Z)

    for(i in 2:iters){    
      var_w<-data.matrix(solve(XTXg + diag((1/sigsq_W[g_ind]),length(g_ind)) ))
      mean_w<-var_w%*%((t(Xg)%*%mu_y[i-1,]) + Mu_W[g_ind]/sigsq_W[g_ind])
      wf[i,g_ind]<-rmvnorm(1,mean_w,sigsq_e[i-1]*var_w)
      Xw<-Xg%*%wf[i,g_ind]
      
      # mu_y
      var_mu_y<-(1/(rep((1/sigsq_e[i-1]),NM) + (1/sigsq_y[i-1])*MTMc))
      mean_mu_y<-var_mu_y * (Xw/sigsq_e[i-1] + (1/sigsq_y[i-1])*c(matrix(t(M)%*%(y-Z%*%b[i-1,]))))
      mu_y[i,]<-rnorm(NM,mean_mu_y,sqrt(var_mu_y))
      
      # sigsq_y
      beta_sigsq_y<-0.5*sum((y - c(mu_y[i,Mvec]) - c(matrix(Z%*%b[i-1,])))^2) + beta_y
      sigsq_y[i]<-rigamma(1,alpha_sigsq_y,beta_sigsq_y)
      
      # b
      var_b<-data.matrix(solve(ZTZ/sigsq_y[i-1] + diag(1/sigsq_B)))
      mu_b<-var_b%*%(t(Z)%*%(y -  mu_y[i,Mvec]))/sigsq_y[i-1]
      b[i,]<-rmvnorm(1,mu_b,var_b)
      
      # sigsq_b
      for(k in 1:numcolZ){
        b_b<-beta_b[k] + 0.5*sum((b[i,((sumprevcolZ[k]+1):(sumprevcolZ[k]+levelcolZ[k]))])^2)
        sigsq_b[i,k]<-rigamma(1,a_b[k],b_b)
      }
      sigsq_B<-rep(sigsq_b[i,],levelcolZ)
      
      # mu_w
      MW<-solve(diag(1,NM) + Xg%*%diag(sigsq_W[g_ind],length(sigsq_W[g_ind]))%*%t(Xg))
      var_mu_w<-(1/((1/var_mu_w0) + t(rowSums(data.matrix(Xg[,-1])))%*%MW%*%rowSums(data.matrix(Xg[,-1]))))
      mean_mu_w<-var_mu_w*(mu_w0/var_mu_w0 + t(rowSums(data.matrix(Xg[,-1])))%*%MW%*%(mu_y[i,] - mu_intercept0))
      mu_w[i]<-rnorm(1,mean_mu_w,sqrt(sigsq_e[i-1]*var_mu_w)) # mean((REGRES_PAIRS[,which(REGRES_PAIRS[,,12]!=0),12])[-1]) # 
      Mu_W<-c(mu_intercept0,rep(mu_w[i],J))
      
      #sigsq_w
      if(sum(gamma[i-1,-1])>1){
        alpha_sigsq_w<-sum(gamma[i-1,-1])/2 + alpha_w
        beta_sigsq_w<-(1/(2*sigsq_e[i-1]))*sum((wf[i,g_ind[-1]]-mu_w[i])^2) + beta_w
        sigsq_w[i]<-rigamma(1,alpha_sigsq_w,beta_sigsq_w)
      }
      else{
        sigsq_w[i]<-0.001
      }
      sigsq_W<-c(sigsq_w0[1],rep(sigsq_w[i],J))
      
      # pi
      alpha_pi<-sum(gamma[i-1,-1])+pi_a
      beta_pi<-sum(gamma[i-1,-1]==0)+pi_b
      pi[i]<-rbeta(1,alpha_pi,beta_pi)
      
      #gamma
      if(method=="fixed"){
        gamma[i,]<-gamma[i-1,]
        m<-gamma_density(gamma[i,],pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)
        XTXg<-m$XTXg
        Xg<-m$Xg
        ySy<-m$YSY
        INV.SIGMA<-m$inv.sigma
      }
      else if(method=="Bern"){
        gamma[i,1]<-gamma[i-1,1]
        samp<-sample((2:NX),(NX-1))
        gam1 <- gamma_curr
        gam1[samp[1]]<-1
        gam0 <- gamma_curr
        gam0[samp[1]] <-0
        
        P1<-gamma_density(gam1,pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)$logmarg
        P0<-gamma_density(gam0,pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)$logmarg
        prob<-(1/(1 + exp(P0-P1)))
        
        if(is.na(prob>-1)=="FALSE"){
          gamma_curr[samp[1]]<-rbern(1,prob)
          gamma[i,samp[1]]<-gamma_curr[samp[1]]
          if(gamma_curr[samp[1]]==1){
            PP<-P1
          }
          else{
            PP<-P0
          }
        }
        else{
          gamma[i,samp[1]]<-gamma[i-1,samp[1]]
        }
        
        for(j in samp[2:(NX-1)]){
          gam1 <- gamma_curr
          gam1[j]<-1
          gam0 <- gamma_curr
          gam0[j] <-0
          if(all(gamma_curr==gam1)){
            P1<-PP
            P0<-gamma_density(gam0,pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)$logmarg
          }
          else{
            P1<-gamma_density(gam1,pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)$logmarg
            P0<-PP
          }
          prob<-(1/(1 + exp(P0-P1)))
          if(is.na(prob>-1)=="FALSE"){
            gamma_curr[j]<-rbern(1,prob)
            gamma[i,j]<-gamma_curr[j]
            if(gamma_curr[j]==1){
              PP<-P1
            }
            else{
              PP<-P0
            }
          }
          else{
            gamma[i,j]<-gamma[i-1,j]
          }
        }
        gind<-(1:NX)*gamma[i,]
        g_ind<-gind[gind!=0]
        m<-gamma_density(gamma[i,],pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)
        XTXg<-m$XTXg
        Xg<-m$Xg
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
        
        q_curr<-sum(dbern(curr,PIP2,log=TRUE))    
        post_prop<-gamma_density(gamma_prop,pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)$logmarg
        q_prop<-sum(dbern(prop,PIP1,log=TRUE))
        post_curr<-gamma_density(gamma_curr,pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)$logmarg
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
          post_prop<-gamma_density(gamma_prop,pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)$logmarg
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
        gamma[i,]<-gamma_curr
        gind<-(1:NX)*gamma[i,]
        g_ind<-gind[gind!=0]
        m<-gamma_density(gamma[i,],pi_a,pi_b,m0,V0,mu_y[i,],X,NM,sigsq_W,XTX,alpha_e,beta_e)
        XTXg<-m$XTXg
        Xg<-m$Xg
        ySy<-m$YSY
      }
      else{
        stop("Invalid method selected")
      }
      
      # sigsq_e
      b_e<-beta_e + 0.5*ySy 
      sigsq_e[i]<-rigamma(1,a_e,b_e)
      
      # Iteration Counter
      if(i/it_count==floor(i/it_count)){
        print(i)
      }
      if(time.check=="TRUE"){
        t<-c(t,difftime(Sys.time(),t0,u="secs")) 
      }
    }
    
    return(list(t=t,Zname=Zname,X=X,y=y,M=M,levelcolZ=levelcolZ,Virus1=Virus1,Virus2=Virus2,Z=Z,Mvec=Mvec,gamma=gamma,wf=wf,mu_w=mu_w,mu_y=mu_y,b=b,sigsq_e=sigsq_e,sigsq_b=sigsq_b,sigsq_w=sigsq_w,sigsq_y=sigsq_y,pi=pi))
  }
}