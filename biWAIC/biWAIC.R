library("mvtnorm")
library("Rlab")
library("coda")
library("msm")
library("pscl")
library("lme4")
library("dummies")
library("doBy")
library("matrixStats")

WAIC_post<-function(log_lik){
  offset <- log_lik[cbind(max.col(abs(t(log_lik))), 1:ncol(log_lik))]/2
  lpd <- log(1/nrow(log_lik)) + log(colSums(exp(sweep(log_lik, 2, offset))))+offset
  var<-colVars(log_lik)
  WAIC<-(-2*sum(lpd -var))
  return(WAIC=WAIC)
}

mvdnorm2 <- function(x, mu, sigma) {
  if (is.vector(x)) x <- t(x)
  x.minus.mu <- t(sweep(x,2,mu,'-'))
  sigma.chol <- chol(sigma)
  sqrt.det <- prod(diag(sigma.chol))
  exp.arg <- -0.5 * colSums(x.minus.mu * backsolve(sigma.chol,forwardsolve(sigma.chol,x.minus.mu,upper.tri=TRUE,transpose=TRUE)))
  drop( -log(sqrt.det) + exp.arg - (ncol(x)/2)*log(2*pi) )
}

WAIC_eSABRE<-function(model,it_seq){
  if(sum(model$b)!=0){
    B<-model$b[it_seq,]
    SIGSQ_Y<-model$sigsq_y[it_seq]
    SIGSQ_E<-model$sigsq_e[it_seq]
    GAM<-model$gamma[it_seq-1,]
    W<-model$wf[it_seq,]
    MU_Y<-model$mu_y[it_seq,]
    nLIKE<-c()
    biLIKE<-matrix(0,ncol=max(model$Mvec),nrow=length(it_seq))
    for(i in 1:length(it_seq)){
      l1<-dnorm(c(model$y), (MU_Y[i,])[model$Mvec]+ c(matrix(model$Z%*%B[i,])),sqrt(SIGSQ_Y[i]),log=TRUE)
      nLIKE<-rbind(nLIKE,l1)
      for(j in 1:max(model$Mvec)){
        ch<-which(model$Mvec==j)
        biLIKE[i,j]<-mvdnorm2(c(model$y[ch]), rep(model$X[j,]%*%W[i,],length(ch)) + c(matrix(model$Z[ch,]%*%B[i,])), diag(SIGSQ_Y[i],length(ch)) + matrix(SIGSQ_E[i],ncol=length(ch),nrow=length(ch)))
      }
      if(i/100==floor(i/100)){
        print(i)
      }
    }
    nWAIC<-WAIC_post(nLIKE)
    biWAIC<-WAIC_post(biLIKE)
    return(list(biWAIC=biWAIC,nWAIC=nWAIC,biLIKE=biLIKE,nLIKE=nLIKE))
  }
  else{
    SIGSQ_Y<-model$sigsq_y[it_seq]
    SIGSQ_E<-model$sigsq_e[it_seq]
    GAM<-model$gamma[it_seq-1,]
    W<-model$wf[it_seq,]
    MU_Y<-model$mu_y[it_seq,]
    nLIKE<-c()
    biLIKE<-matrix(0,ncol=max(model$Mvec),nrow=length(it_seq))
    for(i in 1:length(it_seq)){
      l1<-dnorm(c(model$y), (MU_Y[i,])[model$Mvec],sqrt(SIGSQ_Y[i]),log=TRUE)
      nLIKE<-rbind(nLIKE,l1)
      for(j in 1:max(model$Mvec)){
        ch<-which(model$Mvec==j)
        biLIKE[i,j]<-mvdnorm2(c(model$y[ch]), rep(model$X[j,]%*%W[i,],length(ch)), diag(SIGSQ_Y[i],length(ch)) + matrix(SIGSQ_E[i],ncol=length(ch),nrow=length(ch)))
      }
      if(i/100==floor(i/100)){
        print(i)
      }
    }
    nWAIC<-WAIC_post(nLIKE)
    biWAIC<-WAIC_post(biLIKE)
    return(list(biWAIC=biWAIC,nWAIC=nWAIC,biLIKE=biLIKE,nLIKE=nLIKE))
  }
}
