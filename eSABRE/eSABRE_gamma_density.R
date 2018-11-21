gamma_density<-function(GAMMA,pi_a,pi_b,m0,V0,mu_y,X,NM,sigsq_W,XTX,alpha_e,beta_e){
  g_ind<-which(GAMMA==1)
  Xg<-data.matrix(X[,g_ind])
  sigsq_W<-sigsq_W[g_ind]
  XTXg<-XTX[g_ind,g_ind]
  m0g<-m0[g_ind]
  muy.minus.mu0<-mu_y-Xg%*%m0g
  V0g<-V0[g_ind,g_ind]
  Sigma1<-diag(sigsq_W,length(sigsq_W)) + V0g
  Sigma1.chol<-chol(Sigma1)
  sigma<-chol2inv(Sigma1.chol) + XTXg
  sigma.chol <- chol(sigma)
  log.sqrt.det<-sum(log(diag(sigma.chol))) + sum(log(diag(Sigma1.chol)))
  inv.sigma<-chol2inv(sigma.chol)
  inv<-diag(1,NM) - Xg%*%inv.sigma%*%t(Xg)
  YSY<-t(muy.minus.mu0)%*%inv%*%muy.minus.mu0
  logmarg<-lgamma((sum(GAMMA[2:length(GAMMA)]) + pi_a)) + lgamma((length(GAMMA)-1-sum(GAMMA[2:length(GAMMA)]) + pi_b)) - lgamma((length(GAMMA)-1 + pi_a + pi_b)) + drop( -log.sqrt.det) -(NM/2 + alpha_e)*log(beta_e + 0.5*YSY)
  return(list(logmarg=logmarg,Xg=Xg,XTXg=XTXg,inv.sigma=inv.sigma,YSY=YSY))  
}