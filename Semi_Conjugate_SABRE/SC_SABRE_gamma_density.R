SSdensity<-function(N,X,XTX,GAMMA,sigsq_e,yZB,Mu_w,sigsq_W,pi_a,pi_b){
  g_ind<-which(GAMMA==1)
  Xg<-data.matrix(X[,g_ind])
  sigsq_W<-sigsq_W[g_ind]
  Mu_W<-Mu_w[g_ind]
  XTXg<-XTX[g_ind,g_ind]
  y.minus.mu<-yZB-Xg%*%Mu_W
  sigma<-diag((sigsq_W)^(-1),length(g_ind)) + ((sigsq_e)^(-1))*XTXg
  sigma.chol <- chol(sigma)
  log.sqrt.det<-sum(log(diag(sigma.chol))) + 0.5*sum(log(sigsq_W))
  inv<-diag(sigsq_e^(-1),N) - (sigsq_e^(-2))*Xg%*%chol2inv(sigma.chol)%*%t(Xg)
  logmarg<-lgamma((sum(GAMMA[2:length(GAMMA)]) + pi_a)) + lgamma((length(GAMMA)-1-sum(GAMMA[2:length(GAMMA)]) + pi_b)) - lgamma((length(GAMMA)-1 + pi_a + pi_b)) + drop( -log.sqrt.det -0.5*t(y.minus.mu)%*%inv%*%y.minus.mu)
  return(list(logmarg=logmarg,Xg=Xg,XTXg=XTXg))
}
