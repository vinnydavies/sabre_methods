SSdensity.BM<-function(N,X,XTX,GAMMA,alpha_e,beta_e,yZB,sigsq_W,pi_a,pi_b,M0,V0){
  XG<-data.matrix(sweep(X,2,GAMMA,"*"))
  XTXg<-XTX*(GAMMA%*%t(GAMMA))
  y.minus.mu<-yZB-XG%*%M0
  Sigma1<-diag(sigsq_W) + V0
  sigma<-solve(Sigma1) + XTXg
  log.sqrt.det<-0.5*determinant(sigma,logarithm=TRUE)$modulus[1] + 0.5*determinant(Sigma1,logarithm=TRUE)$modulus[1]
  inv.sigma<-solve(sigma)
  inv<-diag(1,N) - XG%*%inv.sigma%*%t(XG)
  YSY<-t(y.minus.mu)%*%inv%*%y.minus.mu
  logmarg<-lgamma((sum(GAMMA[2:length(GAMMA)]) + pi_a)) + lgamma((length(GAMMA)-1-sum(GAMMA[2:length(GAMMA)]) + pi_b)) - lgamma((length(GAMMA)-1 + pi_a + pi_b)) + drop( -log.sqrt.det) -(N/2 + alpha_e)*log(beta_e + 0.5*YSY)
  return(list(logmarg=logmarg,XG=XG,XTXg=XTXg,inv.sigma=inv.sigma,YSY=YSY))
}
