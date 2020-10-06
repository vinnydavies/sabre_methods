# Adjust below line to set your working directory
setwd("C:/Users/Vinny/work/sabre_methods")

# install any necessary packages here
# required packages are: mvtnorm, Rlab, coda, msm, pscl, lme4, dummies, doBy, matrixStats

###### load code for eSABRE and SABRE methods #######

source("./original_methods/eSABRE/eSABRE.R")
source("./original_methods/eSABRE/eSABRE_TryCatch.R")
source("./original_methods/Semi_Conjugate_SABRE/SC_SABRE.R")
source("./original_methods/Binary_Mask_SABRE/BM_SABRE.R")
source("./original_methods/Conjugate_SABRE/C_SABRE.R")
source("./original_methods/biWAIC/biWAIC.R")

#### end ####

######## Generate the simulated datasets ######

N<-1000 # length of y
Num_Datasets<-20 # number of datasets
J<-50 # number of regression coefficients in each dataset
n_chal<-10 # number of challenge strains in each dataset
n_prot<-10 # number of protective strains in each dataset

sigsq_e = 0.1
sigsq_y = 0.1

challenge<-array(as.character(sample(1:n_chal,N*Num_Datasets,replace=TRUE)),dim=c(1,N,Num_Datasets))
protective<-array(as.character(sample(1:n_prot,N*Num_Datasets,replace=TRUE)),dim=c(1,N,Num_Datasets))
z1<-array(sample(c("a","b","c","d","e","f","g","h"),N,replace=TRUE),dim=c(1,N,Num_Datasets))
z2<-array(sample(c("L1","L2","L3","L4","L5","L6","L7","L8"),N,replace=TRUE),dim=c(1,N,Num_Datasets))

pi_prob<-runif(Num_Datasets,0.2,0.4)
regres<-array(NA,dim=c(1,J+1,Num_Datasets))
regres[,1,]<-10
for(i in 1:Num_Datasets){
  regres[,-1,i]<-runif(J,-0.4,-0.2)*rbinom(J,1,pi_prob[i])
}

X_base<-array(rbinom(n_chal*n_prot*J*Num_Datasets,1,0.3),dim=c(n_chal*n_prot,J,Num_Datasets))
mu_y_base<-array(NA,dim=c(1,n_chal*n_prot,Num_Datasets)) 
for(i in 1:Num_Datasets){
  mu_y_base[,,i]<-rnorm(n_chal*n_prot,c(cbind(1,X_base[,,i])%*%regres[,,i]),sd=sqrt(sigsq_e))
}

# z_select must be updated if a diffent selection of random effects is required
# Similarly, Z must be updated in the methods
z_select<-array(c(1,1,0,0),dim=c(1,4,Num_Datasets))
#z_select<-array(rbinom(4*Num_Datasets,1,0.5),dim=c(1,4,Num_Datasets))
Z<-array(NA,dim=c(N,8+8+n_chal+n_prot,Num_Datasets))
b<-array(NA,dim=c(1,8+8+n_chal+n_prot,Num_Datasets))
for(i in 1:Num_Datasets){
  Z[,,i]<-cbind(dummy(z1[,,i]),dummy(z2[,,i]),dummy(challenge[,,i]),dummy(protective[,,i]))
  b[,,i]<-c(z_select[,1,i]*rnorm(8,0,sd=sqrt(runif(1,0.2,0.4))),z_select[,2,i]*rnorm(8,0,sd=sqrt(runif(1,0.2,0.4))),z_select[,3,i]*rnorm(n_chal,0,sd=sqrt(runif(1,0.2,0.4))),z_select[,4,i]*rnorm(n_prot,0,sd=sqrt(runif(1,0.2,0.4))))
}

y<-array(NA,dim=c(1,N,Num_Datasets))
X<-array(NA,dim=c(N,J,Num_Datasets))
for(j in 1:Num_Datasets){
  for(i in 1:N){
    y[,i,j]<-rnorm(1,mu_y_base[,(as.numeric(challenge[i])-1)*n_prot + as.numeric(protective[i]),j],sqrt(sigsq_y))
    X[i,,j]<-X_base[(as.numeric(challenge[i])-1)*n_prot + as.numeric(protective[i]),,j]
  }
}
  

#### end ####

########## Evalutate the Simulated Data for SABRE methods #########

# we recommend a larger number of iterations in practise
model_conjugate_sabre<-C_SABRE(y[,,1],as.matrix(X[,,1]),cbind(z1[,,1],z2[,,1]),iters=1000,it_count=10)
model_semi_conjugate_sabre<-SC_SABRE(y[,,1],as.matrix(X[,,1]),cbind(z1[,,1],z2[,,1]),iters=1000,it_count=10)
model_BM_conjugate_sabre<-BM_SABRE(y[,,1],as.matrix(X[,,1]),cbind(z1[,,1],z2[,,1]),iters=1000,it_count=10)

#### end ####

###### Evalutate the Simulated Data for the eSABRE methods ####

# we recommend a larger number of iterations in practise
model_esabre_1<-eSABRE(y[,,1],as.matrix(X[,,1]),list(z1[,,1],z2[,,1]),Challenge=challenge[,,1],Protective=protective[,,1],iters=1000,it_count=10)
model_esabre_2<-eSABRE(y[,,1],as.matrix(X[,,1]),NA,Challenge=challenge[,,1],Protective=protective[,,1],iters=1000,it_count=10)

#### end ####

###### Evaluate the Simulated Data for the eSABRE method with Try-Catch ####

model_esabre_TC<-eSABRE_TryCatch(y[,,1],as.matrix(X[,,1]),list(z1[,,1],z2[,,1]),Challenge=challenge[,,1],Protective=protective[,,1],iters=1000,it_count=10)

#### end ####

###### Calculate biWAIC and nWAIC ########

# we recommend a larger number of iterations in practise

model_esabre_waic_1<-WAIC_eSABRE(model_esabre_1,seq(501,1000,10))
nWAIC_1<-model_esabre_waic_1$nWAIC
biWAIC_1<-model_esabre_waic_1$biWAIC

model_esabre_waic_2<-WAIC_eSABRE(model_esabre_2,seq(501,1000,10))
nWAIC_2<-model_esabre_waic_2$nWAIC
biWAIC_2<-model_esabre_waic_2$biWAIC

#### end ###
