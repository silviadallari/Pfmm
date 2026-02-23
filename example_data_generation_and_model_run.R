#Illustrative example of data generation and application of the proposed model:


# Data Generation illustrated in Section 4.1 of the paper: ---------------------------------------------------------

library(mvtnorm)
library(MASS)
library(combinat)

nrep=100
set.seed(123)
numobs=300
r=2
k=3
p=50
m=1

w=matrix(0,k)
w[1]=0.3
w[2]=0.3
w[3]=0.4


# Mean and variance generation: -------------------------------------------

mu<-matrix(0,k,r)

mu[2,1]=-1.19
mu[2,2]=0.77
mu[1,1]=1.20
mu[1,2]=0.76
mu[3,]=-(w[1]*mu[1,]+w[2]*mu[2,])/w[3]
# w[1]*mu[1,]+w[2]*mu[2,]+w[3]*mu[3,] check mu

sigma<-array(0,c(k,r,r))
sigma[1,1,1]=0.17
sigma[1,1,2]=0.08
sigma[1,2,1]=0.08
sigma[1,2,2]=0.14


sigma[2,1,1]=0.16
sigma[2,1,2]=-0.08
sigma[2,2,1]=-0.08
sigma[2,2,2]=0.12
temp1=w[1]*sigma[1,,]+w[2]*sigma[2,,]
temp2=w[1]*mu[1,]%*%t(mu[1,])+w[2]*mu[2,]%*%t(mu[2,])+w[3]*mu[3,]%*%t(mu[3,])
sigma[3,,]<-(diag(r)-temp1 - temp2)/w[3]

## check 
##temp1+w[3]*sigma[3,,]+temp2
#library(matrixcalc)
#is.positive.semi.definite(sigma[3,,])



# Loadings: ----------------------------------------------------------------

set.seed(3)

alpha.zero=runif(p,-1,1)
alpha.zero[1]=0

alpha.1.1=runif(p/2,0.5,1)
alpha.1.2=runif(p/2,-0.5,0)

alpha.2.1=runif(p/2,-0.5,0)
alpha.2.2=runif(p/2,0.5,1)

alpha=matrix(0,p,r)
alpha[,1]=c(alpha.1.1,alpha.1.2)
alpha[,2]=c(alpha.2.1,alpha.2.2)
alpha[1,2]=0 

alpha=cbind(alpha.zero,alpha)
alpha.true=alpha

mu.true=mu
sigma.true=sigma
w.true=w


# Covariates: --------------------------------------------------------------

mu.x=c(-2,0,4)


# Generation: -------------------------------------------------------------

index.true=list()
y=list()
lambda.true=list()

cov= vector("list", nrep)
cov[]=list(matrix(0, nrow = numobs, ncol = m))

for(jj in 1:nrep){
  set.seed(jj)
  index<-sample(1:k,numobs,prob=w,replace=TRUE)
  z=matrix(0,numobs,r)
  for (i in 1: numobs) {
    z[i,]=rmvnorm(1,mu[index[i],],sigma[index[i],,]) 
    cov[[jj]][i]=rnorm(1, mu.x[index[i]])  
  }
  
  num=matrix(alpha.zero,p,numobs) +alpha[,-1]%*%t(z)
  lambda=exp(num)
  lambda=t(lambda)
  x=matrix(0,numobs,p)
  for (i in 1:numobs) for (j in 1:p) x[i,j]=rpois(1,lambda[i,j])
  y[[jj]]=x
  lambda.true[[jj]]=lambda
  index.true[[jj]]=index
  
}


# Examples: ---------------------------------------------------------------

source(pfmm_functions.R)

# Example with no covariates: ---------------------------------------------

example_nocov=pfmm.em(y[[1]],k,r,seed=2,eps=0.00001)  
library(mclust)
adjustedRandIndex(example_nocov$index,index.true[[1]])



# Example with covariates: -------------------------------------------------

example_cov=pfmm.em.cov(y[[1]],k,r,seed=2,eps=0.00001,x=cov[[1]])  
library(mclust)
adjustedRandIndex(example_cov$index,index.true[[1]])

