# Main function ####
# y: data matrix.
# k: number of groups
# r: number of latent variables
# x: covariate matrix


pfmm.em=function(y,k,r,it=100,eps=0.0001,init=NULL,seed=7,pp=gauher(8),maxstep=100)
{
  
  ######################################################################################
  # y is the data of dimension n x p, with n sample size and p number of variables
  # k is the number of groups
  # r is the number of latent variables/factors
  ######################################################################################
  
  library(combinat)
  ptm = proc.time()
  set.seed(seed)
  
  ### initialization
  
  numobs=nrow(y)
  p=ncol(y)
  lik=-1000000000
  
  if (is.null(init$alpha)) {alpha=factanal(y,r)$loadings[,r:1]
  alpha.0=runif(p,-1,1)
  alpha=cbind(alpha.0,alpha)} else alpha=init$alpha
  
  ## Identifiability conditions:
  if (r>1) for (i in 2:r) alpha[1:(i-1),(i+1)]=0
  alpha[1,1]=0
  
  if (k>1) memb=cutree(hclust(dist(y),method="ward.D"),k) else memb=rep(1,numobs)
  if (is.null(init$w)) {  w=table(memb)/numobs
  w=t(t(w))} else w=init$w
  
  if (is.null(init$mu)) {mu=matrix(runif(k*r,-1,1),k,r)} else mu=init$mu
  
  if (is.null(init$sigma)) { sigma<-array(0,c(k,r,r))
  for (i in 1:k) sigma[i,,]=0.5*diag(r) } else sigma=init$sigma
  
  out=pfmm.em.alg(y,numobs,r,k,p,it,alpha,w,mu,sigma,eps,lik,pp,maxstep)
  
  Ez.y=out$Ez.y
  likelihood=out$likelihood
  ps.y=out$ps.y
  py.s=out$py.s
  alpha<-out$alpha
  w<-out$w
  mu<-out$mu
  sigma<-out$sigma
  likelihood<-matrix(likelihood[!likelihood==0])
  py=out$py
  h=p*(r+1)-(r*(r-1)/2+1)+(r*(r+1)/2)*(k-1)+r*(k-1)+k-1
  lik=likelihood[length(likelihood)]
  bic=-2*lik+h*log(numobs)
  aic=-2*lik+2*h
  EN=entr(ps.y)
  clc=-2*lik+2*EN
  icl.bic<--2*lik+2*EN+h*log(numobs)
  criteri=cbind(bic,aic,clc,icl.bic)
  
  if (k>1)  {  index<-(apply(ps.y, 1, order))[k,]} else index<-rep(k,numobs)
  
  output=list(y=y,py=py,index=index,alpha=alpha,lik=likelihood,w=w,mu=mu,k=k,sigma=sigma,ps.y=ps.y,p=p,numobs=numobs,r=r,py.s=py.s,tempo=proc.time()-ptm,
              aic=aic,bic=bic,clc=clc,Ez.y=Ez.y,pp=pp,criteri=criteri)
  invisible(output)
}


pfmm.em.alg=function(y,numobs,r,k,p,it,alpha,w,mu,sigma,eps,lik,pp,maxstep)
  
  
{
  likelihood=NULL
  hh=0
  ratio=1000
  
  nq=length(pp$x)^r
  py.s=matrix(0,numobs,k)
  ps.y=matrix(0,numobs,k)
  if (r>1) punteggi=sqrt(2)*z_ext(pp$x,r) else punteggi=matrix(sqrt(2)*pp$x,nq)
  if (r>1) pesi=z_ext(pp$w,r)/sqrt(pi) else pesi=matrix(pp$w/sqrt(pi),nq)
  E.z.sy=array(0,c(numobs,k,r))
  E.zz.sy=array(0,c(numobs,k,r,r))
  zex=array(0,c(nq,r,k))
  p.z.ys=array(0,c(nq,numobs,k))
  CC=array(0,c(r,r,numobs))
  
  while ((hh < it) & (ratio > eps )) {
    
    hh<-hh+1
    
    #### E step determine p(y|s=i)
    
    for (i in 1:k) {
      Ci=chol(sigma[i,,])
      zex[,,i]=(punteggi%*%(Ci))+t(matrix(mu[i,],r,nq))
      ogg<- exp(alpha%*%t(cbind(1,zex[,,i])))
      ogg<-ifelse(is.infinite(ogg),1/eps,ogg)
      lambda=t(ogg)
      lambda.g=array(lambda,c(nq,p,numobs))
      yg=aperm(array(y,c(numobs,p,nq)),c(3,2,1))
      
      py.zex=apply(yg*log(lambda.g)-lambda.g- lfactorial(yg),c(1,3),sum)
      py.zex=exp(py.zex)
      
      py.s[,i]=apply(pesi,1,prod)%*%py.zex
      py.s[,i]=ifelse(is.infinite(py.s[,i]),1/eps,py.s[,i])
      py.s[,i]=ifelse(is.na(py.s[,i]),eps,py.s[,i])
      ### E.z.sy
      den=t(matrix(py.s[,i],numobs,nq))
      den=ifelse(den==0,eps,den)
      p.z.ys[,,i]=matrix(apply(pesi,1,prod),nq,numobs)*py.zex/den
      E.z.sy[,i,]=t(t(zex[,,i])%*%p.z.ys[,,i])
      E.z.sy[,i,]=ifelse(is.na(E.z.sy[,i,]),exp(-10),E.z.sy[,i,])
      E.z.sy[,i,]=ifelse(is.infinite(E.z.sy[,i,]),10^10,E.z.sy[,i,])
      ### E.zz.sy
      temp=((matrix(zex[,,i],nrow=nq)))%o%t(matrix(zex[,,i],nrow=nq))
      temp=apply(temp,c(2,3),diag)
      temp=aperm(temp,c(2,3,1))
      for (l in 1:r) CC[l,,]=temp[l,,]%*%p.z.ys[,,i]
      E.zz.sy[,i,,]=aperm(CC,c(3,1,2))
      E.zz.sy[,i,,]=ifelse(is.na(E.zz.sy[,i,,]),eps,E.zz.sy[,i,,])
      E.zz.sy[,i,,]=ifelse(is.infinite(E.zz.sy[,i,,]),1/eps,E.zz.sy[,i,,])
      
    }
    
    
    p.y=py.s%*%w
    p.y=ifelse(p.y==0,eps,p.y)
    Ez.y=matrix(0,numobs,r)
    for (i in 1:k) {ps.y[,i]=w[i]*py.s[,i]/p.y
    Ez.y=Ez.y+matrix(ps.y[,i],numobs,r)*E.z.sy[,i,]
    }
    
    
    
    #### M step estimate w
    w=colMeans(ps.y)
    
    #### M step estimate mu and sigma
    
    temp1=matrix(0,r,r)
    temp2=matrix(0,r,r)
    temp3=matrix(0,r)
    for (i in 1:k) {mu[i,]=colSums(matrix(ps.y[,i],numobs,r)*E.z.sy[,i,])/sum(ps.y[,i])
    sigma[i,,]=apply((array(ps.y[,i],c(numobs,r,r))*(E.zz.sy[,i,,]-aperm(array((mu[i,])%*%t(mu[i,]),c(r,r,numobs)),c(3,1,2)))),c(2,3),sum)/sum(ps.y[,i])
    temp1=temp1+w[i]*sigma[i,,]
    temp2=temp2+w[i]*(mu[i,]%*%t(mu[i,]))
    temp3=temp3+w[i]*mu[i,]}
    
    
    for (j in 1:p) {#opt=optim(alpha[j,],fn=da.max,gr=NULL,method="L-BFGS-B",lower=-25,upper=25,y=y[,j],zex=zex,k=k,ps.y=ps.y,p.z.ys=p.z.ys)
      opt=nlm(f=da.max,p=alpha[j,],stepmax=maxstep,y=y[,j],zex=zex,k=k,ps.y=ps.y,p.z.ys=p.z.ys)
      alpha[j,]=opt$estimate}
    
    
    ## Identifiability corrections:
    
    alpha[1,1]=0
    
    var.z=temp1+temp2-temp3%*%t(temp3)
    A<-(chol(var.z))
    for (i in 1:k) {sigma[i,,]<-t(ginv(A))%*%sigma[i,,]%*%ginv(A)
    mu[i,]<-t(ginv(A))%*%mu[i,]
    }
    
    mu.tot=t(w)%*%mu
    mu=mu-t(matrix(mu.tot,r,k))
    alpha[,-1]=alpha[,-1]%*%t(A)
    
    
    ## Zero constraint:
    if (r>1) for (i in 2:r) alpha[1:(i-1),(i+1)]=0 
    
    
    
    
    ######################################
    temp<- sum(log(p.y))
    likelihood<-c(likelihood,temp)
    ratio<-(temp-lik)/abs(lik)
    if (hh < 6) ratio=2*eps
    lik<-temp
  }
  
  
  out<-list(alpha=alpha,w=w,mu=mu,sigma=sigma,likelihood=likelihood,ps.y=ps.y,py.s=py.s,Ez.y=Ez.y,py=p.y)
  return(out)
}


pfmm.em.cov=function(y,k,r,it=500,eps=0.00001,init=NULL,seed=7,pp=gauher(8),maxstep=100,x=x){
  
  ######################################################################################
  # y is the data of dimension n x p, with n sample size and p number of variables
  # k is the number of groups
  # r is the number of latent variables/factors
  # x is the data of dimension n x m, with n sample size and m number of covariates
  ######################################################################################
  library(psych)
  library(combinat)
  ptm = proc.time()
  set.seed(seed)
  
  ### initialization
  m=ncol(x)
  numobs=nrow(y)
  p=ncol(y)
  lik=-1000000000
  
  if (is.null(init$alpha)) {alpha=fa(r=cor(y),nfactors=r,rotate="none",fm="pa",max.iter = 500)$loadings[,r:1]
  alpha.0=runif(p,-1,1)
  alpha=cbind(alpha.0,alpha)} else alpha=init$alpha
  
  ## identifiability conditions:
  if (r>1) for (i in 2:r) alpha[1:(i-1),(i+1)]=0
  alpha[1,1]=0
  
  if (k>1) memb=cutree(hclust(dist(y),method="ward.D"),k) else memb=rep(1,numobs)
  if (is.null(init$w)) {  w=table(memb)/numobs
  w=t(t(w))} else w=init$w
  
  if (is.null(init$mu)) {mu=matrix(runif(k*r,-1,1),k,r)} else mu=init$mu
  
  if (is.null(init$sigma)) { sigma<-array(0,c(k,r,r))
  for (i in 1:k) sigma[i,,]=0.5*diag(r) } else sigma=init$sigma
  
  out=pfmm.em.alg.cov(y,numobs,r,k,p,it,alpha,w,mu,sigma,eps,lik,pp,maxstep,x)
  
  Ez.y=out$Ez.y
  likelihood=out$likelihood
  ps.y=out$ps.y
  py.s=out$py.s
  alpha<-out$alpha
  beta<-out$beta
  w<-out$w
  mu<-out$mu
  sigma<-out$sigma
  phi=out$phi
  q.w=out$q.w
  likelihood<-matrix(likelihood[!likelihood==0])
  py=out$py
  h=p*(r+1)-(r*(r-1)/2+1)+(r*(r+1)/2)*(k-1)+r*(k-1)+(k-1)*q.w
  lik=likelihood[length(likelihood)]
  bic=-2*lik+h*log(numobs)
  aic=-2*lik+2*h
  EN=entr(ps.y)
  clc=-2*lik+2*EN
  icl.bic<--2*lik+2*EN+h*log(numobs)
  criteri=cbind(bic,aic,clc,icl.bic)
  
  if (k>1)  {  index<-(apply(ps.y, 1, order))[k,]} else index<-rep(k,numobs)
  
  output=list(y=y,py=py,index=index,alpha=alpha,lik=likelihood,w=w,mu=mu,k=k,sigma=sigma,ps.y=ps.y,p=p,numobs=numobs,r=r,py.s=py.s,tempo=proc.time()-ptm,aic=aic,bic=bic,clc=clc,Ez.y=Ez.y,pp=pp,criteri=criteri,
              phi=phi)
  invisible(output)
}


pfmm.em.alg.cov=function(y,numobs,r,k,p,it,alpha,w,mu,sigma,eps,lik,pp,maxstep,x){
  likelihood=NULL
  hh=0
  ratio=1000
  x.w<-as.matrix(cbind(rep(1,numobs),x))
  q.w<-ncol(x.w)
  phi<-matrix(0,k,q.w)
  
  
  
  nq=length(pp$x)^r
  py.s=matrix(0,numobs,k)
  ps.y=matrix(0,numobs,k)
  if (r>1) punteggi=sqrt(2)*z_ext(pp$x,r) else punteggi=matrix(sqrt(2)*pp$x,nq)
  if (r>1) pesi=z_ext(pp$w,r)/sqrt(pi) else pesi=matrix(pp$w/sqrt(pi),nq)
  E.z.sy=array(0,c(numobs,k,r))
  E.zz.sy=array(0,c(numobs,k,r,r))
  zex=array(0,c(nq,r,k))
  p.z.ys=array(0,c(nq,numobs,k))
  CC=array(0,c(r,r,numobs))
  w=matrix(w,nrow=numobs,ncol=k, byrow=T)
  
  while ((hh < it) & (ratio > eps )) {
    
    hh<-hh+1
    
    #### E step determine p(y|s=i)
    
    for (i in 1:k) {
      Ci=chol(sigma[i,,])
      zex[,,i]=(punteggi%*%(Ci))+t(matrix(mu[i,],r,nq))
      ogg<- exp(alpha%*%t(cbind(1,zex[,,i])))
      ogg<-ifelse(is.infinite(ogg),1/eps,ogg)
      lambda=t(ogg)
      lambda.g=array(lambda,c(nq,p,numobs))
      yg=aperm(array(y,c(numobs,p,nq)),c(3,2,1))
      
      py.zex=apply(yg*log(lambda.g)-lambda.g- lfactorial(yg),c(1,3),sum)
      py.zex=exp(py.zex)
      
      py.s[,i]=apply(pesi,1,prod)%*%py.zex
      py.s[,i]=ifelse(is.infinite(py.s[,i]),1/eps,py.s[,i])
      py.s[,i]=ifelse(is.na(py.s[,i]),eps,py.s[,i])
      ### E.z.sy
      den=t(matrix(py.s[,i],numobs,nq))
      den=ifelse(den==0,eps,den)
      p.z.ys[,,i]=matrix(apply(pesi,1,prod),nq,numobs)*py.zex/den
      E.z.sy[,i,]=t(t(zex[,,i])%*%p.z.ys[,,i])
      E.z.sy[,i,]=ifelse(is.na(E.z.sy[,i,]),exp(-10),E.z.sy[,i,])
      E.z.sy[,i,]=ifelse(is.infinite(E.z.sy[,i,]),10^10,E.z.sy[,i,])
      ### E.zz.sy
      temp=((matrix(zex[,,i],nrow=nq)))%o%t(matrix(zex[,,i],nrow=nq))
      temp=apply(temp,c(2,3),diag)
      temp=aperm(temp,c(2,3,1))
      for (l in 1:r) CC[l,,]=temp[l,,]%*%p.z.ys[,,i]
      E.zz.sy[,i,,]=aperm(CC,c(3,1,2))
      E.zz.sy[,i,,]=ifelse(is.na(E.zz.sy[,i,,]),eps,E.zz.sy[,i,,])
      E.zz.sy[,i,,]=ifelse(is.infinite(E.zz.sy[,i,,]),1/eps,E.zz.sy[,i,,])
      
    }
    
    
    p.y=(apply(py.s*w,1,sum))
    p.y=ifelse(p.y==0,eps,p.y)
    Ez.y=matrix(0,numobs,r)
    for (i in 1:k) {ps.y[,i]=w[,i]*py.s[,i]/p.y
    Ez.y=Ez.y+matrix(ps.y[,i],numobs,r)*E.z.sy[,i,]
    }
    
    
    
    #### M step estimate w
    
    if (q.w > 1) {
      if (k>1) for (i in 2:k) {
        
        a<-phi[i,]
        for (g in 1:20) {
          A<-pi.greco.hess(a,i,ps.y,x.w,phi)
          b<-pi.greco.grad(a,i,ps.y,x.w,phi)
          
          library(matrixcalc)
          if(is.singular.matrix(A)){
            I=diag(q.w)
            AA=A+eps*I
            solve_A=solve(AA)
          }else(solve_A=solve(A))
          
          a<-a-solve_A%*%t(b)
          a<-c(a)
          
        }
        
        phi[i,]<-a
      }
      
      num=phi%*%t(x.w)
      library(matrixStats)
      
      den=apply(num,2,logSumExp)
      w=apply(num,1, function(x) exp(x-den))
    }
    
    #### M step estimate mu and sigma
    
    temp1=matrix(0,r,r)
    temp2=matrix(0,r,r)
    temp3=matrix(0,r)
    for (i in 1:k) {mu[i,]=colSums(matrix(ps.y[,i],numobs,r)*E.z.sy[,i,])/sum(ps.y[,i])
    sigma[i,,]=apply((array(ps.y[,i],c(numobs,r,r))*(E.zz.sy[,i,,]-aperm(array((mu[i,])%*%t(mu[i,]),c(r,r,numobs)),c(3,1,2)))),c(2,3),sum)/sum(ps.y[,i])
    temp1=temp1+mean(w[,i])*sigma[i,,]
    temp2=temp2+mean(w[,i])*(mu[i,]%*%t(mu[i,]))
    temp3=temp3+mean(w[,i])*mu[i,]}
    
    
    for (j in 1:p) {#opt=optim(alpha[j,],fn=da.max,gr=NULL,method="L-BFGS-B",lower=-25,upper=25,y=y[,j],zex=zex,k=k,ps.y=ps.y,p.z.ys=p.z.ys)
      opt=nlm(f=da.max,p=alpha[j,],stepmax=maxstep,y=y[,j],zex=zex,k=k,ps.y=ps.y,p.z.ys=p.z.ys)
      alpha[j,]=opt$estimate}
    
    ## Identifiability corrections:
    alpha[1,1]=0
    
    var.z=temp1+temp2-temp3%*%t(temp3)
    A<-(chol(var.z))
    for (i in 1:k) {sigma[i,,]<-t(ginv(A))%*%sigma[i,,]%*%ginv(A)
    mu[i,]<-t(ginv(A))%*%mu[i,]
    }
    
    mu.tot=t(colMeans(w))%*%mu
    mu=mu-t(matrix(mu.tot,r,k))
    alpha[,-1]=alpha[,-1]%*%t(A)
    
    
    
    ## Zero constraint:
    if (r>1) for (i in 2:r) alpha[1:(i-1),(i+1)]=0
    
    
    ######################################
    temp<- sum(log(p.y))
    likelihood<-c(likelihood,temp)
    ratio<-(temp-lik)/abs(lik)
    if (hh < 6) ratio=2*eps
    
    lik<-temp
    print(hh)
  }
  
  
  
  
  
  out<-list(alpha=alpha,w=w,mu=mu,sigma=sigma,likelihood=likelihood,ps.y=ps.y,py.s=py.s,Ez.y=Ez.y,py=p.y,phi=phi,q.w=q.w)
  return(out)
}


da.max=function(alpha,y,zex,k,ps.y,p.z.ys)
{
  nq=nrow(zex)
  numobs=length(y)
  temp=matrix(0,numobs)
  for (i in 1:k) {
    t1=t(t(t(y))%*%alpha%*%t(cbind(1,zex[,,i])))
    t2=exp(alpha%*%t(cbind(1,zex[,,i])))
    t2=(matrix(t2,nq,numobs))
    t3=lfactorial(y)
    t3=matrix(t3,nq,numobs,byrow=TRUE)
    log.p.y.z=t1-t2-t3
    temp=temp+ps.y[,i]*colSums(p.z.ys[,,i]*log.p.y.z)
  }
  return(-sum(temp))
}


pi.greco.grad<-function(phi,i,ph.y,x.w,phi.tutte) 
{
  phi.tutte[i,]<-phi
  numobs<-nrow(x.w)
  q.w<-ncol(x.w)
  num1=exp(phi%*%t(x.w))
  num1=ifelse(is.infinite(num1),1/(1e-7),num1)
  den=colSums(exp(phi.tutte%*%t(x.w)))
  out<-colSums(ph.y[,i]*x.w)-(num1/den)%*%x.w
  return(out/numobs)
}

pi.greco.hess<- function(phi,i,ph.y,x.w,phi.tutte) 
{
  phi.tutte[i,]<-phi
  numobs<-nrow(x.w)
  num1=exp(phi%*%t(x.w))
  num1=ifelse(is.infinite(num1),1/(1e-7),num1)
  den=colSums(exp(phi.tutte%*%t(x.w)))
  a<- num1/den
  out<-(t(x.w)%*%x.w)%o%(a^2-a)
  out<-apply(out[,,1,],1,rowMeans)
  return(out/numobs)
}



gauher <- function(n) {# Gauss-Hermite:  returns x,w so that
  #\int_-\infty^\infty exp(-x^2) f(x) dx \doteq \sum w_i f(x_i)
  EPS <- 3.e-14
  PIM4 <- .7511255444649425
  MAXIT <- 10
  m <- trunc((n+1)/2)
  x <- w <- rep(-1,n)
  for (i in 1:m) {
    if (i==1) {
      z <- sqrt(2*n+1)-1.85575*(2*n+1)^(-.16667)
    } else if(i==2) {
      z <- z-1.14*n^.426/z
    } else if (i==3) {
      z <- 1.86*z-.86*x[1]
    } else if (i==4) {
      z <- 1.91*z-.91*x[2]
    } else {
      z <- 2.*z-x[i-2]
    }
    for (its in 1:MAXIT) {
      p1 <- PIM4
      p2 <- 0
      for (j in 1:n) {
        p3 <- p2
        p2 <- p1
        p1 <- z*sqrt(2/j)*p2-sqrt((j-1)/j)*p3
      }
      pp <- sqrt(2*n)*p2
      z1 <- z
      z <- z1-p1/pp
      if(abs(z-z1) <= EPS) break
    }
    x[i] <- z
    x[n+1-i] <- -z
    w[i] <- 2/(pp*pp)
    w[n+1-i] <- w[i]
  }
  list(x=x,w=w)
}


entr<-function(z)
{
  numobs<-nrow(z)
  numg<-ncol(z)
  temp<-0
  z<-ifelse(z==0,z+0.000000000000000000000001,z)
  for (i in 1:numg) for (j in 1:numobs) temp<-temp+(z[j,i]*log(z[j,i]))
  return(-temp)
}

misc=function(classification, truth)
{
  q <- function(map, len, x) {
    x <- as.character(x)
    map <- lapply(map, as.character)
    y <- sapply(map, function(x) x[1])
    best <- y != x
    if (all(len) == 1)
      return(best)
    errmin <- sum(as.numeric(best))
    z <- sapply(map, function(x) x[length(x)])
    mask <- len != 1
    counter <- rep(0, length(len))
    k <- sum(as.numeric(mask))
    j <- 0
    while (y != z) {
      i <- k - j
      m <- mask[i]
      counter[m] <- (counter[m]%%len[m]) + 1
      y[x == name(map)[m]] <- map[[m]][counter[m]]
      temp <- y != x
      err <- sum(as.numeric(temp))
      if (err < errmin) {
        errmin <- err
        best <- temp
      }
      j <- (j + 1)%%k
    }
    best
  }
  if (any(isNA <- is.na(classification))) {
    classification <- as.character(classification)
    nachar <- paste(unique(classification[!isNA]), collapse = "")
    classification[isNA] <- nachar
  }
  MAP <- mapClass(classification, truth)
  len <- sapply(MAP[[1]], length)
  if (all(len) == 1) {
    CtoT <- unlist(MAP[[1]])
    I <- match(as.character(classification), names(CtoT),
               nomatch = 0)
    one <- CtoT[I] != truth
  }
  else {
    one <- q(MAP[[1]], len, truth)
  }
  len <- sapply(MAP[[2]], length)
  if (all(len) == 1) {
    TtoC <- unlist(MAP[[2]])
    I <- match(as.character(truth), names(TtoC), nomatch = 0)
    two <- TtoC[I] != classification
  }
  else {
    two <- q(MAP[[2]], len, classification)
  }
  err <- if (sum(as.numeric(one)) > sum(as.numeric(two)))
    as.vector(one)
  else as.vector(two)
  bad <- seq(along = classification)[err]
  errorRate = length(bad)/length(truth)
  return(errorRate)
}

mapClass=function (a, b)
{
  l <- length(a)
  x <- y <- rep(NA, l)
  if (l != length(b)) {
    warning("unequal lengths")
    return(x)
  }
  aChar <- as.character(a)
  bChar <- as.character(b)
  Tab <- table(a, b)
  Ua <- dimnames(Tab)[[1]]
  Ub <- dimnames(Tab)[[2]]
  aTOb <- rep(list(Ub), length(Ua))
  names(aTOb) <- Ua
  bTOa <- rep(list(Ua), length(Ub))
  names(bTOa) <- Ub
  k <- nrow(Tab)
  Map <- rep(0, k)
  Max <- apply(Tab, 1, max)
  for (i in 1:k) {
    I <- match(Max[i], Tab[i, ], nomatch = 0)
    aTOb[[i]] <- Ub[I]
  }
  if (is.numeric(b))
    aTOb <- lapply(aTOb, as.numeric)
  k <- ncol(Tab)
  Map <- rep(0, k)
  Max <- apply(Tab, 2, max)
  for (j in (1:k)) {
    J <- match(Max[j], Tab[, j])
    bTOa[[j]] <- Ua[J]
  }
  if (is.numeric(a))
    bTOa <- lapply(bTOa, as.numeric)
  list(aTOb = aTOb, bTOa = bTOa)
}


z_ext <-function(x,nfac){
  
  
  ## x    : quadrature points or weights
  ## nfac : number of factors
  
  nq <- length(x)                           # number of quadrature points
  zx <- hcube(rep(nq,nfac))                 # Compute all the possible permutations with repetition of 8 elements taken 3 at a time
  zx <- zx[,dim(zx)[2]:1]                   # zx contains all the "reversed" permutations
  z2 <- matrix(x[zx],dim(zx)[1],dim(zx)[2]) # For each position, I select the corresponding node. 
  #In this way, I obtain all possible combinations of quadrature points for each latent factor
  return(z2)
}


