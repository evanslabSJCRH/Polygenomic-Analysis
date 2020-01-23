#This is the core code to compute the TAP test P value for one gene.

#Suppose M genomic features are annotated to the gene, the association of each feature with the phenotype of interested has been tested by a proper statistical model, producing a length-M vector of P values, ‘pvals’. Suppose the association tests are also performed on permuted data in B round of permutations, producing a M-by-B matrix of P values, ‘pvals.perm.mat’. The sequence of calling the functions in the core code is as follows:

#Tobs = TP.Fisher.stat(pvals,delta=0.05)
#T.perm = rep(NA,B)
#for(j in 1:B) T.perm[j] = TP.Fisher.stat(pvals.perm.mat[,j],delta=0.05)
## The TAP test P value:#This is the core code to compute the TAP test P value for one gene.

#Suppose M genomic features are annotated to the gene, the association of each feature with the phenotype of interested has been tested by a proper statistical model, producing a length-M vector of P values, ‘pvals’. Suppose the association tests are also performed on permuted data in B round of permutations, producing a M-by-B matrix of P values, ‘pvals.perm.mat’. The sequence of calling the functions in the core code is as follows:

#Tobs = TP.Fisher.stat(pvals,delta=0.05)
#T.perm = rep(NA,B)
#for(j in 1:B) T.perm[j] = TP.Fisher.stat(pvals.perm.mat[,j],delta=0.05)
## The TAP test P value:
#P.TAP = 1- TPcdf0.Fisher(Tobs,T.perm,taup=0.9,M,delta=0.05,pii0=NA,tol=0.02,maxiter=100,plt=F,prt=F)



#===========================================================
# The truncated P values test statistic with Fisher's transformation
#=============================================================
TP.Fisher.stat <- function(Pvs,delta=0.2) {
  PP <- Pvs[Pvs<=delta]
  PP[PP<1e-78]=1e-78
  TP <- 0
  if(length(PP)>0) TP <- sum(-log(PP))
  return(TP)
}

#==============================================================
# Estimator of null cdf from nonparametric modification of the
# template cdf by the estimated copmarison distribution (Fhat.VD)
# Fisher's transformation
#=================================================================
TPcdf0.Fisher <- function(tgrid,T0.dat,taup=0.9,M,delta=0.2,pii0=NA,tol=0.02,maxiter=100,plt=F,prt=F) {
  T0.dat <- T0.dat[!is.na(T0.dat)]
  if(length(T0.dat)<30) {
    cat("** Too few non-missing values for estiamting the null cdf\n")
    return(NULL)
  }
  gm <- 1
  Qhat.taup <- K <- NA
  nt <- length(tgrid)
  F0 <- rep(0,nt)
  indx0 <- (1:nt)[tgrid<=0]
  indx1 <- (1:nt)[tgrid>0]
  if(length(indx0)>0) F0[indx0] <- 0
  if(length(indx1)>0)
  {
    tgrid1 <- tgrid[indx1]
    if(!is.na(taup)) {
      # To circumvent the possibility that the right tail of the template cdf does not match 
      # well with the null-hypothesis test stat EDF, find a scalar 'gm' that matches the 
      # upper 'taup'-th quantile (taup=0.8 by defualt):
      # Bernstein polynomial estimator of the tau.p-th quantile:
      #K <- as.integer(length(T0.dat)^(1/3)/sqrt(dnorm(qnorm(taup))))+1
      K <- length(T0.dat)
      Qhat.taup <- sum(dbinom((0:K),K,taup)*EQF(T0.dat,(0:K)/K)) 
      #Qhat.taup <- EQF(T0.dat,taup)
      # Binary search to find the taup-th quantile of the template distribution
      gotQ <- F
      T0lo <- 0; T0hi <- max(T0.dat)
      T0mid <- NA
      if(TPcdf0.Fisher.temp(T0hi,gm=1,M,delta,pii0=pii0)<taup) {T0mid <- T0hi; gotQ <- T}
      if(TPcdf0.Fisher.temp(0,gm=1,M,delta,pii0=pii0)>taup) {T0mid <- Qhat.taup; gotQ <- T}
      if(prt) cat("M =",M,"| initial T0Lo, T0hi =",T0lo,T0hi,"initial T0mid =",T0mid,"\n")
      iter <- 0
      while(!gotQ) {
        T0mid <- (T0lo+T0hi)*0.5
        prmid <- TPcdf0.Fisher.temp(T0mid,gm=1,M,delta,pii0=pii0)
        if(prt) cat("T0lo =",T0lo,"T0hi =",T0hi,"T0mid =",T0mid,"M =",M,"delta =",delta," prmid =",prmid,"\n")
        if(abs(prmid-taup)<=tol) gotQ <- T
        if(prmid>taup) T0hi <- T0mid
        if(prmid<taup) T0lo <- T0mid
        iter <- iter+1
        if(iter>maxiter) break
      }
      gm <- max(Qhat.taup/T0mid,1)
    } 
    if(prt) cat("** gamma =",gm,"taup =",taup,"**\n")
    if (is.na(gm)) gm <- 1  #set gm=1 when gm is NaN
    W <- TPcdf0.Fisher.temp(T0.dat,gm=gm,M,delta,pii0=pii0)
    # Append 0 and 1:
    W <- c(0,W,1)
    F00 <- TPcdf0.Fisher.temp(tgrid1,gm=gm,M,delta,pii0=pii0,echo=prt)
    if(prt) cat("in TPcdf0.Fisher, pii0 =",pii0,"\n")
    iknots <- sort(c(seq(from=1e-8,to=0.1,length=10),
                     EQF(W,(1:9)/10),
                     seq(from=0.8,to=1-1e-8,length=20)))
    if(prt) {cat("in TPcdf0.Fisher, summary of W","\n"); print(summary(W))}
    if(prt) {cat("in TPcdf0.Fisher, summary of F00","\n"); 
      if(length(F00)==1) print(F00) else print(summary(F00))
    }
    F0[indx1] <- TP.Compcdf.VD(F00,W,3,iknots,plt,prt)
  }
  return(list(FTP0=F0,tau=taup,gamma=gm,Q.tau=Qhat.taup,K=K))
}

#======================================================================
# Parametric template for the null cdf of TP (the truncated P values statistic)
# Fisher's transformation
#==========================================================================
TPcdf0.Fisher.temp <- function(tgrid,gm=1,M,delta=0.2,pii0=NA,echo=F) {
  q.delta <- -log(delta)
  F0 <- dbinom(0,M,delta)*as.integer(tgrid/gm>=0)
  weight <- 1
  if(!is.na(pii0)) {
    F0 <- pii0
    weight <- (1-pii0)/(1-(1-delta)^M)
  }
  if(echo) {
    cat("** In TPcdf0.Fisher.temp --------****\n")
    cat("  gm =",gm,"delta =",delta,"M =",M,"pii0 =",pii0,"\n")
    cat("  init F0 =",F0,"weight =",weight,"\n")
    cat("-----------------------------------\n")
  }
  F0a <- 0
  for(j in 1:M) {
    #pj <- min(pgamma(j*q.delta/gm,j,1),1-1e-10)
    pj <- min(pgamma(j*q.delta,j,1),1-1e-10)
    F0a <- F0a+dbinom(j,M,delta)*as.integer(tgrid/gm>=j*q.delta)*(pgamma(tgrid/gm,j,1)-pj)/(1-pj)
  }
  F0 <- F0+weight*F0a
  #return(F0)
  return(pmin(F0,1-1e-16))
}

#============================================================
# VD-spline estimator of the copmarison distribution
#==============================================================
TP.Compcdf.VD <- function(w,W.dat,degree=3,iknots,plt=T,prt=T)
{
  nW <- length(W.dat)
  k <- degree+1 #order of spline = degree + 1
  knotsW <- sort(unique(c(0,iknots,1)))
  nknotsW <- length(knotsW)
# To be conformal with de Boor's notation:
  n <- nknotsW+k-2  # so n = number of unique interior knotsW + order of spline
  knotsW.ext <- c(rep(0,k-1),knotsW,rep(1,k-1))
  tsW <- apply(matrix(1:n,n,1),1,
           function(i,tk,k) {
             return(sum(tk[(i+1):(i+k-1)])/(k-1))
           },knotsW.ext,k
         )
  edfW <- EDF(W.dat,tsW)
  #edfW <- pmax(EDF(W.dat,tsW),tsW)
  FW.VD <- Bspline(w,edfW,knotsW[2:(nknotsW-1)],degree,a=0,b=1)
  #FW.VD[w>1-1e-8] <- 1
  if(prt) {
    cat("nknotsW =",nknotsW,"\n")
    cat("knotsW =",knotsW,"\n")
    cat("tsW =",round(tsW,4),"\n")
    cat("edfW =",round(edfW,4),"\n")
    #print(cbind(w,FW.VD))
  }
  if(plt) {
    X11(); par(mfrow=c(2,1),lab=c(10,10,7))
    u <- (0:nW)/nW
    FW.est <- Bspline(u,edfW,knotsW[2:(nknotsW-1)],degree,a=0,b=1)
    plot(u,FW.est,type="l",xlab="w",ylab="F~ and F^ of W",xlim=c(0,1),ylim=c(0,1),lwd=2)
    lines(u,EDF(W.dat,u),lty=2,col=3)
    lines(u,u,lty=3,col=2)
    u <- seq(from=0.99999,to=1,length=10000)
    FW.est <- Bspline(u,edfW,knotsW[2:(nknotsW-1)],degree,a=0,b=1)
    plot(u,FW.est,type="l",xlab="w",ylab="F~ and F^ of W",xlim=c(0.99999,1),ylim=c(min(c(0.99999,FW.est)),1),lwd=2)
    lines(u,EDF(W.dat,u),lty=2,col=3)
    lines(u,u,lty=3,col=2)
  }
  if(prt) {
    cat("in TP.Compcdf.VD, FW.VD =\n"); print(cbind(w,FW.VD))
  }
  return(FW.VD)
}

#=============================================================
# Low-level functions:
#=============================================================
#==========================================================
# The empirical quantile function
#==========================================================
EQF <- function(X,u) {
  Xodr <- sort(X)
  n <- length(X)
  EQF <- Xodr[ceiling(n*u)]
  i0 <- (1:length(u))[u==0]
  if(length(i0)>0) EQF <- c(rep(Xodr[1],length(i0)),EQF)
  return(EQF)
}
#========================================================
# The eimprical distribution function
# Input: X -- Data (missing values are deleted)
#        x1 -- grid over which to evaluate the EDF
# Output: vector of EDF values over x1
#========================================================
EDF <- function(X,x1) {
  Xodr <- sort(X)
  n <- length(Xodr)
  if(n<=1) return(NA)
  z <- cumsum(table(Xodr))/n
  uniqueXodr <- as.double(names(z))
  LB <- uniqueXodr[1]-(uniqueXodr[2]-uniqueXodr[1])
  y <- approx(c(LB,uniqueXodr),c(0,z),x1,method="constant",
         rule=2,f=0)
  return(y[[2]])
}
#==================================================================
# Compute a B-spline approximation of a continuous function over a
# closed interval [a,b].
# Input:
#   x -- a grid over which to evaluate the approximation
#   fv -- a vector of sampled function values
#   knots.int -- the interior knot sequence.
#            m=length(knots).
#   degree -- spline degree (default 3). degree=order-1
#   a -- left boundary of interval (default=0)
#   b -- right boundary of interval (default=1)
#   The length of fv must be equal to m+degree+1, i.e., the number of
#   interior knots plus order of the spline; otherwise NULL is returned.
# Output: a vector of the spline values over x.
#====================================================================
Bspline <- function(x,fv,knots.int,degree=3,a=0,b=1) {
  m <- length(knots.int)
  nf <- length(fv)
  if(nf!=m+degree+1) return(NULL)
  Smat <- bsmat(x,knots.int=knots.int,degree=degree,
            a=a,b=b)
  func.val <- Smat%*%fv
  return(func.val)
}
#=====================================================================
#  Compute B-spline matrix evaluated overa grid of values on a closed
#  interval [a,b].
# Input:
#   x -- grid
#   knots.int -- the interior knot sequance
#   degree -- spline degree (default 3). degree=order-1
#   a -- left boundary of interval (default=0)
#   b -- right boundary of interval (default=1)
#======================================================================
bsmat <- function(x,knots.int,degree=3,a=0,b=1) {
# To be comformal with de Boor's notation:
  k <- degree+1  # order of spline
# Number of interior knots:
  m <- length(knots.int)
# Number of Splines:
  n <- m+k
  d <- length(x)
  smat1 <- smat2 <- matrix(0,d,n)
# The initial (i.e., k=1) B spline matrix
  kt <- unique(sort(c(a,knots.int,b))) # initial extended knot seq
  el <- 1; nel <- m+el
  indx <- numeric(0)
  for(j in 1:nel) {
    indx <- (1:d)[x>=kt[j]&x<kt[j+1]]
    smat1[indx,j] <- 1
  }
# Make the last piece (order-1 B spline) left-continuous at b:
  indx <- (1:d)[x==b]
  smat1[indx,nel] <- 1
  if(k==1) return(smat1)
# Start the iteration, done if el=k:
  while(el<k) {
    el <- el+1
    nel <- m+el
    kt <- c(a,kt,b)
    for(j in 2:(nel-1)) {
      smat2[,j] <- ((x-kt[j])/(kt[j+el-1]-kt[j]))*smat1[,j-1]+
                   ((kt[j+el]-x)/(kt[j+el]-kt[j+1]))*smat1[,j]
    }
    smat2[,1] <- ((kt[1+el]-x)/(kt[1+el]-kt[2]))*smat1[,1]
    smat2[,nel] <- ((x-kt[nel])/(kt[nel+el-1]-kt[nel]))*smat1[,nel-1]
    smat1 <- smat2
  }
  return(smat1)
}

#P.TAP = 1- TPcdf0.Fisher(Tobs,T.perm,taup=0.9,M,delta=0.05,pii0=NA,tol=0.02,maxiter=100,plt=F,prt=F)
