
# (copied from feature 1.2.13 to avoid tcltk dependency)
#####################################################################
## Matt Wand's version of binned kernel density derivative estimation
##
## Computes the mth derivative of a binned
## d-variate kernel density estimate based
## on grid counts.
#############################################################
#' @importFrom ks binning
drvkde <- function(x,drv,bandwidth,gridsize,range.x,binned=FALSE,se=TRUE, w)
{  
   d <- length(drv)
   if (d==1) x <- as.matrix(x)

   ## Rename common variables
   h <- bandwidth
   tau <- 4 + max(drv)    
   if (length(h)==1) h <- rep(h,d)

   if (missing(gridsize))
     if (!binned)   ## changes 16/02/2009
     {  
       if (d==1) gridsize <- 401
       else if (d==2) gridsize <- rep(151,d)
       else if (d==3) gridsize <- rep(51, d)
       else if (d==4) gridsize <- rep(21, d)
     }
     else
     {
       if (d==1) gridsize <- dim(x)[1]
       else gridsize <- dim(x)
     }

   if(missing(w)) w <- rep(1,nrow(x))
   ## Bin the data if not already binned
  
   if (missing(range.x)) 
   {
     range.x <- list()
     for (id in 1:d)
       range.x[[id]] <- c(min(x[,id])-tau*h[id],max(x[,id])+tau*h[id])  
   }
   
   a <- unlist(lapply(range.x,min))
   b <- unlist(lapply(range.x,max))
   
   M <- gridsize
   gpoints <- list()

   for (id in 1:d)
     gpoints[[id]] <- seq(a[id],b[id],length=M[id])

   if (binned==FALSE)
   {
     if (d==1) gcounts <- binning(x=x, bgridsize=gridsize, h=h, xmin=a, xmax=b, w=w)$counts
     else if (d>1) gcounts <- binning(x=x, bgridsize=gridsize, H=diag(h^2), xmin=a, xmax=b, w=w)$counts
   }
   else
     gcounts <- x

   n <- sum(gcounts)

   kapmid <- list()
   for (id in (1:d))
   {
     ## changes to Lid 13/02/2009
     Lid <- max(min(floor(tau*h[id]*(M[id]-1)/(b[id]-a[id])),M[id]),d)
     lvecid <- (0:Lid)
     facid  <- (b[id]-a[id])/(h[id]*(M[id]-1))
     argid <- lvecid*facid
     kapmid[[id]] <- dnorm(argid)/(h[id]^(drv[id]+1))
     hmold0 <- 1
     hmold1 <- argid
     if (drv[id]==0) hmnew <- 1
     if (drv[id]==1) hmnew <- argid
     if (drv[id] >= 2) 
       for (ihm in (2:drv[id])) 
       {
         hmnew <- argid*hmold1 - (ihm-1)*hmold0
         hmold0 <- hmold1   # Compute drv[id] degree Hermite polynomial
         hmold1 <- hmnew    # by recurrence.
       }
     kapmid[[id]] <- hmnew*kapmid[[id]]*(-1)^drv[id]
   }
  
   if (d==1) kappam <- kapmid[[1]]/n
   if (d==2) kappam <- outer(kapmid[[1]],kapmid[[2]])/n
   if (d==3) kappam <- outer(kapmid[[1]],outer(kapmid[[2]],kapmid[[3]]))/n
   if (d==4) kappam <- outer(kapmid[[1]],outer(kapmid[[2]],outer(kapmid[[3]],kapmid[[4]])))/n

   if (!any(c(d==1,d==2,d==3,d==4))) stop("only for d=1,2,3,4")

   if (d==1) 
   {
     est <- symconv.ks(kappam,gcounts,skewflag=(-1)^drv)
     if (se) est.var <- ((symconv.ks((n*kappam)^2,gcounts)/n) - est^2)/(n-1) 
   }

   if (d==2) 
   {
     ##est <- ks:::symconv.nd(kappam,gcounts,d=d)
     ##if (se) est.var <- ((ks:::symconv.nd((n*kappam)^2,gcounts,d=d)/n) - est^2)/(n-1)
     est <- symconv2D.ks(kappam,gcounts,skewflag=(-1)^drv)
     if (se) est.var <- ((symconv2D.ks((n*kappam)^2,gcounts)/n) - est^2)/(n-1)
   }
     
   if (d==3)
   {
     est <- symconv3D.ks(kappam,gcounts,skewflag=(-1)^drv) 
     if (se) est.var <- ((symconv3D.ks((n*kappam)^2,gcounts)/n) - est^2)/(n-1)
   }
     
   if (d==4)
   {
     est <- symconv4D.ks(kappam,gcounts,skewflag=(-1)^drv) 
     if (se) est.var <- ((symconv4D.ks((n*kappam)^2,gcounts)/n) - est^2)/(n-1) 
   }
   
   if (se)
   {
     est.var[est.var<0] <- 0
     return(list(x.grid=gpoints,est=est,se=sqrt(est.var)))
   }
   else if (!se)
     return(list(x.grid=gpoints,est=est))
}





########################################################################
## Discrete convolution
########################################################################


## Computes the discrete convolution of
## a symmetric or skew-symmetric response 
## vector r and a data vector s.
## If r is symmetric then "skewflag"=1.
## If r is skew-symmetric then "skewflag"=-1.

 
symconv.ks <- function (rr,ss,skewflag = 1) 
{
  L <- length(rr) - 1
  M <- length(ss)
  P <- 2^(ceiling(log(M + L)/log(2)))
  rp <- rep(0,P)
  rp[1:(L+1)] <- rr
  if (L>0) rp[(P-L+1):P] <- skewflag*rr[(L+1):2]
  sp <- rep(0,P) 
  sp[1:M] <- ss
  R <- fft(rp)
  S <- fft(sp)
  t <- fft(R * S, TRUE)
  return((Re(t)/P)[1:M])
}


symconv2D.ks <- function(rr, ss, skewflag=rep(1,2))
{  
  L <- dim(rr)-1
  M <- dim(ss) 
  L1 <- L[1]
  L2 <- L[2]               # find dimensions of r,s
  M1 <- M[1]
  M2 <- M[2]
  P1 <- 2^(ceiling(log(M1+L1)/log(2))) # smallest power of 2 >= M1+L1         
  P2 <- 2^(ceiling(log(M2+L2)/log(2))) # smallest power of 2 >= M2+L2         

  rp <- matrix(0,P1,P2)
  rp[1:(L1+1),1:(L2+1)] <- rr
  if (L1>0)
    rp[(P1-L1+1):P1,1:(L2+1)] <- skewflag[1]*rr[(L1+1):2,]
  if (L2>0)
    rp[1:(L1+1),(P2-L2+1):P2] <- skewflag[2]*rr[,(L2+1):2]
  if (L1 > 0 & L2 > 0)
    rp[(P1-L1+1):P1,(P2-L2+1):P2] <- prod(skewflag)*rr[(L1+1):2,(L2+1):2]   
                                      # wrap around version of rr
  sp <- matrix(0,P1,P2)
  sp[1:M1,1:M2] <- ss                 # zero-padded version of ss

  RR <- fft(rp)        # Obtain FFT's of rr and ss  
  SS <- fft(sp) 
  tt <- fft(RR*SS,TRUE)               # invert element-wise product of FFT's 
  return((Re(tt)/(P1*P2))[1:M1,1:M2]) # return normalized truncated tt
}


symconv3D.ks <- function(rr, ss, skewflag=rep(1,3))
{  
   L <- dim(rr) - 1
   M <- dim(ss) 
   P <- 2^(ceiling(log(M+L)/log(2))) # smallest powers of 2 >= M+L
   L1 <- L[1] ; L2 <- L[2] ; L3 <- L[3]
   M1 <- M[1] ; M2 <- M[2] ; M3 <- M[3]               
   P1 <- P[1] ; P2 <- P[2] ; P3 <- P[3]
   sf <- skewflag

   rp <- array(0,P) 
   rp[1:(L1+1),1:(L2+1),1:(L3+1)] <- rr
   if (L1>0)
     rp[(P1-L1+1):P1,1:(L2+1),1:(L3+1)] <- sf[1]*rr[(L1+1):2,1:(L2+1),1:(L3+1)]
   if (L2>0)
     rp[1:(L1+1),(P2-L2+1):P2,1:(L3+1)] <- sf[2]*rr[1:(L1+1),(L2+1):2,1:(L3+1)]
   if (L3>0)
     rp[1:(L1+1),1:(L2+1),(P3-L3+1):P3] <- sf[3]*rr[1:(L1+1),1:(L2+1),(L3+1):2]
   if (L1>0 & L2>0)
     rp[(P1-L1+1):P1,(P2-L2+1):P2,1:(L3+1)] <- sf[1]*sf[2]*rr[(L1+1):2,(L2+1):2,1:(L3+1)]
   if (L2>0 & L3>0)
     rp[1:(L1+1),(P2-L2+1):P2,(P3-L3+1):P3] <- sf[2]*sf[3]*rr[1:(L1+1),(L2+1):2,(L3+1):2]
   if (L1>0 & L3>0)
     rp[(P1-L1+1):P1,1:(L2+1),(P3-L3+1):P3] <- sf[1]*sf[3]*rr[(L1+1):2,1:(L2+1),(L3+1):2]
   if (L1>0 & L2>0 & L3>0)
     rp[(P1-L1+1):P1,(P2-L2+1):P2,(P3-L3+1):P3] <- sf[1]*sf[2]*sf[3]*rr[(L1+1):2,(L2+1):2,(L3+1):2]

   sp <- array(0,P)
   sp[1:M1,1:M2,1:M3] <- ss            # zero-padded version of ss

   RR <- fft(rp)                       # Obtain FFT's of rr and ss  
   SS <- fft(sp)   
   tt <- fft(RR*SS,TRUE)               # invert element-wise product of FFT's 
   return((Re(tt)/(P1*P2*P3))[1:M1,1:M2,1:M3]) # return normalized truncated tt
}

 
symconv4D.ks <- function(rr, ss, skewflag=rep(1,4) , fftflag=rep(TRUE,2))
{  
   L <- dim(rr) - 1
   M <- dim(ss) 
   P <- 2^(ceiling(log(M+L)/log(2))) # smallest powers of 2 >= M+L
   L1 <- L[1] ; L2 <- L[2] ; L3 <- L[3] ; L4 <- L[4]
   M1 <- M[1] ; M2 <- M[2] ; M3 <- M[3] ; M4 <- M[4]               
   P1 <- P[1] ; P2 <- P[2] ; P3 <- P[3] ; P4 <- P[4] 
   sf <- skewflag

   rp <- array(0,P) 
   rp[1:(L1+1),1:(L2+1),1:(L3+1),1:(L4+1)] <- rr

   if (L1>0)
     rp[(P1-L1+1):P1,1:(L2+1),1:(L3+1),1:(L4+1)] <- sf[1]*rr[(L1+1):2,1:(L2+1),1:(L3+1),1:(L4+1)]
   if (L2>0)
     rp[1:(L1+1),(P2-L2+1):P2,1:(L3+1),1:(L4+1)] <- sf[2]*rr[1:(L1+1),(L2+1):2,1:(L3+1),1:(L4+1)]
   if (L3>0)
     rp[1:(L1+1),1:(L2+1),(P3-L3+1):P3,1:(L4+1)] <- sf[3]*rr[1:(L1+1),1:(L2+1),(L3+1):2,1:(L4+1)]
   if (L4>0)
     rp[1:(L1+1),1:(L2+1),1:(L3+1),(P4-L4+1):P4] <- sf[4]*rr[1:(L1+1),1:(L2+1),1:(L3+1),(L4+1):2]

   if (L1>0 & L2 >0)
     rp[(P1-L1+1):P1,(P2-L2+1):P2,1:(L3+1),1:(L4+1)] <- sf[1]*sf[2]*rr[(L1+1):2,(L2+1):2,1:(L3+1),1:(L4+1)]
   if (L2>0 & L3>0)
     rp[1:(L1+1),(P2-L2+1):P2,(P3-L3+1):P3,1:(L4+1)] <- sf[2]*sf[3]*rr[1:(L1+1),(L2+1):2,(L3+1):2,1:(L4+1)]
   if (L3>0 & L4>0)
     rp[1:(L1+1),1:(L2+1),(P3-L3+1):P3,(P4-L4+1):P4] <- sf[3]*sf[4]*rr[1:(L1+1),1:(L2+1),(L3+1):2,(L4+1):2]
   if (L1>0 & L3>0)
     rp[(P1-L1+1):P1,1:(L2+1),(P3-L3+1):P3,1:(L4+1)] <- sf[1]*sf[3]*rr[(L1+1):2,1:(L2+1),(L3+1):2,1:(L4+1)]
   if (L2>0 & L4>0)
     rp[1:(L1+1),(P2-L2+1):P2,1:(L3+1),(P4-L4+1):P4] <- sf[2]*sf[4]*rr[1:(L1+1),(L2+1):2,1:(L3+1),(L4+1):2]
   if (L1>0 & L4>0)
     rp[(P1-L1+1):P1,1:(L2+1),1:(L3+1),(P4-L4+1):P4] <- sf[1]*sf[4]*rr[(L1+1):2,1:(L2+1),1:(L3+1),(L4+1):2]
   
   if (L1>0 & L2>0 & L3>0)
     rp[(P1-L1+1):P1,(P2-L2+1):P2,(P3-L3+1):P3,1:(L4+1)] <- sf[1]*sf[2]*sf[3]*rr[(L1+1):2,(L2+1):2,(L3+1):2,1:(L4+1)]
   if (L1>0 & L2>0 & L4>0)
     rp[(P1-L1+1):P1,(P2-L2+1):P2,1:(L3+1),(P4-L4+1):P4] <- sf[1]*sf[2]*sf[4]*rr[(L1+1):2,(L2+1):2,1:(L3+1),(L4+1):2]
   if (L2>0 & L3>0 & L4>0)
     rp[1:(L1+1),(P2-L2+1):P2,(P3-L3+1):P3,(P4-L4+1):P4] <- sf[2]*sf[3]*sf[4]*rr[1:(L1+1),(L2+1):2,(L3+1):2,(L4+1):2]

   if (L1>0 & L2>0 & L3>0 & L4>0)
     rp[(P1-L1+1):P1,(P2-L2+1):P2,(P3-L3+1):P3,(P4-L4+1):P4] <- sf[1]*sf[2]*sf[3]*sf[4]*rr[(L1+1):2,(L2+1):2,(L3+1):2,(L4+1):2]
   
   sp <- array(0,P)
   sp[1:M1,1:M2,1:M3,1:M4] <- ss            # zero-padded version of ss

   RR <- fft(rp)                       # Obtain FFT's of rr and ss  
   SS <- fft(sp)   
   tt <- fft(RR*SS,TRUE)               # invert element-wise product of FFT's 
   return((Re(tt)/(P1*P2*P3*P4))[1:M1,1:M2,1:M3,1:M4]) # return normalized truncated tt
}


