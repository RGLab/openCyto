
# (copied from feature 1.2.13 to avoid tcltk dependency)
#####################################################################
## Matt Wand's version of binned kernel density derivative estimation
##
## Computes the mth derivative of a binned
## d-variate kernel density estimate based
## on grid counts.
#############################################################
#' @importFrom ks binning
drvkde <- function(x,drv,bandwidth,gridsize,range.x,se=TRUE, w)
{  
   d <- length(drv)
   stopifnot (d==1) 
   x <- as.matrix(x)

   ## Rename common variables
   h <- bandwidth
   tau <- 4 + drv    
   if (length(h)==1) h <- rep(h,d)

   if (missing(gridsize))
    gridsize <- 401
     
     

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

   gcounts <- binning(x=x, bgridsize=gridsize, h=h, xmin=a, xmax=b, w=w)$counts
   
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
  
    kappam <- kapmid[[1]]/n
   
 
   
  est <- symconv.ks(kappam,gcounts,skewflag=(-1)^drv)
  if (se) est.var <- ((symconv.ks((n*kappam)^2,gcounts)/n) - est^2)/(n-1) 
   if (se)
   {
     est.var[est.var<0] <- 0
     return(list(x.grid=gpoints,est=est,se=sqrt(est.var)))
   }
   
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

