#' a function to find interesting features in density (most notably the minimum)
#' as precursor to creating a 1D gate
#' 
#' @param D data to operate on
#' @param adjust smoothing for building the density
#' @param range data range to search for the minimum, if NA all data is used
#' @author Greg Finak, Phu T. Van
.improvedMindensity <- function(D,adjust=2,range=NA, plot = FALSE, ...){
  # construct the density from data and adjust params we were given
  dens <- density(D,adjust=adjust)
  
  # restrict data to range
  if(!is.na(range)&length(range)==2 & range[1]<range[2]){
    filter <- dens$x>range[1]&dens$x<range[2]
    dens$x <- dens$x[filter]
    dens$y <- dens$y[filter]
  }else{
    #no range provided, do nothing
  }
  
  sp <- smooth.spline(dens$x,dens$y)
  pred <- predict(sp)
  
  d1 <- predict(sp,deriv = 1)
  d2 <- predict(sp,deriv = 2)
  d3 <- predict(sp,deriv = 3)
  
  # find features
  inf1 <- sign(d2$y[-length(d2$y)])>0&sign(d2$y[-1L])<0
  inf2 <- sign(d2$y[-length(d2$y)])<0&sign(d2$y[-1L])>0
  minima <- sign(d1$y[-1])>0&sign(d1$y[-length(d1$y)])<0
  maxima <- sign(d1$y[-1])<0&sign(d1$y[-length(d1$y)])>0
  shoulders <- sign(d3$y[-1])<0&sign(d3$y[-length(d3$y)])>0
  
  
  if (length(which(minima == TRUE)) == 0){ # no minima found, look through shoulders
    pt <- sp$x[median(which(shoulders))] # pick the median shoulder
    
  } else if (length(which(minima == TRUE)) > 1) { # multiple minima, filter them
    pt <- sp$x[min(which(minima))]
    
  } else if (length(which(minima == TRUE)) == 1){ # only 1 minima, use it as cut point
    pt <- sp$x[which(minima)]
    
  }
  
  .plots = function(){
    abline(v = d2$x[which(inf1)],col="red")
    abline(v = d2$x[which(inf2)],col="green")
    abline(v = d2$x[which(maxima)],col="blue")
    abline(v = d2$x[which(minima)],col="orange")
    abline(v = d2$x[which(shoulders)],col="pink")
  }
  if(plot){
    par(mfrow=c(2,1))
    plot(sp,type="l",main="features")
    .plots()
    plot(sp,type="l",main="final_cut")
    abline(v = pt, col="black", lwd=2)  
  }    
  

 return(list(density = sp,
              inf_rising = sp$x[which(inf1)],
              inf_falling = sp$x[which(inf2)],
              maxima = sp$x[which(maxima)],
              minima = sp$x[which(minima)],
              shoulders = sp$x[which(shoulders)], 
              final_cut = pt))
  
}

#' An improved version of mindensity used to determines a cutpoint as the minimum point 
#' of a kernel density estimate between two peaks.
#' 
#' Analogous to the original openCyto::mindensity(), mindensity2 operates on a standard flowFrame. Its behavior is closely modeled on 
#' the original mindensity() whenever possible. However, the underlying peak-finding algorithm 
#' (improvedMindensity) behaves significantly differently. 
#' 
#' @author Greg Finak, Phu T. Van
#' 
#' @param fr a \code{flowFrame} object
#' @param channel the channel to operate on
#' @param filter_id a name to refer to this filter
#' @param positive If \code{TRUE}, then the gate consists of the entire real
#' line to the right of the cutpoint. Otherwise, the gate is the entire real
#' line to the left of the cutpoint. (Default: \code{TRUE})
#' @param pivot logical value. If \code{TRUE}, we choose as the two peaks the
#' largest peak and its neighboring peak. See details.
#' @param gate_range numeric vector of length 2. If given, this sets the bounds
#' on the gate applied. If no gate is found within this range, we set the gate to
#' the minimum value within this range if \code{positive} is \code{TRUE} and the
#' maximum value of the range otherwise.
#' @param min a numeric value that sets the lower boundary for data filtering
#' @param max a numeric value that sets the upper boundary for data filtering
#' @param peaks \code{numeric} vector. If not given , then perform peak detection first by .find_peaks
#' @param ... Additional arguments for peak detection.
#' @return a \code{rectangleGate} object based on the minimum density cutpoint
#' @export
#' @examples
#' \dontrun{
#'  gate <- mindensity2(fr, channel = "APC-A") # fr is a flowFrame
#' }
#' 
mindensity2 <- function(fr, channel, filter_id = "", positive = TRUE, pivot = FALSE, 
                         gate_range = NULL, min = NULL, max = NULL, peaks = NULL, 
                         ...) {
  if (missing(channel) || length(channel) != 1) {
    stop("A single channel must be specified.")
  }
  # if there is 
  if (!(is.null(min) && is.null(max))) {
    fr <- .truncate_flowframe(fr, channels = channel, min = min, 
                              max = max)
  }
  x <- exprs(fr[,channel])
  # this line does the actual searching of the cut point
  g <- .improvedMindensity(D=x,...)
  
  # get the cut point
  if (positive) {
    coords <- list(c(g$final_cut, Inf))
  }
  else {
    coords <- list(c(-Inf, g$final_cut))
  }
  # name the gate and return it
  names(coords) <- channel
  return(rectangleGate(coords, filterId = filter_id)) 
}



#' wrapper for mindensity2
#' 
#' It does some parameter preprocessing before calling the mindensity
#' 

#' @param ... arguments to be passed to \link{mindensity}
#' @inheritParams .flowClust.1d 
#' 
#' @return a \code{filter} object
.mindensity2 <- function(fr, pp_res, channels, positive = TRUE, ...) {
  
  if(length(channels) != 1)
    stop("invalid number of channels for mindensity!")
  gate <- mindensity2(fr, channel = channels, positive = positive, ...)
  gate
}