#' a function to find interesting features in density (most notably the minimum)
#' as precursor to creating a 1D gate
#' 
#' @author Greg Finak, Phu T. Van
#' @param D a \code{density} containing the data to operate on
#' @param gate_range a \code{character} specifying the data range to operate on
#' @param adjust a \code{numeric} specifying the amount of smoothing to be used
#' @param plot a \code{boolean} specifying whether to output a plot
#' @noRd 
.improvedMindensity <- function(D,adjust=2,gate_range=NA, plot = FALSE, ...){
  # construct the density from data and adjust params we were given
  dens <- density(D,adjust=adjust)
  
  # restrict data to range
  if (!is.null(gate_range)) {
    if (length(gate_range) == 2 & gate_range[1] < gate_range[2]) {
      filter <- dens$x > gate_range[1] & dens$x < gate_range[2]
      dens$x <- dens$x[filter]
      dens$y <- dens$y[filter]
    }else{
      #no range provided, do nothing
    }
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
  
  minima_xcoords <- sp$x[which(minima)] # x-coords of minima
  maxima_xcoords <- sp$x[which(maxima)] # x-coords of maxima
  
  minima_ycoords <- sp$y[which(sp$x %in% minima_xcoords)] # y-coords of minima
  maxima_ycoords <- sp$y[which(sp$x %in% maxima_xcoords)] # y-coords of maxima
  
  #### TO-DO: determine a range in X-coord where the peaks are tiny (< some data-derived threshold)
  # exclude features from this range when considering where to cut. 
  # this is motivated by cyTOF data, where there is a long tail of extreme values that 
  # cannot be eliminated by higher `adjust` values, and shouldn't be eliminated by gate_range()
  # 
  
  if (length(which(minima == TRUE)) == 0){ # no minima found, look through shoulders
    # if there is a peak, pick first shoulder to the right of peak  
    if (length(which(maxima == TRUE)) == 1  ){ 
      pkidx <- which(maxima == TRUE)
      pt <- sp$x[median(which(shoulders)[which(shoulders) > pkidx])]
      if (is.na(pt)){
        # there is only one peak, and there is no shoulder to the right of it
        # pick an inflection point as a last resort
        pt <- sp$x[median(which(inf2)[which(inf2) > pkidx])]
      }
      if( is.na(pt)){
        # There also is no inflection point to the right of it. Resort
        # to the overall median shoulder as the cutpoint like below
        pt <- sp$x[median(which(shoulders))]
      }
      
    }  else {
      # no peak (or multiple peaks), select overall median shoulder as cutpoint (for now...)
      pt <- sp$x[median(which(shoulders))]    
    }
    
    # This is an absolute last resort to make sure this block never returns an NA
    # due to a lack of shoulders and inflection points. In this case just return the
    # simple min
    if( is.na(pt)){
      pt <- sp$x[which.min(sp$y)]
    }
    
  } else if (length(which(minima == TRUE)) > 1) { # multiple minima
    
    m <- min(sp$y[which(sp$x %in% minima_xcoords)])  # pick the minima with lowest y
    pt <- sp$x[which(sp$y == m)]
    
    
  } else if (length(which(minima == TRUE)) == 1){ # only 1 minima, use it as cut point
    pt <- sp$x[which(minima)]
    
  }
  
  .plots = function(){
    abline(v = d2$x[which(inf1)],col="red")       # red    == inflection point
    abline(v = d2$x[which(inf2)],col="green")     # green  == inflection point
    abline(v = d2$x[which(maxima)],col="blue")    # blue   == maxima
    abline(v = d2$x[which(minima)],col="orange")  # orange == minima
    abline(v = d2$x[which(shoulders)],col="pink") # pink   == shoulders
  }
  if(plot){
    par(mfrow=c(2,1))
    plot(sp,type="l",main="features")
    .plots()
    plot(sp,type="l",main="final_cut")
    abline(v = pt, col="black", lwd=2)
    par(mfrow=c(1,1)) # be nice, restore plot settings
  }    
  
  return(list(density = sp,
              inf_rising = sp$x[which(inf1)],
              inf_falling = sp$x[which(inf2)],
              maxima = maxima_xcoords,
              minima = minima_xcoords,
              maxima_heights = maxima_ycoords,
              minima_heights = minima_ycoords,
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
#' @name gate_mindensity2
#' @aliases mindensity2
#' @param fr a \code{flowFrame} object
#' @param channel the channel to operate on
#' @param filterId a name to refer to this filter
#' @param gate_range numeric vector of length 2. If given, this sets the bounds
#' on the gate applied. 
#' @param min a numeric value that sets the lower boundary for data filtering
#' @param max a numeric value that sets the upper boundary for data filtering
#' @param peaks \code{numeric} vector. If not given , then perform peak detection first by .find_peaks
#' @param ... Additional arguments for peak detection.
#' @return a \code{rectangleGate} object based on the minimum density cutpoint
#' @examples
#' \dontrun{
#'  gate <- gate_mindensity2(fr, channel = "APC-A") # fr is a flowFrame
#' }
#' @export
gate_mindensity2 <- function(fr, channel, filterId = "", 
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
  g <- .improvedMindensity(D=x,gate_range=gate_range,...)
  
  coords <- list(c(g$final_cut, Inf))
  # name the gate and return it
  names(coords) <- channel
  return(rectangleGate(coords, filterId = filterId)) 
}

#' @export
mindensity2 <- function(fr, channel, filterId = "", 
                        gate_range = NULL, min = NULL, max = NULL, peaks = NULL, 
                        ...){
  .Deprecated("gate_mindensity2")
  gate_mindensity2(fr, channel, filterId, gate_range, min, max, peaks, ...)
}

#' wrapper for mindensity2
#' 
#' It does some parameter preprocessing before calling the mindensity
#' 

#' @param ... arguments to be passed to \link{mindensity}
#' @inheritParams .gate_flowclust_1d 
#' 
#' @return a \code{filter} object
#' @noRd 
.gate_mindensity2 <- function(fr, pp_res, channels, ...) {
  
  if(length(channels) != 1)
    stop("invalid number of channels for mindensity!")
  gate <- gate_mindensity2(fr, channel = channels, ...)
  gate
}

.mindensity2 <- .gate_mindensity2
