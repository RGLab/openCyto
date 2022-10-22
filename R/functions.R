#' reading flowFrame from csv spreadsheet.It is for internal usage.
#' @param file \code{character} csv file name that contains the fluorescence intensities
#' @param stains \code{character} a vector specifying the stains for channels. Default is \code{NA}
#' @return a \code{flowFrame} object
#' @importFrom flowCore parameters<-
#' @noRd 
.read.FCS.csv <- function(file, stains = NA) {
  mat <- as.matrix(read.csv(file, check.names = FALSE))
  
  fr <- new("flowFrame", exprs = mat)
  
  pd <- pData(parameters(fr))
  pd$desc <- as.character(pd$desc)
  pd$name <- as.character(pd$name)

  ## update the desc with marker name
  if (!is.na(stains)) {
    ind <- match(names(stains), pd$name)
    pd[ind, ]$desc <- as.character(stains)

    ## update SSC and FSC description with NA
    ind <- grepl("[F|S]SC", pd$desc)
    pd[ind, ]$desc <- NA
  } else pd$desc <- NA
  
  ## update minRange with -111 for proper display of the data
  pd$minRange[pd$minRange < (-111)] <- -111
  pData(parameters(fr)) <- pd
  fr
}
#' reading flowSet from multiple csv spreadsheets.It is for internal usage.
#' @param files \code{character} csv file names
#' @param ... arguments passed to \link{read.FCS.csv}
#' @return a \code{flowSet} object
#' @noRd 
.read.flowSet.csv <- function(files, ...) {
  fs <- flowSet(lapply(files, .read.FCS.csv, ...))
  sampleNames(fs) <- basename(files)
  fs
}

#' Creates a matrix of points on an ellipse from a fitted flowClust model.
#'
#' The ellipse is constructed from a contour from the fitted flowClust model.
#'
#' By default, the contour level is extracted from the \code{ruleOutliers} slot
#' in the \code{filter}. The user can override this level by specifying value
#' between 0 and 1 in \code{quantile}.
#'
#' @param filter object containing the fitted \code{flowClust} model.
#' @param include the mixture component in the fitted \code{flowClust} model for
#' which the contour (ellipse) is returned
#' @param ecol \code{numeric} to be documented
#' @param elty \code{numeric} to be documented
#' @param quantile the contour level of the ellipse. See details.
#' @param npoints the number of points on the ellipse
#' @param subset the dimensions of the mixture component to return
#' @param \code{...} additional parameters
#' @return matrix containing the points of the ellipse from the flowClust contour
#' @importFrom flowClust rbox
#' @noRd 
.getEllipse <- function(filter = NULL, include = seq_len(filter@K), ecol = 1, elty = 1, 
  quantile = NULL, npoints = 50, subset = c(1, 2),...) {
  .Defunct(".getEllipseGate")
  # Sets the quantile of the ellipse.
  if (is.null(quantile)) {
    quantile <- filter@ruleOutliers[2]
  } else {
    if (!is.numeric(quantile) || quantile < 0 || quantile > 1) {
      stop("The 'quantile' must be a numeric value between 0 and 1.")
    }
  }
  
  # py is the degrees of freedom?
  py <- 2
  ecol <- matrix(ecol, length(include))
  elty <- matrix(elty, length(include))
  
  if (all(filter@nu != Inf)) {
    if (filter@ruleOutliers[1] == 0) {
      # 0 means quantile
      cc <- py * qf(p = quantile, py, filter@nu)
    } else {
      # 1 means u.cutoff
      cc <- ((filter@nu + py)/quantile - filter@nu)
    }
  } else {
    cc <- qchisq(p = quantile, py)
  }
  
  j <- 0
  
  #Does trans exist in the extra parameter list?
  #If not, set it to true by default
  ellipsis<-as.environment(list(...))
  if(exists("trans",envir=ellipsis)){
    trans<-get("trans",ellipsis)
  }else{
    trans<-1
  }

  #Test for trans==0 when lambda is defined to get around the off 
  #by one bug due to the reverse box-cox transformation
  if ((length(filter@lambda) > 0)&&trans==0) {
    lambda <- rep(filter@lambda, length.out = filter@K)
  } else {
    lambda <- numeric(0)
  }
  cc <- rep(cc, length.out = filter@K)
  for (i in include) {
    eigenPair <- eigen(filter@sigma[i, subset, subset])
    l1 <- sqrt(eigenPair$values[1]) * sqrt(cc)
    l2 <- sqrt(eigenPair$values[2]) * sqrt(cc)
    angle <- atan(eigenPair$vectors[2, 1]/eigenPair$vectors[1, 1]) * 180/pi
    
    if ((length(lambda) > 0)&trans==1) {
      res <- rbox(flowClust:::.ellipsePoints(a = l1[i], b = l2[i], alpha = angle, 
        loc = filter@mu[i, subset], n = npoints), lambda[i])
    } else {
      res <- flowClust:::.ellipsePoints(a = l1[i], b = l2[i], alpha = angle, 
        loc = filter@mu[i, subset], n = npoints)
    }
  }
  res
}
#' revised based on .getEllipse.
#' try to construct ellipsoidGate directly from tmixFilter result instead of polygon fitting
#' when trans = 0  
#' @noRd 
.getEllipseGate <- function(filter = NULL, include = seq_len(filter@K), ecol = 1, elty = 1, 
    quantile = NULL, npoints = 50, subset = c(1, 2),...) {
  
  # Sets the quantile of the ellipse.
  if (is.null(quantile)) {
    quantile <- filter@ruleOutliers[2]
  } else {
    if (!is.numeric(quantile) || quantile < 0 || quantile > 1) {
      stop("The 'quantile' must be a numeric value between 0 and 1.")
    }
  }
  
  # py is the degrees of freedom?
  py <- 2
  ecol <- matrix(ecol, length(include))
  elty <- matrix(elty, length(include))
  
  if (all(filter@nu != Inf)) {
    if (filter@ruleOutliers[1] == 0) {
      # 0 means quantile
      cc <- py * qf(p = quantile, py, filter@nu)
    } else {
      # 1 means u.cutoff
      cc <- ((filter@nu + py)/quantile - filter@nu)
    }
  } else {
    cc <- qchisq(p = quantile, py)
  }
  
  j <- 0
  
  #Does trans exist in the extra parameter list?
  #If not, set it to true by default
  ellipsis<-as.environment(list(...))
  if(exists("trans",envir=ellipsis)){
    trans<-get("trans",ellipsis)
  }else{
    trans<-1
  }
  
  #Test for trans==0 when lambda is defined to get around the off 
  #by one bug due to the reverse box-cox transformation
  if ((length(filter@lambda) > 0)&&trans==0) {
    lambda <- rep(filter@lambda, length.out = filter@K)
  } else {
    lambda <- numeric(0)
  }
  cc <- rep(cc, length.out = filter@K)
  cc <- sqrt(cc)
  coln <- as.vector(filter@varNames)
  for (i in include) {
    cov.mat <- filter@sigma[i,subset, subset]
    
    #fit polygon points 
    if ((length(lambda) > 0)&trans==1) {
      eigenPair <- eigen(cov.mat)
      l1 <- sqrt(eigenPair$values[1]) * cc 
      l2 <- sqrt(eigenPair$values[2]) * cc
      angle <- atan(eigenPair$vectors[2, 1]/eigenPair$vectors[1, 1]) * 180/pi
      
      contour_ellipse <- rbox(flowClust:::.ellipsePoints(a = l1[i], b = l2[i], alpha = angle, 
              loc = filter@mu[i, subset], n = npoints), lambda[i])
      
      res <- polygonGate(.gate = matrix(contour_ellipse, ncol = 2,
                        dimnames = list(NULL, coln)),
                      filterId = filter@filterId)
    } else {
      #construct ellipsoidGate directly from covaraince matrix and mu
      dimnames(cov.mat) <- list(coln, coln)
      res <- ellipsoidGate(cov.mat, mean = filter@mu[i, subset], distance = cc[1])
      
    }
  }
  res
}

#' Removes any observation from the given flowFrame object that has values
#' outside the given range for the specified channels
#'
#' The minimum/maximum values are ignored if \code{NULL}.
#'
#' @param flow_frame a \code{flowFrame} object
#' @param min a numeric vector that sets the lower bounds for data filtering
#' @param max a numeric vector that sets the upper bounds for data filtering
#' @param channels \code{character} specifying which channel to operate on
#' @return a \code{flowFrame} object
#' @examples
#' \dontrun{
#'  library(flowClust)
#'  data(rituximab)
#'  # Consider the range of values for FSC.H and SSC.H
#'  summary(rituximab)
#' 
#'  # Truncates any observations with FSC.H outside [100, 950]
#'  rituximab2 <- .truncate_flowframe(rituximab, channels = "FSC.H", min = 100, max = 950)
#'  summary(rituximab2)
#'  # Next, truncates any observations with SSC.H outside [50, 1000]
#'  rituximab3 <- .truncate_flowframe(rituximab2, channels = "SSC.H", min = 50, max = 1000)
#'  summary(rituximab3)
#'
#'  # Instead, truncates both channels at the same time
#'  rituximab4 <- .truncate_flowframe(rituximab, channels = c("FSC.H", "SSC.H"),
#'  min = c(100, 50), max = c(950, 1000))
#'  summary(rituximab4)
#' }
#' 
#' @noRd 
.truncate_flowframe <- function(flow_frame, channels, min = NULL, max = NULL) {
  channels <- as.character(channels)
  num_channels <- length(channels)

  # For comparison purposes, we update the min and max values to -Inf and Inf,
  # respectively, if NULL.
  if (is.null(min)) {
    min <- rep(-Inf, num_channels)
  }
  if (is.null(max)) {
    max <- rep(Inf, num_channels)
  }

  if (!(num_channels == length(min) && num_channels == length(max))) {
    stop("The lengths of 'min' and 'max' must match the number of 'channels' given.")
  }
 
  gate_coordinates <- lapply(seq_len(num_channels), function(i) {
    c(min[i], max[i])
  })
  names(gate_coordinates) <- channels

  truncate_filter <- rectangleGate(gate_coordinates)
  Subset(flow_frame, truncate_filter)
}

#' Removes any observation from each flowFrame object within the flowSet that has
#' values outside the given range for the specified channels
#'
#' The minimum/maximum values are ignored if \code{NULL}.
#'
#' @param flow_set a \code{flowSet} object
#' @inheritParams .truncate_flowframe
#' @return a \code{flowSet} object
#' @noRd 
.truncate_flowset <- function(flow_set, channels, min = NULL, max = NULL) {
  channels <- as.character(channels)
  num_channels <- length(channels)

  # For comparison purposes, we update the min and max values to -Inf and Inf,
  # respectively, if NULL.
  if (is.null(min)) {
    min <- rep(-Inf, num_channels)
  }
  if (is.null(max)) {
    max <- rep(Inf, num_channels)
  }

  if (!(num_channels == length(min) && num_channels == length(max))) {
    stop("The lengths of 'min' and 'max' must match the number of 'channels' given.")
  }
 
  gate_coordinates <- lapply(seq_len(num_channels), function(i) {
    c(min[i], max[i])
  })
  names(gate_coordinates) <- channels

  truncate_filter <- rectangleGate(gate_coordinates)
  Subset(flow_set, truncate_filter)
}


#' Computes the quantile from flowClust for a given vector of probabilties
#'
#' We estimate the quantile from a \code{flowClust} fit with a combination of
#' numerical integration and a root-finding method. We are effectively
#' estimating the cumulative distribution function (CDF) of the mixture density
#' estimated by \code{flowClust}.
#'
#' Because we are using numerical methods, we also need an \code{interval} of
#' values in which we will attempt to find the specified quantile.
#'
#' @param p vector of probabilities
#' @param object an object containing the \code{flowClust} fit
#' @param interval a vector of length 2 containing the end-points of the interval
#' of values to find the quantile
#' @param ... Additional arguments that are passed to \code{uniroot} to find the
#' quantile.
#' @return the quantile corresponding to the specified probabilities
#' @noRd 
.quantile_flowClust <- function(p, object, interval, ...) {
  cdf_target <- function(x, p, object) {
    cdf_values <- sapply(seq_len(object@K), function(k) {
      nu <- ifelse(length(object@nu) == 1, object@nu, object@nu[k])
      lambda <- ifelse(length(object@lambda) == 1, object@lambda, object@lambda[k])
      
      # TODO: Incorporate the Box-Cox transformation (i.e., box(qt(...), lambda =
      # lambda)) into quantile The case of 'lambda = 1' may be not be trivial -- this
      # case is largely ignored in flowClust.
      pt((x - object@mu[k])/sqrt(object@sigma[k]), df = nu)
    })
    weighted.mean(cdf_values, w = object@w) - p
  }
  
  uniroot(cdf_target, interval = interval, p = p, object = object, ...)$root
}

#' Extracts the quadrants of a quadGate as a list of rectangleGates
#'
#' The quadrants are numbered in a clockwise manner with the top-left quadrant
#' numbered 1, the top-right quadrant numbered 2, and so on.
#'
#' @param quad_gate a \code{quadGate} object
#' @param markers character vector of the marker names for the x- and y-axes
#' @param channels character vector of the channel names for the x- and y-axes
#' @param quadrants a vector indicating the quadrants to extract
#' @return a \code{filters} object containing a list of the rectangle gates
#' @noRd 
.quadGate2rectangleGates <- function(quad_gate, markers, channels, quadrants = 1:4) {
  x_gate <- quad_gate@boundary[1]
  y_gate <- quad_gate@boundary[2]
  
  gates_list <- list()
  
  # Top-left quadrant
  gates <- list(c(-Inf, x_gate), c(y_gate, Inf))
  names(gates) <- channels
  gates_list[[paste0(markers[1], "-", markers[2], "+")]] <- rectangleGate(gates)
  
  # Top-right quadrant
  gates <- list(c(x_gate, Inf), c(y_gate, Inf))
  names(gates) <- channels
  gates_list[[paste0(markers[1], "+", markers[2], "+")]] <- rectangleGate(gates)
  
  # Lower-right quadrant
  gates <- list(c(x_gate, Inf), c(-Inf, y_gate))
  names(gates) <- channels
  gates_list[[paste0(markers[1], "+", markers[2], "-")]] <- rectangleGate(gates)
  
  # Lower-left quadrant
  gates <- list(c(-Inf, x_gate), c(-Inf, y_gate))
  names(gates) <- channels
  gates_list[[paste0(markers[1], "-", markers[2], "-")]] <- rectangleGate(gates)
  
  filters(gates_list[quadrants])
}

#' Constructs a vector of all the combinations of A & B & C
#'
#' The \code{permutations} function is from the \code{gregmisc} package on CRAN.
#' @param markers character vector of marker names
#' @return vector containing all combinations of the markers
#' @examples
#' .polyfunction_nodes(c('IFNg', 'IL2', 'TNFa', 'GzB', 'CD57'))
#' @noRd 
.polyfunction_nodes <- function(markers) {
  
  markers <- paste0(markers, "+")
  num_markers <- length(markers)
  and_list <- as.vector(permutations(n = 1, r = num_markers - 1, c("&"), repeats.allowed = TRUE))
  isnot_list <- permutations(n = 2, r = num_markers, c("!", ""), repeats.allowed = TRUE)
  apply(isnot_list, 1, function(isnot_row) {
    isnot_row[-1] <- paste0(and_list, isnot_row[-1])
    paste(paste0(isnot_row, markers), collapse = "")
  })
} 

#' Finds values of a vector in an interval
#'
#' @param x numeric vector
#' @param interval numeric vector of length 2
#' @return numeric vector containing the values of \code{x} between
#' \code{interval}. If no values are found, \code{NA} is returned.
#' @examples
#' z <- seq.int(1, 9, by = 2)
#' .between_interval(z, interval = c(2, 8))
#' @noRd 
.between_interval <- function(x, interval) {
  x <- x[findInterval(x, interval) == 1]
  if (length(x) == 0) {
    x <- NA
  }
  x
}

#' Constructs a 2x2 rotation matrix for a given angle
#' @param theta \code{numeric} the degree of rotation that ensures the angle between the x-axis and the eigenvector is between 0 and pi
#' @noRd 
.rotation_matrix <- function(theta) {
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
}

#' Validate the arguments to flowClust and tell the user early if there are problems
#' For now we make sure that the K argument to flowClust and to prior_flowClust are present and agree.
#' @param df a \code{data.frame} of the template
#' @return a \code{data.frame}
#' @noRd 
.validateFlowClustArgs <- function(df){
  inds<-which(df$gating_method%like%"flowClust")
  for(i in inds){
    fc_args<-.argParser(df[i,"gating_args"],TRUE)
    prepro_args<-.argParser(df[i,"preprocessing_args"],TRUE)
    if(length(fc_args)==0){
      #default K
      fc_args$K<-2
    }
    
    if(exists("K",fc_args)){
      K<-fc_args$K
    }

    if(length(prepro_args)==0){
      prepro_args$K<-K 
    }else if(exists("K",prepro_args)){
      if(!(prepro_args$K==K)){
        warning("K in preprocessing_args doesn't match K in gating_args for flowClust.\nOverridiing with K from flowClust.")
      }
      prepro_args$K<-K
    }
    df[i,"preprocessing_args"] <- gsub("\\)$","",gsub("list\\(","",deparse(prepro_args,control=c())))
    df[i,"gating_args"] <- gsub("\\)$","",gsub("list\\(","",deparse(fc_args,control=c())))
  }
  df
}

#' Centers a vector of data using the mode of the kernel density estimate
#'
#' @param x numeric vector
#' @param ... additional arguments passed to \code{\link{density}}
#' @return numeric vector containing the centered data
#' @noRd 
.center_mode <- function(x, ...) {
  x <- as.vector(x)
  density_x <- density(x, ...)
  mode <- density_x$x[which.max(density_x$y)]
  
  x <- as.vector(scale(x, center = mode, scale = FALSE))
  attributes(x) <- list(`mode` = mode)
  x
}

#' Scales a vector of data using the Huber robust estimator for mean and
#' standard deviation
#'
#' This function is an analog to \code{\link{scale}} but using Huber robust
#' estimators instead of the usual sample mean and standard deviation.
#'
#' @param x numeric vector
#' @param center logical value. Should \code{x} be centered?
#' @param scale logical value. Should \code{x} be scaled?
#' @return numeric vector containing the scaled data
#' @noRd 
.scale_huber <- function(x, center = TRUE, scale = TRUE) {
  
  
  x <- as.vector(x)
  sd <- mad(x)
  huber_x <- robust_m_estimator(x, sd)
  
  # If 'center' is set to TRUE, we center 'x' by the Huber robust location
  # estimator.
  center_x <- FALSE
  if (center) {
    center_x <- huber_x
  }
  
  # If 'scale' is set to TRUE, we scale 'x' by the Huber robust standard
  # deviation estimator.
  scale_x <- FALSE
  if (scale) {
    scale_x <- sd
  }
  
  x <- as.vector(base::scale(x, center = center_x, scale = scale_x))
  
  if (!center) {
    center_x <- NULL
  }
  if (!scale) {
    scale_x <- NULL
  }
  attributes(x) <- list(center = center_x, scale = scale_x)
  x
}

#' Standardizes a channel within a \code{flowSet} object using the mode of the
#' kernel density estimate and the Huber estimator of the standard deviation
#' 
#' @param fs a \code{flowSet} object
#' @param channel the channel to standardize
#' @param data \code{logical} indicating whether to return the transformed flow data.
#' 
#' @return the \code{transformation} list, center, scale and optionally the transformed \code{flowFrame}
#' @noRd 
.standardize_flowFrame <- function(fr, channel, data = TRUE) {
  
  x <- exprs(fr)[, channel]
  
  if (length(x) >= 2) {
    # First, centers the values by the mode of the kernel density estimate.
    x <- .center_mode(x)
    mode <- attr(x, "mode")
    
    # Scales the marker cells by the Huber estimator of the standard deviation.
    x <- .scale_huber(x, center = FALSE)
    sd_huber <- attr(x, "scale")
    
    exprs(fr)[, channel] <- x
  } else {
    mode <- NA
    sd_huber <- NA
  }
  res <- list(center = mode, scale = sd_huber)
  
  if(data)
    res$flow_frame <- fr
  
  attr(res, "openCyto_preprocessing") <- "standardize"  
  
  
  res
  
}
#' convert gate to a filterResult
#' 
#' used for computing the gate indices on the fly
#' 
#' @param fr flowFrame
#' @param channel \code{character} channel to used for gate
#' @param gate \code{filter} object
#' @param positive \code{logical}
#' @noRd 
.gateToFilterResult <- function(fr, channel, gate, positive){

  x <- exprs(fr)[, channel]
  gate_coordinates <- c(gate@min, gate@max)
  cutpoint <- gate_coordinates[!is.infinite(gate_coordinates)]
  nCount <- length(x)
  #deal with dummy gate where both boundaries are Inf
  if(length(cutpoint) == 0){
    if(gate@min < gate@max){
      ind <- rep(positive, nCount)  
    }else{
      ind <- rep(!positive, nCount)
    }
  }
  if(positive)
    ind <- x >= cutpoint
  else
    ind <- x < cutpoint
#  fres <- as(ind, "filterResult") 
#  filterDetails(fres, identifier(gate)) <- gate
  #pack into bit vector to ease the traffic
  
  fres <- as(gate, "ocRectangleGate")
  fres@ind <- ncdfFlow:::toBitVec(ind) 
  fres  
}


