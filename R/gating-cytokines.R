#' Constructs cytokine gates from the derivative of a kernel density estimate
#' after standardizing and collapsing flowFrames
#'
#' @param fr a \code{flowFrame} object
#' @param channel the channel from which the cytokine gate is constructed
#' @param filter_id the name of the filter
#' @param num_peaks the number of peaks expected to see. This effectively removes
#' any peaks that are artifacts of smoothing
#' @param ref_peak After \code{num_peaks} are found, this argument provides the
#' index of the reference population from which a gate will be obtained.
#' @param tol the tolerance value used to construct the cytokine gate from the
#' derivative of the kernel density estimate
#' @param positive If \code{TRUE}, then the gate consists of the entire real
#' line to the right of the cutpoint. Otherwise, the gate is the entire real
#' line to the left of the cutpoint. (Default: \code{TRUE})
#' @param ... additional arguments passed to \code{\link{cytokine_cutpoint}}
#' @return a \code{filterList} containing the gates (cutpoints) for each sample
#' @export
cytokine <- function(fr, channel, filter_id = "", num_peaks = 1,
                     ref_peak = 1, tol = 1e-2, positive = TRUE, ...) {
  # Standardizes the flowFrame's for a given channel using the mode of the kernel
  # density estimate and the Huber estimator of the standard deviation
  standardize_out <- standardize_flowset(as(fr,"flowSet"), channel = channel)

  # Coerces the standardized flowSet into a single flowFrame from which a single
  # cutpoint is calculated using the first derivative of the kernel density
  # estimate. 
  cutpoint <- cytokine_cutpoint(flow_frame = as(standardize_out$flow_set, "flowFrame"),
                                channel = channel, num_peaks = num_peaks,
                                ref_peak = ref_peak, tol = tol, ...)

  # Backtransforms the cutpoint with respect to each sample to place
  # them on the scales of the original samples
  cutpoints <- lapply(standardize_out$transformation, function(transf_i) {
    with(transf_i, center + scale * cutpoint)
  })

  # If a sample has no more than 1 observation when the 'cytokine' gate is
  # attempted, the 'center' and/or 'scale' will result be NA, in which case we
  # replace the resulting NA cutpoints with the average of the remaining
  # cutpoints. If all of the cutpoints are NA, we set the mean to 0, so that
  # all of the cutpoints are 0.
  cutpoints_unlisted <- unlist(cutpoints)
  if (sum(!is.na(cutpoints_unlisted)) > 0) {
    mean_cutpoints <- mean(cutpoints_unlisted, na.rm = TRUE)
  } else {
    mean_cutpoints <- 0
  }
  cutpoints <- as.list(replace(cutpoints_unlisted, is.na(cutpoints_unlisted),
                               mean_cutpoints))

  # Creates a list of filters for each set of cutpoints.
  # Note that the gate consists of the entire real line to the right of the
  # cutpoint.
  cytokine_gates <- lapply(cutpoints, function(cutpoint) {
    # After the 1D cutpoint is set, we set the gate coordinates used in the
    # rectangleGate that is returned. If the `positive` argument is set to TRUE,
    # then the gate consists of the entire real line to the right of the cut point.
    # Otherwise, the gate is the entire real line to the left of the cut point.
    if (positive) {
      gate_coordinates <- list(c(cutpoint, Inf))
    } else {
      gate_coordinates <- list(c(-Inf, cutpoint))
    }
    names(gate_coordinates) <- channel
    rectangleGate(gate_coordinates, filterId = filter_id)
  })
  
  cytokine_gates[[1]]
}

#' Standardizes a channel within a \code{flowSet} object using the mode of the
#' kernel density estimate and the Huber estimator of the standard deviation
#' 
#' @param flow_set a \code{flowSet} object
#' @param channel the channel to standardize
#' @return list containing the transformed \code{flowSet} object along with the
#' \code{transformation} list, where each element contains the transformation
#' parameters for each \code{flowFrame}
#' @export
standardize_flowset <- function(flow_set, channel = "FSC-A") {
  transform_out <- fsApply(flow_set, function(flow_frame) {
    x <- exprs(flow_frame)[, channel]

    if (length(x) >= 2) {
      # First, centers the values by the mode of the kernel density estimate.
      x <- center_mode(x)
      mode <- attr(x, "mode")
  
      # Scales the marker cells by the Huber estimator of the standard deviation.
      x <- scale_huber(x, center = FALSE)
      sd_huber <- attr(x, "scale")

      exprs(flow_frame)[, channel] <- x
    } else {
      mode <- NA
      sd_huber <- NA
    }

    list(flow_frame = flow_frame, center = mode, scale = sd_huber)
  })

  # Creates a flowSet object from the transformed flowFrame objects
  flow_set <- flowSet(lapply(transform_out, function(x) x$flow_frame))

  # Extracts the transformation parameters
  transformation <- lapply(transform_out, function(x) {
    x$flow_frame <- NULL
    x
  })

  list(flow_set = flow_set, transformation = transformation)
}

#' Constructs a cutpoint for a flowFrame by using a derivative of the kernel
#' density estimate
#'
#' We determine a gating cutpoint using either the first or second derivative of
#' the kernel density estimate (KDE) of the \code{channel} specified within the
#' \code{flow_frame}.
#'
#' By default, we compute the first derivative of the kernel density estimate
#' from the channel specified within the given \code{flow_frame}. Next, we
#' determine the lowest valley from the derivative, which corresponds to the
#' density's mode for cytokines. We then contruct a gating cutpoint as the value
#' less than the tolerance value \code{tol} in magnitude and is also greater
#' than the lowest valley.
#'
#' Alternatively, if the \code{method} is selected as \code{second_deriv}, we
#' select a cutpoint from the second derivative of the KDE. Specifically, we
#' choose the cutpoint as the largest peak of the second derivative of the KDE
#' density which is greater than the reference peak.
#'
#' @param flow_frame a \code{flowFrame} object
#' @param channel the channel name
#' @param num_peaks the number of peaks expected to see. This effectively removes
#' any peaks that are artifacts of smoothing
#' @param ref_peak After \code{num_peaks} are found, this argument provides the
#' index of the reference population from which a gate will be obtained. By
#' default, the peak farthest to the left is used.
#' @param method the method used to select the cutpoint. See details.
#' @param tol the tolerance value
#' @param adjust the scaling adjustment applied to the bandwidth used in the
#' first derivative of the kernel density estimate
#' @param ... additional arguments passed to \code{\link{deriv_density}}
#' @return the cutpoint along the x-axis
#' @export
cytokine_cutpoint <- function(flow_frame, channel, num_peaks = 1, ref_peak = 1,
                              method = c("first_deriv", "second_deriv"),
                              tol = 1e-2, adjust = 1, ...) {

  method <- match.arg(method)

  x <- as.vector(exprs(flow_frame)[, channel])
  peaks <- sort(find_peaks(x, num_peaks = num_peaks, adjust = adjust))

  if (ref_peak > num_peaks) {
    warning("The reference peak is larger than the number of peaks found.",
            "Setting the reference peak to 'num_peaks'...",
            call. = FALSE)
    ref_peak <- num_peaks
  }

  # TODO: Double-check that a cutpoint minimum found via 'first_deriv'
  # passes the second-derivative test.

  if (method == "first_deriv") {
    # Finds the deepest valleys from the kernel density and sorts them.
    # The number of valleys identified is determined by 'num_peaks'
    deriv_out <- deriv_density(x = x, adjust = adjust, deriv = 1, ...)

    deriv_valleys <- with(deriv_out, find_valleys(x = x, y = y, adjust = adjust))
    deriv_valleys <- deriv_valleys[deriv_valleys > peaks[ref_peak]]
    deriv_valleys <- sort(deriv_valleys)[1]
    
    cutpoint <- with(deriv_out, x[x > deriv_valleys & abs(y) < tol])
    cutpoint <- cutpoint[1]
  } else {
    # The cutpoint is selected as the first peak from the second derivative
    # density which is to the right of the reference peak.
    deriv_out <- deriv_density(x = x, adjust = adjust, deriv = 2, ...)
    deriv_peaks <- with(deriv_out, find_peaks(x, y, adjust = adjust))
    deriv_peaks <- deriv_peaks[deriv_peaks > peaks[ref_peak]]
    cutpoint <- sort(deriv_peaks)[1]
  }
  
  cutpoint
}

#' Constructs the derivative specified of the kernel density estimate of a
#' numeric vector
#'
#' The derivative is computed with \code{\link[ks:drvkde]{drvkde}}.
#'
#' For guidance on selecting the bandwidth, see this CrossValidated post:
#' \url{http://bit.ly/12LkJWz}
#' 
#' @param x numeric vector
#' @param deriv a numeric value specifying which derivative should be calculated.
#' By default, the first derivative is computed.
#' @param bandwidth the bandwidth to use in the kernel density estimate. If
#' \code{NULL} (default), the bandwidth is estimated using the plug-in estimate
#' from \code{\link[ks]{hpi}}.
#' @param adjust a numeric weight on the automatic bandwidth, analogous to the
#' \code{adjust} parameter in \code{\link{density}}
#' @param num_points the length of the derivative of the kernel density estimate
#' @param ... additional arguments passed to \code{\link[ks:drvkde]{drvkde}}
#' @return list containing the derivative of the kernel density estimate
#' @export
#' @importFrom ks hpi,drvkde
deriv_density <- function(x, deriv = 1, bandwidth = NULL, adjust = 1,
                          num_points = 10000, ...) {
#  require('feature')
#  require('ks')
  if (is.null(bandwidth)) {
    bandwidth <- hpi(x, deriv.order = deriv)
  }
  deriv_x <- drvkde(x = x, drv = deriv, bandwidth = adjust * bandwidth,
                    gridsize = num_points, ...)
  list(x = deriv_x$x.grid[[1]], y = deriv_x$est)
}


#' Centers a vector of data using the mode of the kernel density estimate
#'
#' @param x numeric vector
#' @param ... additional arguments passed to \code{\link{density}}
#' @return numeric vector containing the centered data
#' @export
center_mode <- function(x, ...) {
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
#' @export
scale_huber <- function(x, center = TRUE, scale = TRUE) {
  require('MASS')

  x <- as.vector(x)
  huber_x <- huber(x)

  # If 'center' is set to TRUE, we center 'x' by the Huber robust location
  # estimator.
  center_x <- FALSE
  if (center) {
    center_x <- huber_x$mu
  }

  # If 'scale' is set to TRUE, we scale 'x' by the Huber robust standard
  # deviation estimator.
  scale_x <- FALSE
  if (scale) {
    scale_x <- huber_x$s
  }

  x <- as.vector(base:::scale(x, center = center_x, scale = scale_x))

  if (!center) {
    center_x <- NULL
  }
  if (!scale) {
    scale_x <- NULL
  }
  attributes(x) <- list(center = center_x, scale = scale_x)
  x
}
