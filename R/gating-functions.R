#' Applies flowClust to 1 feature to determine a cutpoint between the minimum
#' cluster and all other clusters.
#'
#' We cluster the observations in \code{fr} into \code{K} clusters.
#'
#' By default, the cutpoint is chosen to be the boundary of the first two
#' clusters. That is, between the first two cluster centroids, we find the
#' midpoint between the largest observation from the first cluster and the
#' smallest observations from the second cluster. Alternatively, if the
#' \code{cutpoint_method} is \code{min_density}, then the cutpoint is the point
#' at which the density between the first and second smallest cluster centroids
#' is minimum.
#'
#' @param fr a \code{flowFrame} object
#' @param params TODO
#' @param filterId TODO
#' @param K the number of clusters to find
#' @param usePrior Should we use the Bayesian version of \code{flowClust}?
#' Answers are "yes", "no", or "vague". The answer is passed along to
#' \code{flowClust}.
#' @param prior_list list of prior parameters for the Bayesian \code{flowClust}.
#' If \code{usePrior} is set to 'no', then the list is unused.
#' @param trans numeric indicating whether the Box-Cox transformation parameter
#' is estimated from the data. May take 0 (no estimation), 1 (estimation) or 2
#' (cluster-speciﬁc estimation). NOTE: For the Bayesian version of
#' \code{flowClust}, this value cannot be 2.
#' @param positive If \code{TRUE}, then the gate consists of the entire real
#' line to the right of the cutpoint. Otherwise, the gate is the entire real
#' line to the left of the cutpoint. (Default: \code{TRUE})
#' @param cutpoint_method How should the cutpoint be chosen from the fitted
#' \code{flowClust} model? See Details.
#' @param neg_cluster integer. The index of the negative cluster. The cutpoint
#' is computed between clusters \code{neg_cluster} and \code{neg_cluster + 1}.
#' @param truncate_min Truncate observations less than this minimum value. By
#' default, this value is \code{NULL} and is ignored.
#' @param truncate_max Truncate observations greater than this maximum value. By
#' default, this value is \code{NULL} and is ignored.
#' @param quantile the quantile for which we will find the cutpoint using
#' the quantile \code{cutpoint_method}. If the \code{cutpoint_method} is not set
#' to \code{quantile}, this argument is ignored.
#' @param quantile_interval a vector of length 2 containing the end-points of the
#' interval of values to find the quantile cutpoint. If the
#' \code{cutpoint_method} is not setto \code{quantile}, this argument is ignored.
#' @param ... additional arguments that are passed to \code{flowClust}
#' @return a \code{rectangleGate} object ranging consisting of all values greater
#' than the cutpoint determined
flowClust.1d <- function(fr, params, filterId = "", K = 2,  adjust = 1, trans = 0,
                         positive = TRUE, plot = FALSE, usePrior = 'no', prior = NULL,
                         cutpoint_method = c("boundary", "min_density", "quantile", "posterior_mean"),
                         neg_cluster = 1, truncate_min = NULL, truncate_max = NULL,
                         quantile = 0.99, quantile_interval = c(0, 10),...) {

  cutpoint_method <- match.arg(cutpoint_method)

  # If appropriate, we generate prior parameters for the Bayesian version of flowClust.
  if (usePrior == "yes" && is.null(prior)) {
    prior <- prior_flowClust1d(fr = fr, channel = params[1], K = K, adjust = adjust)
  }

  # HACK: Circumvents a bug in flowClust.
  if (missing(prior) || is.null(prior)) {
    prior <- list(NA)
  }

  # If a truncation value is specified, we remove all observations less than
  # (greater than) this value for the marker specified to construct the gate.
  # NOTE: These observations are removed from the 'flowFrame' locally and are gated
  # out only for the determining the gate.
  if (!is.null(truncate_min) || !is.null(truncate_max)) {
    fr <- truncate_flowframe(flow_frame = fr, channel = params[1],
                             min = truncate_min, max = truncate_max)
  }

  # Applies `flowClust` to the feature specified in the `params` argument using
  # the data given in `fr`. We use priors with hyperparameters given by the
  # elements in the list `prior_list`.
  filter1 <- tmixFilter(filterId, params[1], K = K, trans = trans,
                        usePrior = usePrior, prior = prior, ...)
  tmixRes1 <- filter(fr, filter1)
  
  # To determine the cutpoint, we first sort the centroids so that we can determine
  # the second largest centroid.
  centroids_sorted <- sort(getEstimates(tmixRes1)$locations)

  # Also, because the cluster labels are arbitrary, we determine the cluster
  # the cluster labels, sorted by the ordering of the cluster means.
  labels_sorted <- order(getEstimates(tmixRes1)$locations)

  # Grabs the data matrix that is being gated.
  x <- exprs(fr)[, params[1]]

  # Determines the cutpoint between clusters 1 and 2.
  if (cutpoint_method == "boundary") {
    # Choose the cutpoint as the boundary between the first two clusters.
    # First, we sort the data.
    order_x <- order(x)
    x_sorted <- x[order_x]
    labels <- flowClust::Map(tmixRes1, rm.outliers = FALSE)[order_x]

    # Determine which observations are between the first two centroids and their
    # corresponding cluster labels.
    which_between <- which(centroids_sorted[neg_cluster] < x_sorted & x_sorted < centroids_sorted[neg_cluster + 1])
    x_between <- x_sorted[which_between]
    labels_between <- labels[which_between]

    # For the observations between the first two centroids, we find the last
    # observation that belongs to the first cluster and the first observation
    # that belongs to the second cluster. In the rare occurrence that no
    # observations from one of the clusters is between the labels, we set the
    # max/min index as NA. This results in the cutpoint being set to the max/min
    # observation on the boundary.
    which_cluster1 <- which(labels_between == labels_sorted[neg_cluster])
    which_cluster2 <- which(labels_between == labels_sorted[neg_cluster + 1])
    max_obs_cluster1 <- ifelse(length(which_cluster1) == 0, NA, max(which_cluster1))
    min_obs_cluster2 <- ifelse(length(which_cluster2) == 0, NA, min(which_cluster2))

    # We define the cutpoint to be the midpoint between the two clusters.
    cutpoint <- mean(x_between[c(max_obs_cluster1, min_obs_cluster2)], na.rm = TRUE)
    if (is.nan(cutpoint)) {
      cutpoint <- centroids_sorted[neg_cluster]
    }
  } else if (cutpoint_method == "min_density") {
    # Determine the minimum density value of the observations between clusters 1 and 2.
    # Sets the cutpoint at the observation that attained this minimum density.
    x_between <- x[centroids_sorted[neg_cluster] < x & x < centroids_sorted[neg_cluster + 1]]
    x_dens <- dmvtmix(x = as.matrix(x_between), object = tmixRes1)
    cutpoint <- x_between[which.min(x_dens)]
  } else if (cutpoint_method == "quantile") {
    cutpoint <- quantile_flowClust(p = quantile, object = tmixRes1, interval = quantile_interval)
  } else { # cutpoint_method == "posterior_mean"
    cutpoint <- centroids_sorted[neg_cluster]
  }

  
  # After the 1D cutpoint is set, we set the gate coordinates used in the
  # rectangleGate that is returned. If the `positive` argument is set to TRUE,
  # then the gate consists of the entire real line to the right of the cut point.
  # Otherwise, the gate is the entire real line to the left of the cut point.
  if (positive) {
    gate_coordinates <- list(c(cutpoint, Inf))
  } else {
    gate_coordinates <- list(c(-Inf, cutpoint))
  }

  names(gate_coordinates) <- params
  
  fres <- rectangleGate(gate_coordinates, filterId = filterId)

  # Saves posterior point estimates
  postList <- list()
  posteriors <- list(mu = tmixRes1@mu, lamdda = tmixRes1@lambda,
                     sigma = tmixRes1@sigma, nu = tmixRes1@nu, min = min(x),
                     max = max(x))
  postList[[params[1]]] <- posteriors

  # Saves prior point estimates
  priorList <- list()
  priorList[[params[1]]] <- prior
  
  if (plot) {
    gate_pct <- round(100 * mean(x > cutpoint), 3)
    plot_title <- paste0(filterId, " (", gate_pct, "%)")
    plot(fr, tmixRes1, main = plot_title, labels = FALSE)
    abline(v = centroids_sorted, col = rainbow(K))
    abline(v = cutpoint, col = "black", lwd = 3, lty = 2)

    if (!is.null(prior)) {
      x_dens <- seq(min(x), max(x), length = 1000)

      for(k in seq_len(K)) {
        prior_density <- with(prior, dnorm(x_dens, mean = Mu0[k], sd = sqrt(Omega0[k])))

        # Grab posterior estimates for the degrees of freedom (nu) and the
        # transformation parameters. Because these can either be of length 1 or
        # of length K, we grab the appropriate posterior estimates for the
        # current mixture component.
        nu <- ifelse(length(tmixRes1@nu) == 1, tmixRes1@nu, tmixRes1@nu[k])
        lambda <- ifelse(length(tmixRes1@lambda) == 1, tmixRes1@lambda, tmixRes1@lambda[k])        
        posterior_density <- dmvt(x_dens, mu = tmixRes1@mu[k,],
                                  sigma = tmixRes1@sigma[k,,], nu = nu,
                                  lambda = lambda)$value
        lines(x_dens, prior_density, col = rainbow(K)[k], lty = 2, lwd = 1)
        lines(x_dens, posterior_density, col = rainbow(K)[k], lwd = 1)
      }
    }
  }

  fcRectangleGate(fres, priorList, postList)
}

#' Applies flowClust to two features in a flowFrame
#'
#' We cluster the observations in \code{fr} into \code{K} clusters. We set the
#' cutpoint to be the point at which the density between the first and second
#' smallest cluster centroids is minimum.
#'
#' The cluster for the population of interest is selected as the one with cluster
#' centroid nearest the \code{target} in Euclidean distance.
#'
#' When constructing the axis gate, we translate the gate from the cluster
#' centroid of along the eigenvector with positive slope from the bottom-left
#' corner to top-right corner. The \code{axis_translation} argument scales the
#' translation shift as a function of the appropriate chi-squared coefficient.
#' The larger \code{axis_translation} is, the more gate is shifted in a positive
#' direction.
#'
#' @param fr a \code{flowFrame} object
#' @param xChannel TODO
#' @param yChannel TODO
#' @param filterId TODO
#' @param K the number of clusters to find
#' @param usePrior Should we use the Bayesian version of \code{flowClust}?
#' Answers are "yes", "no", or "vague". The answer is passed along to
#' \code{flowClust}.
#' @param prior list of prior parameters for the Bayesian \code{flowClust}.
#' If \code{usePrior} is set to 'no', then the list is unused.
#' @param trans numeric indicating whether the Box-Cox transformation parameter
#' is estimated from the data. May take 0 (no estimation), 1 (estimation) or 2
#' (cluster-speciﬁc estimation). NOTE: For the Bayesian version of
#' \code{flowClust}, this value cannot be 2.
#' @param plot a logical value indicating if the fitted mixture model should be
#' plotted. By default, no.
#' @param target a numeric vector of length \code{K} containing the location of
#' the cluster of interest. See details.
#' @param gate_type character value specifying the type of gate to construct from
#' the \code{flowClust} fit
#' @param quantile the contour level of the target cluster from the
#' \code{flowClust} fit to construct the gate
#' @param axis_translation a numeric value between 0 and 1 used to position the
#' gate if \code{gate_type} is selected as \code{"axis"}. See details.
#' @param truncate_min A vector of length 2. Truncate observations less than this
#' minimum value. The first value truncates the \code{xChannel}, and the second
#' value truncates the \code{yChannel}. By default, this vector is \code{NULL}
#' and is ignored.
#' @param truncate_max A vector of length 2. Truncate observations greater than
#' this maximum value. The first value truncates the \code{xChannel}, and the
#' second value truncates the \code{yChannel}. By default, this vector is
#' \code{NULL} and is ignored.
#' @param ... additional arguments that are passed to \code{flowClust}
#' @return a \code{polygonGate} object containing the contour (ellipse) for 2D
#' gating.
flowClust.2d <- function(fr, xChannel, yChannel, filterId = "", K = 2,
                         usePrior = 'no', prior = list(NA), trans = 0,
                         plot = FALSE, target = rep(0, K),
                         gate_type = c("ellipse", "axis"), quantile = 0.9,
                         axis_translation = 0.25, truncate_min = NULL,
                         truncate_max = NULL, ...) {

  if (length(target) != 2) {
    stop("The 'target' location must be a numeric vector of length 2.")
  }
  gate_type <- match.arg(gate_type)

  # If a truncation value is specified, we remove all observations less than this
  # value for the marker specified to construct the gate.
  # NOTE: These observations are removed from the 'flowFrame' locally and are gated
  # out only for the determining the gate.
  if (!is.null(truncate_min)) {
    exprs(fr) <- exprs(fr)[exprs(fr)[, xChannel] >= truncate_min[1], ]
    exprs(fr) <- exprs(fr)[exprs(fr)[, yChannel] >= truncate_min[2], ]
  }

  if (!is.null(truncate_max)) {
    exprs(fr) <- exprs(fr)[exprs(fr)[, xChannel] <= truncate_max[1], ]
    exprs(fr) <- exprs(fr)[exprs(fr)[, yChannel] <= truncate_max[2], ]    
  }

  x <- exprs(fr)[, xChannel]
  y <- exprs(fr)[, yChannel]

  # If appropriate, we generate prior parameters for the Bayesian version of flowClust.
  if (usePrior == "yes" && identical(prior, list(NA))) {
    prior <- prior_flowClust(fr = fr, channels = c(xChannel, yChannel), K = K)
  }

  # Applies `flowClust` to the feature specified in the `params` argument using
  # the data given in `fr`. We use priors with hyperparameters given by the
  # elements in the list `prior`.
  filter1 <- tmixFilter(filterId, c(xChannel, yChannel), K = K, trans = trans,
                        usePrior = usePrior, prior = prior, ...)
  tmixRes1 <- filter(fr, filter1)

  # Converts the tmixFilterResult object to a polygonGate.
  # We select the cluster with the minimum 'yChannel' to be the subpopulation from
  # which we obtain the contour (ellipse) to generate the polygon gate.
  fitted_means <- getEstimates(tmixRes1)$locations
  target_dist <- apply(fitted_means, 1, function(x) {
    dist(rbind(x, target))
  })
  cluster_selected <- which.min(target_dist)

  if (gate_type == "ellipse") {
    contour_ellipse <- .getEllipse(filter = tmixRes1, include = cluster_selected,
                                   quantile = quantile)
    flowClust_gate <- polygonGate(.gate = matrix(contour_ellipse, ncol = 2,
                                    dimnames = list(NULL, tmixRes1@varNames)),
                                  filterId = filterId)
  } else if (gate_type == "axis") {
    chisq_quantile <- qchisq(quantile, df = 2)

    xbar <- tmixRes1@mu[cluster_selected, ]
    Sigma <- tmixRes1@sigma[cluster_selected, , ]

    Sigma_eigen <- eigen(Sigma, symmetric = TRUE)
    u1 <- Sigma_eigen$vectors[, 1]
    u2 <- Sigma_eigen$vectors[, 2]
    lambda1 <- Sigma_eigen$values[1]
    lambda2 <- Sigma_eigen$values[2]

    # Determines which eigenvector points towards the first quadrant (has
    # positive slope). Because both eigenvectors can potentially point in the
    # negative direction, we also check to see the negated eigenvectors point
    # towards the first quadrant.
    if (all(u1 >= 0)) {
      axis <- sqrt(lambda1 * chisq_quantile) * u1
      if (u2[1] > 0) {
        axis_perp <- sqrt(lambda2 * chisq_quantile) * u2    
      } else {
        axis_perp <- -sqrt(lambda2 * chisq_quantile) * u2    
      }
    } else if (all(-u1 >= 0)) {
      axis <- -sqrt(lambda1 * chisq_quantile) * u1
      if (u2[1] > 0) {
        axis_perp <- sqrt(lambda2 * chisq_quantile) * u2    
      } else {
        axis_perp <- -sqrt(lambda2 * chisq_quantile) * u2    
      }
    } else if (all(u2 >= 0)) {
      axis <- sqrt(lambda2 * chisq_quantile) * u2
      if (u1[1] > 0) {
        axis_perp <- sqrt(lambda1 * chisq_quantile) * u1   
      } else {
        axis_perp <- -sqrt(lambda1 * chisq_quantile) * u1    
      }
    } else if (all(-u2 >= 0)) {
      axis <- -sqrt(lambda2 * chisq_quantile) * u2
      if (u1[1] > 0) {
        axis_perp <- sqrt(lambda1 * chisq_quantile) * u1   
      } else {
        axis_perp <- -sqrt(lambda1 * chisq_quantile) * u1    
      }
    }

    # The gate location is the frame of reference for the gate. If it is xbar,
    # then the frame of reference is the cross-section along the eigenvector
    # from top-left to bottom-right. We translate this reference as a function
    # of the appropriate chi-squared coefficient.
    gate_location <- xbar + axis_translation * axis

    # To construct the gate, we have effectively shifted the eigenvector with the
    # negative slope. We then extend the gate both horizontally and vertically
    # to the maximum observed values in the horizontal and vertical directions.
    # NOTE: We extend the gate one standard deviation beyond the maximum values
    # observed to mimic a gate that extends without limit in the positive
    # directions. However, because flowCore cannot handle such a shape, we force
    # the gate to have the same shape for the observed data.
    x_max <- max(x) + sd(x)
    y_max <- max(y) + sd(y)
    gate_x <- c(x_max, (gate_location + axis_perp)[2])
    gate_y <- c((gate_location - axis_perp)[1], y_max)

    # We extend the gate to the min and max values
    polygon_gate <- rbind(gate_location + axis_perp,
                          gate_x,
                          c(x_max, y_max),
                          gate_y,
                          gate_location - axis_perp)

    colnames(polygon_gate) <- c(xChannel, yChannel)
    flowClust_gate <- polygonGate(filterId = filterId, boundaries = polygon_gate)
  }
  
  # TODO: need to verify if this is the right way to grab posteriors from 2-D tmixRes
  posteriors <- list(mu = tmixRes1@mu, lamdda = tmixRes1@lambda,
                     sigma = tmixRes1@sigma, nu = tmixRes1@nu)

  if (plot) {
    plot(fr, tmixRes1, main = filterId)

    if (gate_type == "axis") {
      # The major and minor axes (eigenvectors) scaled by their respective
      # eigenvalues and the chi-squared quantile.
      lines(rbind(xbar, xbar + sqrt(lambda1 * chisq_quantile) * u1), col = "darkgreen")
      lines(rbind(xbar, xbar - sqrt(lambda1 * chisq_quantile) * u1), col = "darkgreen")
      lines(rbind(xbar, xbar + sqrt(lambda2 * chisq_quantile) * u2), col = "darkgreen")
      lines(rbind(xbar, xbar - sqrt(lambda2 * chisq_quantile) * u2), col = "darkgreen")

      # Draws the polygon gate.
      lines(polygon_gate, col = "red")
      lines(rbind(gate_location - axis_perp, gate_location + axis_perp), col = "red")

      # Also, draws points at the vertices of the polygon gate.
      points(polygon_gate, col = "red", pch = 16)
    }
  }

  fcPolygonGate(flowClust_gate, prior, posteriors)
}


#' Finds the local maxima (peaks) in the given vector after smoothing the data
#' with a kernel density estimator.
#'
#' First, we smooth the data using kernel density estimation (KDE) with the
#' \code{\link{density}} function. Then, we find every local maxima such that the
#' density is concave (downward).
#'
#' Effectively, we find the local maxima with a discrete analogue to a second
#' derivative applied to the KDE. For details, see this StackOverflow post:
#' \url{http://bit.ly/Zbl7LV}.
#'
#' @param x numeric vector
#' @param order_peaks logical value. If \code{TRUE}, the peaks will be sorted by
#' the height of the peaks in decreasing order. (Default: \code{FALSE})
#' @param ... additional arguments passed to the \code{\link{density}} function
#' @return the values where the peaks are attained.
#' @examples
#' library(flowClust)
#' set.seed(42)
#' # 2 peaks with a minor hump
#' y <- SimulateMixture(10000, c(.5, .3, .2), c(2, 5, 7), c(1, 1, 1), nu = 10)
#' plot(density(y))
#' peaks <- find_peaks(y, adjust = 1.5)
#' abline(v = peaks, col = "red")
find_peaks <- function(x, order_peaks = FALSE, ...) {
  x <- as.vector(x)
  dens <- density(x, ...)

  # Discrete analogue to a second derivative applied to the KDE. See details.
  second_deriv <- diff(sign(diff(dens$y)))
  which_maxima <- which(second_deriv == -2) + 1

  if (order_peaks) {
    which_maxima <- which_maxima[order(dens$y[which_maxima], decreasing = TRUE)]
  }

  dens$x[which_maxima]
}

#' Finds the local minima (valleys) in the given vector after smoothing the data
#' with a kernel density estimator.
#'
#' First, we smooth the data using kernel density estimation (KDE) with the
#' \code{\link{density}} function. Then, we find every local minima such that the
#' density is concave (downward).
#'
#' Effectively, we find the local minima with a discrete analogue to a second
#' derivative applied to the KDE. For details, see this StackOverflow post:
#' \url{http://bit.ly/Zbl7LV}.
#'
#' @param x numeric vector
#' @param order_valleys logical value. If \code{TRUE}, the valleys will be sorted by
#' the height of the valleys in increasing order. (Default: \code{FALSE})
#' @param ... additional arguments passed to the \code{\link{density}} function
#' @return the values where the valleys are attained.
#' @examples
#' library(flowClust)
#' set.seed(42)
#' # 3 peaks and 2 valleys
#' y <- SimulateMixture(10000, c(.25, .5, .25), c(1, 5, 9), c(1, 1, 1), nu = 10)
#' plot(density(y))
#' valleys <- find_valleys(y, adjust = 1.5)
#' abline(v = valleys, col = "red")
find_valleys <- function(x, order_valleys = FALSE, ...) {
  x <- as.vector(x)
  dens <- density(x, ...)

  # Discrete analogue to a second derivative applied to the KDE. See details.
  second_deriv <- diff(sign(diff(dens$y)))
  which_minima <- which(second_deriv == 2) + 1

  if (order_valleys) {
    which_minima <- which_minima[order(dens$y[which_minima], decreasing = FALSE)]
  }

  dens$x[which_minima]
}

quantileGate <- function(fr, probs, stain, plot = FALSE, positive = TRUE,
                         filterId = "", ...) {
  x <- exprs(fr[, stain])
  cutpoint <- quantile(x, probs = probs)

  if (positive) {
    gate_coordinates <- list(c(cutpoint, Inf))
  } else {
    gate_coordinates <- list(c(-Inf, cutpoint))
  }
  names(gate_coordinates) <- stain
  
  if (plot) {
    hist(x)
    plot(density(x))
    abline(v = cutpoint, col = "red")
    text(y = 0.5, x = cutpoint, labels = paste("quantile:", probs))
  }
  rectangleGate(gate_coordinates, filterId = filterId)
}

#' Determines a cutpoint as the minimum point of a kernel density estimate
#' between two peaks
#'
#' We fit a kernel density estimator to the cells in the \code{flowFrame} and
#' identify the two largest peaks using the \code{find_peaks} function. We then
#' select as the cutpoint the value at which the minimum density is attained.
#'
#' @param flow_frame a \code{flowFrame} object
#' @param channel TODO
#' @param filter_id TODO
#' @param positive If \code{TRUE}, then the gate consists of the entire real
#' line to the right of the cutpoint. Otherwise, the gate is the entire real
#' line to the left of the cutpoint. (Default: \code{TRUE})
#' @param gate_min TODO: For now, this is the minimum value allowed for the
#' gate. If the constructed gate is less than this value, the minimum value
#' (or maybe min(x)) is returned.
#' @param ... Additional arguments pased on to the \code{find_peaks} function.
#' @return a \code{rectangleGate} object based on the minimum density cutpoint
mindensity <- function(flow_frame, channel, filter_id = "", positive = TRUE,
                       gate_min = NULL, ...) {
  
  if (missing(channel) || length(channel) != 1) {
    stop("A single channel must be specified.")
  }
  # Grabs the data matrix that is being gated.
  x <- exprs(flow_frame)[, channel]

  # Find the minimum density between the two highest peaks and set the cutpoint
  # there.
  largest_peaks <- find_peaks(x, order_peaks = TRUE, ...)[1:2]
  x_between <- x[findInterval(x, sort(largest_peaks)) == 1]

  # The cutpoint is the deepest valley between the two largest peaks
  cutpoint <- find_valleys(x_between, order_valleys = TRUE)[1]

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
}
