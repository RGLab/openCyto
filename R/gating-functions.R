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
#' @param trans numeric indicating whether the Box-Cox transformation parameter
#' is estimated from the data. May take 0 (no estimation), 1 (estimation) or 2
#' (cluster-speciﬁc estimation). NOTE: For the Bayesian version of
#' \code{\link{flowClust}}, this value cannot be 2.
#' @param positive If \code{TRUE}, then the gate consists of the entire real
#' line to the right of the cutpoint. Otherwise, the gate is the entire real
#' line to the left of the cutpoint. (Default: \code{TRUE})
#' @param prior list of prior parameters for the Bayesian
#' \code{\link{flowClust}}. If \code{NULL}, no prior is used.
#' @param criterion a character string stating the criterion used to choose the
#' best model. May take either "BIC" or "ICL". This argument is only relevant
#' when \code{K} is \code{NULL} or if \code{length(K) > 1}. The value selected
#' is passed to \code{\link{flowClust}}.
#' @param cutpoint_method How should the cutpoint be chosen from the fitted
#' \code{\link{flowClust}} model? See Details.
#' @param neg_cluster integer. The index of the negative cluster. The cutpoint
#' is computed between clusters \code{neg_cluster} and \code{neg_cluster + 1}.
#' @param cutpoint_min numeric value that sets a minimum thresold for the
#' cutpoint. If a value is provided, any cutpoint below this value will be set
#' to the given minimum value. If \code{NULL} (default), there is no minimum
#' cutpoint value.
#' @param cutpoint_max numeric value that sets a maximum thresold for the
#' cutpoint. If a value is provided, any cutpoint above this value will be set
#' to the given maximum value. If \code{NULL} (default), there is no maximum
#' cutpoint value.
#' @param min a numeric value that sets the lower boundary for data filtering
#' @param max a numeric value that sets the upper boundary for data filtering
#' @param quantile the quantile for which we will find the cutpoint using
#' the quantile \code{cutpoint_method}. If the \code{cutpoint_method} is not set
#' to \code{quantile}, this argument is ignored.
#' @param quantile_interval a vector of length 2 containing the end-points of
#' the interval of values to find the quantile cutpoint. If the
#' \code{cutpoint_method} is not set to \code{quantile}, this argument is
#' ignored.
#' @param plot logical value indicating that the fitted \code{\link{flowClust}}
#' model should be plotted along with the cutpoint
#' @param ... additional arguments that are passed to \code{\link{flowClust}}
#' @return a \code{rectangleGate} object consisting of all values beyond the
#' cutpoint calculated
flowClust.1d <- function(fr, params, filterId = "", K = NULL, trans = 0,
                         positive = TRUE, prior = NULL,
                         criterion = c("BIC", "ICL"),
                         cutpoint_method = c("boundary", "min_density", "quantile", "posterior_mean"),                         
                         neg_cluster = 1, cutpoint_min = NULL,
                         cutpoint_max = NULL, min = NULL, max = NULL,
                         quantile = 0.99, quantile_interval = c(0, 10),
                         plot = FALSE, ...) {

  cutpoint_method <- match.arg(cutpoint_method)

  # TODO: Determine if Bayesian flowClust works when 'K' is specified and has a
  # different length than the number of components given in 'prior'.
  # If not, add GitHub issue and then fix.

  # If 'K' is NULL, then 'K' is autoselected using either 'BIC' or 'ICL'.
  # The candidate values for 'K' range from 1 to the number of mixture components
  # given in the 'prior' list.
  if (is.null(K)) {
    criterion <- match.arg(criterion)

    if (!is.null(prior)) {
      K <- length(as.vector(prior$Mu0))
    } else {
      stop("Values for 'K' must be provided if no 'prior' is given.")
    }
  }

  usePrior <- ifelse(is.null(prior), "no", "yes")

  # HACK: Circumvents a bug in flowClust.
  # TODO: Add an issue to the Github for flowClust to allow prior to be NULL.
  if (is.null(prior)) {
    prior <- list(NA)
  }

  # Filter out values less than the minimum and above the maximum, if they are
  # given. NOTE: These observations are removed from the 'flowFrame' locally and
  # are gated out only for the determining the gate.
  fr <- truncate_flowframe(fr, channel = params[1], min = min, max = max)

  # Applies `flowClust` to the feature specified in the `params` argument using
  # the data given in `fr`. We use priors with hyperparameters given by the
  # elements in the `prior` list.
  
  tmix_filter <- tmixFilter(filterId, params[1], K = K, trans = trans,
                            usePrior = usePrior, prior = prior,
                            criterion = criterion, ...)
  tmix_results <- filter(fr, tmix_filter)
  
  # To determine the cutpoint, we first sort the centroids so that we can determine
  # the second largest centroid.
  centroids_sorted <- sort(getEstimates(tmix_results)$locations)

  # Also, because the cluster labels are arbitrary, we determine the cluster
  # the cluster labels, sorted by the ordering of the cluster means.
  labels_sorted <- order(getEstimates(tmix_results)$locations)

  # Grabs the data matrix that is being gated.
  x <- exprs(fr)[, params[1]]

  # Determines the cutpoint between clusters 1 and 2.
  if (cutpoint_method == "boundary") {
    # Choose the cutpoint as the boundary between the first two clusters.
    # First, we sort the data.
    order_x <- order(x)
    x_sorted <- x[order_x]
    labels <- flowClust::Map(tmix_results, rm.outliers = FALSE)[order_x]

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
    x_dens <- dmvtmix(x = as.matrix(x_between), object = tmix_results)
    cutpoint <- x_between[which.min(x_dens)]
  } else if (cutpoint_method == "quantile") {
    cutpoint <- quantile_flowClust(p = quantile, object = tmix_results,
                                   interval = quantile_interval)
  } else { # cutpoint_method == "posterior_mean"
    cutpoint <- centroids_sorted[neg_cluster]
  }

  # In some cases, we wish that a gating cutpoint not exceed some threshold, in
  # which case we allow the user to specify a minimum and/or maximum value for
  # the cutpoint. For instance, when constructing a debris gate for samples from
  # various batches, some samples may have had a debris gate applied already
  # while other samples have debris that should be removed. For the former, the
  # user may wish to avoid gating, in which case, we provide the option to
  # threshold the gate cutpoint.
  if (!is.null(cutpoint_min) && cutpoint < cutpoint_min) {
    cutpoint <- cutpoint_min
  }
  if (!is.null(cutpoint_max) && cutpoint > cutpoint_max) {
    cutpoint <- cutpoint_max
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
  posteriors <- list(mu = tmix_results@mu, lamdda = tmix_results@lambda,
                     sigma = tmix_results@sigma, nu = tmix_results@nu, min = min(x),
                     max = max(x))
  postList[[params[1]]] <- posteriors

  # Saves prior point estimates
  priorList <- list()
  priorList[[params[1]]] <- prior
  
  if (plot) {
    gate_pct <- round(100 * mean(x > cutpoint), 3)
    plot_title <- paste0(filterId, " (", gate_pct, "%)")
    plot(fr, tmix_results, main = plot_title, labels = FALSE)
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
        nu <- ifelse(length(tmix_results@nu) == 1, tmix_results@nu, tmix_results@nu[k])
        lambda <- ifelse(length(tmix_results@lambda) == 1, tmix_results@lambda, tmix_results@lambda[k])        
        posterior_density <- dmvt(x_dens, mu = tmix_results@mu[k,],
                                  sigma = tmix_results@sigma[k,,], nu = nu,
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
#' @param usePrior Should we use the Bayesian version of
#' \code{\link{flowClust}}?  Answers are "yes", "no", or "vague". The answer is
#' passed along to \code{\link{flowClust}}.
#' @param prior list of prior parameters for the Bayesian
#' \code{\link{flowClust}}.  If \code{usePrior} is set to 'no', then the list is
#' unused.
#' @param trans numeric indicating whether the Box-Cox transformation parameter
#' is estimated from the data. May take 0 (no estimation), 1 (estimation) or 2
#' (cluster-speciﬁc estimation). NOTE: For the Bayesian version of
#' \code{\link{flowClust}}, this value cannot be 2.
#' @param plot a logical value indicating if the fitted mixture model should be
#' plotted. By default, no.
#' @param target a numeric vector of length \code{K} containing the location of
#' the cluster of interest. See details.
#' @param gate_type character value specifying the type of gate to construct from
#' the \code{\link{flowClust}} fit
#' @param quantile the contour level of the target cluster from the
#' \code{\link{flowClust}} fit to construct the gate
#' @param axis_translation a numeric value between 0 and 1 used to position the
#' gate if \code{gate_type} is selected as \code{"axis"}. See details.
#' @param min A vector of length 2. Truncate observations less than this minimum
#' value. The first value truncates the \code{xChannel}, and the second value
#' truncates the \code{yChannel}. By default, this vector is \code{NULL} and is
#' ignored.
#' @param max A vector of length 2. Truncate observations greater than this
#' maximum value. The first value truncates the \code{xChannel}, and the second
#' value truncates the \code{yChannel}. By default, this vector is \code{NULL}
#' and is ignored.
#' @param ... additional arguments that are passed to \code{\link{flowClust}}
#' @return a \code{polygonGate} object containing the contour (ellipse) for 2D
#' gating.
flowClust.2d <- function(fr, xChannel, yChannel, filterId = "", K = 2,
                         usePrior = 'no', prior = list(NA), trans = 0,
                         plot = FALSE, target = rep(0, K),
                         gate_type = c("ellipse", "axis"), quantile = 0.9,
                         axis_translation = 0.25, min = NULL, max = NULL, ...) {

  if (length(target) != 2) {
    stop("The 'target' location must be a numeric vector of length 2.")
  }
  gate_type <- match.arg(gate_type)

  # If a truncation value is specified, we remove all observations less than
  # this value for the marker specified to construct the gate. NOTE: These
  # observations are removed from the 'flowFrame' locally and are gated out only
  # for the determining the gate.
  fr <- truncate_flowframe(fr, channel = xChannel, min = min[1], max = max[1])
  fr <- truncate_flowframe(fr, channel = yChannel, min = min[2], max = max[2])

  x <- exprs(fr)[, xChannel]
  y <- exprs(fr)[, yChannel]

  # If appropriate, we generate prior parameters for the Bayesian version of flowClust.
  if (usePrior == "yes" && identical(prior, list(NA))) {
    prior <- prior_flowClust(fr = fr, channels = c(xChannel, yChannel), K = K)
  }

  # Applies `flowClust` to the feature specified in the `params` argument using
  # the data given in `fr`. We use priors with hyperparameters given by the
  # elements in the list `prior`.
  tmix_filter <- tmixFilter(filterId, c(xChannel, yChannel), K = K, trans = trans,
                        usePrior = usePrior, prior = prior, ...)
  tmix_results <- filter(fr, tmix_filter)

  # Converts the tmixFilterResult object to a polygonGate.
  # We select the cluster with the minimum 'yChannel' to be the subpopulation from
  # which we obtain the contour (ellipse) to generate the polygon gate.
  fitted_means <- getEstimates(tmix_results)$locations
  target_dist <- apply(fitted_means, 1, function(x) {
    dist(rbind(x, target))
  })
  cluster_selected <- which.min(target_dist)

  if (gate_type == "ellipse") {
    contour_ellipse <- .getEllipse(filter = tmix_results, include = cluster_selected,
                                   quantile = quantile)
    flowClust_gate <- polygonGate(.gate = matrix(contour_ellipse, ncol = 2,
                                    dimnames = list(NULL, tmix_results@varNames)),
                                  filterId = filterId)
  } else if (gate_type == "axis") {
    chisq_quantile <- qchisq(quantile, df = 2)

    xbar <- tmix_results@mu[cluster_selected, ]
    Sigma <- tmix_results@sigma[cluster_selected, , ]

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
  posteriors <- list(mu = tmix_results@mu, lamdda = tmix_results@lambda,
                     sigma = tmix_results@sigma, nu = tmix_results@nu)

  if (plot) {
    plot(fr, tmix_results, main = filterId)

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
#' select as the cutpoint the value at which the minimum density is attained
#' between the two peaks of interest.
#'
#' In the default case, the two peaks of interest are the two largest peaks
#' obtained from the \code{link{\density}} function. However, if \code{pivot} is
#' \code{TRUE}, we choose the largest peak and its neighboring peak as the two
#' peaks of interest. In this case, the neighboring peak is the peak immediately
#' to the left of the largest peak if \code{positive} is \code{TRUE}. Otherwise,
#' the neighboring peak is selected as the peak to the right.
#'
#' In the special case that there is only one peak, we are conservative and set
#' the cutpoint as the \code{min(x)} if \code{positive} is \code{TRUE}, and the
#' \code{max(x)} otherwise.
#'
#' @param flow_frame a \code{flowFrame} object
#' @param channel TODO
#' @param filter_id TODO
#' @param positive If \code{TRUE}, then the gate consists of the entire real
#' line to the right of the cutpoint. Otherwise, the gate is the entire real
#' line to the left of the cutpoint. (Default: \code{TRUE})
#' @param pivot logical value. If \code{TRUE}, we choose as the two peaks the
#' largest peak and its neighboring peak. See details.
#' @param min a numeric value that sets the lower boundary for data filtering
#' @param max a numeric value that sets the upper boundary for data filtering
#' @param ... Additional arguments passed on to the \code{find_peaks} function
#' @return a \code{rectangleGate} object based on the minimum density cutpoint
mindensity <- function(flow_frame, channel, filter_id = "", positive = TRUE,
                       pivot = FALSE, min = NULL, max = NULL, ...) {
  
  if (missing(channel) || length(channel) != 1) {
    stop("A single channel must be specified.")
  }

  # Filter out values less than the minimum and above the maximum, if they are
  # given.
  flow_frame <- truncate_flowframe(flow_frame, channel = channel, min = min,
                                     max = max)

  # Grabs the data matrix that is being gated.
  x <- exprs(flow_frame)[, channel]

  peaks <- find_peaks(x, order_peaks = TRUE, ...)

  if (pivot) {
    # If 'pivot' is selected, we choose the largest peak and its neighbor, which
    # is chosen based on the current value of 'positive'
    largest_peak <- peaks[1]
    peaks <- sort(peaks)
    which_largest <- which(peaks == largest_peak)
    if (positive) {
      peaks <- peaks[c(which_largest - 1, which_largest)]
    } else {
      peaks <- peaks[c(which_largest, which_largest + 1)]
    }
    peaks <- peaks[!is.na(peaks)]
  } else {
    # Otherwise, we choose the two largest peaks and sort them
    peaks <- sort(peaks[1:2])
  }

  # In the special case that there is only one peak, we are conservative and set
  # the cutpoint as the min(x) if 'positive' is TRUE, and the max(x) otherwise.
  # value otherwise.
  if (length(peaks) == 1) {
    cutpoint <- ifelse(positive, min(x), max(x))
  } else {
    # Find the minimum density between the two peaks selected and set the
    # cutpoint there.
    x_between <- x[findInterval(x, peaks) == 1]

    # The cutpoint is the deepest valley between the two peaks selected. In the
    # case that there are no valleys (i.e., if 'x_between' has an insufficient
    # number of observations), we are conservative and set the cutpoint as the
    # minimum value if 'positive' is TRUE, and the maximum value otherwise.
    cutpoint <- try(find_valleys(x_between, order_valleys = TRUE)[1],
                    silent = TRUE)

    if (class(cutpoint) == "try-error" || is.na(cutpoint)) {
      cutpoint <- ifelse(positive, peaks[1], peaks[2])
    }
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

  names(gate_coordinates) <- channel
  
  rectangleGate(gate_coordinates, filterId = filter_id)
}
