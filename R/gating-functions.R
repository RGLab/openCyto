#' @templateVar old gate_flowClust_1d
#' @templateVar new gate_flowclust_1d
#' @template template-depr_pkg
NULL
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
#' @rdname gate_flowclust_1d
#' @aliases gate_flowclust_1d gate_flowClust_1d flowClust.1d
#' @param fr a \code{flowFrame} object
#' @param params \code{character} channel to be gated on
#' @param filterId A \code{character} string that identifies the filter created.
#' @param K the number of clusters to find
#' @param trans,min.count,max.count,nstart some flowClust parameters. see \code{\link{flowClust}}
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
#' @param min a numeric value that sets the lower bound for data filtering. If
#' \code{NULL} (default), no truncation is applied.
#' @param max a numeric value that sets the upper bound for data filtering. If
#' \code{NULL} (default), no truncation is applied.
#' @param quantile the quantile for which we will find the cutpoint using
#' the quantile \code{cutpoint_method}. If the \code{cutpoint_method} is not set
#' to \code{quantile}, this argument is ignored.
#' @param quantile_interval a vector of length 2 containing the end-points of
#' the interval of values to find the quantile cutpoint. If the
#' \code{cutpoint_method} is not set to \code{quantile}, this argument is
#' ignored.
#' @param plot logical value indicating that the fitted \code{\link{flowClust}}
#' model should be plotted along with the cutpoint
#' @param debug \code{logical} indicating whether to carry the prior and posterious with the gate
#'                               for debugging purpose. Default is FALSE.
#' @param ... additional arguments that are passed to \code{\link{flowClust}}
#' @return a \code{rectangleGate} object consisting of all values beyond the
#' cutpoint calculated
#' @export
#' @importFrom flowClust getEstimates tmixFilter dmvtmix dmvt flowClust
#' @examples
#' \dontrun{
#'  gate <- gate_flowclust_1d(fr, params = "APC-A", K =2) # fr is a flowFrame
#' }
gate_flowclust_1d <- function(fr, params, filterId = "", K = NULL 
                              , trans = 0 #no box transform
                              , min.count = -1, max.count = -1 #no flowClust data filering
                              , nstart = 1 #change kmeans nstart
                              , prior = NULL,
                              criterion = c("BIC", "ICL"),
                              cutpoint_method = c("boundary", "min_density",
                                                  "quantile", "posterior_mean", "prior_density"),
                              neg_cluster = 1, cutpoint_min = NULL,
                              cutpoint_max = NULL, min = NULL, max = NULL,
                              quantile = 0.99, quantile_interval = c(0, 10),
                              plot = FALSE, debug = FALSE, ...) {
  options("cores" = 1L) #suppress parallelism for 1d gating   
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
  
  #  L<-list(...)
  #  if("useprior"%in%names(L)){
  #	  usePrior <- L["useprior"]
  #  	L$useprior<-NULL
  #  }
  
  # HACK: Circumvents a bug in flowClust.
  # TODO: Add an issue to the Github for flowClust to allow prior to be NULL.
  if (is.null(prior)) {
    prior <- list(NA)
  }
  
  # Filter out values less than the minimum and above the maximum, if they are
  # given. NOTE: These observations are removed from the 'flowFrame' locally and
  # are gated out only for the determining the gate.
  if (!(is.null(min) && is.null(max))) {
    fr <- .truncate_flowframe(fr, channels = params[1], min = min, max = max)
  }
  if (nrow(fr) < 2) {
    warning("Less than two observations are present in the given flowFrame.",
            "Constructing gate from prior...")
  }
  
  # Applies `flowClust` to the feature specified in the `params` argument using
  # the data given in `fr`. We use priors with hyperparameters given by the
  # elements in the `prior` list.
  # call via do.call. L contains the rest of the pairlist, after extracting the passed value of usePrior
  #  	tmix_filter<-do.call(tmixFilter,c(filterId=filterId,params[1],K=K,
  #                                      trans=trans,usePrior=usePrior,
  #                                      prior=list(prior),criterion=list(criterion),
  #                                      L))
  tmix_results <- flowClust(fr, varNames = params[1]
                            , K = K, trans = trans
                            , usePrior = usePrior
                            , prior = prior
                            , min.count = min.count, max.count = max.count
                            , nstart = nstart
                            , criterion = criterion
                            , ...)
  
  # In the case an error occurs when applying 'flowClust', the gate is
  # constructed from the density of the prior distributions. This error
  # typically occurs when there are less than 2 observations in the flow frame.
  if (class(tmix_results) != "try-error") {
    # To determine the cutpoint, we first sort the centroids so that we can
    # determine the second largest centroid.
    centroids_sorted <- sort(getEstimates(tmix_results)$locations)
    
    # Also, because the cluster labels are arbitrary, we determine the cluster
    # the cluster labels, sorted by the ordering of the cluster means.
    labels_sorted <- order(getEstimates(tmix_results)$locations)
    
    # Grabs the data matrix that is being gated.
    x <- exprs(fr)[, params[1]]
  } else {
    cutpoint_method <- "prior_density"
  }
  
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
    cutpoint <- .quantile_flowClust(p = quantile, object = tmix_results,
                                    interval = quantile_interval)
  } else if (cutpoint_method == "posterior_mean") {
    cutpoint <- centroids_sorted[neg_cluster]
  } else { # cutpoint_method == "prior_density"
    # The prior_density cutpoint is determined as the point at which the density
    # of the negative cluster becomes smaller than the density of its adjacent
    # cluster.
    prior_x <- as.matrix(seq(prior$Mu0[neg_cluster], prior$Mu0[neg_cluster + 1], length = 1000))
    
    prior_proportions <- with(prior, w0 / sum(w0))
    prior_y <- lapply(c(neg_cluster, neg_cluster + 1), function(k) {
      prior_density <- dmvt(x = prior_x, mu = prior$Mu0[k],
                            sigma = prior$Omega0[k], nu = 4)$value
      prior_proportions[k] * prior_density
    })
    prior_y <- do.call(cbind, prior_y)
    
    # For the two prior densities, we determine the cutpoint as the first value
    # where the second density is larger than the first.
    diff_densities <- apply(prior_y, 1, diff)
    cutpoint <- prior_x[which(diff_densities > 0)[1]]
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
  
  
  gate_coordinates <- list(c(cutpoint, Inf))
  
  names(gate_coordinates) <- params
  
  fres <- rectangleGate(gate_coordinates, filterId = filterId)
  if(debug){
    # Saves posterior point estimates
    # In the case that an error is thrown, the posterior is set to the prior
    # because no prior updating was performed.
    postList <- list()
    if (class(tmix_results) != "try-error") {
      posteriors <- list(mu = tmix_results@mu, lambda = tmix_results@lambda,
                         sigma = tmix_results@sigma, nu = tmix_results@nu, min = min(x)
                         ,w = tmix_results@w
                         , max = max(x))
    } else {
      posteriors <- prior
      posteriors$min <- NA
      posteriors$max <- NA
    }
    postList[[params[1]]] <- posteriors
    
    # Saves prior point estimates
    priorList <- list()
    priorList[[params[1]]] <- prior
    fres <- fcRectangleGate(fres, priorList, postList)
  }
  
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
  fres
  
}

#' @export
gate_flowClust_1d <- function(fr, params, filterId = "", K = NULL 
                              , trans = 0 #no box transform
                              , min.count = -1, max.count = -1 #no flowClust data filering
                              , nstart = 1 #change kmeans nstart
                              , prior = NULL,
                              criterion = c("BIC", "ICL"),
                              cutpoint_method = c("boundary", "min_density",
                                                  "quantile", "posterior_mean", "prior_density"),
                              neg_cluster = 1, cutpoint_min = NULL,
                              cutpoint_max = NULL, min = NULL, max = NULL,
                              quantile = 0.99, quantile_interval = c(0, 10),
                              plot = FALSE, debug = FALSE, ...) {
  .Deprecated("gate_flowclust_1d")
  gate_flowclust_1d(fr, params, filterId, K , trans, min.count, max.count
                    , nstart, prior, criterion, cutpoint_method, neg_cluster 
                    , cutpoint_min, cutpoint_max, min, max, quantile, quantile_interval
                    , plot, debug, ...)
}

#' @export 
flowClust.1d <- gate_flowClust_1d

#' @templateVar old gate_flowClust_2d
#' @templateVar new gate_flowclust_2d
#' @template template-depr_pkg
NULL
#' Automatic identification of a population of interest via flowClust based on
#' two markers
#'
#' We cluster the observations in \code{fr} into \code{K} clusters. We set the
#' cutpoint to be the point at which the density between the first and second
#' smallest cluster centroids is minimum.
#'
#' The cluster for the population of interest is selected as the one with
#' cluster centroid nearest the \code{target} in Euclidean distance. By default,
#' the largest cluster (i.e., the cluster with the largest proportion of
#' observations) is selected as the population of interest.
#'
#' We also provide the option of constructing a \code{transitional} gate from
#' the selected population of interest. The location of the gate can be
#' controlled with the \code{translation} argument, which translates the gate
#' along the major axis of the targest cluster as a function of the appropriate
#' chi-squared coefficient. The larger \code{translation} is, the more gate is
#' shifted in a positive direction. Furthermore, the width of the
#' \code{transitional} gate can be controlled with the \code{quantile} argument.
#'
#' The direction of the transitional gate can be controlled with the
#' \code{transitional_angle} argument. By default, it is \code{NULL}, and we use
#' the eigenvector of the \code{target} cluster that points towards the first
#' quadrant (has positive slope). If \code{transitional_angle} is specified, we
#' rotate the eigenvectors so that the angle between the x-axis (with the cluster
#' centroid as the origin) and the major eigenvector (i.e., the eigenvector with
#' the larger eigenvalue) is \code{transitional_angle}. 
#' So based on range that the angle falls in, the final rectangleGate will be constructed 
#' at the corresponding quadrant. i.e. Clockwise, [0,pi/2] UR, (pi/2, pi] LR,
#' (pi, 3/2 * pi] LL, (3/2 * pi, 2 * pi] UL
#'
#' @name gate_flowclust_2d
#' @aliases gate_flowclust_2d gate_flowClust_2d flowClust.2d
#' @param fr a \code{flowFrame} object
#' @param xChannel,yChannel \code{character} specifying channels to be gated on
#' @param filterId A \code{character} string that identifies the filter created.
#' @param K the number of clusters to find
#' @param usePrior Should we use the Bayesian version of \code{\link{flowClust}}?
#' Answers are "yes", "no", or "vague". The answer is passed along to
#' \code{\link{flowClust}}.
#' @param prior list of prior parameters for the Bayesian version of
#' \code{\link{flowClust}}. If \code{usePrior} is set to \code{no}, then the
#' list is unused.
#' @param trans,min.count,max.count,nstart some flowClust parameters. see \code{\link{flowClust}}
#' @param plot a logical value indicating if the fitted mixture model should be
#' plotted. By default, no.
#' @param target a numeric vector of length \code{2} (number of dimensions) containing the location of
#' the cluster of interest. See details.
#' @param transitional logical value indicating if a transitional gate should be
#' constructed from the target \code{\link{flowClust}} cluster. By default, no.
#' @param quantile the contour level of the target cluster from the
#' \code{\link{flowClust}} fit to construct the gate
#' @param translation a numeric value between 0 and 1 used to position a
#' transitional gate if \code{transitional = TRUE}. This argument is ignored if
#' \code{transitional = FALSE}. See details
#' @param transitional_angle the angle (in radians) of the transitional
#' gate. It is also used to determine which quadrant the final gate resides in.
#' See details. Ignored if \code{transitional = FALSE}. 
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
#' @export
#' @examples
#' \dontrun{
#'  gate <- gate_flowclust_2d(fr, xChannel = "FSC-A", xChannel = "SSC-A", K = 3) # fr is a flowFrame
#' }
gate_flowclust_2d <- function(fr, xChannel, yChannel, filterId = "", K = 2,
                         usePrior = 'no', prior = list(NA)
                         , trans = 0 #no box transform
                         , min.count = -1, max.count = -1 #no flowClust data filering
                         , nstart = 1 #change kmeans nstart
                         , plot = FALSE, target = NULL, transitional = FALSE,
                         quantile = 0.9, translation = 0.25, transitional_angle = NULL,
                         min = NULL, max = NULL, ...) {
  options("cores" = 1L) ##suppress parallelism since it may lock the process due to the lack of resource when parallel gating was already using up the cores. 
  
  if (!is.null(target)) {
    target <- as.numeric(target)
    if (length(target) != 2) {
      warning("The 'target' location must be a numeric vector of length 2.
               Using largest cluster instead...")
      target <- NULL
    }    
  }

  # If specified, truncates all observations outside the 'min' and 'max' values.
  # NOTE: These observations are removed from the 'flowFrame' locally and are
  # gated out only for the determining the gate.
  if (!(is.null(min) && is.null(max))) {
    fr <- .truncate_flowframe(fr, channels = c(xChannel, yChannel), min = min,
                             max = max)
  }

  # If appropriate, we generate prior parameters for the Bayesian version of flowClust.
  if (usePrior == "yes" && identical(prior, list(NA))) {
    prior <- prior_flowclust(fr = fr, channels = c(xChannel, yChannel), K = K)
  }

  # Applies `flowClust` to the feature specified in the `params` argument using
  # the data given in `fr`. We use priors with hyperparameters given by the
  # elements in the list `prior`.
  
  tmix_results <- flowClust(fr, varNames = c(xChannel, yChannel)
                                , K = K, trans = trans
                                , usePrior = usePrior
                                , prior = prior
                                , min.count = min.count, max.count = max.count
                                , nstart = nstart
                                , ...)
                    
  # tmix_results <- filter(fr, tmix_filter)
  # In the case an error occurs when applying 'flowClust', the gate is
  # constructed from the prior distributions. Errors typically occur when there
  # are less than 2 observations in the flow frame.
  if (class(tmix_results) == "try-error") {
    tmix_results <- new("flowClust", varNames = c(xChannel, yChannel), K = K,
                        w = prior$w0, mu = prior$Mu0, sigma = prior$Lambda0,
                        nu = 4, prior = prior, ruleOutliers = c(0, quantile, quantile))
  }

  fitted_means <- getEstimates(tmix_results)$locations

  # By default, the cluster with the largest number of observations is
  # selected. Otherwise, the cluster centroid nearest to the 'target' is
  # selected.
  if (is.null(target)) {
    cluster_selected <- which.max(tmix_results@w)
  } else {
    target_dist <- as.matrix(dist(rbind(fitted_means, target)))
    target_dist <- tail(target_dist, n = 1)[seq_len(K)]
    cluster_selected <- which.min(target_dist)
  }

  if (!transitional) {
    flowClust_gate <- .getEllipseGate(filter = tmix_results, include = cluster_selected,
                                   quantile = quantile,trans=trans)
    
  } else {
    chisq_quantile <- qchisq(quantile, df = 2)
    tol <- sqrt(.Machine$double.eps)

    xbar <- tmix_results@mu[cluster_selected, ]
    Sigma <- tmix_results@sigma[cluster_selected, , ]

    Sigma_eigen <- eigen(Sigma, symmetric = TRUE)
    u1 <- Sigma_eigen$vectors[, 1]
    u2 <- Sigma_eigen$vectors[, 2]
    lambda1 <- Sigma_eigen$values[1]
    lambda2 <- Sigma_eigen$values[2]

    # Computes the angles of each eigenvector with the x-axis in terms of polar
    # coordinates. Note that each vector has magnitude 1, so the radius is 1.
    u1_angle <- atan2(u1[2], u1[1])
    if (u1_angle < 0) {
      u1_angle <- u1_angle + 2 * pi
    }
    u2_angle <- atan2(u2[2], u2[1])
    if (u2_angle < 0) {
      u2_angle <- u2_angle + 2 * pi
    }

    # We ensure that each eigenvector is pointing vertically. That is, we ensure
    # the angle between the x-axis and the eigenvector is between 0 and pi. If
    # If they are not, we rotate them by pi (i.e., 180 degrees).
    if (u1_angle > pi) {
      R <- .rotation_matrix(pi)
      u1 <- as.vector(R %*% u1)
    }
    if (u2_angle > pi) {
      R <- .rotation_matrix(pi)
      u2 <- as.vector(R %*% u2)
    }

    u1_angle <- atan2(u1[2], u1[1])
    if (u1_angle < 0) {
      u1_angle <- u1_angle + 2 * pi
    }
    u2_angle <- atan2(u2[2], u2[1])
    if (u2_angle < 0) {
      u2_angle <- u2_angle + 2 * pi
    }
    eigen_angles <- c(u1_angle, u2_angle)

    # If the transitional angle is not provided, we set it as the angle between
    # the x-axis as the eigenvector pointing towards the postive quadrant.
    # If the transitional angle is provided, we first calculate the angle between
    # the major eigenvector and the x-axis (with xbar as the origin). We then
    # construct a rotation matrix, where the angle applied is the difference
    # between the angle provided and the angle calculated. We then rotate the
    # eigenvectors u1 and u2 accordingly
    #
    # In both cases, we compute the axis from which the transitional gate is
    # constructed as well as the perpendicular axis. We ensure that 'axis_perp'
    # is pi/2 radians counterclockwise from 'axis', following the right-hand
    # rule.
    if (is.null(transitional_angle)) {

      # Determines which eigenvector points towards the positive quadrant
      which_pos_quadrant <- which(0 < eigen_angles & eigen_angles < pi/2)

      # Calculates the angle of the transitional gate as the angle of the
      # eigenvector pointing toward the positive quadrant
      transitional_angle <- eigen_angles[which_pos_quadrant]

      if (which_pos_quadrant == 1) {
        axis <- sqrt(lambda1 * chisq_quantile) * u1
        axis_perp <- sqrt(lambda2 * chisq_quantile) * u2
      } else {
        axis <- sqrt(lambda2 * chisq_quantile) * u2
        axis_perp <- sqrt(lambda1 * chisq_quantile) * u1
      }
    } else {
      # Rotation angle
      theta_u1 <- transitional_angle - eigen_angles[1]
      theta_u2 <- transitional_angle - (pi/2) - eigen_angles[2]

      # Rotation matrix
      R1 <- .rotation_matrix(theta_u1)
      R2 <- .rotation_matrix(theta_u2)

      # Rotates the eigenvectors
      u1 <- as.vector(R1 %*% u1)
      u2 <- as.vector(R2 %*% u2)

      axis <- sqrt(lambda1 * chisq_quantile) * u1
      axis_perp <- -sqrt(lambda2 * chisq_quantile) * u2
    }

    # The gate location is the frame of reference for the gate. If it is xbar,
    # then the frame of reference is the cross-section along the eigenvector
    # from top-left to bottom-right. We translate this reference as a function
    # of the appropriate chi-squared coefficient.
    gate_location <- xbar + translation * axis

    # To construct the gate, we have effectively shifted the eigenvector with the
    # negative slope. We then extend the gate both horizontally and vertically
    # to the maximum observed values in the horizontal and vertical directions.
    # NOTE: We extend the gate one standard deviation beyond the maximum values
    # observed to mimic a gate that extends without limit in the positive
    # directions. However, because flowCore cannot handle such a shape, we force
    # the gate to have the same shape for the observed data.
    x <- exprs(fr)[, xChannel]
    y <- exprs(fr)[, yChannel]

    x_min <- min(x) - sd(x)
    y_min <- min(y) - sd(y)
    x_max <- max(x) + sd(x)
    y_max <- max(y) + sd(y)

    # We construct the gate in clockwise fashion. The vertices of the gate
    # depend on the quadrant towards which the 'transitional_angle' is pointed.
    # No matter what the transitional gate's angle, the first and last vertices
    # will be the same
    first_vertex <- gate_location + axis_perp
    fifth_vertex <- gate_location - axis_perp
    
    if (0 <= transitional_angle && transitional_angle <= pi/2) {
      # First quadrant
      second_vertex <- c(first_vertex[1], y_max)
      third_vertex <- c(x_max, y_max)
      fourth_vertex <- c(x_max, fifth_vertex[2])
    } else if (pi/2 < transitional_angle && transitional_angle <= pi) {
      # Second quadrant
      second_vertex <- c(x_min, first_vertex[2])
      third_vertex <- c(x_min, y_max)
      fourth_vertex <- c(fifth_vertex[1], y_max)
    } else if (pi < transitional_angle && transitional_angle <= 3*pi/2) {
      # Third quadrant
      second_vertex <- c(first_vertex[1], y_min)
      third_vertex <- c(x_min, y_min)
      fourth_vertex <- c(x_min, fifth_vertex[2])
    } else {
      # Fourth quadrant
      second_vertex <- c(x_max, first_vertex[2])
      third_vertex <- c(x_max, y_min)
      fourth_vertex <- c(fifth_vertex[1], y_min)
    }
    polygon_gate <- rbind(first_vertex,
                          second_vertex,
                          third_vertex,
                          fourth_vertex,
                          fifth_vertex,
                          first_vertex)
    colnames(polygon_gate) <- c(xChannel, yChannel)
    flowClust_gate <- polygonGate(filterId = filterId, .gate = polygon_gate)
  }
  
  # List of posterior point estimates
  posteriors <- list(mu = tmix_results@mu, lambda = tmix_results@lambda,
                     sigma = tmix_results@sigma, nu = tmix_results@nu)

  if (plot) {
    plot(data = fr, tmix_results, main = filterId)

    if (transitional) {
      # The major and minor axes (eigenvectors) scaled by their respective
      # eigenvalues and the chi-squared quantile.    
      lines(rbind(xbar - axis, xbar + axis), col = "darkgreen")
      lines(rbind(xbar - axis_perp, xbar + axis_perp), col = "darkgreen")

      # Draws the polygon gate.
      lines(polygon_gate, col = "red")
      lines(rbind(gate_location - axis_perp, gate_location + axis_perp), col = "red")

      # Also, draws points at the vertices of the polygon gate.
      points(polygon_gate, col = "red", pch = 16)
    }
  }
  if(class(flowClust_gate) == "polygonGate")
    fcPolygonGate(flowClust_gate, prior, posteriors)
  else
    fcEllipsoidGate(flowClust_gate, prior, posteriors)
    
}

#' @export
gate_flowClust_2d <- function(fr, xChannel, yChannel, filterId = "", K = 2,
                              usePrior = 'no', prior = list(NA)
                              , trans = 0 #no box transform
                              , min.count = -1, max.count = -1 #no flowClust data filering
                              , nstart = 1 #change kmeans nstart
                              , plot = FALSE, target = NULL, transitional = FALSE,
                              quantile = 0.9, translation = 0.25, transitional_angle = NULL,
                              min = NULL, max = NULL, ...) {
  .Deprecated("gate_flowclust_2d")
  gate_flowclust_2d(fr, xChannel, yChannel, filterId, K, usePrior, prior, trans
                    , min.count, max.count, nstart, plot, target, transitional
                    , quantile, translation, transitional_angle, min, max, ...)
}


#' @export 
flowClust.2d <- gate_flowClust_2d

#' Determine the cutpoint by the events quantile.
#' 
#' It is possible that the cutpoint calculated by quantile function may not 
#' produce the exact the probability set by 'probs' argument if there are not enough
#' cell events to reach that precision. Sometime the difference could be significant.
#' 
#' @name gate_quantile
#' @aliases quantileGate
#' @param fr a \code{flowFrame} object
#' @param channel the channel from which the cytokine gate is constructed
#' @param filterId the name of the filter
#' @param plot whether to plot the gate result
#' @param probs probabilities passed to 'stats::quantile' function.
#' @param min a numeric value that sets the lower boundary for data filtering
#' @param max a numeric value that sets the upper boundary for data filtering
#' @param ... additional arguments passed to 'stats::quantile' function.
#' @return a \code{rectangleGate} 
#' @export
#' @examples
#' \dontrun{
#'  gate <- gate_quantile(fr, Channel = "APC-A", probs = 0.995) # fr is a flowFrame
#' }
gate_quantile <- function(fr, channel, probs = 0.999, plot = FALSE,
                         filterId = "", min = NULL, max = NULL, ...) {
   if (missing(channel) || length(channel) != 1) {
     stop("A single channel must be specified.")
   }
   
   # Filter out values less than the minimum and above the maximum, if they are
   # given.
   if (!(is.null(min) && is.null(max))) {
     fr <- .truncate_flowframe(fr, channels = channel, min = min,
         max = max)
   }
                       
  x <- exprs(fr)[, channel]
  cutpoint <- quantile(x, probs = probs, ...)

  gate_coordinates <- list(c(cutpoint, Inf))
  names(gate_coordinates) <- channel
  
  if (plot) {
    hist(x)
    plot(density(x))
    abline(v = cutpoint, col = "red")
    actual_probs <- sum(x>cutpoint)/length(x)
    text(y = 0.5, x = cutpoint, labels = paste("%:", actual_probs))
  }
  rectangleGate(gate_coordinates, filterId = filterId)
}
#' @templateVar old quantileGate
#' @templateVar new gate_quantile
#' @template template-depr_pkg
NULL
#' @export
quantileGate <- function(fr, channel, probs = 0.999, plot = FALSE,
                        filterId = "", min = NULL, max = NULL, ...){
  .Deprecated("gate_quantile")
  gate_quantile(fr, channel, probs, plot, filterId, min, max, ...)
}

#' Determines a cutpoint as the minimum point of a kernel density estimate
#' between two peaks
#'
#' We fit a kernel density estimator to the cells in the \code{flowFrame} and
#' identify the two largest peaks. We then
#' select as the cutpoint the value at which the minimum density is attained
#' between the two peaks of interest.
#'
#' In the default case, the two peaks of interest are the two largest peaks
#' obtained from the \code{link{density}} function. 
#'
#' In the special case that there is only one peak, we are conservative and set
#' the cutpoint as the \code{min(x)} if \code{positive} is \code{TRUE}, and the
#' \code{max(x)} otherwise.
#'
#' @name gate_mindensity
#' @aliases mindensity
#' @param fr a \code{flowFrame} object
#' @param channel TODO
#' @param filterId TODO
#' @param positive If \code{TRUE}, then the gate consists of the entire real
#' line to the right of the cutpoint. Otherwise, the gate is the entire real
#' line to the left of the cutpoint. (Default: \code{TRUE})
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
#'  gate <- gate_mindensity(fr, channel = "APC-A") # fr is a flowFrame
#' }
#' @importFrom flowCore %in% identifier filterDetails<-
gate_mindensity <- function(fr, channel, filterId = "", positive = TRUE,
                       gate_range = NULL, min = NULL, max = NULL,
                       peaks = NULL, ...) {
  
  if (missing(channel) || length(channel) != 1) {
    stop("A single channel must be specified.")
  }

  # Filter out values less than the minimum and above the maximum, if they are
  # given.
  if (!(is.null(min) && is.null(max))) {
    fr <- .truncate_flowframe(fr, channels = channel, min = min,
                                     max = max)
  }
  # Grabs the data matrix that is being gated.
  x <- exprs(fr)[, channel]
  
  if(is.null(peaks))
    peaks <- .find_peaks(x, ...)[, "x"]
  
  if (is.null(gate_range)) {
    gate_range <- c(min(x), max(x))
  } else {
    gate_range <- sort(gate_range)
  }

  
  
  # In the special case that there is only one peak, we are conservative and set
  # the cutpoint as min(x) if 'positive' is TRUE, and max(x) otherwise.
  if (length(peaks) == 1) {
    cutpoint <- ifelse(positive, gate_range[1], gate_range[2])
  } else {
    # The cutpoint is the deepest valley between the two peaks selected. In the
    # case that there are no valleys (i.e., if 'x_between' has an insufficient
    # number of observations), we are conservative and set the cutpoint as the
    # minimum value if 'positive' is TRUE, and the maximum value otherwise.
    valleys <- try(.find_valleys(x, ...), silent = TRUE)
    valleys <- .between_interval(x = valleys, interval = gate_range)

    if (any(is.na(valleys))) {
    #FIXME:currently it is still returning the first peak,
    #we want to pass density instead of x_between to 'min'
    #because x_between is the signal values
      cutpoint <- ifelse(positive, gate_range[1], gate_range[2])
    } else if (length(valleys) == 1) {
      cutpoint <- as.vector(valleys)
    } else if (length(valleys) > 1) {
      # If there are multiple valleys, we determine the deepest valley between
      # the two largest peaks.
      peaks <- sort(peaks[1:2])
      cutpoint <- .between_interval(valleys, peaks)[1]

      # If none of the valleys detected are between the two largest peaks, we
      # select the deepest valley.
      if (is.na(cutpoint)) {
        cutpoint <- valleys[1]
      }      
    }
  }
  gate_coordinates <- ifelse(positive, list(c(cutpoint, Inf)), list(c(-Inf, cutpoint)))
  
  names(gate_coordinates) <- channel
  
  rectangleGate(gate_coordinates, filterId = filterId)
  
}
#' @export
mindensity <- gate_mindensity


#' sequential quadrant gating function
#' 
#' The order of 1d-gating is determined so that the gates better capture the 
#' distributions of flow data.
#' 
#' @name gate_quad_sequential
#' @aliases quadGate.seq
#' @param fr \code{flowFrame}
#' @param channels \code{character} two channels used for gating
#' @param gFunc the name of the 1d-gating function to be used for either dimension
#' @param min a numeric vector that sets the lower bounds for data filtering
#' @param max a numeric vector that sets the upper bounds for data filtering
#' @param ... other arguments passed to \code{.find_peak} (e.g. 'num_peaks' and 'adjust'). see \link{tailgate}
#' @return a \code{filters} that contains four rectangleGates
#' @export 
gate_quad_sequential <- function(fr, channels, gFunc, min = NULL, max = NULL, ...){
  if (missing(channels) || length(channels) != 2) {
    stop("two channels must be specified.")
  }
  
  # Filter out values less than the minimum and above the maximum, if they are
  # given.
  if (!(is.null(min) && is.null(max))) {
    fr <- .truncate_flowframe(fr, channels = channels, min = min,max = max)
  }
  
  res <- sapply(channels, function(channel){
        
        x <- exprs(fr)[, channel]
        
        peaks <-.find_peaks(x,..., num_peaks = 2)
        #get peak and scores
        if(nrow(peaks)<2)
          score <- 0
        else
          score <- abs(diff(peaks[, "x"]))
        list(peaks = peaks[, "x"], score = score)
      }, simplify = FALSE)
  
  #order channel by scores
  channels.ordered <- names(sort(sapply(res, `[[`, "score"), decreasing = T))
  x.first <- all(channels == channels.ordered)
    
  #get first cut
  chnl <- channels.ordered[1]
  thisFunc <- function(fr, chnl, ...){
    thisCall <- substitute(f(data, channel
                          , peaks = p
                        , ...)
                    ,list(f=as.symbol(gFunc), data = fr
                        , channel = chnl
                      , p = res[[chnl]][["peaks"]]
                    )
    )  
  
    eval(thisCall)
  }
  
  g1 <- thisFunc(fr, chnl, ...)
  
  coord1 <- c(g1@min, g1@max)
  cut.1 <- coord1[!is.infinite(coord1)]
#  if(length(cut.1) == 0)
#    cut.1 <- Inf
  
  
  
  #get ind
  fres <- filter(fr, g1)
  ind <- as(fres, "logical")
  
  #gate on the second channel
  chnl <- channels.ordered[2]
  
  #+
  g2 <- thisFunc(fr[ind, ], chnl, ...)
  coord2 <- c(g2@min, g2@max)
  cut.2.pos <- coord2[!is.infinite(coord2)]
  
#  if(length(cut.2.pos) == 0)
#    cut.2.pos<- Inf
  #-
  g3 <- thisFunc(fr[!ind, ], chnl, ...)
  coord3 <- c(g3@min, g3@max)
  cut.2.neg <- coord3[!is.infinite(coord3)]
  
  if(x.first)
    coords <-list(q1 = list(c(-Inf, cut.1), c(cut.2.neg, Inf))
                  , q2 = list(c(cut.1, Inf), c(cut.2.pos, Inf))
                  , q3 = list(c(cut.1, Inf), c(-Inf, cut.2.pos))
                  , q4 = list(c(-Inf, cut.1), c(-Inf, cut.2.neg))
                  )
  else
    coords <-list(q1 = list(c(-Inf, cut.2.pos), c(cut.1, Inf))
                  , q2 = list(c(cut.2.pos, Inf), c(cut.1, Inf))
                  , q3 = list(c(cut.2.neg, Inf), c(-Inf, cut.1))
                  , q4 = list(c(-Inf, cut.2.neg), c(-Inf, cut.1))
                )
    
  
  gates <- lapply(coords, function(coord){
            
      names(coord) <- as.character(channels)
      rectangleGate(coord)
    })
  
  filters(gates)
}
#' @templateVar old quadGate.seq
#' @templateVar new gate_quad_sequential
#' @template template-depr_pkg
NULL
#' @export 
quadGate.seq <- function(fr, channels, gFunc, min = NULL, max = NULL, ...){
  .Deprecated("gate_quad_sequential")
  gate_quad_sequential(fr, channels, gFunc, min, max, ...)
}

##TODO: come up a more general design so that polygonGates can be derived from any two quadrants
#' quadGate based on flowClust::tmixFiler
#' 
#' This gating method identifies two quadrants (first, and third quadrants) by fitting the data with tmixture model.
#' It is particually useful when the two markers are not well resolved thus the regular quadGate method
#' based on 1d gating will not find the perfect cut points on both dimensions.
#' 
#' @name gate_quad_tmix
#' @aliases quadGate.tmix
#' @param fr \code{flowFrame}
#' @param channels \code{character} vector specifies two channels
#' @param usePrior see \link{gate_flowclust_2d}
#' @param K see \link{gate_flowclust_2d}
#' @param prior see \link{gate_flowclust_2d}
#' @param trans see \link{gate_flowclust_2d}
#' @param plot logical whether to plot flowClust clustering results
#' @param quantile1 \code{numeric} specifies the  quantile level(see 'level' in \link{flowClust}) for the first quadrant (x-y+)
#' @param quantile3 \code{numeric} specifies the  quantile level see 'level' in \link{flowClust} for third quadrant (x+y-)
#' @param ... other arguments passed to \link{flowClust}
#' @return a \code{filters} object that contains four \code{polygonGate}s following the order of (-+,++,+-,--)
#' @export 
gate_quad_tmix <- function(fr, channels, K, usePrior = "no", prior = list(NA)
    , quantile1 = 0.8, quantile3 = 0.8
    , trans = 0
    , plot = FALSE
    , ...){
  
  max.v <- .Machine$integer.max #.Machine$double.xmax double max no longer to be valid coordinates in the plot
  min.v <- - max.v
# browser()  
  # fit tmix models 
  tmix_filter <- tmixFilter(parameters = channels
      , K = K
      , usePrior = usePrior, prior = prior, trans = trans
      , ...)

  tmix_results <- try(filter(fr, tmix_filter), silent = F)
  if(plot)
    plot(tmix_results, data = fr, level = 0.85)
  #find the 1st and 3rd quads
  fitted_means <- getEstimates(tmix_results)$locations
  q1.ind <- which.max(fitted_means[,2])
  q3.ind <- which.max(fitted_means[,1])
  
  #construct polygon gates 
  q1.gate <- as(.getEllipseGate(filter = tmix_results
          , include = q1.ind,
          quantile = quantile1
          ,trans = 0), "polygonGate")
  q3.gate <- as(.getEllipseGate(filter = tmix_results
          , include = q3.ind,
          quantile = quantile3
          ,trans = 0), "polygonGate")
  #find intersection of 1,3 quad gates
  q1.bottom <- min(q1.gate@boundaries[, 2])
  q1.right <- max(q1.gate@boundaries[, 1])
  
  q3.top <- max(q3.gate@boundaries[, 2])
  q3.left <- min(q3.gate@boundaries[, 1])
  
  #two gates overlap on both dimensions
  p.BL <- c(q3.left, q1.bottom)
  p.TR <- c(q1.right, q3.top)
  if(q1.bottom >= q3.top){
    thisY <- mean(c(q1.bottom , q3.top))
    p.BL[2] <- p.TR[2] <- thisY
  }
  
  if(q3.left >= q1.right){
    thisX <- mean(c(q1.right , q3.left))
    p.BL[1] <- p.TR[1] <- thisX
  }
  
  
  #construct four quadGates based on these two points
  q1.coord <- matrix(c(min.v, max.v
          , min.v, p.BL[2]
          , p.BL
          , p.TR
          , p.TR[1], max.v
      )
      , 5, 2
      , byrow = TRUE
  )  
  colnames(q1.coord) <- channels    
  q1.g <- polygonGate(q1.coord, filterId = paste(paste0(channels, c("-", "+")), collapse = ""))
  
  
  q2.coord <- matrix(c(p.TR
          , max.v, p.TR[2]
          , max.v, max.v
          , p.TR[1], max.v
      )
      , 4, 2
      , byrow = TRUE
  )  
  colnames(q2.coord) <- channels    
  q2.g <- polygonGate(q2.coord, filterId = paste(paste0(channels, c("+", "+")), collapse = ""))
  
  q3.coord <- matrix(c(p.TR
          , max.v, p.TR[2]
          , max.v, min.v
          , p.BL[1], min.v
          ,p.BL
      )
      , 5, 2
      , byrow = TRUE
  )
  colnames(q3.coord) <- channels    
  q3.g <- polygonGate(q3.coord, filterId = paste(paste0(channels, c("+", "-")), collapse = ""))
  
  q4.coord <- matrix(c(min.v, p.BL[2]
          , p.BL
          , p.BL[1], min.v
          , min.v, min.v
      )
      , 4, 2
      , byrow = TRUE
  )
  colnames(q4.coord) <- channels    
  q4.g <- polygonGate(q4.coord, filterId = paste(paste0(channels, c("-", "-")), collapse = ""))
  
  filters(list(q1.g, q2.g, q3.g,q4.g))
}
#' @templateVar old quadGate.tmix
#' @templateVar new gate_quad_tmix
#' @template template-depr_pkg
NULL
#' @export 
quadGate.tmix <- function(fr, channels, K, usePrior = "no", prior = list(NA)
                          , quantile1 = 0.8, quantile3 = 0.8
                          , trans = 0
                          , plot = FALSE
                          , ...){
  .Deprecated("gate_quad_tmix")
  quad_gate_tmix(fr, channels, K, usePrior, prior, quantile1, quantile3, trans, plot, ...)
}
