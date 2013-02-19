#' Elicits data-driven priors from a flowSet object for specified channels
#'
#' We elicit data-driven prior parameters from a \code{flowSet} object for
#' specified channels. For each sample in the \code{flowSet} object, we apply the
#' given \code{method} to elicit the priors parameters.
#' 
#' Currently, we have implemented only two methods. In the case that one channel
#' is given, we use the kernel-density estimator (KDE) approach for each sample
#' to obtain \code{K} peaks from which we elicit prior parameters. Otherwise,
#' if more than one channel is specified, we apply K-Means to each of the samples
#' in the \code{flowSet} and aggregate the clusters to elicit the prior
#' parameters.
#'
#' @param flow_set a \code{flowSet} object
#' @param channels a character vector containing the channels in the
#' \code{flowSet} from which we elicit the prior parameters for the Student's t
#' mixture
#' @param K the number of mixture components to identify
#' @param nu0 prior degrees of freedom of the Student's t mixture components.
#' @param w0 the number of prior pseudocounts of the Student's t mixture components.
#' @param ... Additional arguments passed to the prior elicitation method selected
#' @return list of the necessary prior parameters
prior_flowClust <- function(flow_set, channels, method = c("kmeans"), K = 2,
                            nu0 = 4, w0 = 10, ...) {

  if (length(channels) == 1) {
    prior_list <- prior_flowClust1d(flow_set = flow_set, channel = channels,
                                    K = K, nu0 = nu0, w0 = w0, ...)
  } else {
    method <- match.arg(method)
    if (method == "kmeans") {
      prior_list <- prior_kmeans(flow_set = flow_set, channels = channels, K = K,
                                 nu0 = nu0, w0 = w0, ...)
    }
  }

  prior_list
}


#' Elicits data-driven priors from a flowSet object for a specified channel
#'
#' We elicit data-driven prior parameters from a \code{flowSet} object for a
#' specified channel. For each sample in the \code{flowSet} object, we apply a
#' kernel-density estimator (KDE) and obtain \code{K} peaks. We then aggregate
#' these peaks to elicit a prior mean and variance for each of the \code{K}
#' mixture components. We elicit a vague prior for the variance of each mixture
#' component by computing the sample variance of each sample and then averaging
#' these values.
#'
#' Some samples in the \code{flowSet} object may not have \code{K} distinct peaks
#' apparent in the KDEs, in which case we will find less than the specified
#' number with the remaining values being set to \code{NA}. For any sample that
#' has less than \code{K} peaks, we align its peaks with those of the first
#' sample that has \code{K} peaks.
#'
#' We recommend that the Huber estimators be used to aggregate the values. This
#' option can be set in the \code{estimator} argument. Alternatively, the maximum
#' likelihood estimators (MLE) under normality can be used.
#'
#' In the case that the majority of the samples have a peak with NA (i.e., K is
#' overspecified for these samples), the 'huber' function can throw the following
#' error. "Cannot estimate scale: MAD is zero for this sample." In this case, we
#' throw a warning and use vague priors for the problematic peaks.
#'
#' To ensure that the KDEs are smooth, we recommend that the bandwidth set in the
#' \code{adjust} argument be sufficiently large. We have defaulted this value to
#' 1.25. If the bandwidth is not large enough, the KDE may contain numerous
#' bumps, resulting in erroneous peaks.
#'
#' @param flow_set a \code{flowSet} object
#' @param channel the channel in the \code{flowSet} from which we elicit the
#' prior parameters for the Student's t mixture
#' @param K the number of mixture components to identify
#' @param nu0 prior degrees of freedom of the Student's t mixture components.
#' @param w0 the number of prior pseudocounts of the Student's t mixture components.
#' @param adjust the bandwidth to use in the kernel density estimation. See
#' \code{\link{density}} for more information.
#' @param estimator the method to aggregate the samples in the \code{flow_set}.
#' See details.
#' @return list of the necessary prior parameters
prior_flowClust1d <- function(flow_set, channel, K = 2, nu0 = 4, w0 = 10,
                              adjust = 1, estimator = c("huber", "mle")) {

  estimator <- match.arg(estimator)

  # For each sample in the flow_set, we find the K peaks after smoothing the
  # data. We also compute the variance of the 
  estimates <- lapply(seq_along(flow_set), function(i) {
    # Grabs the channel data for the ith sample
    x <- exprs(flow_set[[i]])[, channel]
    peaks <- sort(find_peaks(x, peaks = K, adjust = adjust), na.last = TRUE)

    if (estimator == "huber") {
      # In the case that the MAD is 0, an error results, in which case we use the
      # MLE instead.
      variance <- try(huber(x)$s^2, silent = TRUE)
      if (class(variance) == "try-error") {
        variance <- var(x)
      }
    } else {
      variance <- var(x)
    }
    
    list(peaks = peaks, variance = variance)
  })

  peaks <- do.call(rbind, lapply(estimates, function(x) x$peaks))
  variances <- sapply(estimates, function(x) x$variance)

  # For each sample that has less than K peaks, we align its peaks with the peaks
  # of the first sample having K peaks.
  peaks <- align_peaks(peaks)

  if (estimator == "huber") {
    # In the case that the majority of the samples have a peak with NA (i.e., K
    # is overspecified for these samples), the 'huber' function can throw the
    # following error. "Cannot estimate scale: MAD is zero for this sample."
    # In this case, we throw a warning and use vague priors for the problematic peaks.
    huber_out <- try(apply(peaks, 2, huber), silent = TRUE)
    if (class(huber_out) == "try-error") {
      warning("The number of peaks is overspecified. Using vague priors.")

      huber_out <- lapply(seq_len(K), function(k) {
        huber_peak <- try(huber(peaks[, k]), silent = TRUE)
        if (class(huber_peak) == "try-error") {
          huber_peak <- list(mu = mean(peaks[, k], na.rm = TRUE),
               s = mean(variances))
        }
        huber_peak
      })
    }
    Mu0 <- sapply(huber_out, function(x) x$mu)
    Omega0 <- sapply(huber_out, function(x) x$s^2)
    Lambda0 <- huber(variances)$mu
  } else {
    Mu0 <- colMeans(peaks, na.rm = TRUE)
    Omega0 <- apply(peaks, 2, var, na.rm = TRUE)
    Lambda0 <- mean(variances, na.rm = TRUE)
    # If any of the variances in Omega0 are NA, we keep them vague by setting them to the value in Lambda0
    Omega0[is.na(Omega0)] <- Lambda0
  }

  # Mu0 dimensions: K x p (p is the number of features. Here, p = 1 for 1D)
  Mu0 <- matrix(Mu0, K, 1)

  # Lambda0 dimensions: K x p x p
  Lambda0 <- array(Lambda0, c(K, 1, 1))

  # Omega0 dimensions: K x p x p
  Omega0 <- array(Omega0, c(K, 1, 1))

  # We assume that the degrees of freedom is the same for each mixture component.
  nu0 <- rep(nu0, K)

  # We assume that the prior probability of mixture component membership is equal
  # across all mixture components and use 10 pseudocounts.
  w0 <- rep(w0, K)

  list(Mu0 = Mu0, Lambda0 = Lambda0, Omega0 = Omega0, nu0 = nu0, w0 = w0)
}

#' Elicits data-driven priors from a flowSet object for specified channels using
#' the K-Means clustering algorithm
#'
#' We elicit data-driven prior parameters from a \code{flowSet} object for
#' specified channels. For each sample in the \code{flowSet} object, we apply
#' \code{kmeans} to obtain \code{K} clusters. From each cluster, we determine its
#' centroid and the sample covariance matrix. We then aggregate these two sample
#' moments across all samples for each cluster. 
#'
#' Because the cluster labels returned from \code{kmeans} are arbitrary, we align
#' the clusters based on the centroids that are closest to a randomly selected
#' reference sample. We apply the Hungarian algorithm implemented using the
#' \code{solve_LSAP} function from the \code{clue} package to assist with the
#' alignment.
#'
#' If each frame within \code{flow_set} has a large number of cells, the
#' computational costs of \code{kmeans} can be a burden. We provide the option
#' to randomly select \code{pct}, a percentage of the cells from each flow frame
#' to which \code{kmeans} is applied.
#'
#' @param flow_set a \code{flowSet} object
#' @param channels a character vector containing the channels in the
#' \code{flowSet} from which we elicit the prior parameters for the Student's t
#' mixture
#' @param K the number of mixture components to identify
#' @param nu0 prior degrees of freedom of the Student's t mixture components.
#' @param w0 the number of prior pseudocounts of the Student's t mixture
#' components.
#' @param nstart number of random starts used by \code{kmeans} algorithm
#' @param pct percentage of randomly selected cells in each \code{flowFrame}
#' that is used to elicit the prior parameters. The value should must be greater
#' than 0 and less than or equal to 1.
#' @param ... Additional arguments passed to \code{kmeans}
#' @return list of \code{flowClust} prior parameters
prior_kmeans <- function(flow_set, channels, K, nu0 = 4, w0 = 10, nstart = 1,
                         pct = 0.1, ...) {
  require('clue')

  # For each randomly selected sample in the flow_set, we apply K-means with to
  # find K clusters and retain additional summary statistics to elicit the prior
  # parameters for flowClust.
  num_samples <- length(flow_set)
  
  kmeans_summary <- lapply(seq_len(num_samples), function(i) {
    # Grabs the channel data for the ith sample and randomly selects the
    # specified percentage of cells to which we apply 'kmeans'
    x <- exprs(flow_set[[i]])[, channels]
    x <- x[sample(seq_len(nrow(x)), pct * nrow(x)), ]
    kmeans_out <- kmeans(x = x, centers = K, nstart = nstart) #, ...)

    cluster_centroids <- kmeans_out$centers
    cluster_sizes <- kmeans_out$size

    cluster_covs <- tapply(seq_len(nrow(x)), kmeans_out$cluster, function(i) {
      cov(x[i, ])
    })

    list(centroids = cluster_centroids, sizes = cluster_sizes, covs = cluster_covs)
  })

  # We randomly select one of the flowFrame's within the flow_set as a reference
  # sample.
  ref_sample <- sample.int(length(kmeans_summary), 1)
  kmeans_ref_sample <- kmeans_summary[[ref_sample]]
  kmeans_ref_sample$covs <- unname(kmeans_ref_sample$covs)
  kmeans_ref_sample$covs <- aperm(simplify2array(kmeans_ref_sample$covs), c(3, 1, 2))
  ref_centroids <- kmeans_ref_sample$centroids

  # Because the cluster labels returned from 'kmeans' are arbitrary, we align the
  # clusters based on the centroids that are closest to the reference sample.
  kmeans_summary <- kmeans_summary[-ref_sample]
  kmeans_centroids <- lapply(kmeans_summary, function(x) x$centroids)
  
  # For each sample, we obtain the alignment indices.
  ref_indices <- lapply(kmeans_centroids, function(centroids) {
    dist_centroids <- as.matrix(dist(rbind(ref_centroids, centroids)))
    dist_centroids <- dist_centroids[seq_len(K), seq_len(K) + K]
    as.vector(solve_LSAP(dist_centroids))
  })

  # Now, we align the 'kmeans' summary information using the alignment indices.
  kmeans_summary <- mapply(function(kmeans_sample, align_idx) {
    kmeans_sample$centroids <- kmeans_sample$centroids[align_idx, ]
    kmeans_sample$sizes <- kmeans_sample$sizes[align_idx]
    kmeans_sample$covs <- unname(kmeans_sample$covs[align_idx])

    # Converts covariance matrices to a 3D array (K x p x p)
    kmeans_sample$covs <- aperm(simplify2array(kmeans_sample$covs), c(3, 1, 2))

    kmeans_sample
  }, kmeans_summary, ref_indices, SIMPLIFY = FALSE)

  # Appends the reference sample to the 'kmeans' summary list.
  kmeans_summary <- c(kmeans_summary, list(kmeans_ref_sample))

  # Mu0 dimensions: K x p
  # Prior Mean
  # p is the number of features, corresponding to the number of 'channels'
  # We average each of the cluster centroids across all samples to elicit Mu0.
  p <- length(channels)
  Mu0 <- Reduce("+", lapply(kmeans_summary, function(x) x$centroids)) / num_samples
  Mu0 <- matrix(Mu0, K, p)

  # Lambda0 dimensions: K x p x p
  # Prior Covariance Matrix
  # For each cluster, we pool the covariance matrices from all samples to elicit
  # Lambda0.
  kmeans_sizes <- do.call(rbind,lapply(kmeans_summary, function(x) x$sizes))
  kmeans_sizes_sum <- colSums(kmeans_sizes)
  kmeans_scatter <- lapply(kmeans_summary, function(x) x$sizes * x$covs)
  Lambda0 <- Reduce("+", kmeans_scatter) / kmeans_sizes_sum
  Lambda0 <- array(Lambda0, c(K, p, p))

  # Omega0 dimensions: K x p x p
  # Hyperprior covariance matrix for the prior mean Mu0
  # For each cluster we calculate the covariance matrix of the cluster centroids
  # across all samples to elicit Omega0.
  kmeans_centroids <- lapply(kmeans_summary, function(x) x$centroids)
  centroids_split <- lapply(kmeans_centroids, split, f = seq_len(K))
  Omega0 <- lapply(seq_len(K), function(k) {
    cov(do.call(rbind, lapply(centroids_split, function(x) x[[k]])))
  })
  Omega0 <- array(do.call(rbind, lapply(Omega0, as.vector)), dim = c(K, p, p))

  # We assume that the degrees of freedom is the same for each mixture component.
  nu0 <- rep(nu0, K)

  # We assume that the prior probability of mixture component membership is equal
  # across all mixture components and use 10 pseudocounts.
  w0 <- rep(w0, K)

  list(Mu0 = Mu0, Lambda0 = Lambda0, Omega0 = Omega0, nu0 = nu0, w0 = w0)
}

#' Peak alignment for prior elicitation
#'
#' We elicit data-driven priors by applying a kernel-density estimator (KDE) to
#' obtain \code{K} peaks for several samples. Some samples may not have \code{K}
#' distinct peaks apparent in the KDEs, in which case we will find less than the
#' specified number with the remaining values being set to \code{NA}. For any
#' sample that has less than \code{K} peaks, we align its peaks with those of the
#' first sample that has \code{K} peaks.
#'
#' We apply the Hungarian algorithm implemented using the \code{solve_LSAP}
#' function from the \code{clue} package.
#'
#' @param peaks matrix. The rows corresponds to the samples, and the columns
#' correspond to the peaks. A value of \code{NA} is used if no peak is present.
#' @return matrix where the peaks are aligned
align_peaks <- function(peaks) {
  require('clue')
  K <- ncol(peaks)
  
  # For each sample that has less than K peaks, we align its peaks with the peaks
  # of the first sample having K peaks.
  num_peaks <- apply(peaks, 1, function(x) sum(!is.na(x)))
  peaks_align <- peaks[min(which(num_peaks == K)), ]

  aligned_peaks <- lapply(which(num_peaks < K), function(i) {
    x <- peaks[i, ]
    x <- x[!is.na(x)]

    dist_align <- as.matrix(dist(c(x, peaks_align)))
    dist_align[is.na(dist_align)] <- 0
    dist_align <- dist_align[seq_along(x), seq_len(K) + length(x)]

    # The 'solve_LSAP' function expects a matrix. In the case that only one peak
    # is present, 'dist_align' is a vector. We coerce it to a matrix to prevent
    # an error with 'solve_LSAP'.
    if (is.vector(dist_align)) {
      dist_align <- matrix(dist_align, nrow = 1)
    }

    # We extract the indices of alignment and then return a vector with the
    # aligned peaks
    align_idx <- as.vector(solve_LSAP(dist_align))
    replace(rep(NA, K), align_idx, x)
  })
  aligned_peaks <- do.call(rbind, aligned_peaks)
  peaks[num_peaks < K] <- aligned_peaks
  peaks
}
    
