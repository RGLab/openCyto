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

    # We extract the indices of alignment and then return a vector with the
    # aligned peaks
    align_idx <- as.vector(solve_LSAP(dist_align))
    replace(rep(NA, K), align_idx, x)
  })
  aligned_peaks <- do.call(rbind, aligned_peaks)
  peaks[num_peaks < K] <- aligned_peaks
  peaks
}


#' Generates data-driven priors for the Bayesian version of flowClust for 2-D
#' gating.
#'
#' We use vague priors (i.e., scaled identity matrices), where the scales are
#' determined by the covariance matrix of the overall data without regard to the
#' mixture components.
#'
#' To construct a prior for the \code{K} means (i.e., Mu0), we use a common prior
#' mean along the x-channel to be the highest density point. For the y-channel
#' prior mean, we find \code{K} peaks along the y-channel using a kernel-density
#' estimator.
#'
#' @param fr a \code{flowFrame} object
#' @param xChannel TODO
#' @param yChannel TODO
#' @param nu0 prior degrees of freedom of the Student's t mixture components.
#' @param w0 the number of prior pseudocounts.
#' @return list of the necessary prior parameters
prior_flowClust2d <- function(fr, xChannel, yChannel, K = 2, nu0 = 4, w0 = 10, adjust = 2) {

  # From the given flowFrame, grabs the data.
  x <- exprs(fr)[, c(xChannel, yChannel)]
  cov_x <- cov(x)

  # Mu0 dimensions: K x p (p is the number of features. Here, p = 2 for 2D)
  x_peak <- find_peaks(x[, xChannel], peaks = K, adjust = adjust)
  y_peaks <- sort(find_peaks(x[, yChannel], peaks = K, adjust = adjust), na.last = TRUE)

  # HACK: In the event that the number of peaks along the y-axis is overspecified,
  # then the remaining 'y_peaks' will be NA. In this case, we replace the NA peaks
  # with the first peak. Effectively, this sets the prior means as equal.
  # TODO: Devise a better prior elicitation in the future.
  if (any(is.na(y_peaks))) {
    y_peaks[is.na(y_peaks)] <- y_peaks[1]
  }
  
  Mu0 <- cbind(x_peak, y_peaks)

  # We assume that the prior probability of mixture component membership is equal
  # across all mixture components and use 10 pseudocounts.
  w0 <- rep(w0, K)

  # Lambda0 dimensions: K x p x p
  Lambda0 <- aperm(replicate(K, cov_x), perm = c(3, 2, 1))

  # Omega0 dimensions: K x p x p
  Omega0 <- aperm(replicate(K, cov_x), perm = c(3, 2, 1))

  list(Mu0 = Mu0, nu0 = nu0, w0 = w0, Lambda0 = Lambda0, Omega0 = Omega0)
}

#' Prior means elicitation for robust 1D-gating with the Bayesian flowClust.
#'
#' We elicit prior means for 1D-gating using the Bayesian version of
#' \code{flowClust}. For each sample, we apply a kernel-density estimator to each
#' of the specified channels and obtain two peaks. The minimum of these peaks is
#' labeled the 'Negative' peak, while the maximum is labeled the 'Positive' peak.
#' We pool each of the 'Negative' and 'Positive' peaks across the specified
#' channels for each of the given samples to obtain a prior mean for the
#' 'Negative' and 'Positive' peaks, respectively.
#'
#' Some channels will not have two distinct peaks in the kernel-density estimate.
#' By default, we do not consider a channel in the prior-means elicitation if
#' this channel has only one clear peak in any of the samples. If
#' \code{remove_missing_peaks} is set to \code{FALSE}, then the default behavior
#' is altered and all of the channels for which we have two distinct peaks will
#' be included in the prior-means elicitation.
#'
#' To ensure that the kernel-density estimates are smooth, we recommend that the
#' bandwidth set in the \code{adjust} argument be sufficiently large. We have
#' defaulted this to 3. If the bandwidth is not large enough, the density
#' estimate may contain numerous bumps, resulting in erroneous peaks.
#'
#' @param wf the workflow that contains the samples from which the prior means
#' should be calculated
#' @param view_name the workflow view used to elicit the prior means. By default,
#' the prior means are computed conditional on the cellular debris being gated
#' out.
#' @param channels the channels that are used to elicit the prior means
#' @param remove_missing_peaks Should we remove channels that do not have two
#' clearly defined peaks? By default, yes. See Details.
#' @param adjust the bandwidth to use in the kernel density estimation. See
#' \code{\link{density}} for more information.
#' @return list containing the elicited prior means, the kernel-density peaks
#' computed for each marker by sample, and a vector of the markers included
#' in the elication of the prior means.
#' @examples \dontrun{
#' # For a given workflow and view, we use all of the channels in the workflow
#' # that do not include 'Time' and forward- and side-scatter channels.
#' view_name <- 'nonDebris+'
#' wf_colnames <- colnames(Data(wf[[view_name]]))
#' channels <- wf_colnames[!(wf_colnames %in% c("Time", "FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H"))]
#' prior_means_1d(wf = wf, view_name = view_name, channels = channels)
#' }
prior_means_flowClust1d <- function(wf, view_name = 'nonDebris+', channels,
                                    remove_missing_peaks = TRUE, adjust = 3) {
  require('plyr')
  require('reshape2')

  # For each of the samples, grabs the markers
  markers_samples <- lapply(seq_along(Data(wf[[view_name]])), function(i) {
    # For the current view, we grab the ith sample data matrix for the channels
    # that correspond to markers.
    markers_x <- exprs(Data(wf[[view_name]])[[i]])[, channels]

    # Updates the channel names to the corresponding markers.
    colnames(markers_x) <- channels2markers(wf = wf, channels = channels)

    cbind.data.frame(Sample = paste("Sample", i), markers_x)
  })
  markers_samples <- do.call(rbind, markers_samples)

  # The markers are the column names from the above data.frame except for the
  # first column, which denotes the "Sample" number.
  markers <- colnames(markers_samples)[-1]

  # Melts the markers into a ggplot2-friendly data.frame
  melt_markers <- melt(markers_samples, id = "Sample", variable.name = "Marker")

  # For each sample and for each of its markers, we find the highest 2 peaks from
  # a kernel-density estimator.
  marker_peaks <- ddply(melt_markers, .(Sample, Marker), function(x) {
    find_peaks(x$value, peaks = 2, adjust = adjust)
  })
  colnames(marker_peaks) <- c("Sample", "Marker", "Negative", "Positive")

  # Melts the data.frame of peaks for each marker by sample.
  peaks_by_sample <- melt(marker_peaks, id = c("Sample", "Marker"), variable = "Peak")

  if (remove_missing_peaks) {
    # Determines which markers have two clear peaks.
    markers_two_peaks <- ddply(peaks_by_sample, .(Marker), summarize, two_peaks = !any(is.na(value)))
    markers_two_peaks <- as.character(subset(markers_two_peaks, two_peaks == TRUE)$Marker)
    markers <- markers_two_peaks

    # Subsets the peaks for each marker by sample.
    peaks_by_sample <- subset(peaks_by_sample, Marker %in% markers_two_peaks)
  }

  # We calculate the data-driven prior mean as the average of the two peaks across all samples
  # and across all markers for which we were able to find two clear peaks.
  prior_means <- sort(ddply(peaks_by_sample, .(Peak), summarize, avg = mean(value, na.rm = TRUE))$avg)

  list(prior_means = prior_means, peaks_by_sample = peaks_by_sample, markers = markers)
}

