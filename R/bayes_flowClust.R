#' Generates data-driven priors for the Bayesian version of flowClust for 1-D
#' gating.
#'
#' We use vague priors (i.e., scaled identity matrices), where the scale is
#' determined by the variance of the overall data without regard to the mixture
#' components.
#'
#' To construct a prior for the \code{K} means (i.e., Mu0), we find \code{K}
#' peaks along the channel specified using a kernel-density estimator.
#'
#' @param fr a \code{flowFrame} object
#' @param channel TODO
#' @param nu0 prior degrees of freedom of the Student's t mixture components.
#' @param adjust the bandwidth to use in the kernel density estimation. See
#' \code{\link{density}} for more information.
#' 
#' @return list of the necessary prior parameters
prior_flowClust1d <- function(fr, channel, K = 2, nu0 = 4, adjust = 3) {

  # From the given flowFrame, grabs the data.
  x <- exprs(fr)[, channel]

  var_x <- var(x)

  # Mu0 dimensions: K x p (p is the number of features. Here, p = 1 for 1D)
  Mu0 <- matrix(sort(find_peaks(x, peaks = K, adjust = adjust)), K, 1)

  # We assume that the prior probability of mixture component membership is equal
  # across all mixture components.
  w0 <- rep(1, K) / K

  # Lambda0 dimensions: K x p x p
  Lambda0 <- array(rep(var_x, K), c(K, 1, 1))

  # Omega0 dimensions: K x p x p
  Omega0 <- array(rep(var_x, K), c(K, 1, 1))

  list(Mu0 = Mu0, nu0 = nu0, w0 = w0, Lambda0 = Lambda0, Omega0 = Omega0)
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
#' @return list of the necessary prior parameters
prior_flowClust2d <- function(fr, xChannel, yChannel, K = 2, nu0 = 4) {

  # From the given flowFrame, grabs the data.
  x <- exprs(fr)[, c(xChannel, yChannel)]
  cov_x <- cov(x)

  # Mu0 dimensions: K x p (p is the number of features. Here, p = 2 for 2D)
  x_peak <- find_peaks(x[, xChannel], peaks = 1)
  y_peaks <- sort(find_peaks(x[, yChannel], peaks = K), na.last = TRUE)

  # HACK: In the event that the number of peaks along the y-axis is overspecified,
  # then the remaining 'y_peaks' will be NA. In this case, we replace the NA peaks
  # with the first peak. Effectively, this sets the prior means as equal.
  # TODO: Devise a better prior elicitation in the future.
  if (any(is.na(y_peaks))) {
    y_peaks[is.na(y_peaks)] <- y_peaks[1]
  }
  
  Mu0 <- cbind(x_peak, y_peaks)

  # We assume that the prior probability of mixture component membership is equal
  # across all mixture components.
  w0 <- rep(1, K) / K

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
