#' @templateVar old prior_flowClust
#' @templateVar new prior_flowclust
#' @template template-depr_pkg
NULL
#' Elicits data-driven priors from a flowSet object for specified channels
#'
#' We elicit data-driven prior parameters from a \code{flowSet} object for
#' specified channels. For each sample in the \code{flowSet} object, we apply the
#' given \code{prior_method} to elicit the priors parameters.
#' 
#' Currently, we have implemented only two methods. In the case that one channel
#' is given, we use the kernel-density estimator (KDE) approach for each sample
#' to obtain \code{K} peaks from which we elicit prior parameters. Otherwise,
#' if more than one channel is specified, we apply K-Means to each of the samples
#' in the \code{flowSet} and aggregate the clusters to elicit the prior
#' parameters.
#'
#' In the rare case that a prior covariance matrix is singular, we shrink the
#' eigenvalues of the matrix slightly to ensure that it is positive definite. For
#' instance, if the \code{flow_set} has two samples, this case can occur. The
#' amount of shrinkage is controlled in \code{shrink}.
#'
#' @name prior_flowclust
#' @aliases prior_flowClust
#' @param flow_set a \code{flowSet} object
#' @param channels a character vector containing the channels in the
#' \code{flowSet} from which we elicit the prior parameters for the Student's t
#' mixture
#' @param prior_method the method to elicit the prior parameters
#' @param K the number of mixture components to identify
#' @param nu0 prior degrees of freedom of the Student's t mixture components.
#' @param w0 the number of prior pseudocounts of the Student's t mixture components. (only the first element is used and the rest is ignored at the moment)
#' @param shrink the amount of eigenvalue shrinkage to add in the case the prior
#' covariance matrices are singular. See details.
#' @param ... Additional arguments passed to the prior elicitation method selected
#' @return list of the necessary prior parameters
#' @export 
#' @examples 
#' \dontrun{
#' library(flowCore)
#' data(GvHD)
#' prior_flowclust(GvHD[1:3], c("FSC-H", "SSC-H"))
#' }
prior_flowclust <- function(flow_set, channels, prior_method = c("kmeans"),
                            K = 2, nu0 = 4, w0 = c(10,10), shrink = 1e-6, ...) {
  #pass only the first element of w0 since it will be replicated ..
  #that said, the second element never gets used at this moment
  if (length(channels) == 1) {
    prior_list <- .prior_flowClust1d(flow_set = flow_set, channel = channels,
                                    K = K, nu0 = nu0, w0 = w0[1], ...)
  } else {
    prior_method <- match.arg(prior_method)
    if (prior_method == "kmeans") {
      prior_list <- .prior_kmeans(flow_set = flow_set, channels = channels, K = K,
                                 nu0 = nu0, w0 = w0[1], ...)  
    }
    # In the rare case a covariance matrix is singular, we shrink the eigenvalues
    # of the matrix. The amount of shrinkage is controlled in 'shrink'.

    prior_list$Lambda0 <- .aaply(prior_list$Lambda0, 1, function(cov_mat) {
      if (is(try(solve(cov_mat), silent = TRUE), "try-error")) {
        cov_mat <- cov_mat + shrink * diag(nrow(cov_mat))
      }
      cov_mat
    })
    prior_list$Lambda0 <- unname(prior_list$Lambda0)
    
    #Lambda0 should be Kxpxp.. but it's reduced to pxp when K=1
    #correct the dimensions here
    L0<-prior_list$Lambda0
    if(length(dim(L0))==2){
      dim(L0)<-c(1,dim(L0))
      prior_list$Lambda0<-L0
    }
    
    prior_list$Omega0 <- .aaply(prior_list$Omega0, 1, function(cov_mat) {
      if (class(try(solve(cov_mat), silent = TRUE)) == "try-error") {
        cov_mat <- cov_mat + shrink * diag(nrow(cov_mat))
      }
      cov_mat
    })
    prior_list$Omega0 <- unname(prior_list$Omega0)
  }
 
  prior_list
   
}

#' @export
prior_flowClust <- function(flow_set, channels, prior_method = c("kmeans"),
                            K = 2, nu0 = 4, w0 = c(10,10), shrink = 1e-6, ...){
  .Deprecated("prior_flowclust")
  prior_flowclust(flow_set, channels, prior_method, K, nu0, w0, shrink, ...)
}

#' Elicits data-driven priors from a flowSet object for a specified channel
#'
#' We elicit data-driven prior parameters from a \code{flowSet} object for a
#' specified channel. For each sample in the \code{flowSet} object, we apply a
#' kernel-density estimator (KDE) and identify its local maxima (peaks). 
#' We then aggregate these peaks to elicit a prior
#' parameters for each of \code{K} mixture components.
#'
#' Here, we outline the approach used for prior elicitation. First, we apply a
#' KDE to each sample and extract all of its peaks (local maxima). It is
#' important to note that different samples may have a different number of
#' peaks. Our goal then is to align the peaks before aggregating the information
#' across all samples.  To do this, we utilize a technique similar to the peak
#' probability contrasts (PPC) method from Tibshirani et al (2004). Effectively,
#' we apply hierarchical clustering to the peaks from all samples to find
#' clusters of peaks. We compute the sample mean and variance of the peaks
#' within each cluster to elicit the prior means and its hyperprior variance,
#' respectively, for a \code{\link{flowClust}} mixture component. We elicit the
#' prior variance for each mixture component by first assigning the observations
#' within each sample to the nearest prior mean. Then, we compute the variance
#' of the observations within each cluster. Finally, we average the variances
#' corresponding to each mixture component across all samples in the
#' \code{flowSet} object.
#'
#' Following Tibshirani et al. (2004), we cluster the peaks from each sample
#' using complete-linkage hierarchical clustering. The linkage type can be
#' changed via the \code{hclust_method} argument. This argument is passed
#' directly to \code{\link{hclust}}.
#'
#' To cluster the peaks, we must cut the hierarchical tree by selecting either a
#' value for \code{K} or by providing a height of the tree to cut. By default,
#' we cut the tree using as the height the median of the distances between
#' adjacent peaks within each sample. This value can be changed via the
#' \code{hclust_height} argument and, if provided, will be passed to
#' \code{\link{cutree}}. Also, by default, the number of mixture components
#' \code{K} is \code{NULL} and is ignored.  However, if \code{K} is provided,
#' then it has priority over \code{hclust_height} and is passed instead directly
#' to \code{\link{cutree}}.
#'
#' To ensure that the KDEs are smooth, we recommend that the bandwidth set in
#' the \code{adjust} argument be sufficiently large. We have defaulted this
#' value to 2. If the bandwidth is not large enough, the KDE may contain
#' numerous bumps, resulting in erroneous peaks.
#'
#' @references Tibshirani, R et al. (2004), "Sample classification from protein
#' mass spectrometry, by 'peak probability contrasts'," Bioinformatics, 20, 17,
#' 3034-3044. \url{http://bioinformatics.oxfordjournals.org/content/20/17/3034}.
#' @param flow_set a \code{flowSet} object
#' @param channel the channel in the \code{flowSet} from which we elicit the
#' prior parameters for the Student's t mixture
#' @param K the number of mixture components to identify. By default, this value
#' is \code{NULL} and determined automatically
#' @param hclust_height the height of the \code{\link{hclust}} tree of peaks,
#' where the should be cut By default, we use the median of the distances
#' between adjacent peaks. If a value is specified, we pass it directly to
#' \code{\link{cutree}}.
#' @param clust_method the method used to cluster peaks together when for prior
#' elicitation. By default, \code{kmeans} is used. However, if \code{K} is not
#' specified, \code{hclust} will be used instead.
#' @param hclust_method the agglomeration method used in the hierarchical
#' clustering. This value is passed directly to \code{\link{hclust}}. Default is
#' complete linkage.
#' @param artificial a numeric vector containing prior means for artificial
#' mixture components. The remaining prior parameters for the artificial
#' components are copied directly from the most informative prior component
#' elicited. If \code{NULL} (default), no artificial prior components are added.
#' @param nu0 prior degrees of freedom of the Student's t mixture components.
#' @param w0 the number of prior pseudocounts of the Student's t mixture components.
#' @param adjust the bandwidth to use in the kernel density estimation. See
#' \code{\link{density}} for more information.
#' @param min a numeric value that sets the lower bound for data filtering. If
#' \code{NULL} (default), no truncation is applied.
#' @param max a numeric value that sets the upper bound for data filtering. If
#' \code{NULL} (default), no truncation is applied.
#' @param vague \code{logical} Whether to elicit a vague prior. If \code{TRUE}, we first calculate the median of standard
#'                              deviations from all flowFrames. Then, we divide the overall standard
#'                              deviation by the number of groups to the scale the standard deviation.
#' @return list of prior parameters
#' @rdname prior_flowClust1d
#' @noRd 
.prior_flowClust1d <- function(flow_set, channel, K = NULL, hclust_height = NULL,
                              clust_method = c("kmeans", "hclust"),
                              hclust_method = "complete", artificial = NULL,
                              nu0 = 4, w0 = 10, adjust = 2, min = -200,
                              max = NULL, vague = TRUE) {

  channel <- as.character(channel)
  clust_method <- match.arg(clust_method)
  if (length(channel) != 1) {
    stop("There can be only 1...channel.")
  }

  if (!(is.null(min) && is.null(max))) {
    flow_set <- .truncate_flowset(flow_set, channels = channel, min = min, max = max)
  }
  # For each sample in 'flow_set', we identify the peaks after smoothing.
  peaks <- fsApply(flow_set, function(flow_frame, adjust) {
    x <- exprs(flow_frame)[, channel]
    peaks_found <- .find_peaks(x, adjust = adjust)[, "x"]
    
    # If K is specified and is smaller than the number of peaks found,
    # we keep only the K largest peaks from the sample.
    if (!is.null(K) && length(peaks_found) > K) {
      peaks_found <- peaks_found[seq_len(K)]
    }
    peaks_found
  }, adjust = adjust, simplify = FALSE)
  peaks_collapsed <- as.vector(do.call(c, peaks))

  # Remove any NAs that may have resulted in the peaks.
  peaks_collapsed <- peaks_collapsed[!is.na(peaks_collapsed)]

  # For each sample, we sort the peaks and find the distance between the
  # adjacent peaks.
  peaks_dist <- lapply(peaks, function(x) {
    diff(sort(x))
  })
  peaks_dist <- do.call(c, peaks_dist)

  if (is.null(K) || clust_method == "hclust") {
    # Applies hierarchical clustering to the peaks.
    hclust_out <- hclust(dist(peaks_collapsed), method = hclust_method)

    # To cluster the peaks, we must cut the hierarchical tree. By default, we
    # use the median of the distances between adjacent peaks.
    # By default, 'K' is 'NULL' and is ignored by 'cutree'. If 'K' is provided,
    # then it has priority over 'hclust_height'.
    if (is.null(hclust_height)) {
      hclust_height <- median(peaks_dist)
    }
    clust_labels <- factor(cutree(hclust_out, k = K, h = hclust_height))
    prior_means <- as.numeric(tapply(peaks_collapsed, clust_labels, mean))
  } else if (clust_method == "kmeans") {
    # Applies kmeans clustering to the peaks.
    kmeans_out <- kmeans(x = peaks_collapsed, centers = K, nstart = 10,algorithm="MacQueen")
    clust_labels <- factor(kmeans_out$cluster)
    prior_means <- as.vector(kmeans_out$centers)
  }
  K <- nlevels(clust_labels)

  if (vague) {
    prior_means <- sort(prior_means)

    # To elicit a vague prior, we first calculate the median of standard
    # deviations from all flowFrames. Then, we divide the overall standard
    # deviation by the number of groups to the scale the standard deviation. 
    sd_x <- fsApply(flow_set, function(flow_frame) {
      sd(exprs(flow_frame)[, channel])
    })
    sd_x <- median(sd_x, na.rm = TRUE) / K
    var_x <- sd_x^2

    hyperprior_vars <- rep.int(var_x, K)
    prior_vars <- rep.int(var_x, K)
                   
  } else {
    # For each cluster of peaks, we elicit the hyperprior variance (Omega0) for
    # the prior mean Mu0 by computing the variance of each cluster.
    hyperprior_vars <- as.numeric(tapply(peaks_collapsed, clust_labels, var))

    # Because the cluster labels assigned by the hierarchical clustering algorithm
    # are arbitrary, the 'prior_means' are not necessarily sorted. Here, we sort
    # them and then apply the new odering to the other prior parameters.
    prior_order <- order(prior_means)
    prior_means <- prior_means[prior_order]
    hyperprior_vars <- hyperprior_vars[prior_order]

    # Elicitation of prior variances for each mixture component:
    # For each flowFrame in the original flowSet object, we cluster its
    # observations by assigning them to the nearest mean in 'prior_means.' Then,
    # we compute the variance of the observations within each cluster. Finally, we
    # aggregate the variances across all of the flowFrame objects.
    prior_vars <- fsApply(flow_set, function(flow_frame, prior_means) {
      x <- exprs(flow_frame)[, channel]

      # To determine the nearest mean, we find the midpoints of the peaks and then
      # 'cut' the observations into the regions.
      peak_midpoints <- rowMeans(embed(prior_means, 2))
      labels <- cut(x, breaks = c(-Inf, peak_midpoints, Inf))
      tapply(x, labels, var)
    }, prior_means = prior_means)
    prior_vars <- as.numeric(colMeans(prior_vars, na.rm = TRUE))
  }

  # Here, we add any 'artificial' prior components, if provided, and increment
  # 'K' accordingly.
  if (!is.null(artificial)) {
    prior_df <- cbind.data.frame(Mu0 = prior_means, Omega0 = hyperprior_vars,
                                 Lambda0 = prior_vars)

    # For the remaining prior parameters for each artificial component, we select
    # minimum values of Omega0 and Lambda0.
    min_Omega0 <- min(prior_df$Omega0)
    min_Lambda0 <- min(prior_df$Lambda0)
    prior_df <- rbind(prior_df, cbind(Mu0 = artificial, Omega0 = min_Omega0,
                                      Lambda0 = min_Lambda0))

    prior_df <- prior_df[order(prior_df$Mu0), ]

    prior_means <- prior_df$Mu0
    hyperprior_vars <- prior_df$Omega0
    prior_vars <- prior_df$Lambda0

    K <- length(prior_means)
  }

  # Mu0 dimensions: K x p (p is the number of features. Here, p = 1 for 1D)
  Mu0 <- matrix(prior_means, K, 1)

  # Lambda0 dimensions: K x p x p
  Lambda0 <- array(prior_vars, c(K, 1, 1))

  # Omega0 dimensions: K x p x p
  Omega0 <- array(hyperprior_vars, c(K, 1, 1))

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
#' @param min a numeric vector that sets the lower bounds for data filtering. If
#' \code{NULL} (default), no truncation is applied.
#' @param max a numeric vector that sets the upper bounds for data filtering. If
#' \code{NULL} (default), no truncation is applied.
#' @param ... Additional arguments passed to \code{kmeans}
#' @return list of \code{flowClust} prior parameters
#' @rdname prior_kmeans
#' @noRd 
.prior_kmeans <- function(flow_set, channels, K, nu0 = 4, w0 = 10, nstart = 10,
                         pct = 0.1, min = NULL, max = NULL, ...) {

  channels <- as.character(channels)
  
  # Truncates flow_set before eliciting priors when necessary 
  if (!(is.null(min) && is.null(max))) {
    flow_set <- .truncate_flowset(flow_set, channels = channels, min = min,
                                 max = max)
  }
  # For each randomly selected sample in the flow_set, we apply K-means with to
  # find K clusters and retain additional summary statistics to elicit the prior
  # parameters for flowClust.
  num_samples <- length(flow_set)

  # The number of features, corresponding to the number of 'channels'
  p <- length(channels)
  
  kmeans_summary <- lapply(seq_len(num_samples), function(i) {
    # Grabs the channel data for the ith sample and randomly selects the
    # specified percentage of cells to which we apply 'kmeans'
    x <- as.matrix(exprs(flow_set[[i]])[, channels])
    x <- as.matrix(x[sample(seq_len(nrow(x)), pct * nrow(x)), ])

    # If there are too few rows in 'x', 'kmeans' will throw an error.
    # Instead, we do not elicit summary statistics from this sample and return NA.
    if (K >= nrow(x)) {
      warning("The number of clusters 'K' exceeds the number of subsampled rows.",
              "Try increasing the subsample 'pct'.")
      return(NA)
    }
    
    kmeans_out <- kmeans(x = x, centers = K, nstart = nstart, algorithm="MacQueen",...)

    cluster_centroids <- kmeans_out$centers
    cluster_sizes <- kmeans_out$size

    cluster_covs <- tapply(seq_len(nrow(x)), kmeans_out$cluster, function(i) {
      # In the case that a singleton cluster results, a vague covariance matrix
      # is reported: the covariance matrix of the entire data set.
      if (length(i) >= 2) {
        cov(as.matrix(x[i, ]))
      } else {
        warning("A singleton cluster resulted. Using vague covariance matrix.",
                "Try increasing the subsample 'pct'.")
        cov(x)
      }
    })

    list(centroids = cluster_centroids, sizes = cluster_sizes, covs = cluster_covs)
  })

  # Removes samples marked with NA that had too few observations.
  kmeans_summary <- kmeans_summary[!is.na(kmeans_summary)]

  # We select the flowFrame with the largest number of observations as the
  # reference sample.
  sample_sizes <- sapply(kmeans_summary, function(x) sum(x$sizes))
  ref_sample <- which.max(sample_sizes)
  kmeans_ref_sample <- kmeans_summary[[ref_sample]]
  kmeans_ref_sample$covs <- simplify2array(unname(kmeans_ref_sample$covs))

  # Converts covariance matrices to a 3D array (K x p x p)
  # We must treat the case of 1 channel as a special case to ensure the
  # dimensions of the 3-dimensional array is correct.
  if (p > 1) {
    kmeans_ref_sample$covs <- aperm(kmeans_ref_sample$covs, c(3, 1, 2))
  } else {
    kmeans_ref_sample$covs <- array(kmeans_ref_sample$covs, dim = c(p, 1, 1))
  }
  ref_centroids <- kmeans_ref_sample$centroids

  # Because the cluster labels returned from 'kmeans' are arbitrary, we align the
  # clusters based on the centroids that are closest to the reference sample.
  kmeans_summary <- kmeans_summary[-ref_sample]
  kmeans_centroids <- lapply(kmeans_summary, "[[", "centroids")
  
  # For each sample, we obtain the alignment indices.
    ref_indices <- lapply(kmeans_centroids, function(centroids) {
      dist_centroids <- as.matrix(dist(rbind(ref_centroids, centroids)))
      dist_centroids <- dist_centroids[seq_len(K), seq_len(K) + K,drop=FALSE]
      as.vector(solve_LSAP(dist_centroids))
    })
  # Now, we align the 'kmeans' summary information using the alignment indices.
  kmeans_summary <- mapply(function(kmeans_sample, align_idx) {
    kmeans_sample$centroids <- kmeans_sample$centroids[align_idx, ]
    kmeans_sample$sizes <- kmeans_sample$sizes[align_idx]
    kmeans_sample$covs <- unname(kmeans_sample$covs[align_idx])

    # Converts covariance matrices to a 3D array (K x p x p) As above,
    # we must treat the case of 1 channel as a special case to ensure
    # the dimensions of the 3-dimensional array is correct.
    kmeans_sample$covs <- simplify2array(kmeans_sample$covs)
    if (p > 1) {
      kmeans_sample$covs <- aperm(kmeans_sample$covs, c(3, 1, 2))
    } else {
      kmeans_sample$covs <- array(kmeans_sample$covs, dim = c(p, 1, 1))
    }

    kmeans_sample
  }, kmeans_summary, ref_indices, SIMPLIFY = FALSE)

  # Appends the reference sample to the 'kmeans' summary list.
  kmeans_summary <- c(kmeans_summary, list(kmeans_ref_sample))

  # Calculates the total number of observations assigned to each cluster across
  # all samples.
  kmeans_sizes <- do.call(rbind, lapply(kmeans_summary, "[[", "sizes"))
  kmeans_sizes_sum <- colSums(kmeans_sizes)

  # For each cluster after alignment, we extract a matrix of the centroids from
  # each sample.
  kmeans_centroids <- lapply(kmeans_summary, "[[", "centroids")
  centroids_split <- lapply(kmeans_centroids, split, f = seq_len(K))
  centroids_split <- lapply(seq_len(K), function(k) {
    do.call(rbind, lapply(centroids_split, function(x) x[[k]]))
  })

  # Next, for each cluster, we compute summary statistics of the centroids from
  # each sample. The statistics are weighted by the number of samples. This is
  # especially imporant when the number of samples is small and some of them
  # have few observations.
  weighted_stats <- lapply(centroids_split, cov.wt, wt = sample_sizes,
                           method = "ML")

  # Mu0 dimensions: K x p
  # Prior Mean
  # We aggregate the cluster centroids across all samples weighted by the
  # cluster sample sizes to elicit Mu0.
  centroids_aggregate <- lapply(weighted_stats, "[[", "center")
  Mu0 <- do.call(rbind, centroids_aggregate)
  Mu0 <- matrix(Mu0, K, p)

  # Omega0 dimensions: K x p x p
  # Hyperprior covariance matrix for the prior mean Mu0
  # For each cluster we calculate the weighted covariance matrix of the cluster
  # centroids across all samples to elicit Omega0.
  Omega0 <- lapply(weighted_stats, "[[", "cov")
  Omega0 <- array(do.call(rbind, lapply(Omega0, as.vector)), dim = c(K, p, p))

  # Lambda0 dimensions: K x p x p
  # Prior Covariance Matrix
  # For each cluster, we pool the covariance matrices from all samples to elicit
  # Lambda0.
  kmeans_scatter <- lapply(kmeans_summary, function(x) x$sizes * x$covs)
  Lambda0 <- Reduce("+", kmeans_scatter) / kmeans_sizes_sum
  Lambda0 <- array(Lambda0, c(K, p, p))

  # We assume that the degrees of freedom is the same for each mixture component.
  nu0 <- rep(nu0, K)

  # We assume that the prior probability of mixture component membership is equal
  # across all mixture components and use 10 pseudocounts.
  w0 <- rep(w0, K)

  list(Mu0 = Mu0, Lambda0 = Lambda0, Omega0 = Omega0, nu0 = nu0, w0 = w0)
}

#' Finds the local maxima (peaks) in the given vector after smoothing the data
#' with a kernel density estimator.
#'
#' First, we smooth the data using kernel density estimation (KDE) with the
#' \code{\link{density}} function. Then, we find every local maxima such that the
#' density is concave (downward).
#'
#' The \code{num_peaks} argument returns only the largest peaks. If
#' \code{num_peaks} is greater than the number of peaks found, then all the peaks
#' are returned, and a warning is issued.
#'
#' Effectively, we find the local maxima with a discrete analogue to a second
#' derivative applied to the KDE. For details, see this StackOverflow post:
#' \url{http://bit.ly/Zbl7LV}.
#'
#' @param x numeric vector
#' @param y numeric vector. If given, it is treated as the density values for
#' \code{x}. The length of \code{y} must equal that of \code{x}.
#' @param num_peaks the number of peaks to find. By default, all peaks are
#' returned. See details.
#' @param adjust the bandwidth to use in the kernel density estimation. See
#' \code{\link{density}} for more information.
#' @param ... additional arguments passed to the \code{\link{density}} function
#' @return a \code{data.frame} that contains the peaks(and their density heights) attained. The peaks are sorted in
#' descending order based on the density heights.
#' @examples
#' library(flowClust)
#' set.seed(42)
#' # 2 peaks with a minor hump
#' y <- flowStats::SimulateMixture(10000, c(.5, .3, .2), c(2, 5, 7), c(1, 1, 1), nu = 10)
#' plot(density(y))
#' peaks <- .find_peaks(y)
#' abline(v = peaks[, "x"], col = "red")
#' @noRd 
.find_peaks <- function(x, y = NULL, num_peaks = NULL, adjust = 2, plot = FALSE, ...) {
  x <- as.vector(x)

  if (length(x) < 2) {
    warning("At least 2 observations must be given in 'x' to find peaks.")
    return(NA)
  }

  if (is.null(y)) {
    dens <- density(x, adjust = adjust, ...)
  } else {
    y <- as.vector(y)
    if (length(x) != length(y)) {
      stop("The lengths of 'x' and 'y' must be equal.")
    }
    dens <- list(x = x, y = y)
  }

  # Discrete analogue to a second derivative applied to the KDE. See details.
  second_deriv <- diff(sign(diff(dens$y)))
  which_maxima <- which(second_deriv == -2) + 1

  # The 'density' function can consider observations outside the observed range.
  # In rare cases, this can actually yield peaks outside this range.  We remove
  # any such peaks.
  which_maxima <- which_maxima[findInterval(dens$x[which_maxima], range(x)) == 1]

  # Next, we sort the peaks in descending order based on the density heights.
  which_maxima <- which_maxima[order(dens$y[which_maxima], decreasing = TRUE)]
  
  # Returns the local maxima. If there are none, we return 'NA' instead.
  if (length(which_maxima) > 0) {
    peaks <- dens$x[which_maxima]
    if (is.null(num_peaks) || num_peaks > length(peaks)) {
      num_peaks <- length(peaks)
    }
    peaks <- peaks[seq_len(num_peaks)]
  } else {
    peaks <- NA
  }
  
  peaks <- data.frame(x = peaks, y = dens$y[which_maxima][seq_len(num_peaks)])
  if(plot){
    plot(dens, main = paste("adjust =" ,  adjust))
    points(peaks, ,col = "red")  
  }
  
  peaks  
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
#' @param y numeric vector. If given, it is treated as the density values for
#' \code{x}. The length of \code{y} must equal that of \code{x}.
#' @param num_valleys the number of valleys to find. By default, all valleys are
#' returned. See details.
#' @param adjust the bandwidth to use in the kernel density estimation. See
#' \code{\link{density}} for more information.
#' @param ... additional arguments passed to the \code{\link{density}} function
#' @return the values where the valleys are attained.
#' @examples
#' library(flowClust)
#' set.seed(42)
#' # 3 peaks and 2 valleys
#' y <- flowStats::SimulateMixture(10000, c(.25, .5, .25), c(1, 5, 9), c(1, 1, 1), nu = 10)
#' plot(density(y))
#' valleys <- .find_valleys(y)
#' abline(v = valleys, col = "red")
#' @noRd 
.find_valleys <- function(x, y = NULL, num_valleys = NULL, adjust = 2, ...) {

  x <- as.vector(x)

  if (length(x) < 2) {
    warning("At least 2 observations must be given in 'x' to find valleys.")
    return(NA)
  }
  
  if (is.null(y)) {
    dens <- density(x, adjust = adjust, ...)
  } else {
    y <- as.vector(y)
    if (length(x) != length(y)) {
      stop("The lengths of 'x' and 'y' must be equal.")
    }
    dens <- list(x = x, y = y)
  }

  # Discrete analogue to a second derivative applied to the KDE. See details.
  second_deriv <- diff(sign(diff(dens$y)))
  which_minima <- which(second_deriv == 2) + 1

  # The 'density' function can consider observations outside the observed range.
  # In rare cases, this can actually yield valleys outside this range. We remove
  # any such valleys.
  which_minima <- which_minima[findInterval(dens$x[which_minima], range(x)) == 1]

  # Next, we sort the valleys in descending order based on the density heights.
  which_minima <- which_minima[order(dens$y[which_minima], decreasing = FALSE)]

  # Returns the local minima. If there are none, we return 'NA' instead.
  if (length(which_minima) > 0) {
    valleys <- dens$x[which_minima]
    if (is.null(num_valleys) || num_valleys > length(valleys)) {
      num_valleys <- length(valleys)
    }
    valleys <- valleys[seq_len(num_valleys)]
  } else {
    valleys <- NA
  }
  valleys
}
