#' Estimates a common logicle transformation for a flowSet.
#'
#' Of the negative values for each channel specified, the median of the specified
#' quantiles are used.
#'
#' @param flow_set object of class 'flowSet'
#' @param channels character vector of channels to transform
#' @param m TODO -- default value from .lgclTrans
#' @param q quantile
#' @return TODO
estimateMedianLogicle <- function(flow_set, channels, m = 4.5, q = 0.05) {
  if (!is(flow_set, "flowSet")) {
    stop("flow_set has to be an object of class 'flowSet'")
  }
  if (missing(channels)) {
    stop("Please specify the channels to be logicle transformed")
  }
  indx <- channels %in% unname(colnames(exprs(flow_set[[1]])))
  if (!all(indx)) {
    stop(paste("Channels", channels[!indx], "were not found in flow_set "))
  }

  neg_marker_quantiles <- fsApply(flow_set, function(sample) {
    apply(exprs(sample), 2, function(markers) {
      quantile(markers[markers < 0], probs = q)
    })
  })
  # Replaces 'r' in flowCore:::.lgclTrans
  neg_marker_quantiles <- apply(neg_marker_quantiles, 2,
                                median, na.rm = TRUE)[channels]

  # In the case that no negative markers are present, we set this quantile to the
  # default value of 1/2.
  neg_marker_quantiles <- replace(neg_marker_quantiles,
                                  is.na(neg_marker_quantiles), 0.5)

  # Replaces 't' in flowCore:::.lgclTrans
  max_range <- do.call(rbind, lapply(fsApply(flow_set, range), function(x) {
    x[2, channels]
  }))
  max_range <- apply(max_range, 2, max)

  # Replaces 'w' in flowCore:::.lgclTrans
  w <- (m - log10(max_range / abs(neg_marker_quantiles))) / 2

  transformation <- lapply(channels, function(channel) {
    transId <- paste(channel, "medianLogicleTransform", sep = "_")

    logicleTransform(transformationId = transId, w = w[channel],
                     t = max_range[channel], m = m, a = 0)
  })

  transformList(channels, transformation,
                transformationId = "medianLogicleTransform")
}
