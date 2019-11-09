#' @templateVar old tautString
#' @templateVar new tautstring
#' @template template-depr_pkg
NULL

#' @rdname gate_tautstring
#' @aliases tautString
#' @title Taut String Density Estimator Gating
#' @param sorted_vector \code{numeric} vector of single cell expresion from a single cytometric channel.
#' @param modeprior \code{numeric} scalar specifying the expected number of modes. Default 0 (autodetect). Rarely should this be
#' set by the user.
#' @export
tautstring <- function(sorted_vector, modeprior = 0){
  tsGates(xVec = sorted_vector, modePrior = modeprior )
}
#' @export 
tautString <- function(sorted_vector, modeprior = 0){
  .Deprecated("tautstring")
  tautstring(sorted_vector, modeprior)
}

#'@name gate_tautstring
#'@aliases tautStringGate
#'@description The taut string density estimator gating returns 0, 1, or more gates, depending on how many modes it identifies in the data.
#'@param fr a flowFrame object
#'@param channel The channel to gate.
#'@param gate_range The range to look for a gate, no truncation occurs.
#'@param min The min range of the data to truncate the flowFrame
#'@param max The max range of the data to truncate the flowFrame
#'@param filterId The id / name of the gate.
#'@export
gate_tautstring <- function(fr, channel, gate_range = NULL, min = NULL, max = NULL, filterId = "") {
  if (missing(channel) || length(channel) != 1) {
    stop("One channel must be specified for a tautString gate.")
  }
  if (!(is.null(min) && is.null(max))) {
    fr <- .truncate_flowframe(fr, channels = channel, min = min,
                              max = max)
  }
  sorted_vector <- sort(exprs(fr[,channel])[,1])
  if (singleDip(sorted_vector) < 0.25) {
      gate_locations <- tautString(sorted_vector)
  } else {
    gate_locations <- tautString(sorted_vector, modeprior = 1)
  }
  l <- length(gate_locations)
  if (l == 2) {
    # We did not find any gates.
    gate_def <- list(c(-Inf,Inf))
    names(gate_def) <- channel
    filters <- list(rectangleGate(gate_def, filterId = filterId))
  } else if (l > 2) {
    gate_locations <- gate_locations[-c(1, l)]
    gate_locations <- c(-Inf,gate_locations,Inf)
    filters <- list()
    for (i in seq_along(gate_locations)[-l]) {
      for (j in i + 1) {
        gate_def <- list(c(gate_locations[i], gate_locations[j]))
        names(gate_def) <- channel
        filters <- c(filters, rectangleGate(gate_def, filterId = filterId))
      }
    }
  } else {
    filters <- list()
  }
  return(filters)
}
#' @templateVar old tautStringGate
#' @templateVar new gate_tautstring
#' @template template-depr_pkg
NULL

#'@export
tautStringGate <- function(fr, channel, gate_range = NULL, min = NULL, max = NULL, filterId = ""){
  .Deprecated("gate_tautstring")
  gate_tautstring(fr, channel, gate_range, min, max, filterId)
}