setMethod("Gating2D", signature = c("flowSet"), definition = function(obj, ...) {
  .Gating2D.flowSet(obj, ...)
})

#setMethod("Gating2D", signature = c("workFlow"), definition = function(obj, parent, ...) {
#  .Gating2D.flowSet(Data(get(parent, obj)), ...)
#})


#setMethod("Gating2D", signature = c("list"), definition = function(obj, ...) {
#  .Gating2D.list(obj, ...)
#})


##by default the flowset is split into groups by FCS filename indicating sample specific gating
##if groupBy is set as subject id (like "pub_id"),then samples are collapsed and generate subject speicific gating
##in general, gating can be performed at any arbitery granule level based on groupBy arugment
.Gating2D.flowSet <- function(obj, x, y, groupBy = "name", type = "multicore",
                              nslaves = NULL, collapse = FALSE, ...) {
  # guess the x channel name
  xChnl <- getChannelMarker(obj[[1]], x)
  stainX <- xChnl$name

  # guess the y channel name
  yChnl <- getChannelMarker(obj[[1]], y)
	stainY <- yChnl$name

  f1 <- eval(substitute(pData(obj)$f, list(f = groupBy)))
  if (!collapse) {
    if(class(obj) == "ncdfFlowSet") {
      fslist <- split(obj, f1)@datalist
    } else {
      fslist <- split(obj, f1)
    }

    if (any(grepl("parallel", loadedNamespaces()))
        && (is.null(nslaves) || nslaves > 1)) {
      if (is.null(nslaves)) {
        nslaves <- min(length(fslist), parallel::detectCores())
      }
      message("Running in parallel mode with ", nslaves, " nodes.")
      if (type == "multicore") {
        fres.list <- mclapply(fslist, .Gating2D, xChannel = stainX, yChannel = stainY,
                              mc.cores = nslaves, ...)
      } else {
        cl <- parallel::makeCluster(nslaves, type = "SOCK")
        clusterExport(cl, "flowClust.2d")
        fres.list <- parallel::parLapply(cl, fslist, .Gating2D, xChannel = stainX,
                                       yChannel = stainY, ...)
        parallel::stopCluster(cl)
      }
	  } else {
      fres.list <- lapply(fslist, .Gating2D, xChannel = stainX, yChannel = stainY, ...)
    }
    res <- filterList(fres.list)
  } else {
    ## one common gate from all samples
    res <- .Gating2D(obj, xChannel = stainX, yChannel = stainY, ...)
  }

  res
}

#.Gating2D.list <- function(obj, xChannel, yChannel, alpha = "min", groupBy = NULL,
#                           isNcdf = FALSE, plot = FALSE, plotType = "xyplot") {
#
#  stop("This function is stubbed and not yet implemented.")
#  NULL
#}


.Gating2D <- function(nc, xChannel, yChannel, isflowClust = FALSE, absolute = FALSE, ...) {
  if (!require('flowStats', quietly = TRUE)) {
    stop("The package 'flowStats' is needed but not installed.")
  }
  if (isflowClust) {
    if (!require('flowClust', quietly = TRUE)) {
      stop("The package 'flowClust' is needed but not installed.")
    }
    flowClust.2d(fr = nc[[1]], xChannel = xChannel, yChannel = yChannel, ...)
  } else {
    quadrantGate(nc, stain = stains, absolute = FALSE, inBetween = TRUE, ...)
  }
}
