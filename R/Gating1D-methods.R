setMethod("Gating1D", signature = c("flowSet"), definition = function(obj, ...) {
  .Gating1D.flowSet(obj, ...)
})

##by default the flowset is split into groups by FCS filename indicating sample specific gating
##if groupBy is set as subject id (like "pub_id"),then samples are collapsed and generate subject speicific gating
##in general, gating can be performed at any arbitery granule level based on groupBy arugment
.Gating1D.flowSet <- function(obj, x = NULL, y, groupBy, ref, nslaves = NULL,
                              collapse = FALSE, tol = 1e-3, filterId = "",
                              type = "multicore", ...) {

  # guess the y channel name
  yChnl <- getChannelMarker(obj[[1]], y)
	stainY <- yChnl$name

  # guess the x channel name
  if (!missing(x) && !is.null(x)) {
    xChnl <- getChannelMarker(obj[[1]], x)
    stainX <- xChnl$name
  } else {
    stainX <- NULL
  }

  if (missing(groupBy)) {
    f1 <- eval(substitute(name), pData(obj))
  } else {
    f1 <- eval(substitute(groupBy), pData(obj))
  }

  if (!collapse) {
    fslist <- split(obj, f1)

    if (any(grepl("parallel", loadedNamespaces()))
        && (is.null(nslaves) || nslaves > 1)) {
			if (is.null(nslaves)) {
        nslaves <- min(length(fslist), parallel::detectCores())
      }
      message("Running in parallel mode with ", nslaves, " nodes.")
      
      if (type == "multicore") {
        gating_list <- mclapply(fslist, .Gating1D, c(stainX, stainY),
                                mc.cores = nslaves, tol = tol,
                                filterId = filterId, ...)
      } else {
        cl <- parallel::makeCluster(nslaves, type = "SOCK")
        clusterExport(cl, "flowClust.1d")
        gating_list <- parallel::parLapply(cl, fslist, .Gating1D, c(stainX, stainY),
                                           tol = tol, filterId = filterId, ...)
        parallel::stopCluster(cl)
      }  
    } else {
      if (missing(ref)) {
        gating_list <- lapply(fslist, .Gating1D, c(stainX, stainY), tol = tol,
                              filterId = filterId, ...)
      } else {
        ref <- substitute(ref)
        gating_list <- lapply(fslist, function(fs) {
          curPd <- pData(fs)
          # use one particular reference sample within each group to produce gate
          selected <- curPd[eval(ref, curPd), "name"]
          if (length(selected) != 1) {
            stop("invalid number of reference sample: ", selected)
          }
          g <- .Gating1D(fs[selected], c(stainX, stainY), tol = tol,
                         filterId = filterId, ...)

          # replicate the single gate for all flowFrames within the group
          glist <- replicate(nrow(curPd), g)
          names(glist) <- sampleNames(fs)
          glist
        })
        gating_list <- unlist(unname(gating_list))
      }
    }
	
    gating_results <- filterList(gating_list)
  } else {
    # We __collapse__ all samples (i.e., 'flowFrame' objects) in the 'flowSet'
    # object to a single 'flowFrame' object from which we construct a common gate
    # for all samples. We coerce the collapsed 'flowFrame' back to a 'flowSet'
    # for usage with 'flowClust'.
    obj <- flowSet(as(obj, "flowFrame"))

    gating_results <- .Gating1D(obj, c(stainX, stainY), tol = tol,
                                filterId = filterId, ...)

    # Next, we replicate the single gate for all flowFrames into a list and
    # combine them.
    gating_list <- replicate(length(f1), gating_results)
    names(gating_list) <- f1
    gating_results <- filterList(gating_list)
  }

  gating_results
}

.Gating1D <- function(nc, stains, split = TRUE, absolute = FALSE, filterId = "",
                      method = c("flowClust", "flowStats", "density1d", "quadseq", "quantile"),
                      tol, prior = NULL, probs = 0.999, ...) {

	method <- match.arg(method)
  
  if (length(stains) == 1) {
    if (method == "flowClust") {
      if (!require('flowClust', quietly = TRUE)) {
        stop("The package 'flowClust' is needed but not installed.")
      }
	  
      flowClust.1d(fr = nc[[1]], params = stains, tol = tol, filterId = filterId,
                   prior = prior, ...)
    } else if (method == "flowStats") {
      if (!require('flowStats', quietly = TRUE)) {
        stop("The package 'flowStats' is needed but not installed.")
      }
      #it is important to set inBetween as T and abs as FALSE
      #otherwise the anchor based on 25% and 50% of quantiles will be used
      #for determining pos and neg and theoretical range will be used instead of
      #the real data range
      rangeGate(nc, stain = stains, inBetween = TRUE, absolute = absolute,
                filterId = filterId, ...)
    } else if (method == "density1d") {
      if (!require('flowStats', quietly = TRUE)) {
        stop("The package 'flowStats' is needed but not installed.")
      }

      rangeGate(nc, stain = stains, inBetween = TRUE, absolute = absolute,
                filterId = filterId, simple = TRUE, ...)
    } else if (method == "quantile") {
      quantileGate(fr = nc[[1]], probs = probs, stain = stains, filterId = filterId, ...)
    }
  } else {

    if (method == "flowClust") {
      if (!require('flowClust', quietly = TRUE)) {
        stop("The package 'flowClust' is needed but not installed.")
      }

      # If two stains are given, we apply flowClust to each marker and then
      # combine the two gates into a quadrant gate.
      gate_x <- flowClust.1d(fr = nc[[1]], params = stains[1], tol = tol,
                             filterId = as.character(getChannelMarker(nc[[1]], stains[1])$desc),
                             prior = prior$xChannel, ...)
      gate_y <- flowClust.1d(fr = nc[[1]], params = stains[2], tol = tol,
                             filterId = as.character(getChannelMarker(nc[[1]], stains[2])$desc),
                             prior = prior$yChannel, ...)

      gate_x <- ifelse(is.finite(gate_x@min), gate_x@min, gate_x@max)
      gate_y <- ifelse(is.finite(gate_y@min), gate_y@min, gate_y@max)

      flowClust_quadGate <- list(gate_x, gate_y)
      names(flowClust_quadGate) <- stains
      quadGate(filterId = filterId, flowClust_quadGate)
    } else if (method == "flowStats") {
      quadrantGate(nc, stain = stains, absolute = FALSE, inBetween = TRUE, ...)
    } else {
      stop("Invalid method for 2-parameter gating.")
    }
  }
}


