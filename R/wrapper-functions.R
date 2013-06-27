# This file contains all wrapper methods for dispatching data and arguments to
# gating algorithms.

# wrapper for prior_flowClust to return sample-specific priors
.prior_flowClust <- function(flow_set, prior_group = NULL, ...) {
  if (is.null(prior_group)) {
    prior_list <- prior_flowClust(flow_set, ...)
    # replicate the prior for each sample
    sapply(sampleNames(flow_set), function(i) prior_list, simplify = FALSE)
  } else {
    splitBy1 <- factor(pData(flow_set)[, prior_group])
    fs_list <- split(flow_set, splitBy1)
    names(fs_list) <- NULL  #strip group name
    plist <- lapply(fs_list, function(fs) {
      prior_list <- prior_flowClust(fs, ...)

      # replicate the prior for each sample within this prior group
      sapply(sampleNames(fs), function(i) prior_list, simplify = FALSE)
    })
    unlist(plist, recursive = FALSE)
  }
}

## wrappers for the different gating routines
.singletGate <- function(fs, xChannel = "FSC-A", yChannel = "FSC-H",
                         prediction_level = 0.99, ...) {
  require(openCyto)
  fs <- fs[, c(xChannel,yChannel)]
  # Creates a list of polygon gates based on the prediction bands at the minimum
  # and maximum x_channel observation using a robust linear model trained by
  # flowStats.
  singletGate(fs[[1]], area = xChannel, height = yChannel,
              prediction_level = prediction_level)
}

.flowClust.1d <- function(fs, xChannel = NA, yChannel, tol = 1e-5, prior = NULL,
                          filterId = "", split = TRUE, ...) {
  require(openCyto)
  fs <- fs[, yChannel]
  sname <- sampleNames(fs)
  fr <- fs[[sname]]
  priorList <- list()
  priorList[[yChannel]] <- prior[[yChannel]][[sname]]
  
  if (is.na(xChannel)) {
    # 1d gate
    flowClust.1d(fr = fr, params = yChannel, tol = tol, filterId = filterId, 
      prior = priorList[[yChannel]], ...)
  } else {
    stop("flowClust1d does not support 2d gate!")
  }
}

.cytokine <- function(fs, xChannel = NA, yChannel = "FSC-A", filterId = "",
                      ...) {
  require(openCyto)
  fs <- fs[, yChannel]
  cytokine(flow_set = fs, channel = yChannel, filter_id = filterId, ...)
}

.mindensity <- function(fs, yChannel = "FSC-A", filterId = "", ...) {
  require(openCyto)
  fs <- fs[, yChannel]
  # TODO: Iterate through the flowFrames within 'fs', given that 'fs' may
  # contain more than one flowFrame if 'split' is specified in the CSV file.
  mindensity(flow_frame = fs[[1]], channel = yChannel, filter_id = filterId, ...)
}

.flowClust.2d <- function(fs, xChannel, yChannel, usePrior = "yes", prior = NULL,
                          ...) {
  require(openCyto)
  fs <- fs[,c(xChannel, yChannel)]
  sname <- sampleNames(fs)
  #collapse if necessary
  if(length(sname)>1){
    fr <- as(fs,"flowFrame")
  }else{
    fr <- fs[[sname]]  
  }
  
  flowClust.2d(fr = fr, xChannel = xChannel, yChannel = yChannel, usePrior = usePrior
#               ,prior = prior[[sname]]
                ,prior = prior[[1]] #TODO:this is a hack to get collapsed gating work.
                , ...)
}

.rangeGate <- function(fs, xChannel = NA, yChannel, absolute = FALSE, filterId = "", 
                       ...) {
  require(openCyto)
  fs <- fs[, yChannel]
  # TODO: Iterate through the flowFrames within 'fs', given that 'fs' may
  # contain more than one flowFrame if 'split' is specified in the CSV file.
  fr <- fs[[1]]
  rangeGate(x = fr, stain = yChannel, inBetween = TRUE, absolute = absolute,
            filterId = filterId, ...)
}

.quantileGate <- function(fs, xChannel = NA, yChannel, probs = 0.999, filterId = "",
                          ...) {
  require(openCyto)
  fs <- fs[, yChannel]
  # TODO: Iterate through the flowFrames within 'fs', given that 'fs' may
  # contain more than one flowFrame if 'split' is specified in the CSV file.
  fr <- fs[[1]]
  quantileGate(fr = fr, probs = probs, stain = yChannel, filterId = filterId, ...)
}

.quadrantGate <- function(fs, xChannel = NA, yChannel, ...) {
  require(openCyto)
  fs <- fs[,c(xChannel, yChannel)]
  # TODO: Iterate through the flowFrames within 'fs', given that 'fs' may
  # contain more than one flowFrame if 'split' is specified in the CSV file.
  fr <- fs[[1]]
  
  qfilter <- quadrantGate(fr, stain = c(xChannel, yChannel), absolute = FALSE, 
    inBetween = TRUE, ...)
  
  ###############################################################     
  #construct rectangleGates based on the cuts and popNames,clock-wise
  ###############################################################
  cut.x <- qfilter@boundary[xChannel]
  cut.y <- qfilter@boundary[yChannel]
  gateList <- new("filters")
  
  chnls <- c(xChannel, yChannel)
  markers <- chnls
  
  coord <- list(c(-Inf, cut.x), c(cut.y, Inf))
  names(coord) <- as.character(chnls)
  gateList[[paste(paste0(markers, c("-", "+")), collapse = "")]] <- rectangleGate(coord)
  
  coord <- list(c(cut.x, Inf), c(cut.y, Inf))
  names(coord) <- as.character(chnls)
  gateList[[paste(paste0(markers, c("+", "+")), collapse = "")]] <- rectangleGate(coord)
  
  coord <- list(c(cut.x, Inf), c(-Inf, cut.y))
  names(coord) <- as.character(chnls)
  gateList[[paste(paste0(markers, c("+", "-")), collapse = "")]] <- rectangleGate(coord)
  
  coord <- list(c(-Inf, cut.x), c(-Inf, cut.y))
  names(coord) <- as.character(chnls)
  gateList[[paste(paste0(markers, c("-", "-")), collapse = "")]] <- rectangleGate(coord)
  
  gateList
}
 
