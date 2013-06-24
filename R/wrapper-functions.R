# This file contains all wrapper methods for dispatching data and arguments to
# gating algorithms.

#'  wrapper for prior_flowClust

.prior_flowClust <- function(fs, gs, gm, xChannel, yChannel
                                , prior_source = NULL
                                , k, neg, pos
                                , min, max, ...){
    prior_list <- list()
  
  # prior estimation is done separately from flowClust routine because
  # prior_flowClust1d requires the entire parent flowSet yet flowClust only
  # takes one flowFrame
    
    if (is.null(prior_source)) {
      prior_data <- fs
    } else {
      prior_data <- getData(gs, prior_source)
    }
#    browser()
    if (names(gm) == "flowClust.1d") {
      # If any of 'K, 'neg' or 'pos' are given, we extract them from the
      # 'args' and coerce them as integers. Otherwise, they are defaulted to
      # NULL.
      # NOTE: the named elements within 'args' are coerced to lowercase.
      K <- neg_cluster <- pos_cluster <- NULL
      
      if(!missing(k)) {
        K <- as.integer(k)
      }
      if(!missing(neg)) {
        neg_cluster <- as.integer(neg)
      }
      if(!missing(pos)) {
        pos_cluster <- as.integer(pos)
      }
      
      # Examines the various cases of pos, neg and K.
      # In the case that K is not given, we automatically select it.
      # Exception: If both 'pos' and 'neg' are given, we set: K <- pos + neg
      if (is.null(K)) {
        if (!is.null(neg_cluster) && !is.null(pos_cluster)) {
          K <- neg_cluster + pos_cluster
        }
      } else {
        # If pos and neg are given, throw a warning and set: K <- pos + neg
        # In the case that one of pos and neg are given and either exceeds K,
        # we throw an error.         
        if (!is.null(neg_cluster) && !is.null(pos_cluster)) {
          warning("Values given for 'K', 'neg', and 'pos'. Setting K = neg + pos")
          K <- neg_cluster + pos_cluster
        } else if (!is.null(pos_cluster) && K < pos_cluster) {
          stop("The number of positive clusters exceeds 'K'.")
        } else if (!is.null(neg_cluster) && K < neg_cluster) {
          stop("The number of negative clusters exceeds 'K'.")
        }
      }
      
      # If 'min' and/or 'max' are given, we pass this value along to the
      # prior-elicitation method as well as flowClust. Otherwise, these values
      # are set to NULL.
      min_values <- -Inf
      max_values <- Inf
      if (!missing(min)) {
        min_values <- as.numeric(min)
      }
      if (!missing(max)) {
        max_values <- as.numeric(max)
      }                          
      
      if (!is.na(xChannel)) {
        prior_list[[xChannel]] <- prior_flowClust(flow_set = prior_data,
            channels = xChannel, K = K,
            min = min_values, max = max_values, ...)
      }
      
      prior_list[[yChannel]] <- prior_flowClust(flow_set = prior_data,
          channels = yChannel, K = K,
          min = min_values, max = max_values, ...)
      
    } else {
      # get the value of neg and pos
      if (!missing(k)) {
        K <- as.integer(k)
      } else {
        message("'K' argument is missing! Using default setting: K = 2")
        K <- 2
      }
      
      prior_list <- prior_flowClust(flow_set = prior_data, channels = c(xChannel, 
              yChannel), K = K,  ...)
    }
    
    prior_list
}


## wrappers for the different gating routines
.singletGate <- function(fs, xChannel = "FSC-A", yChannel = "FSC-H",
                         prediction_level = 0.99, pp_res = NULL, ...) {
  require(openCyto)
  # Creates a list of polygon gates based on the prediction bands at the minimum
  # and maximum x_channel observation using a robust linear model trained by
  # flowStats.
  singletGate(fs[[1]], area = xChannel, height = yChannel,
              prediction_level = prediction_level)
}

.flowClust.1d <- function(fs, pp_res, xChannel = NA, yChannel, tol = 1e-5,
                          filterId = "", ...) {
  require(openCyto)
  
  prior <- pp_res
  
  sname <- sampleNames(fs)
  fr <- fs[[sname]]
  priorList <- list()
  if(!is.null(prior))
    prior <- prior[[yChannel]]
  
  if (is.na(xChannel)) {
    # 1d gate
    flowClust.1d(fr = fr, params = yChannel, tol = tol, filterId = filterId, 
      prior = prior, ...)
  } else {
    stop("flowClust1d does not support 2d gate!")
  }
}

.cytokine <- function(fs, pp_res, xChannel = NA, yChannel = "FSC-A", filterId = "", 
                      ...) {
  require(openCyto)
  cytokine(flow_set = fs, channel = yChannel, filter_id = filterId, ...)
}

.mindensity <- function(fs, pp_res, yChannel = "FSC-A", filterId = "", ...) {
  require(openCyto)
  # TODO: Iterate through the flowFrames within 'fs', given that 'fs' may
  # contain more than one flowFrame if 'split' is specified in the CSV file.
  mindensity(flow_frame = fs[[1]], channel = yChannel, filter_id = filterId, ...)
}

.flowClust.2d <- function(fs, pp_res, xChannel, yChannel, usePrior = "yes", 
                          ...) {
  require(openCyto)
  prior <- pp_res
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

.rangeGate <- function(fs, pp_res, xChannel = NA, yChannel, absolute = FALSE, filterId = "",
                       ...) {
  require(openCyto)
  # TODO: Iterate through the flowFrames within 'fs', given that 'fs' may
  # contain more than one flowFrame if 'split' is specified in the CSV file.
  fr <- fs[[1]]
  rangeGate(x = fr, stain = yChannel, inBetween = TRUE, absolute = absolute,
            filterId = filterId, ...)
}

.quantileGate <- function(fs, pp_res, xChannel = NA, yChannel, probs = 0.999, filterId = "",
                          ...) {
  require(openCyto)
  # TODO: Iterate through the flowFrames within 'fs', given that 'fs' may
  # contain more than one flowFrame if 'split' is specified in the CSV file.
  fr <- fs[[1]]
  quantileGate(fr = fr, probs = probs, stain = yChannel, filterId = filterId, ...)
}

.quadrantGate <- function(fs, pp_res, xChannel = NA, yChannel, ...) {
  require(openCyto)
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
 
