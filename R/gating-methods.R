#################################
## gatingTemplate
#------------------------------
#basic steps applied to all gating strategies:
#1.comp 
#2.trans
#3.boundary filter
#################################
setMethod("gating",signature=c("gatingTemplate","GatingSet"),definition=function(x,wf,plot=F,batch=F,xbin=128,isCompensated=FALSE,...)
    {
#     browser()
      nodeNames<-getNodes(wf[[1]])
      ##################################################################
      ##compensation 
      #------------------------------------------------------------------
      #(one comp across all samples)workflow can't handle sample specific comp 
      ##################################################################
      
#     comp<-getCompensationMatrices(x)
#     marker<-parameters(comp)
#
#     if(!isCompensated&&!"comp"%in%nodeNames)
#     {
#       message("compensating...")
#       
#       
#       add(wf,comp,name="comp")
#       pview<-"comp"
#     }else
#       pview<-"base view"
      
      ##################################################################
      ##transformation
      #------------------------------------------------------------------
      #currently use one sample to estimate logicle trans and apply to all
      #because wf doesn't support multiple trans list
      ##################################################################
#     if(!"trans"%in%nodeNames)
#     {
#       message("transforming...")
#       
#       trans<-getTransformations(x)
##      browser()
#       if(is.null(trans)||length(trans)==0)
#       {
#         message("transformation not defined,using estimated logicle trans...")
#         trans <- estimateLogicle(Data(wf[[pview]])[[1]], channels = marker)
#       }
#       
#       
#       add(wf,trans,parent=pview,name="trans")
#     }
#   browser()
      if(plot)
      {
#     trans<-get(action(wf[["trans"]])@transform@ID,wf)
        #check log scaled 1-d hist to verify compensation is done properly
#     grid.arrange(
#         xyplot(`PerCP-Cy5-5-A`~`APC-H7-A`,transform(Data(wf[["base view"]])[[1]],trans),smooth=F,main="raw",margin=F)
#         ,xyplot(`PerCP-Cy5-5-A`~`APC-H7-A`,Data(wf[["trans"]])[[1]],smooth=F,main="compensated",margin=F)
#         ,densityplot(~`APC-H7-A`,transform(Data(wf[["base view"]])[[1]],trans),main="raw")
#         ,densityplot(~`APC-H7-A`,Data(wf[["trans"]])[[1]],main="compensated")
#         ,ncol=2,nrow=2)
#     if(!batch)
#     {
#       retval<-readline("Enter to resume,c to exit:")
#       dev.off()
#       if(retval=="c")
#       {
#         opt.outer <- options(show.error.messages=FALSE)
#         on.exit(options(opt.outer))
#         
#         stop()
#       }
#     }
      }
#   nodeNames<-ls(alias(wf))
#   
      if(!"boundary+"%in%nodeNames)
      {
        gating.boundary(x,wf,pViewName="trans",...)
      }
#   
    })

gating.boundary<-function(x,wf,pViewName,step=100,...){
  if(step>=1)
  { message("boundary gating...")
    fid<-"boundary"
    
    bf<-boundaryFilter(c('FSC-A','SSC-A'),filterId=fid)
    add(wf,bf,parent=pViewName)
    Rm(paste(fid,"-",sep=""),wf)
    message("done.")
  }
  
}

#################################
## viable
#################################
gating.viable <- function(x, wf, pViewName, plot = FALSE, batch = FALSE,
                          xbin = 128, step = 100, cond = "name", nslaves = 0,
                          ...) {
  
  tView <- "viable"
  nodeNames <- getNodes(wf[[1]])

  if (step >= 2 && !tView %in% nodeNames) {
	 curData <- getData(wf, pViewName)
    message("viable gating...")
    fid <- "viable"

    yChannel <- vMarkers(x)

    viable.filterList <- Gating1D(curData, y = yChannel, filterId = fid,
                                  nslaves = nslaves, method="flowClust",
                                  positive = FALSE)

    nodeID <- add(wf, viable.filterList, parent = pViewName)
    recompute(wf, nodeID)
    
    message("done.")
    if (plot) {
      print(plotGate(wf,nodeID, main = tView, xbin = xbin))
      
      if (!batch) {
        retval <- readline("Enter to resume, c to exit:")
        dev.off()
        if (retval == "c") {
          opt.outer <- options(show.error.messages = FALSE)
          on.exit(options(opt.outer))
          
          stop()
        }
      }
    }
  }
}

gating.singlet.w <- function(x, wf, pViewName, plot = FALSE, batch = FALSE,
                             xbin = 128, step = 100, cond = "name",
                             gating_method = "flowStats", ...) {

  tView <- "singlet"
  nodeNames <-  getNodes(wf[[1]])
  curData <- getData(wf, pViewName)
  if (step >= 3 && !tView %in% nodeNames) {
    message("singlet gating...")
    fid <- "singlet"
    yChannel <- "SSC-W"
    singlet.filterList <- Gating1D(curData, y = yChannel, bwFac = 3, sd = 4,
                                   filterId = fid, peakNr = 2, positive = FALSE,
                                   nslaves = 0, method = gating_method, ...)
    nodeID <- add(wf, singlet.filterList, parent = pViewName)
    recompute(wf, nodeID)

    message("done.")

    if (plot) {
      print(plotGate(wf,nodeID, xbin = xbin, pos = c(0.5,0.8)))

      if (!batch) {
        retval <- readline("Enter to resume, c to exit:")
        dev.off()
        if (retval == "c") {
          opt.outer <- options(show.error.messages = FALSE)
          on.exit(options(opt.outer))
          stop()
        }
      }
    }
  }
}

gating.singlet.h <- function(x, wf, pViewName, plot = FALSE, batch = FALSE,
                             xbin = 128, step = 100, cond = "name",
                             xChannel = "FSC-A", yChannel = "FSC-H",
                             sidescatter = "SSC-A", prediction_level = 0.99,
                             nslaves = 0, type = "multicore", name="singlet",...) {

  chunksize <- 10
  lwf <- length(wf)
  nchunks <- lwf %/% chunksize
  if (nchunks * chunksize < lwf) {
    nchunks <- nchunks + 1
  }
  splitvar <- gl(nchunks, chunksize)
  chunks <- split(seq_len(lwf), splitvar)
  tView <- "singlet"
  nodeNames <- getNodes(wf[[1]])
  
  if(step >= 1 && !any(grepl(tView, nodeNames))) {
    singlet_filter_list_chunks<-vector('list',nchunks)
    for(j in seq_along(chunks)) {
      
      curData <- getData(wf[chunks[[j]]],pViewName)
      message("singlet gating...")
      fid <- name
    
      # Splits the flow set into a list, where each element in the list is a
      # flowSet containg one flow frame.
      # Below, we access manually the individual flow frame with current_frame[[1]].
      flowset_list <- split(curData, sampleNames(curData))

    	# Creates a list of polygon gates based on the prediction bands at the minimum and maximum
    	# xChannel observation using a robust linear model trained by flowStats.
      if (any(grepl("parallel", loadedNamespaces()))&& (is.null(nslaves) || nslaves > 1)) {
        if (is.null(nslaves)) 
          nslaves <- min(length(fslist), parallel::detectCores())
		
        message("Running in parallel mode with ", nslaves, " nodes.")
        if (type == "multicore") {
          singlet.filterList <- mclapply(flowset_list, mc.cores=nslaves, function(current_frame) {
            flowStats:::singletGate(current_frame[[1]], area = xChannel, height = yChannel,
                                    prediction_level = prediction_level, filter_id = fid)$gate
            message(".")
          })
        } else {
          cl <- parallel::makeCluster(nslaves, type = "SOCK")

          singlet.filterList <- parallel::parLapply(cl,flowset_list, function(current_frame) {
						flowStats:::singletGate(current_frame[[1]], area = xChannel, height = yChannel,
                                    prediction_level = prediction_level, filter_id = fid)$gate
					})
          parallel::stopCluster(cl)
        }
      } else {
        singlet.filterList <- lapply(flowset_list, function(current_frame) {
					flowStats:::singletGate(current_frame[[1]], area = xChannel, height = yChannel,
                                  prediction_level = prediction_level, filter_id = fid)$gate
        })	
      }
      singlet_filter_list_chunks[[j]] <- singlet.filterList
    }

    singlet.filterList <- do.call(c, singlet_filter_list_chunks)

    # Adds the list of singlet polygon gates to the workflow.
    singlet.filterResList <- filterList(singlet.filterList)
    nodeID <- add(wf, singlet.filterResList, parent = pViewName)
    recompute(wf, nodeID)
    
    message("done.")
    if (plot) {
      print(plotGate(wf,nodeID, xbin = xbin, pos = c(0.5,0.8)))
      if (!batch) {
        retval <- readline("Enter to resume, c to exit:")
        dev.off()
        if(retval == "c") {
          opt.outer <- options(show.error.messages = FALSE)
          on.exit(options(opt.outer))
          stop()
        }
      }
    }
  }
}

#################################
## lymph +
#################################

gating.lymph <- function(x, wf, pViewName, plot = FALSE, batch = FALSE,
                         xbin = 128, step = 100, nslaves = 0, cond = "name",
                         filterId = "Lymph", isflowClust = TRUE, gate_1D = FALSE,
                         usePrior = "yes", which_gate = c("right", "left", "top", "bottom"), ...) {

  nodeNames <- getNodes(wf[[1]])
  tView <- "Lymph"

  which_gate <- match.arg(which_gate)

  if (step >= 3 && !any(grepl(tView , nodeNames))) {
    message("lymph gating...")
    curData <- getData(wf, pViewName)
    if (gate_1D) {
      yChannel <- "SSC-A"
      positive <- which_gate == "right"

      if (class(x) == "BCell") {
        prior_mean <- c(55328.58, 200000)
        prior_omega <- c(6267.43, 41956.38)^2
      } else {
        prior_mean <- NULL
        prior_omega <- NULL
      }
                      
      lymph.filterList <- Gating1D(curData, y = yChannel, filterId = filterId,
                                   positive = positive, nslaves = nslaves,
								   method="flowClust",
                                   usePrior = usePrior, prior_mean = prior_mean,
                                   prior_omega = prior_omega, ...)
    } else {
      xChannel <- "FSC-A"
      yChannel <- "SSC-A"

      lymph.filterList <- Gating2D(curData, x = xChannel, y = yChannel,
                                   filterId = filterId, nslaves = nslaves,
                                   isflowClust = isflowClust, usePrior = usePrior,
                                   trans = 0, ...)
    }
    nodeID <- add(wf, lymph.filterList, parent = pViewName)
    recompute(wf, tView)

    message("done.")
    if (plot) {
      print(plotGate(wf, nodeID, main = tView, xbin = xbin))
      
      if (!batch) {
        retval <- readline("Enter to resume, c to exit:")
        dev.off()
        if (retval == "c") {
          opt.outer <- options(show.error.messages = FALSE)
          on.exit(options(opt.outer))
          stop()
        }
      }
    }
  }
}

gating.nonDebris <- function(x, wf, pViewName, plot = FALSE, batch = FALSE,
                             xbin = 128, step = 100, nslaves = 0, cond = "name",
                             isflowClust = TRUE, gate_1D = TRUE, usePrior = "yes",
                             which_gate = c("right", "left", "top", "bottom"), ...) {
  
  nodeNames <- getNodes(wf[[1]])
  tView <- "nonDebris"
  curData <- getData(wf, pViewName)
  which_gate <- match.arg(which_gate)
  
  if (step >= 4 && !any(grepl(tView , nodeNames))) {
    message("not Debris gating...")

    if (gate_1D) {
      yChannel <- "FSC-A"
      positive <- which_gate == "right"
        

      if (class(x) == "BCell" || class(x) == "TCell") {
        Mu0 <- c(19776.55, 106196.40, 200000)
        K <- length(Mu0)
        var_x <- (sd(exprs(as(curData, "flowFrame"))[, "FSC-A"]) / K)^2

        Mu0 <- matrix(Mu0, ncol = 1)
        sigma1 <- sqrt(191511992.15)
        sigma2 <- sqrt(81905737.51)
        sigma3 <- (Mu0[3] - Mu0[2]) / 3 - sigma2
        Omega0 <- array(c(sigma1, sigma2, sigma3)^2, c(K, 1, 1))
        Lambda0 <- array(rep(var_x, K), c(K, 1, 1))

        # We assume that the prior probability of mixture component membership
        # is equal across all mixture components.
        w0 <- rep(5, K)

        nu0 <- 4

        prior <- list(Mu0 = Mu0, nu0 = nu0, w0 = w0, Lambda0 = Lambda0, Omega0 = Omega0)
      } else if(class(x) %in% c("ICS")) {
        prior_mean<-c(0, 30196.40, 200000)
        Mu0 <- matrix(prior_mean,ncol=1)
        K <- length(Mu0)
        var_x <- (sd(exprs(as(curData, "flowFrame"))[, "FSC-A"]) / K)^2
		  
        sigma2 <- sqrt(81905737.51)
        sigma3 <- (prior_mean[3] - prior_mean[2]) / 3 - sigma2
        prior_omega <- c(191511992.15/100000, sigma2^2/10, sigma3^2/10)
        Omega0 <- array(prior_omega^2, c(K, 1, 1))
        Lambda0 <- array(rep(var_x, K), c(K, 1, 1))
        w0 <- rep(1, K) / K
		  
        nu0 <- 4
		  
        prior <- list(Mu0 = Mu0, nu0 = nu0, w0 = w0, Lambda0 = Lambda0, Omega0 = Omega0)
      } else {
        prior <- NULL
      }

      nonDebris.filterList <- Gating1D(curData, y = yChannel, filterId = "nonDebris",
                                       positive = positive, nslaves = nslaves, K = 3,
                                       method = "flowClust", usePrior = usePrior,
                                       prior = prior, truncate_min = 0, ...)
    } else {
      xChannel <- "FSC-A"
      yChannel <- "SSC-A"
      nonDebris.filterList <- Gating2D(curData, x = xChannel, y = yChannel,
                                       filterId = "nonDebris", nslaves = nslaves,
                                       isflowClust = isflowClust, usePrior = usePrior,
                                       which_gate = which_gate, ...)
    }
    nodeID <- add(wf, nonDebris.filterList, parent = pViewName)
    recompute(wf, tView)

    message("done.")
    if (plot) {
      print(plotGate(wf, nodeID, main = tView, xbin = xbin))
      
      if (!batch) {
        retval <- readline("Enter to resume, c to exit:")
        dev.off()
        if (retval == "c") {
          opt.outer <- options(show.error.messages = FALSE)
          on.exit(options(opt.outer))
          stop()
        }
      }
    }
  }
}

#################################
## CD3 
#################################
gating.cd3 <- function(x, wf, pViewName, plot = FALSE, batch = FALSE, xbin = 128,
                       step = 100, cond = "name", nslaves = 0, xChannel = "FITC-A",
                       trans = 1, gating_method = "flowClust", usePrior = "yes",
                       prior = NULL, ...) {

  # Check if tcell is already gated
  nodeNames <-  getNodes(wf[[1]])
  
  marker <- x@tMarkers

  if (class(x) == "BCell") {
    positive <- FALSE
  } else if (class(x) == "TCell") {
    positive <- TRUE
  } else if (class(x) == "ICS") {
    positive <- TRUE
  } else {
    stop("unknown gating template type!")
  }
  
  ## when positive==FALSE, CD3+ is actually CD3-
  ## but workflow currently automatically assign
  ## the label which is not consistent with rangeGate definition.
  tView <- marker

  if (step >= 4 && !any(grepl(tView, nodeNames))) {
    message("CD3 gating...")
    curData <- getData(wf, pViewName)
    yChannel <- marker
  
    cd3.filterList <- Gating1D(curData, y = yChannel, filterId = marker,
                               nslaves = nslaves, method = gating_method,
                               positive = positive, trans = trans,
                               usePrior = usePrior, prior = prior, ...)

    cd3.filterList@filterId <- "CD3"
    nodeID <- add(wf, cd3.filterList, parent = pViewName)
    recompute(wf, nodeID)
    

    message("done.")
    if (plot) {
      print(plotGate(wf, nodeID, main = tView, xbin = xbin, pos = c(0.5,0.8)))
      
      if (!batch) {
        retval <- readline("Enter to resume, c to exit:")
        dev.off()
        if (retval == "c") {
          opt.outer <- options(show.error.messages = FALSE)
          on.exit(options(opt.outer))
          stop()
        }
      }
    }
  }
}

#################################
## NK, TH (CD4+, CD8+) 
#################################
quadGate.sequential <- function(fr, markers, nslaves = 0, split = TRUE,
                                gating_method = "flowClust", usePrior = "yes",
                                prior_x = NULL, prior_y = NULL, ...) {
  xParam <- getChannelMarker(fr, markers[1])$name
  yParam <- getChannelMarker(fr, markers[2])$name

  
  ################################      
  #flowClust on xParam
  ################################
  filter1 <- Gating1D(flowSet(fr), y = xParam, nslaves = nslaves,
                      method = gating_method, usePrior = usePrior,
                      filterId = markers[1], prior = prior_x, ...)
  cut.x <- filter1[[1]]@min
  if (split) {
	  ################################################################
	  ## split the data for the further gating
	  ################################################################
	  newFrList <- flowCore::split(fr, filter1[[1]])
	  
	  negFr <- newFrList[[paste0(markers[1], "-")]]
	  posFr <- newFrList[[paste0(markers[1], "+")]]
	  
	  ################################################################
	  ## gate on - and +
	  ################################################################
    # Find the minimum density between the two highest peaks and set the
    # cutpoint there.
    peaks_negFr <- find_peaks(exprs(negFr)[, yParam], peaks = 2, adjust = 3)
    density_y <- density(exprs(negFr)[, yParam], adjust = 3)

    x <- density_y$x[peaks_negFr[1] <= density_y$x & density_y$x <= peaks_negFr[2]]
    y <- density_y$y[peaks_negFr[1] <= density_y$x & density_y$x <= peaks_negFr[2]]    

    cutpoint <- x[which.min(y)]
    # NOTE: The intended gate should work with Gating1d(..., "density1d"), but
    # this did not work for the CIMR T-Cell panel. The gate was about 5.5, but the
    # largest observation was about 4.3.

    gate_coordinates <- list(c(cutpoint, Inf))
    names(gate_coordinates) <- yParam
    filter2 <- list(rectangleGate(gate_coordinates, filterId = markers[1]))
    names(filter2) <- markers[1]
    filter2 <- filterList(filter2)

	  cut.y.l <- filter2[[1]]@min

    prior_y$Mu0 <- matrix(c(find_peaks(exprs(posFr)[, yParam], 1), 2, 4), nrow = 3)    
    prior_y$nu0 <- rep(30, 3)
    prior_y$w0 <- rep(50, 3)
    prior_y$Omega0 <- array(c(0.05, 0.05, 0.05)^2, dim = c(3, 1, 1))
    prior_y$Lambda0 <- array(c(0.5, 0.5, 0.5)^2, dim = c(3, 1, 1))    

	  filter3 <- Gating1D(flowSet(posFr), y = yParam, nslaves = nslaves,
                        method = gating_method, usePrior = usePrior,
                        filterId = markers[2], prior = prior_y, neg_cluster = 2,
                        K = 3, nu.est = 2, ...)

	  cut.y.r <- filter3[[1]]@min
  } else {
    filter2 <- Gating1D(flowSet(fr), y = yParam, method = "density1d", nslaves = nslaves, ...)
    cut.y.l <- cut.y.r <- filter2[[1]]@min
  }
  
  ###############################################################     
  #construct rectangleGates based on the cuts and popNames,clock-wise
  ###############################################################
  gateList <- new("filters")
  
  coord <- list(c(-Inf, cut.x), c(cut.y.l, Inf))
  names(coord) <- as.character(c(xParam, yParam))
  gateList[[paste(paste0(markers, c("-", "+")), collapse="")]] <- rectangleGate(coord)
  
  coord <- list(c(cut.x, Inf), c(cut.y.r, Inf))
  names(coord) <- as.character(c(xParam, yParam))
  gateList[[paste(paste0(markers, c("+", "+")), collapse="")]] <- rectangleGate(coord)
  
  coord <- list(c(cut.x, Inf), c(-Inf, cut.y.r))
  names(coord) <- as.character(c(xParam, yParam))
  gateList[[paste(paste0(markers, c("+", "-")), collapse="")]] <- rectangleGate(coord)
  
  coord <- list(c(-Inf, cut.x), c(-Inf, cut.y.l))
  names(coord) <- as.character(c(xParam, yParam))
  gateList[[paste(paste0(markers, c("-", "-")), collapse="")]] <- rectangleGate(coord)

  gateList
}

#' Applies flowClust to 1 feature to determine a cutpoint between the minimum
#' cluster and all other clusters.
#'
#' We cluster the observations in \code{fr} into \code{K} clusters.
#'
#' By default, the cutpoint is chosen to be the boundary of the first two
#' clusters. That is, between the first two cluster centroids, we find the
#' midpoint between the largest observation from the first cluster and the
#' smallest observations from the second cluster. Alternatively, if the
#' \code{cutpoint_method} is \code{min_density}, then the cutpoint is the point
#' at which the density between the first and second smallest cluster centroids
#' is minimum.
#'
#' @param fr a \code{flowFrame} object
#' @param params TODO
#' @param filterId TODO
#' @param K the number of clusters to find
#' @param usePrior Should we use the Bayesian version of \code{flowClust}?
#' Answers are "yes", "no", or "vague". The answer is passed along to
#' \code{flowClust}.
#' @param prior_list list of prior parameters for the Bayesian \code{flowClust}.
#' If \code{usePrior} is set to 'no', then the list is unused.
#' @param prior_mean a vector of length \code{K} that provides the prior means
#' for each peak. By default, this is \code{NULL}, in which case 
#' @param trans numeric indicating whether the Box-Cox transformation parameter
#' is estimated from the data. May take 0 (no estimation), 1 (estimation) or 2
#' (cluster-speciﬁc estimation). NOTE: For the Bayesian version of
#' \code{flowClust}, this value cannot be 2.
#' @param positive TODO
#' @param cutpoint_method How should the cutpoint be chosen from the fitted
#' \code{flowClust} model? See Details.
#' @param neg_cluster integer. The index of the negative cluster. The cutpoint
#' is computed between clusters \code{neg_cluster} and \code{neg_cluster + 1}.
#' @param truncate_min Truncate observations less than this minimum value. By
#' default, this value is \code{NULL} and is ignored.
#' @param truncate_max Truncate observations greater than this maximum value. By
#' default, this value is \code{NULL} and is ignored.
#' @param percentile the percentile for which we will find the cutpoint using
#' the quantile \code{cutpoint_method}. If the \code{cutpoint_method} is not set
#' to \code{quantile}, this argument is ignored.
#' @param quantile_interval a vector of length 2 containing the end-points of the
#' interval of values to find the quantile cutpoint. If the
#' \code{cutpoint_method} is not setto \code{quantile}, this argument is ignored.
#' @param ... additional arguments that are passed to \code{flowClust}
#' @return a \code{rectangleGate} object ranging consisting of all values greater
#' than the cutpoint determined
flowClust.1d <- function(fr, params, filterId = "", K = 2,  adjust = 1, trans = 0,
                         positive = TRUE, plot = FALSE, usePrior = 'no', prior = NULL,
                         cutpoint_method = c("boundary", "min_density", "quantile", "posterior_mean"),
                         neg_cluster = 1, truncate_min = NULL, truncate_max = NULL,
                         percentile = 0.99, quantile_interval = c(0, 10), ...) {

  cutpoint_method <- match.arg(cutpoint_method)

  if (neg_cluster + 1 > K) {
    stop("The value for K specified is larger than the index of the positive cluster.")
  }

  # If appropriate, we generate prior parameters for the Bayesian version of flowClust.
  if (usePrior == "yes" && is.null(prior)) {
    prior <- prior_flowClust1d(fr = fr, channel = params[1], K = K, adjust = adjust)
  }

  # HACK: Circumvents a bug in flowClust.
  if (missing(prior) || is.null(prior)) {
    prior <- list(NA)
  }

  # If a truncation value is specified, we remove all observations less than
  # (greater than) this value for the marker specified to construct the gate.
  # NOTE: These observations are removed from the 'flowFrame' locally and are gated
  # out only for the determining the gate.
  if (!is.null(truncate_min) || !is.null(truncate_max)) {
    fr <- truncate_flowframe(flow_frame = fr, channel = params[1],
                             min = truncate_min, max = truncate_max)
  }

  # Applies `flowClust` to the feature specified in the `params` argument using
  # the data given in `fr`. We use priors with hyperparameters given by the
  # elements in the list `prior_list`.
  filter1 <- tmixFilter(filterId, params[1], K = K, trans = trans,
                        usePrior = usePrior, prior = prior, ...)
  tmixRes1 <- filter(fr, filter1)

  # To determine the cutpoint, we first sort the centroids so that we can determine
  # the second largest centroid.
  centroids_sorted <- sort(getEstimates(tmixRes1)$locations)

  # Also, because the cluster labels are arbitrary, we determine the cluster
  # the cluster labels, sorted by the ordering of the cluster means.
  labels_sorted <- order(getEstimates(tmixRes1)$locations)

  # Grabs the data matrix that is being gated.
  x <- exprs(fr)[, params[1]]

  # Determines the cutpoint between clusters 1 and 2.
  if (cutpoint_method == "boundary") {
    # Choose the cutpoint as the boundary between the first two clusters.
    # First, we sort the data.
    order_x <- order(x)
    x_sorted <- x[order_x]
    labels <- Map(tmixRes1, rm.outliers = FALSE)[order_x]

    # Determine which observations are between the first two centroids and their
    # corresponding cluster labels.
    which_between <- which(centroids_sorted[neg_cluster] < x_sorted & x_sorted < centroids_sorted[neg_cluster + 1])
    x_between <- x_sorted[which_between]
    labels_between <- labels[which_between]

    # For the observations between the first two centroids, we find the last
    # observation that belongs to the first cluster and the first observation
    # that belongs to the second cluster. In the rare occurrence that no
    # observations from one of the clusters is between the labels, we set the
    # max/min index as NA. This results in the cutpoint being set to the max/min
    # observation on the boundary.
    which_cluster1 <- which(labels_between == labels_sorted[neg_cluster])
    which_cluster2 <- which(labels_between == labels_sorted[neg_cluster + 1])
    max_obs_cluster1 <- ifelse(length(which_cluster1) == 0, NA, max(which_cluster1))
    min_obs_cluster2 <- ifelse(length(which_cluster2) == 0, NA, min(which_cluster2))

    # We define the cutpoint to be the midpoint between the two clusters.
    cutpoint <- mean(x_between[c(max_obs_cluster1, min_obs_cluster2)], na.rm = TRUE)
    if (is.nan(cutpoint)) {
      cutpoint <- centroids_sorted[neg_cluster]
    }
  } else if (cutpoint_method == "min_density") {
    # Determine the minimum density value of the observations between clusters 1 and 2.
    # Sets the cutpoint at the observation that attained this minimum density.
    x_between <- x[centroids_sorted[neg_cluster] < x & x < centroids_sorted[neg_cluster + 1]]
    x_dens <- dmvtmix(x = as.matrix(x_between), object = tmixRes1)
    cutpoint <- x_between[which.min(x_dens)]
  } else if (cutpoint_method == "quantile") {
    cutpoint <- quantile_flowClust(p = percentile, object = tmixRes1, interval = quantile_interval)
  } else { # cutpoint_method == "posterior_mean"
    cutpoint <- centroids_sorted[neg_cluster]
  }

  # After the 1D cutpoint is set, we set the gate coordinates used in the
  # rectangleGate that is returned. If the `positive` argument is set to TRUE,
  # then the gate consists of the entire real line to the right of the cut point.
  # Otherwise, the gate is the entire real line to the left of the cut point.
  if (positive) {
    gate_coordinates <- list(c(cutpoint, Inf))
  } else {
    gate_coordinates <- list(c(-Inf, cutpoint))
  }

  names(gate_coordinates) <- params
  
  fres <- rectangleGate(gate_coordinates, filterId = filterId)
  
  if (plot) {
    gate_pct <- round(100 * mean(x > cutpoint), 3)
    plot_title <- paste0(filterId, " (", gate_pct, "%)")
    plot(fr, tmixRes1, main = plot_title, labels = FALSE)
    abline(v = centroids_sorted, col = rainbow(K))
    abline(v = cutpoint, col = "black", lwd = 3, lty = 2)

    if (!is.null(prior)) {
      x_dens <- seq(min(x), max(x), length = 1000)

      for(k in seq_len(K)) {
        prior_density <- with(prior, dnorm(x_dens, mean = Mu0[k], sd = sqrt(Omega0[k])))

        # Grab posterior estimates for the degrees of freedom (nu) and the
        # transformation parameters. Because these can either be of length 1 or
        # of length K, we grab the appropriate posterior estimates for the
        # current mixture component.
        nu <- ifelse(length(tmixRes1@nu) == 1, tmixRes1@nu, tmixRes1@nu[k])
        lambda <- ifelse(length(tmixRes1@lambda) == 1, tmixRes1@lambda, tmixRes1@lambda[k])        
        posterior_density <- dmvt(x_dens, mu = tmixRes1@mu[k,],
                                  sigma = tmixRes1@sigma[k,,], nu = nu,
                                  lambda = lambda)$value
        lines(x_dens, prior_density, col = rainbow(K)[k], lty = 2, lwd = 1)
        lines(x_dens, posterior_density, col = rainbow(K)[k], lwd = 1)
      }
    }
  }

  fres
}

#' Applies flowClust to two features in a flowFrame to construct an elliptical
#' gate.
#'
#' We cluster the observations in \code{fr} into \code{K} clusters. We set the
#' cutpoint to be the point at which the density between the first and second
#' smallest cluster centroids is minimum.
#'
#' @param fr a \code{flowFrame} object
#' @param xChannel TODO
#' @param yChannel TODO
#' @param filterId TODO
#' @param K the number of clusters to find
#' @param usePrior Should we use the Bayesian version of \code{flowClust}?
#' Answers are "yes", "no", or "vague". The answer is passed along to
#' \code{flowClust}.
#' @param prior_list list of prior parameters for the Bayesian \code{flowClust}.
#' If \code{usePrior} is set to 'no', then the list is unused.
#' @param trans numeric indicating whether the Box-Cox transformation parameter
#' is estimated from the data. May take 0 (no estimation), 1 (estimation) or 2
#' (cluster-speciﬁc estimation). NOTE: For the Bayesian version of
#' \code{flowClust}, this value cannot be 2.
#' @param plot a logical value indicating if the fitted mixture model should be
#' plotted. By default, no.
#' @param truncate_min Truncate observations less than this minimum value. By
#' default, this value is \code{NULL} and is ignored.
#' @param truncate_max Truncate observations greater than this maximum value. By
#' default, this value is \code{NULL} and is ignored.
#' @param ... additional arguments that are passed to \code{flowClust}
#' @return a \code{polygonGate} object containing the contour (ellipse) for 2D
#' gating.
flowClust.2d <- function(fr, xChannel, yChannel, filterId = "", K = 2,
                         usePrior = 'no', prior_list = list(NA), trans = 0,
                         plot = FALSE, which_gate = c("bottom", "top", "left", "right"),
                         gate_type = c("ellipse", "axis"), quantile = 0.995,
                         truncate_min = NULL, truncate_max = NULL, ...) {

  which_gate <- match.arg(which_gate)
  gate_type <- match.arg(gate_type)

  # If appropriate, we generate prior parameters for the Bayesian version of flowClust.
  if (usePrior == "yes" && all.equal(prior_list, list(NA))) {
    prior_list <- prior_flowClust2d(fr = fr, xChannel = xChannel, yChannel = yChannel, K = K)
  }

  # If a truncation value is specified, we remove all observations less than this
  # value for the marker specified to construct the gate.
  # NOTE: These observations are removed from the 'flowFrame' locally and are gated
  # out only for the determining the gate.
  if (!is.null(truncate_min)) {
    exprs(fr) <- exprs(fr)[exprs(fr)[, xChannel] >= truncate_min, ]
    exprs(fr) <- exprs(fr)[exprs(fr)[, yChannel] >= truncate_min, ]    
  }

  x <- exprs(fr)[, xChannel]
  y <- exprs(fr)[, yChannel]

  # Applies `flowClust` to the feature specified in the `params` argument using
  # the data given in `fr`. We use priors with hyperparameters given by the
  # elements in the list `prior_list`.
  filter1 <- tmixFilter(filterId, c(xChannel, yChannel), K = K, trans = trans,
                        usePrior = usePrior, prior = prior_list, ...)
  tmixRes1 <- filter(fr, filter1)

  # Converts the tmixFilterResult object to a polygonGate.
  # We select the cluster with the minimum 'yChannel' to be the subpopulation from
  # which we obtain the contour (ellipse) to generate the polygon gate.
  fitted_means <- getEstimates(tmixRes1)$locations
  cluster_selected <- switch(which_gate,
                             bottom = which.min(fitted_means[, 2]),
                             top = which.max(fitted_means[, 2]),
                             left = which.min(fitted_means[, 1]),
                             right = which.max(fitted_means[, 1]))    

  if (gate_type == "ellipse") {
    contour_ellipse <- .getEllipse(filter = tmixRes1,
                                   include = cluster_selected)
    flowClust_gate <- polygonGate(.gate = matrix(contour_ellipse, ncol = 2,
                                    dimnames = list(NULL, tmixRes1@varNames)),
                                  filterId = filterId)
  } else if (gate_type == "axis") {
    chisq_quantile <- qchisq(quantile, df = 2)

    xbar <- tmixRes1@mu[cluster_selected, ]
    Sigma <- tmixRes1@sigma[cluster_selected, , ]

    Sigma_eigen <- eigen(Sigma, symmetric = TRUE)
    u1 <- Sigma_eigen$vectors[, 1]
    u2 <- Sigma_eigen$vectors[, 2]
    lambda1 <- Sigma_eigen$values[1]
    lambda2 <- Sigma_eigen$values[2]

    # Determines which eigenvector points towards the first quadrant (has
    # positive slope). Because both eigenvectors can potentially point in the
    # negative direction, we also check to see the negated eigenvectors point
    # towards the first quadrant.
    if (all(u1 >= 0)) {
      axis <- sqrt(lambda1 * chisq_quantile) * u1
      axis_perp <- sqrt(lambda2 * chisq_quantile) * u2
    } else if (all(-u1 >= 0)) {
      axis <- -sqrt(lambda1 * chisq_quantile) * u1
      axis_perp <- sqrt(lambda2 * chisq_quantile) * u2
    } else if (all(u2 >= 0)) {
      axis <- sqrt(lambda2 * chisq_quantile) * u2
      axis_perp <- sqrt(lambda1 * chisq_quantile) * u1
    } else if (all(-u2 >= 0)) {
      axis <- -sqrt(lambda2 * chisq_quantile) * u2
      axis_perp <- sqrt(lambda1 * chisq_quantile) * u1
    }

    # The gate location is the frame of reference for the gate. If it is xbar,
    # then the frame of reference is the cross-section along the eigenvector
    # from top-left to bottom-right. We translate this reference as a function
    # of the appropriate chi-squared coefficient.
    gate_location <- xbar + 0.25 * axis

    # To construct the gate, we have effectively shifted the eigenvector with the
    # negative slope. We then extend the gate both horizontally and vertically
    # to the maximum observed values in the horizontal and vertical directions.
    # NOTE: We extend the gate one standard deviation beyond the maximum values
    # observed to mimic a gate that extends without limit in the positive
    # directions. However, because flowCore cannot handle such a shape, we force
    # the gate to have the same shape for the observed data.
    x_max <- max(x) + sd(x)
    y_max <- max(y) + sd(y)
    gate_x <- c(x_max, (gate_location - axis_perp)[2])
    gate_y <- c((gate_location + axis_perp)[1], y_max)

    # We extend the gate to the min and max values
    polygon_gate <- rbind(gate_location + axis_perp,
                          gate_y,
                          c(x_max, y_max),
                          gate_x,
                          gate_location - axis_perp)

    colnames(polygon_gate) <- c(xChannel, yChannel)
    flowClust_gate <- polygonGate(filterId = filterId, boundaries = polygon_gate)
  }

  if (plot) {
    plot(fr, tmixRes1, main = filterId)

    if (gate_type == "axis") {
      # The major and minor axes (eigenvectors) scaled by their respective
      # eigenvalues and the chi-squared quantile.
      lines(rbind(xbar, xbar + sqrt(lambda1 * chisq_quantile) * u1), col = "darkgreen")
      lines(rbind(xbar, xbar - sqrt(lambda1 * chisq_quantile) * u1), col = "darkgreen")
      lines(rbind(xbar, xbar + sqrt(lambda2 * chisq_quantile) * u2), col = "darkgreen")
      lines(rbind(xbar, xbar - sqrt(lambda2 * chisq_quantile) * u2), col = "darkgreen")

      # Draws the polygon gate.
      lines(polygon_gate, col = "red")
      lines(rbind(gate_location - axis_perp, gate_location + axis_perp), col = "red")

      # Also, draws points at the vertices of the polygon gate.
      points(polygon_gate, col = "red", pch = 16)
    }
  }

  flowClust_gate
}

flowClust.act <- function(fr, params, filterId, nslaves = 0,
                          gating_method = "flowClust", usePrior = "yes",
                          prior_x = NULL, prior_y = NULL, ...) {

  # flowClust on xParam
  gate_xParam <- Gating1D(flowSet(fr), y = params[1], positive = TRUE, nslaves = nslaves,
                          method = gating_method, usePrior = usePrior,
                          filterId = channels2markers(fr, params[1]),
                          prior = prior_x, ...)

  # Subsets the data for further gating on the positive population.
  newFr <- Subset(flowSet(fr), gate_xParam)
  
  # flowClust on yParam
  gate_yParam <- Gating1D(newFr, y = params[2], nslaves = nslaves,
                          method = gating_method, usePrior = usePrior,
                          filterId = channels2markers(fr, params[2]),
                          prior = prior_y, K = 3, neg_cluster = 2, ...)

  # Extracts the gate information from each gate.
  gate_xParam <- c(gate_xParam[[1]]@min, gate_xParam[[1]]@max)
  gate_yParam <- c(gate_yParam[[1]]@min, gate_yParam[[1]]@max)

  # Constructs a rectangleGate based on the flowClust cuts.
  activated_gate <- list(gate_xParam, gate_yParam)
  names(activated_gate) <- params
  rectangleGate(activated_gate, filterId = filterId)
}

#' Finds a specified number of peaks in the given vector after smoothing the data
#' with a kernel density estimator.
#'
#' First, we smooth the data using kernel density estimation with the
#' \code{density} function. Then, we find every local maxima such that the
#' density is concave (downward). We choose the \code{peaks} largest local maxima
#' as the peaks.
#'
#' @param x numeric vector
#' @param peaks the number of peaks to find
#' @param ... additional arguments passed to the \code{density} function.
#' @return the values where the peaks are attained.
#' @examples
#' library(flowClust)
#' set.seed(42)
#' # 2 peaks with a minor hump
#' y <- SimulateMixture(10000, c(.5, .3, .2), c(2, 5, 7), c(1, 1, 1), nu = 10)
#' plot(density(y))
#' peaks <- find_peaks(y, peaks = 2)
#' abline(v = peaks, col = "red")
find_peaks <- function(x, peaks = 2, ...) {
  x <- as.vector(x)

  if (peaks > length(x)) {
    warning("There are not enough observations to find the specified number of peaks.")
    largest_peaks <- c(x, rep(NA, peaks - length(x)))
  } else {
    dens <- density(x, ...)

    # The 'embed' function with dim = 3 generated an N x 3 matrix,
    # where the columns correspond to observations x[i-1], x[i], and x[i+1],
    # which allows us to look at the before and after density value to see if a
    # local maxima is achieved.
    local_maxima <- apply(embed(dens$y, dim = 3), 1, function(z) {
      z[2] > z[1] & z[2] > z[3]
    })
    local_maxima <- which(local_maxima) + 1
    peaks_order <- order(dens$y[local_maxima], decreasing = TRUE)[seq_len(peaks)]
    largest_peaks <- local_maxima[peaks_order]
    largest_peaks <- dens$x[largest_peaks]
  }
  sort(largest_peaks, na.last = TRUE)
}

quantileGate <- function(fr, probs, stain, plot = FALSE, positive = TRUE,
                         filterId = "", ...) {
  x <- exprs(fr[, stain])
  cutpoint <- quantile(x, probs = probs)

  if (positive) {
    gate_coordinates <- list(c(cutpoint, Inf))
  } else {
    gate_coordinates <- list(c(-Inf, cutpoint))
  }
  names(gate_coordinates) <- stain
	
  if (plot) {
		hist(x)
		plot(density(x))
		abline(v = cutpoint, col = "red")
		text(y = 0.5, x = cutpoint, labels = paste("quantile:", probs))
	}
	rectangleGate(gate_coordinates, filterId = filterId)
}
