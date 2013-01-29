#################################
## B cell panel gating
#################################
setMethod("gating", signature = c("BCell", "GatingSet"),
          definition = function(x, wf, priors = NULL, ...) {

  lView <- "root"
  gating.nonDebris(x, wf, pViewName = lView, ...)

  lView <- "nonDebris"
  gating.lymph(x, wf, pViewName = lView, xChannel = "FSC-A", yChannel = "SSC-A", which_gate = "left", ...)

  lView <- "Lymph"
  # To select the appropriate singlet gating, we use the following logic:
  # 1. If FSC-H is present in the channels list, apply singlet gating to FSC-A vs FSC-H.
  # 2. If not, then look for SSC-W. If it's present, apply singlet gating to SSC-W vs SSC-A.
  # 3. If not, skip singlet gating so that the parent node of CD3 is Lymph.
  colnames_wf <- colnames(getData(wf[[1]], lView))
  if ("FSC-H" %in% colnames_wf) {
    gating.singlet.h(x, wf, pViewName = lView, ...)
    lView <- "singlet"
  } else if ("SSC-W" %in% colnames_wf) {
    gating.singlet.w(x, wf, pViewName = lView, ...)
    lView <- "singlet"
  }

  gating.cd3(x, wf, pViewName = lView, xChannel = "SSC-A",
             gating_method = "density1d", ...)

  # CD19 1D gating
  lView <- "CD3"
  gating.cd19(x, wf, pViewName = lView, gating_method = "density1d", ...)

  # CD20 1D gating
  gating.cd20(x, wf, pViewName = lView, gating_method = "density1d", ...)

  # Here, we require intersections of the CD19 and CD20 1D-gates to gate further
  # the following two sub-populations:
  # 1. CD19+/CD20+
  # 2. CD19+/CD20-
  #
  # When we plot the workflow gates, the intersected gates are not plotted
  # because the `plotGate()` function skips nodes that are objects of class
  # `intersectFilter`.
  gating.bool_cd19_cd20(wf, pViewName = lView, ...)

  # Gates plasmablasts: CD38 vs CD27
  lView <- "CD19 and not CD20"
  gating.plasmaB(x, wf, pViewName = lView, collapse = TRUE,
                 prior_CD38 = priors[["CD38"]], prior_CD27 = priors[["CD27"]],
                 truncate_min = -1, ...)

  # Gates B-cell subpopulations: IgD vs CD27 with a quad gate
  lView <- "CD19 and CD20"
  priors[["CD27"]]$nu0 <- rep(30, length(priors[["CD27"]]$nu0))
  priors[["IgD"]]$nu0 <- rep(30, length(priors[["IgD"]]$nu0))  
  gating.bsub(x, wf, pViewName = lView, collapse = TRUE,
              prior_IgD = priors[["IgD"]], prior_CD27 = priors[["CD27"]], ...)

  # Gates CD38 vs CD24
  gating.trans(x, wf, pViewName = lView, prior_CD38 = priors[["CD38"]],
               prior_CD24 = priors[["CD24"]], truncate_min = -1, ...)

  message("done.")	
})

# Here, we require intersections of the CD19 and CD20 1D-gates to gate further
# the following two sub-populations:
# 1. CD19+/CD20+
# 2. CD19+/CD20-
gating.bool_cd19_cd20 <- function(wf, pViewName, step = 100, ...) {
	
  nodeNames <-  getNodes(wf[[1]])
  tNodes <- "CD19 and CD20" 
	
	# Constructs intersections of the CD19 and CD20 gates.
  if (step >= 7 && !any(tNodes %in% nodeNames)) {
    message(tNodes, " gating...")

    cd19pos_cd20pos <- booleanFilter(CD19 & CD20, filterId = tNodes)
    nodeID <- add(wf, cd19pos_cd20pos, parent = pViewName)
		recompute(wf,nodeID)

    message("done.")
  }
	
	tNodes <-"CD19 and not CD20"

  # Creates a node for the CD19+/CD20- intersected gate.
  # Then, we remove the superfluous node that is added here
  if (step >= 7 && !any(tNodes %in% nodeNames)) {
    message(tNodes," gating...")
    cd19pos_cd20neg <- booleanFilter(CD19 & !CD20, filterId = tNodes)
    nodeID <- add(wf, cd19pos_cd20neg, parent = pViewName)
    recompute(wf, nodeID)
    message("done.")
  }

}

############################
#quadGate on CD19 x CD20
#############################
gating.bcell <- function(x, wf, pViewName, plot = FALSE, batch = TRUE,
                         xbin = 128, step = 100, cond = "name", ...) {

  nodeNames <- getNodes(wf[[1]])
  curData <- getData(wf, pViewName)

  markers <- bMarkers(x)
  markers.desc <- unlist(lapply(markers, function(y) {
    getChannelMarker(frm = curData[[1]], y)$desc
  }))
  bsubView <- paste0(markers.desc[1], "\\+", markers.desc[2], "\\-")
  vBsub <- paste0(markers.desc[1], "+", markers.desc[2], "-")

  xChannel <- markers[1]
  yChannel <- markers[2]

  if (step >= 6 && !vBsub %in% nodeNames) {
    message("B cells gating...")
    
    bcell.filterList <- Gating1D(curData, y = yChannel, x = xChannel,
                                 filterId = "B Cells", nslaves = 0, ...)
    nodeID <- add(wf, bcell.filterList, parent = pViewName)
    recompute(wf, nodeID)
    message("done.")	

    if (plot) {
      print(plotGate(x = curData, y = gate(action(wf[[vBsub]])), smooth = FALSE,
                     margin = FALSE, xbin = xbin, cond = cond,
                     main = "CD19 vs CD20"))
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

############################
# 1D-gate on CD19
#############################
gating.cd19 <- function(x, wf, pViewName, plot = FALSE, batch = TRUE, xbin = 128,
                        step = 100, cond = "name", gating_method = "flowClust",
                        usePrior = "yes", prior = NULL, ...) {
  
  nodeNames <-  getNodes(wf[[1]])
  curData <- getData(wf, pViewName)

  xChannel<- "CD19"
  yChannel<- "SSC-A"
  vBsub <- "CD19"

  if (step >= 6 && !vBsub %in% nodeNames) {
    message("CD19 gating...")
    cd19.filterList <- Gating1D(curData, y = xChannel, filterId = "CD19",
                                method = gating_method, nslaves = 0,
                                usePrior = usePrior, prior = prior, ...)
    nodeID <- add(wf, cd19.filterList, parent = pViewName)
    recompute(wf, nodeID)

   message("done.")	

   if (plot) {
     print(plotGate(wf, nodeID, xbin = xbin, pos = c(0.5, 0.8)))
      
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

############################
# 1D-gate on CD20
#############################
gating.cd20 <- function(x, wf, pViewName, plot = FALSE, batch = TRUE, xbin = 128,
                        step = 100, cond = "name", gating_method = "flowClust",
                        usePrior = "yes", prior = NULL, ...) {

  nodeNames <- getNodes(wf[[1]])
  curData <- getData(wf, pViewName)

  xChannel<- "CD20"
  yChannel<- "SSC-A"
  vBsub <- "CD20"

  if (step >= 6 && !vBsub %in% nodeNames) {
    message("CD20 gating...")
    cd20.filterList <- Gating1D(curData, y = xChannel, filterId = "CD20",
                                method = gating_method, nslaves = 0,
                                usePrior = "yes", prior = prior, ...)

    nodeID <- add(wf, cd20.filterList, parent = pViewName)
    recompute(wf, nodeID)

    message("done.")
    if (plot) {
      print(plotGate(wf, nodeID, xbin = xbin, pos = c(0.5, 0.8)))
		
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

# CD38 vs CD27
gating.plasmaB <- function(x, wf, pViewName, plot = FALSE, batch = FALSE,
                           xbin = 128, step = 100, cond = "name", collapse = FALSE,
                           gating_method = "flowClust", usePrior = "yes",
                           prior_CD38 = NULL, prior_CD27 = NULL, ...) {

  nodeNames <-  getNodes(wf[[1]])
  curData <- getData(wf,pViewName)
  markers <- pMarkers(x)
  markers.desc <- unlist(lapply(markers, function(y) {
    getChannelMarker(frm = curData[[1]], y)$desc
  }))
  bsubView <- paste0(markers.desc[1], "\\+", markers.desc[2], "\\-")
  vP <- "plasmablast+"
	
  if (step >= 8 && !any(grepl(vP, nodeNames))) {	
    message("plasma blast cell gating...", pViewName)

    xChannel <- markers[1]
    yChannel <- markers[2]

    ## 1-d gate on cd27 first
    cd27.filterList <- Gating1D(curData, y = yChannel, method = gating_method,
                                filterId = "cd27", nslaves = 0, usePrior = usePrior,
                                prior = prior_CD27, collapse = collapse, ...)

    # Then gate CD38 on the cd27+
    cd27Pos <- flowSet(lapply(sampleNames(curData), function(curSample) {
      split(curData[[curSample]], cd27.filterList[[curSample]])[[1]]
    }))

    sampleNames(cd27Pos) <- sampleNames(curData)

    # We take the prior from the positive population to be the largest peak after
    # applying a lot of smoothing.
    xMarker <- markers2channels(flow_frame = cd27Pos[[1]], markers = xChannel)
    CD38_collapsed <- exprs(as(cd27Pos, "flowFrame"))[, xMarker]
    CD38_positive_peak <- find_peaks(CD38_collapsed, peaks = 1, adjust = 3)

    # We set the prior mean of the positive peak to be the maximum density point,
    # and we set the prior mean of the negative peak to be the same point shifted
    # by 2. Although this is a hack, this reflects our prior beliefs about where
    # the negative population should be. Also, we use the variance of the second
    # prior mean for the first prior mean.
    prior_CD38$Mu0[1,] <- CD38_positive_peak - 2
    prior_CD38$Mu0[2,] <- CD38_positive_peak
    prior_CD38$Omega0[1,,] <- prior_CD38$Omega0[2,,]

    cd38.filterList <- Gating1D(cd27Pos, y = xChannel, method = gating_method,
                                filterId = "cd38", nslaves = 0, usePrior = usePrior,
                                prior = prior_CD38, collapse = collapse, ...)

		# Constructs rectangle gate afterwards.
    bsub.filterList <- sapply(names(cd38.filterList), function(curSample) {
      xfilter <- cd38.filterList[[curSample]]
      yfilter <- cd27.filterList[[curSample]]

      coord <- list(c(xfilter@min, xfilter@max), c(yfilter@min, yfilter@max))					
      names(coord) <- c(parameters(xfilter), parameters(yfilter))
      rectangleGate(coord, filterId = "plasmablast")
    })

    nodeID <- add(wf, filterList(bsub.filterList), parent = pViewName)
    recompute(wf,nodeID)
		
		message("done.")
		
		if (plot) {
			print(plotGate(wf, nodeID, xbin = xbin, pos = c(0.5,0.8)))

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

##################################################
## transitional, memory, naive
## Gates CD38 vs CD24
##################################################
gating.trans <- function(x, wf, pViewName, plot = FALSE, batch = FALSE,
                         xbin = 128, step = 100, cond = "name", collapse = FALSE,
                         gating_method = "flowClust", usePrior = "yes",
                         prior_CD38 = NULL, prior_CD24 = NULL,
                         gate_CD38_pos = FALSE, ...) {

  nodeNames <- getNodes(wf[[1]])
  curData <- getData(wf, pViewName)
  markers <- c("CD38", "CD24")
  markers.desc <- unlist(lapply(markers, function(y) {
    getChannelMarker(frm = curData[[1]], y)$desc
  }))
  tView <- paste0(markers.desc[1], "\\+", markers.desc[2], "\\-")
  vBsub <- paste0(markers.desc[1], "+", markers.desc[2], "-")

  xChannel <- markers[1]
  yChannel <- markers[2]

  if(step >= 8 && !any(grepl(tView, nodeNames))) {	
    message("transitional gating...", pViewName)

    gate_CD38_pos <- TRUE

    # If specified, we gate the CD38 marker with a 1D gate first and retain only
    # the CD38+ cells.
    if (gate_CD38_pos) {
      cd38.filterList <- Gating1D(curData, y = xChannel, filterId = "CD38",
                                  method = gating_method, nslaves = 0,
                                  usePrior = usePrior, prior = prior_CD38,
                                  collapse = collapse, nu = 30, ...)
      
      cd38_pos <- flowSet(lapply(sampleNames(curData), function(curSample) {
        split(curData[[curSample]], cd38.filterList[[curSample]])[[1]]
      }))
      sampleNames(cd38_pos) <- sampleNames(curData)
      xParam <- getChannelMarker(cd38_pos[[1]], xChannel)$name
      curData <- cd38_pos
    }

    # We truncate all CD38 values above 4 because these appear to be boundary
    # events.
    truncate_channel <- markers2channels(curData[[1]], xChannel)
    curData <- fsApply(curData, truncate_flowframe, channel = truncate_channel,
                       max = 4)

    prior <- list(xChannel = prior_CD24, yChannel = prior_CD38)

    # Applies 2D flowClust to CD38 vs CD24 markers. We select the cluster on the
    # right and use the 'axis' gate.
    trans.filterList <- Gating2D(curData, x = xChannel, y = yChannel,
                                   filterId = "transitional", nslaves = 0,
                                   isflowClust = TRUE, usePrior = usePrior,
                                   K = 4, gate_type = "axis", which_gate = "right", ...)
								
    nodeID <- add(wf, trans.filterList, parent = pViewName)
    recompute(wf, nodeID)
    message("done.")

    if (plot) {
      print(plotGate(wf, nodeID, xbin = xbin, pos = c(0.5, 0.8)))
    }
  }
}

##################################################
## B subpopulations
## Gates IgD vs CD27
##################################################
gating.bsub <- function(x, wf, pViewName, plot = FALSE, batch = FALSE,
                        xbin = 128, step = 100, cond = "name", collapse = FALSE,
                        gating_method = "flowClust", usePrior = "yes",
                        prior_IgD = NULL, prior_CD27 = NULL, ...) {

	nodeNames <-  getNodes(wf[[1]])
	curData <- getData(wf, pViewName)
  markers <- bsubMarkers(x)
  markers.desc <- unlist(lapply(markers, function(y) {
    getChannelMarker(frm = curData[[1]], y)$desc
  }))
  bsubView <- paste0(markers.desc[1], "\\+", markers.desc[2], "\\-")
  vBsub <- paste0(markers.desc[1], "+", markers.desc[2], "-")

  if(step >= 8 && !any(grepl(bsubView, nodeNames))) {
    message("B Cells sub gating...", pViewName)

    xChannel <- markers[1]
    yChannel <- markers[2]

    ## Apply the CD27 gate first.
    cd27.filterList <- Gating1D(curData, y = yChannel, method = gating_method,
                                nslaves = 0, filterId = "cd27", usePrior = usePrior,
                                prior = prior_CD27, collapse = collapse, nu = 30, truncate_min = -1, ...)

    # Then gate IgD on the cd27+
    cd27_positive <- flowSet(lapply(sampleNames(curData), function(curSample) {
      split(curData[[curSample]], cd27.filterList[[curSample]])[[1]]
    }))

    sampleNames(cd27_positive) <- sampleNames(curData)

    xParam <- getChannelMarker(cd27_positive[[1]], xChannel)$name

    # Here, we collapse all samples. Then, we take the prior of the negative
    # population to be the largest peak after applying a lot of smoothing.
    # Generally, there is only one visible peak anyways.
    IgD_collapsed <- exprs(as(cd27_positive, "flowFrame"))[, xParam]
    IgD_negative_peak <- find_peaks(IgD_collapsed, peaks = 1, adjust = 3)

    # We use three components to model the IgD marker. We use two components for
    # the negative peak and one component for the positive peak.
    huber_s <- huber(IgD_collapsed)$s
    prior_IgD$Mu0 <- matrix(c(IgD_negative_peak - 3.5 * huber_s, IgD_negative_peak, IgD_negative_peak + 2.75 * huber_s), nrow = 3)
    prior_IgD$Omega0 <- array(rep(0.05^2, 3), dim = c(3, 1, 1))
    prior_IgD$Lambda0 <- array(rep(0.5^2, 3), dim = c(3, 1, 1))
    prior_IgD$w0 <- rep(30, 3)
    prior_IgD$nu0 <- rep(30, 3)

    IgD.filterList <- Gating1D(cd27_positive, y = xChannel, method = gating_method,
                               nslaves = 0, filterId = "IgD", usePrior = usePrior,
                               prior = prior_IgD, collapse = collapse, nu.est = 2,
                               K = 3, neg_cluster = 2, truncate_min = 0, truncate_max = 2.5, ...)

    
    # construct quadGate afterwards
    bsub.filterList <- sapply(names(IgD.filterList), function(curSample) {
      xfilter <- IgD.filterList[[curSample]]
      
      # change from neg to pos for quadGate only recognize (xx,Inf)
      yfilter <- cd27.filterList[[curSample]]
      coord <- list(unname(xfilter@min), unname(yfilter@min))
      names(coord) <- c(parameters(xfilter), parameters(yfilter))
      quadGate(coord)
    })

    nodeID <- add(wf, filterList(bsub.filterList), parent = pViewName)
    recompute(wf, nodeID)
    message("done.")

    if (plot) {
      print(plotGate(wf, nodeID, xbin = xbin, pos = c(0.5,0.8)))
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
