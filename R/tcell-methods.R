#################################
## T cell panel gating
#################################
setMethod("gating", signature = c("TCell", "GatingSet"),
          definition = function(x, wf, priors = NULL, ...) {
            
  lView <- "root"
  gating.nonDebris(x, wf, pViewName = lView, ...)

  lView <- "nonDebris"    
  gating.lymph(x, wf, pViewName=lView,...)
  
  lView <- "Lymph"
  gating.cd3(x, wf, pViewName = lView, xChannel = "SSC-A", gating_method = "density1d", ...)

  lView <- "CD3"
  priors[["CD4"]]$nu0 <- rep(30, length(priors[["CD4"]]$nu0))
  priors[["CD8"]]$nu0 <- rep(30, length(priors[["CD8"]]$nu0))
  gating.tsub(x, wf, pViewName = lView, prior_CD4 = priors[["CD4"]],
              prior_CD8 = priors[["CD8"]],...)
  
  tView.cd4 <- "CD4"
  tView.cd8 <- "CD8"
  priors[["CD38"]]$nu0 <- rep(30, length(priors[["CD38"]]$nu0))
  priors[["HLA-DR"]]$nu0 <- rep(30, length(priors[["HLA-DR"]]$nu0))
  gating.act.cd4(x, wf, pViewName = tView.cd4, prior_CD38 = priors[["CD38"]],
                           prior_HLA_DR = priors[["HLA-DR"]],...)
  gating.act.cd8(x, wf, pViewName = tView.cd8, prior_CD38 = priors[["CD38"]],
                           prior_HLA_DR = priors[["HLA-DR"]],...)

  priors[["CCR7"]]$nu0 <- rep(30, length(priors[["CCR7"]]$nu0))
  priors[["CD45RA"]]$nu0 <- rep(30, length(priors[["CD45RA"]]$nu0))
  gating.mem_CD4(x, wf, pViewName = tView.cd4, prior_CCR7 = priors[["CCR7"]],
                 prior_CD45RA = priors[["CD45RA"]],, ...)
  gating.mem_CD8(x, wf, pViewName = tView.cd8, prior_CCR7 = priors[["CCR7"]],
                 prior_CD45RA = priors[["CD45RA"]],, ...)

  message("finished.")
})

gating.tsub <- function(x, wf, pViewName, split = TRUE, plot = FALSE,
                        batch = TRUE, xbin = 128, step = 100, prior_CD4 = NULL,
                        prior_CD8 = NULL, nslaves = 0, ...) {
  
  nodeNames <- getNodes(wf[[1]])

  tView.cd4 <- "CD4"
  tView.cd8 <- "CD8"
  tViews <- c("CD4", "CD8")

  if (step >= 5 && !any(unlist(lapply(tViews, grepl, nodeNames)))) {
    curData <- getData(wf, pViewName)
    message("T-sub gating...")

    markers <- tsubMarkers(x)
    channels <- markers2channels(getData(wf[[1]]), markers)
    
    xChannel <- markers[1]
    yChannel <- markers[2]

    ############################
    # sequential flowClust
    #############################
    if (class(x) == "HVTN080") {
      prior_list <- list(xChannel = prior_CD4, yChannel = prior_CD8)

      cd4cd8.filterList <- Gating1D(curData, x = xChannel, y = yChannel,
                                    trans = 0, usePrior = "yes", prior = prior_list,
                                    neg_cluster = 2, nslaves = nslaves, ...)
      
      cd4cd8.filterList <- lapply(cd4cd8.filterList, quadGate2rectangleGates,
                                  markers = markers, channels = channels, quadrants = c(1, 3))
    } else {
      cd4cd8.filterList <- fsApply(curData, quadGate.sequential, markers = markers,
                                   split = split, trans = 0, prior_x = prior_CD4,
                                   prior_y = prior_CD8, nu = 30, ...)
      ## only keep CD4+, CD8+
      cd4cd8.filterList <- lapply(cd4cd8.filterList, function(curFilters) {
        curFilters[["CD4+CD8+"]] <- NULL
        curFilters[["CD4-CD8-"]] <- NULL
        curFilters
      })
    }

    cd4cd8.fl <- as(cd4cd8.filterList, "filtersList")
    
    ############################
    #fetch 4 filterList from filtersList
    #before adding them to wf since wf doesn't support filtersList yet
    #############################
  
    cd4.list <- lapply(cd4cd8.fl, function(curFilters) {
      curFilters[["CD4+CD8-"]]@filterId <- "CD4"
      curFilters[["CD4+CD8-"]]
    })
    cd4.list <- filterList(cd4.list)
  
    cd8.list <- lapply(cd4cd8.fl, function(curFilters) {
      curFilters[["CD4-CD8+"]]@filterId <- "CD8"
      curFilters[["CD4-CD8+"]]
    })
    cd8.list <- filterList(cd8.list)
  
    #############
    #add to wf
    ##############
    nodeID1 <- add(wf, cd4.list, parent = pViewName)
    recompute(wf, nodeID1)
  
    nodeID2 <- add(wf, cd8.list, parent = pViewName)
    recompute(wf, nodeID2)
    
    message("done.")  
    if (plot) {
      print(plotGate(wf, c(nodeID1, nodeID2), xbin = xbin))
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
## activated (double +)
#################################
gating.act.cd4 <- function(x, wf, pViewName, plot = FALSE, batch = FALSE,
                           xbin = 128, step = 100, prior_CD38 = NULL,
                           prior_HLA_DR = NULL, ...) {

  nodeNames <- getNodes(wf[[1]])
  curData <- getData(wf, pViewName)
  
  tView <- "activated CD4"
  
  # Checks if activated gating already done. If not, proceed with gating.
  if (step >= 7 && !any(grepl(tView, nodeNames))) {
    message("Activate CD4 gating...")
    markers <- actMarkers(x)
    xChannel <- markers[1]
    yChannel <- markers[2]
    
    xParam <- getChannelMarker(curData[[1]], xChannel)$name
    yParam <- getChannelMarker(curData[[1]], yChannel)$name
    
    ##############
    #flowClust  
    ##############
    filterId <- "activated CD4"

    # Pools the samples in the flowSet into a flowFrame from which we obtain the
    # marker-specific priors.
    x <- exprs(as(curData, "flowFrame"))

    # Prior elicitation for CD38
    peaks_CD38 <- find_peaks(x[, xParam], peaks = 2, adjust = 3)

    # If the second peak is smaller than the max peak, we set the first peak
    # equal to the max peak and set the second peak to two standard deviations to
    # the right.
    if (any(is.na(peaks_CD38))) {
      peaks_CD38[2] <- peaks_CD38[1] + 2 * huber(x[, xParam])$s
    }
    
    prior_CD38$Mu0 <- matrix(peaks_CD38, nrow = 2)
    prior_CD38$Omega0[1,,] <- 0.1
    prior_CD38$Omega0[2,,] <- 0.1  
    prior_CD38$Lambda0[1,,] <- 0.5
    prior_CD38$Lambda0[2,,] <- 0.5  
    prior_CD38$w0 <- rep(10, 2)

    # Prior elicitation for HLA-DR
    peak_HLA_DR <- find_peaks(x[, yParam], peaks = 1)
    huber_s <- huber(x[, yParam])$s

    prior_HLA_DR$Mu0 <- matrix(c(peak_HLA_DR - 2.5 * huber_s, peak_HLA_DR, peak_HLA_DR + 2.5 * huber_s), nrow = 3)
    prior_HLA_DR$w0 <- rep(10, 3)
    prior_HLA_DR$nu0 <- rep(30, 3)
    prior_HLA_DR$Omega0 <- array(rep(0.1^2, 3), dim = c(3, 1, 1))
    prior_HLA_DR$Lambda0 <- array(rep(0.5, 3), dim = c(3, 1, 1))

    act.filterList <- fsApply(curData, flowClust.act, params = c(xParam, yParam),
                              filterId = filterId, prior_x = prior_CD38,
                              prior_y = prior_HLA_DR, nu.est = 2, truncate_min = 0, ...)
    act.fl <- filterList(act.filterList)                  
    
    nodeID <- add(wf, act.fl, parent = pViewName)
    recompute(wf, nodeID)
      
    message("done.")
    if (plot) {
      print(plotGate(wf, nodeID, xbin = xbin, main = tView))

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

gating.act.cd8 <- function(x, wf, pViewName, plot = FALSE, batch = FALSE,
                           xbin = 128, step = 100, prior_CD38 = NULL,
                           prior_HLA_DR = NULL, ...) {
  nodeNames <- getNodes(wf[[1]])
  curData <- getData(wf, pViewName)
  tView <- "activated CD8"
    
  if (step >= 7 && !any(grepl(tView, nodeNames))) {
    message("Activate CD8 gating...")
    markers <- actMarkers(x)
    xChannel <- markers[1]
    yChannel <- markers[2]
      
    xParam <- getChannelMarker(curData[[1]], xChannel)$name
    yParam <- getChannelMarker(curData[[1]], yChannel)$name
    
    ##############
    #flowClust  
    ##############
    filterId <- "activated CD8"

    # Pools the samples in the flowSet into a flowFrame from which we obtain the
    # marker-specific priors.
    x <- exprs(as(curData, "flowFrame"))

    # Prior elicitation for HLA-DR
    peak_HLA_DR <- find_peaks(x[, yParam], peaks = 1)
    huber_s <- huber(x[, yParam])$s

    prior_HLA_DR$Mu0 <- matrix(c(peak_HLA_DR - 2.5 * huber_s, peak_HLA_DR, peak_HLA_DR + 2.5 * huber_s), nrow = 3)
    prior_HLA_DR$w0 <- rep(10, 3)
    prior_HLA_DR$nu0 <- rep(30, 3)
    prior_HLA_DR$Omega0 <- array(rep(0.1^2, 3), dim = c(3, 1, 1))
    prior_HLA_DR$Lambda0 <- array(rep(0.5, 3), dim = c(3, 1, 1))

    act.filterList <- fsApply(curData, flowClust.act, params = c(xParam, yParam),
                              filterId = filterId, prior_x = prior_CD38,
                              prior_y = prior_HLA_DR, nu.est = 2, ...)
    act.fl <- filterList(act.filterList)        
  
    nodeID <- add(wf, act.fl, parent = pViewName)
    recompute(wf, nodeID)
    
    message("done.")  

    if (plot) {
      print(plotGate(wf, nodeID, xbin = xbin, main = tView))

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

####################################################
## central memory, naive, effector, effector memory
####################################################
gating.mem_CD4 <- function(x, wf, pViewName, plot = FALSE, batch = FALSE, xbin = 128,
                       step = 100, nslaves = 0, gating_method = "flowClust",
                       usePrior = "yes", prior_CCR7 = NULL, prior_CD45RA = NULL,
                       ...) {
  nodeNames <- getNodes(wf[[1]])
  curData <- getData(wf, pViewName)
  markers <- memMarkers(x)
  memView <- paste0(pViewName, " ", markers[1], "\\+", markers[2], "\\-")
  vMem <- paste0(pViewName, " ", markers[1], "+", markers[2], "-")
  
  if (step >= 8 && !any(grepl(memView, nodeNames))) {
    message("Memory, effector gating...", pViewName)

    xChannel <- markers[1]
    yChannel <- markers[2]

    # Pools the samples in the flowSet into a flowFrame from which we obtain the
    # marker-specific priors.
    x <- exprs(as(curData, "flowFrame"))[, markers2channels(curData[[1]], xChannel)]

    # Prior elicitation for HLA-DR
    peak_CCR7 <- find_peaks(x, peaks = 1)
    huber_s <- huber(x)$s

    prior_CCR7$Mu0 <- matrix(c(peak_CCR7 - 3 * huber_s, peak_CCR7 - huber_s, peak_CCR7), nrow = 3)
    prior_CCR7$w0 <- rep(10, 3)
    prior_CCR7$Omega0 <- array(rep(0.05^2, 3), dim = c(3, 1, 1))
    prior_CCR7$Lambda0 <- array(rep(0.5, 3), dim = c(3, 1, 1))
    prior_CCR7$nu0 <- rep(50, 3)
    CCR7.filterList <- Gating1D(curData, y = xChannel, nslaves = nslaves,
                                filterId = markers[1], method = gating_method,
                                usePrior = usePrior, prior = prior_CCR7,
                                nu.est = 2, K = 3, ...)

    prior_CD45RA$w0 <- rep(10, 2)
    CD45RA.filterList <- Gating1D(curData, y = yChannel, nslaves = nslaves,
                                  filterId = markers[2], method = gating_method,
                                  usePrior = usePrior, prior = prior_CD45RA,
                                  nu.est = 2, ...)

    memory.filterList <- lapply(seq_along(curData), function(i) {
      gate_CCR7 <- ifelse(is.finite(CCR7.filterList[[i]]@min), CCR7.filterList[[i]]@min, CCR7.filterList[[i]]@max)
      gate_CD45RA <- ifelse(is.finite(CD45RA.filterList[[i]]@min), CD45RA.filterList[[i]]@min, CD45RA.filterList[[i]]@max)

      flowClust_quadGate <- list(gate_CCR7, gate_CD45RA)
      names(flowClust_quadGate) <- markers2channels(curData[[1]], markers)
      quadGate(filterId = "", flowClust_quadGate)
    })
    names(memory.filterList) <- sampleNames(getData(wf))
    memory.filterList <- filterList(memory.filterList)

    memory.filterList@filterId <- pViewName

    nodeIDs <- add(wf, memory.filterList, parent = pViewName,
                   names = c("Effector", "Naive", "Central memory", "Effector memory"))
    lapply(nodeIDs, function(nodeID) {
      recompute(wf, nodeID)
    })
    
    message("done.")
    if (plot) {
      print(plotGate(wf, nodeIDs, xbin = xbin, main = pViewName))

      if (!batch) {
        retval <- readline("Enter to resume, c to exit:")
        dev.off()
        if (retval == "c") {
          opt.outer  <-  options(show.error.messages = FALSE)
          on.exit(options(opt.outer))
          stop()
        }
      }
    }
  }
}

####################################################
## central memory, naive, effector, effector memory
####################################################
gating.mem_CD8 <- function(x, wf, pViewName, plot = FALSE, batch = FALSE, xbin = 128,
                       step = 100, nslaves = 0, gating_method = "flowClust",
                       usePrior = "yes", prior_CCR7 = NULL, prior_CD45RA = NULL,
                       ...) {
  nodeNames <- getNodes(wf[[1]])
  curData <- getData(wf, pViewName)
  markers <- memMarkers(x)
  memView <- paste0(pViewName, " ", markers[1], "\\+", markers[2], "\\-")
  vMem <- paste0(pViewName, " ", markers[1], "+", markers[2], "-")
  
  if (step >= 8 && !any(grepl(memView, nodeNames))) {
    message("Memory, effector gating...", pViewName)

    xChannel <- markers[1]
    yChannel <- markers[2]

    prior_CCR7$w0 <- rep(10, 2)
    prior_CD45RA$w0 <- rep(10, 2)    

    CCR7.filterList <- Gating1D(curData, y = xChannel, nslaves = nslaves,
                                filterId = markers[1], method = gating_method,
                                usePrior = usePrior, prior = prior_CCR7,
                                nu.est = 2, ...)

    # Pools the samples in the flowSet into a flowFrame from which we obtain the
    # marker-specific priors.
    x <- exprs(as(curData, "flowFrame"))[, markers2channels(curData[[1]], yChannel)]

    # Prior elicitation for HLA-DR
    peak_CD45RA <- find_peaks(x, peaks = 1)
    huber_s <- huber(x)$s

    prior_CD45RA$Mu0 <- matrix(c(peak_CD45RA - 4 * huber_s, peak_CD45RA - 2.5 * huber_s, peak_CD45RA), nrow = 3)
    prior_CD45RA$w0 <- rep(10, 3)
    prior_CD45RA$nu0 <- rep(10, 3)
    prior_CD45RA$Omega0 <- array(rep(0.03^2, 3), dim = c(3, 1, 1))
    prior_CD45RA$Lambda0 <- array(rep(0.5, 3), dim = c(3, 1, 1))

    CD45RA.filterList <- Gating1D(curData, y = yChannel, nslaves = nslaves,
                                  filterId = markers[2], method = gating_method,
                                  usePrior = usePrior, prior = prior_CD45RA,
                                  nu.est = 2, K = 3, neg_cluster = 2,
                                  cutpoint_method = "posterior_mean", ...)

    memory.filterList <- lapply(seq_along(curData), function(i) {
      gate_CCR7 <- ifelse(is.finite(CCR7.filterList[[i]]@min), CCR7.filterList[[i]]@min, CCR7.filterList[[i]]@max)
      gate_CD45RA <- ifelse(is.finite(CD45RA.filterList[[i]]@min), CD45RA.filterList[[i]]@min, CD45RA.filterList[[i]]@max)

      flowClust_quadGate <- list(gate_CCR7, gate_CD45RA)
      names(flowClust_quadGate) <- markers2channels(curData[[1]], markers)
      quadGate(filterId = "", flowClust_quadGate)
    })
    names(memory.filterList) <- sampleNames(getData(wf))
    memory.filterList <- filterList(memory.filterList)

    memory.filterList@filterId <- pViewName

    nodeIDs <- add(wf, memory.filterList, parent = pViewName,
                   names = c("Effector", "Naive", "Central memory", "Effector memory"))
    lapply(nodeIDs, function(nodeID) {
      recompute(wf, nodeID)
    })
    
    message("done.")
    if (plot) {
      print(plotGate(wf, nodeIDs, xbin = xbin, main = pViewName))

      if (!batch) {
        retval <- readline("Enter to resume, c to exit:")
        dev.off()
        if (retval == "c") {
          opt.outer  <-  options(show.error.messages = FALSE)
          on.exit(options(opt.outer))
          stop()
        }
      }
    }
  }
}

