#################################
## T cell panel gating
#################################
setMethod("gating", signature = c("ICS", "GatingSet"), definition = function(x, wf, ...) {
  lView <- "root"
  gating.singlet.h(x, wf, pViewName = lView, ...)

  lView <- "singlet"
  gating.viable(x, wf, pViewName = lView, ...)

  lView <- "viable"
  gating.nonDebris(x, wf, pViewName = lView, ...)
  
  lView <- "nonDebris"
  gating.lymph(x, wf, pViewName = lView, ...)
  gc()

  lView <- "Lymph"
  gating.cd3(x, wf, pViewName = lView, xChannel = "SSC-A", gating_method="density1d",trans = 0, ...)
  gc()

  lView <- x@tMarkers
  gating.tsub(x, wf, pViewName = lView,split=FALSE, ...)
  gc()

  tView.cd4 <- "CD4"
  gating.cytokine(x, wf, pViewName = tView.cd4, ...)
  gating.polyfunction(x, wf, pViewName = tView.cd4, ...)
  gc()

  tView.cd8 <- "CD8"

  gating.cytokine(x, wf, pViewName = tView.cd8, ...)
  gating.polyfunction(x, wf, pViewName = tView.cd8, ...)
	
	message("finished.")
})

gating.polyfunction <- function(x, wf, pViewName, plot = FALSE, batch = FALSE,
                                xbin = 128, step = 100, cond = "name",
                                nslaves = 0, ...) {
  cks <- x@ckMarkers
  cks <- paste0(cks, "+")
  ckNodes <- paste(pViewName, cks, sep = "_")
  nMarkers <- length(cks)
  nodeNames <- getNodes(wf[[1]])

  if(step >= 7) {
    ## all the comibnations of A & B & C

    opList <- permutations(n = 1, r = nMarkers - 1, c("&"), repeats = TRUE)
    isNotList <- permutations(n = 2, r = nMarkers, c("!", ""), repeats = TRUE)
    polyExprsList <- apply(opList, 1, function(curOps) {
      apply(isNotList, 1, function(curIsNot) {
        polyExprs <- curIsNot
        polyExprs[-1] <- paste0(curOps, curIsNot[-1])

        paste(paste0(polyExprs, ckNodes), collapse = "")
      })
    })
    polyExprsList <- as.vector(polyExprsList)

    #IL2|IFNg
    polyExpr1 <- paste(ckNodes[1:2], collapse = "|")
    polyExprsList <- c(polyExpr1, polyExprsList)

    #IL2|IFNg|TNFa
    polyExpr2 <- paste(ckNodes, collapse = "|")
    polyExprsList <- c(polyExpr2, polyExprsList)
			
    #actual gating
    lapply(polyExprsList, function(polyExpr) {
      polyExpr <- as.symbol(polyExpr)
      tNodes <- gsub(paste0(pViewName, "_"), "", polyExpr)
      tNodes <-paste(pViewName, tNodes, sep = ":")
      # Constructs intersections of the CD19 and CD20 gates.
      if (!any(tNodes %in% nodeNames)) {
        message(tNodes, " gating...")
        bf <- eval(substitute(booleanFilter(x), list(x = polyExpr)))
        bf@filterId <- tNodes
        invisible(nodeID <- add(wf, bf, parent = pViewName))
        invisible(recompute(wf, nodeID))
        message("done.")
      }
    })					
  }
}

gating.cytokine <- function(x, wf, pViewName, plot = FALSE, batch = FALSE,
                            xbin = 128, step = 100, cond = "name", nslaves = 0,
                            quantile_IL2 = 0.995, quantile_IFNg = 0.995,
                            quantile_TNFa = 0.995, ...) {

	nodeNames <- getNodes(wf[[1]])
	curData <- getData(wf,pViewName)
	markers <- x@ckMarkers
	lapply(markers, function(marker) {
    tRegView <- paste0(pViewName, "_", marker, "\\+")
		tView <- paste0(pViewName, "_", marker, "+")
		
		if (step >= 6 && !any(grepl(tRegView, nodeNames))) {
			message(marker, " gating (", pViewName, ") ...")

			yChannel <- marker
			
      # Pools the data from the flowFrames within the flowSet into a single
      # flowFrame, from which a prior is elicited.
      # NOTE: If there is a memory issue with coercing the flowSet into a single
      # flowFrame, perhaps we should randomly select a subset of the flowFrames
      # and then elicit the prior.
      marker_channel <- markers2channels(curData[[1]], yChannel)
      x_cytokines <- exprs(as(curData, "flowFrame"))[, marker_channel]

      prior <- list()
      cytokines_peak <- find_peaks(x_cytokines, peaks = 1, adjust = 3)

      huber_s <- huber(x_cytokines)$s

      prior$Mu0 <- matrix(c(cytokines_peak - 2 * huber_s, cytokines_peak, cytokines_peak + 2 * huber_s), nrow = 3)

      prior$Omega0 <- array(rep(0.05^2, 3), dim = c(3, 1, 1))
      prior$Lambda0 <- array(rep(0.5^2, 3), dim = c(3, 1, 1))
      prior$w0 <- rep(10, 3)
      prior$nu0 <- rep(30, 3)

      if (yChannel == "IL2") {
        percentile <- quantile_IL2
      } else if (yChannel == "IFNg") {
        percentile <- quantile_IFNg
      } else { # yChannel == "TNFa"
        percentile <- quantile_TNFa
      }

      ck.filterList <- Gating1D(curData, y = yChannel, filterId = tView,
                                method = "flowClust", nslaves = nslaves,
                                usePrior = "yes", prior = prior, K = 3, plot = FALSE,
                                truncate_min = -0.5, trans = 0,
                                cutpoint_method = "quantile", percentile = percentile, ...)

			nodeID <- add(wf, ck.filterList, parent = pViewName)
			recompute(wf, nodeID)			
			
			message("done.")

      if (plot) {
        print(plotGate(wf, nodeID, main = tView, xbin = xbin))

        if(!batch) {
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
  })
}
