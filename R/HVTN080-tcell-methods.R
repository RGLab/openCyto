#################################
## HVTN080 gating
#################################
setMethod("gating", signature = c("HVTN080", "GatingSet"), definition = function(x, wf, ...) {
  parent <- "root"

  gate_singlet(gating_template = x, gating_set = wf, parent = parent,
               x_channel = "FSC-H", y_channel = "FSC-A", prediction_level = 0.999, ...)

  parent <- "singlet"

  channel_vMarkers <- markers2channels(getData(wf)[[1]], vMarkers(x))

  # Elicitation of priors for the viable gate
  prior_viable <- list()
  prior_viable$xChannel <- prior_flowClust1d(flow_set = getData(wf, parent),
                                              channel = channel_vMarkers, K = 3,
                                              adjust = 1.5)
  prior_viable$yChannel <- prior_flowClust1d(flow_set = getData(wf, parent),
                                              channel = "FSC-A", K = 3, adjust = 1.5)

  gating.viable(x, wf, pViewName = parent, xChannel = vMarkers(x), yChannel = "FSC-A",
                usePrior = "yes", K = 3, prior = prior_viable, ...)

  parent <- "viable"

  # Elicitation of prior for the Lymphocyte gate
  prior_lymph <- prior_flowClust1d(flow_set = getData(wf, parent),
                                   channel = "SSC-A", K = 3, adjust = 1.5)
  gating.lymph(x, wf, pViewName = parent, xChannel = "SSC-A", K = 3,
               which_gate = "left", prior = prior_lymph, ...)
  gc()

  parent <- "Lymph"

  # Elicitation of prior for the CD3 gate
  channel_CD3 <- markers2channels(getData(wf)[[1]], tMarkers(x))
  prior_CD3 <- prior_flowClust1d(flow_set = getData(wf, parent),
                                   channel = channel_CD3, K = 3, adjust = 1.5)
  gating.cd3(x, wf, pViewName = parent, xChannel = channel_CD3, usePrior = "yes",
             trans = 0, prior = prior_CD3, neg_cluster = 2, K = 3, ...)
  gc()

  parent <- "CD3"

  # Elicitation of priors for the CD4 and CD8 gates
  channel_tsubMarkers <- markers2channels(getData(wf)[[1]], tsubMarkers(x))
  prior_CD4 <- prior_flowClust1d(flow_set = getData(wf, parent),
                                 channel = channel_tsubMarkers[1], K = 3,
                                 adjust = 1.5)
  prior_CD8 <- prior_flowClust1d(flow_set = getData(wf, parent),
                                 channel = channel_tsubMarkers[2], K = 3,
                                 adjust = 1.5)
  gating.HVTN080_tsub(x, wf, pViewName = parent, K = 3, prior_CD4 = prior_CD4,
                      prior_CD8 = prior_CD8, ...)
  gc()

  parent <- "CD4"

  # Elicitation of prior for the CD57 gate
  channel_CD57 <- markers2channels(getData(wf)[[1]], "CD57")
  prior_CD57 <- prior_flowClust1d(flow_set = getData(wf, parent),
                                  channel = channel_CD57, K = 3)
  gating.CD57(x, wf, pViewName = parent, K = 3, prior = prior_CD57, cutpoint_method = "quantile",
              quantile = 0.975, ...)

  # Elicitation of priors for the Perforin and Granzyme B gates
  channel_Perforin <- markers2channels(getData(wf)[[1]], "Perforin")
  channel_GzB <- markers2channels(getData(wf)[[1]], "Granzyme B")
  prior_Perforin <- prior_flowClust1d(flow_set = getData(wf, parent),
                                      channel = channel_Perforin, K = 3)
  prior_GzB <- prior_flowClust1d(flow_set = getData(wf, parent),
                                 channel = channel_GzB, K = 3)

  gating.Perforin_GzB(x, wf, pViewName = parent, prior_Perforin = prior_Perforin,
                      prior_GzB = prior_GzB, cutpoint_method = "quantile", K = 3, quantile = 0.95, ...)

  # Elicitation of priors for the IL2 and IFNg gates
  channel_IL2 <- markers2channels(getData(wf)[[1]], "IL2")
  channel_IFNg <- markers2channels(getData(wf)[[1]], "IFNg")
  prior_IL2 <- prior_flowClust1d(flow_set = getData(wf, parent),
                                      channel = channel_IL2, K = 3)
  prior_IFNg <- prior_flowClust1d(flow_set = getData(wf, parent),
                                 channel = channel_IFNg, K = 3)

  gating.IL2_IFNg(x, wf, pViewName = parent, prior_IL2 = prior_IL2,
                      prior_IFNg = prior_IFNg, cutpoint_method = "quantile", K = 3, ...)
  
  # Elicitation of prior for the TNFa gate
  channel_TNFa <- markers2channels(getData(wf)[[1]], "TNFa")
  prior_TNFa <- prior_flowClust1d(flow_set = getData(wf, parent),
                                  channel = channel_TNFa, K = 3)
  gating.TNFa(x, wf, pViewName = parent, K = 3, prior = prior_TNFa,
              cutpoint_method = "quantile", quantile = 0.999, ...)

  # Constructs polyfunctional gates following manual gating strategy
  gating.HVTN080_polyfunction(x, wf, pViewName = parent, ...)

  gc()

  parent <- "CD8"

  # Elicitation of prior for the CD57 gate
  prior_CD57 <- prior_flowClust1d(flow_set = getData(wf, parent),
                                  channel = channel_CD57, K = 3)
  gating.CD57(x, wf, pViewName = parent, K = 3, prior = prior_CD57, neg_cluster = 2, ...)

  # Elicitation of priors for the Perforin and Granzyme B gates
  channel_Perforin <- markers2channels(getData(wf)[[1]], "Perforin")
  channel_GzB <- markers2channels(getData(wf)[[1]], "Granzyme B")
  prior_Perforin <- prior_flowClust1d(flow_set = getData(wf, parent),
                                      channel = channel_Perforin, K = 3)
  prior_GzB <- prior_flowClust1d(flow_set = getData(wf, parent),
                                 channel = channel_GzB, K = 3)

  gating.Perforin_GzB(x, wf, pViewName = parent, prior_Perforin = prior_Perforin,
                      prior_GzB = prior_GzB, cutpoint_method = "quantile", K = 3, quantile = 0.95, ...)

  # Elicitation of priors for the IL2 and IFNg gates
  channel_IL2 <- markers2channels(getData(wf)[[1]], "IL2")
  channel_IFNg <- markers2channels(getData(wf)[[1]], "IFNg")
  prior_IL2 <- prior_flowClust1d(flow_set = getData(wf, parent),
                                      channel = channel_IL2, K = 3)
  prior_IFNg <- prior_flowClust1d(flow_set = getData(wf, parent),
                                 channel = channel_IFNg, K = 3)

  gating.IL2_IFNg(x, wf, pViewName = parent, prior_IL2 = prior_IL2,
                      prior_IFNg = prior_IFNg, cutpoint_method = "quantile", K = 3, ...)
  
  # Elicitation of prior for the TNFa gate
  channel_TNFa <- markers2channels(getData(wf)[[1]], "TNFa")
  prior_TNFa <- prior_flowClust1d(flow_set = getData(wf, parent),
                                  channel = channel_TNFa, K = 3)
  gating.TNFa(x, wf, pViewName = parent, K = 3, prior = prior_TNFa,
              cutpoint_method = "quantile", quantile = 0.999, ...)

  # Constructs polyfunctional gates following manual gating strategy
  gating.HVTN080_polyfunction(x, wf, pViewName = parent, ...)

  gc()

  parent <- "not CD4"
  prior_IFNg <- prior_flowClust1d(flow_set = getData(wf, parent),
                                 channel = channel_IFNg, K = 3)

  gating.notCD4_IFNg(x, wf, pViewName = parent, prior_IFNg = prior_IFNg,
                     cutpoint_method = "quantile", K = 3, ...)

	message("finished.")
})

gating.HVTN080_tsub <- function(x, wf, pViewName, split = TRUE, plot = FALSE,
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
    prior_list <- list(xChannel = prior_CD4, yChannel = prior_CD8)

    cd4cd8.filterList <- Gating1D(curData, x = xChannel, y = yChannel,
                                  trans = 0, usePrior = "yes", prior = prior_list,
                                  neg_cluster = 2, nslaves = nslaves, ...)
      
    cd4cd8.filterList <- lapply(cd4cd8.filterList, quadGate2rectangleGates,
                                markers = markers, channels = channels,
                                quadrants = c(1, 3, 4))

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
  
    not_cd4.list <- lapply(cd4cd8.fl, function(curFilters) {
      curFilters[["CD4-CD8-"]]@filterId <- "not CD4"
      curFilters[["CD4-CD8-"]]
    })
    not_cd4.list <- filterList(not_cd4.list)

    #############
    #add to wf
    ##############
    nodeID1 <- add(wf, cd4.list, parent = pViewName)
    recompute(wf, nodeID1)
  
    nodeID2 <- add(wf, cd8.list, parent = pViewName)
    recompute(wf, nodeID2)

    nodeID3 <- add(wf, not_cd4.list, parent = pViewName)
    recompute(wf, nodeID3)

    
    message("done.")  
    if (plot) {
      print(plotGate(wf, c(nodeID1, nodeID2, nodeID3), xbin = xbin))
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
## CD57+ and CD57-
#################################
gating.CD57 <- function(x, wf, pViewName, plot = FALSE, batch = FALSE, xbin = 128,
                        step = 100, nslaves = 0, cond = "name", xChannel = "CD57",
                        filterId = "CD57", isflowClust = TRUE, usePrior = "yes",
                        ...) {

  nodeNames <- getNodes(wf[[1]])
  tView <- xChannel

  if (step >= 3) {
    message("CD57 gating...")
    curData <- getData(wf, pViewName)

    # CD57+ gates
    CD57.filterList <- Gating1D(curData, y = xChannel, filterId = "CD57",
                                    positive = TRUE, nslaves = nslaves,
                                    method = "flowClust", usePrior = usePrior, ...)

    nodeID <- add(wf, CD57.filterList, parent = pViewName)
    recompute(wf, nodeID)

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
## Perforin+, GzB+, and GzB-
#################################
gating.Perforin_GzB <- function(x, wf, pViewName, plot = FALSE, batch = FALSE,
                                xbin = 128, step = 100, nslaves = 0, cond = "name",
                                isflowClust = TRUE, usePrior = "yes",
                                xChannel = "Perforin", yChannel = "Granzyme B",
                                prior_Perforin, prior_GzB, ...) {

  nodeNames <- getNodes(wf[[1]])
  tView <- xChannel

  if (step >= 3) {
    message("Perforin and GzB gating...")
    curData <- getData(wf, pViewName)

    prior_list <- list(xChannel = prior_Perforin, yChannel = prior_GzB)
    Perforin_GzB.filterList <- Gating1D(curData, x = xChannel, y = yChannel,
                                        nslaves = nslaves, method = "flowClust",
                                        usePrior = usePrior, prior = prior_list, ...)

    markers <- c(xChannel, yChannel)
    channels <- markers2channels(getData(wf[[1]]), markers)
    Perforin_GzB.filterList <- lapply(Perforin_GzB.filterList, quadGate2rectangleGates,
                                      markers = markers, channels = channels,
                                      quadrants = c(1, 3))
    Perforin_GzB.filterList <- filtersList(Perforin_GzB.filterList)

    Perforin.list <- lapply(Perforin_GzB.filterList, function(filter_i) {
      filter_i[["Perforin+Granzyme B-"]]@filterId <- "Perforin"
      filter_i[["Perforin+Granzyme B-"]]
    })
    Perforin.list <- filterList(Perforin.list)

    GzB.list <- lapply(Perforin_GzB.filterList, function(filter_i) {
      filter_i[["Perforin-Granzyme B+"]]@filterId <- "GzB"
      filter_i[["Perforin-Granzyme B+"]]
    })
    GzB.list <- filterList(GzB.list)

    nodeID <- add(wf, Perforin.list, parent = pViewName)
    recompute(wf, nodeID)

    nodeID <- add(wf, GzB.list, parent = pViewName)
    recompute(wf, nodeID)

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
## IL2+ and IFNg+
#################################
gating.IL2_IFNg <- function(x, wf, pViewName, plot = FALSE, batch = FALSE,
                                xbin = 128, step = 100, nslaves = 0, cond = "name",
                                isflowClust = TRUE, usePrior = "yes",
                                xChannel = "IL2", yChannel = "IFNg",
                                prior_IL2, prior_IFNg, ...) {

  nodeNames <- getNodes(wf[[1]])
  tView <- xChannel

  if (step >= 3) {
    message("IL2 and IFNg gating...")
    curData <- getData(wf, pViewName)

    prior_list <- list(xChannel = prior_IL2, yChannel = prior_IFNg)
    IL2_IFNg.filterList <- Gating1D(curData, x = xChannel, y = yChannel,
                                        nslaves = nslaves, method = "flowClust",
                                        usePrior = usePrior, prior = prior_list, ...)

    markers <- c(xChannel, yChannel)
    channels <- markers2channels(getData(wf[[1]]), markers)
    IL2_IFNg.filterList <- lapply(IL2_IFNg.filterList, quadGate2rectangleGates,
                                      markers = markers, channels = channels,
                                      quadrants = c(1, 3))
    IL2_IFNg.filterList <- filtersList(IL2_IFNg.filterList)

    IL2.list <- lapply(IL2_IFNg.filterList, function(filter_i) {
      filter_i[["IL2+IFNg-"]]@filterId <- "IL2"
      filter_i[["IL2+IFNg-"]]
    })
    IL2.list <- filterList(IL2.list)

    IFNg_pos.list <- lapply(IL2_IFNg.filterList, function(filter_i) {
      filter_i[["IL2-IFNg+"]]@filterId <- "IFNg"
      filter_i[["IL2-IFNg+"]]
    })
    IFNg_pos.list <- filterList(IFNg_pos.list)

    nodeID <- add(wf, IL2.list, parent = pViewName)
    recompute(wf, nodeID)

    nodeID <- add(wf, IFNg_pos.list, parent = pViewName)
    recompute(wf, nodeID)

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
## TNFa+
#################################
gating.TNFa <- function(x, wf, pViewName, plot = FALSE, batch = FALSE, xbin = 128,
                        step = 100, nslaves = 0, cond = "name", xChannel = "TNFa",
                        filterId = "TNFa", isflowClust = TRUE, usePrior = "yes",
                        ...) {

  nodeNames <- getNodes(wf[[1]])
  tView <- xChannel

  if (step >= 3) {
    message("TNFa gating...")
    curData <- getData(wf, pViewName)

    # TNFa+ gates
    TNFa_pos.filterList <- Gating1D(curData, y = xChannel, filterId = filterId,
                                    positive = TRUE, nslaves = nslaves,
                                    method = "flowClust", usePrior = usePrior,
                                    ...)
    
    nodeID <- add(wf, TNFa_pos.filterList, parent = pViewName)
    recompute(wf, nodeID)

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
## !CD4+/IFNg+
#################################
gating.notCD4_IFNg <- function(x, wf, pViewName, plot = FALSE, batch = FALSE,
                               xbin = 128, step = 100, nslaves = 0, cond = "name",
                               isflowClust = TRUE, usePrior = "yes",
                               xChannel = "IFNg", prior_IFNg, ...) {

  nodeNames <- getNodes(wf[[1]])
  tView <- xChannel

  if (step >= 3) {
    message("IFNg gating...")
    curData <- getData(wf, pViewName)

    IFNg.filterList <- Gating1D(curData, y = xChannel, filterId = "IFNg",
                                nslaves = nslaves, method = "flowClust",
                                usePrior = usePrior, prior = prior_IFNg, ...)

    nodeID <- add(wf, IFNg.filterList, parent = pViewName)
    recompute(wf, nodeID)

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

gating.HVTN080_polyfunction <- function(x, wf, pViewName, plot = FALSE,
                                        batch = FALSE, xbin = 128, step = 100,
                                        cond = "name", nslaves = 0, ...) {
  nodeNames <- getNodes(wf[[1]])
  if(step >= 7) {
    # Cytokines named properly
    polyExprsList <- polyfunction_nodes(c("IFNg", "IL2", "TNFa", "GzB", "CD57"))
    polyExprsList <- c(polyExprsList, polyfunction_nodes(c("IFNg", "IL2", "TNFa", "GzB")))
    polyExprsList <- c(polyExprsList, polyfunction_nodes(c("IFNg", "IL2")))

    # actual gating
    lapply(polyExprsList, function(polyExpr) {
      polyExpr <- as.symbol(polyExpr)
      tNodes <- paste(pViewName, polyExpr, sep = ":")
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
