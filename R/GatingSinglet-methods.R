gate_singlet <- function(gating_template, gating_set, parent,
                         filter_id = "singlet", x_channel = "FSC-A",
                         y_channel = "FSC-H", prediction_level = 0.99,
                         num_nodes = 1, parallel_type = c("multicore", "sock"),
                         plot = FALSE, xbin = 128, ...) {

  require('flowStats')
  require('parallel')

  parallel_type <- match.arg(parallel_type)
  nodes <- getNodes(gating_set[[1]])
  
  if (!any(grepl(filter_id, nodes))) {
    node_data <- getData(gating_set, parent)

    message("Constructing ", filter_id, " gates...")
    
    # Splits the flow set into a list, where each element in the list is a
    # flowSet containg one flow frame.
    # Below, we access manually the individual flow frame with current_frame[[1]].
    flowset_list <- split(node_data, sampleNames(node_data))

    if (num_nodes > 1) {
      message("Running in parallel mode with ", num_nodes, " nodes.")
    }

    # Creates a list of polygon gates based on the prediction bands at the minimum and maximum
    # x_channel observation using a robust linear model trained by flowStats.
    if (parallel_type == "multicore") {
      singlet_list <- mclapply(flowset_list, mc.cores = num_nodes, function(flow_frame) {
        singletGate(flow_frame[[1]], area = x_channel, height = y_channel,
                    prediction_level = prediction_level, filter_id = filter_id)$gate
      })
    } else {
      cl <- makeCluster(num_nodes, type = "SOCK")
      singlet_list <- parLapply(cl, flowset_list, function(flow_frame) {
        singletGate(flow_frame[[1]], area = x_channel, height = y_channel,
                    prediction_level = prediction_level, filter_id = filter_id)$gate
      })
      stopCluster(cl)
    }

    # Adds the list of singlet polygon gates to the workflow.
    singlet_list <- filterList(singlet_list)
    node_id <- add(gating_set, singlet_list, parent = parent)
    recompute(gating_set, node_id)
    message("done.")
  }
  
  if (plot) {
    print(plotGate(gating_set, node_id, xbin = xbin, pos = c(0.5, 0.8)))
  }
}

