#' Pre-processing plugin for manual/templated gating
#'
#' @param fs a \code{flowSet} object containing the samples within a group to be
#'   gated.
#' @param gs the original \code{GatingSet} containing all samples with attached
#'   compensation, transformations and gates.
#' @param gm the name of the gating method.
#' @param channels markers or channels to use for gating.
#' @param groupBy a colon separated character string specifying the names of the
#'   variables to use to group the samples prior to gating, set to NA by default
#'   to group all samples. 
#' @param isCollapse logical indicating whether the data for the current group
#'   is collapsed prior to gating, set to NA by default.
#' @param ... not in use.
#'
#' @return \code{pp_gate_template()} returns the name of the current group so it
#'   can be used to extract the appropriate gate(s) in \code{gate_template()}.
#'
#' @importFrom flowWorkspace sampleNames pData
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@ozette.com}
#'
#' @noRd
.pp_gate_template <- function(fs,
                              gs,
                              gm,
                              channels,
                              groupBy = NA,
                              isCollapse = NA,
                              ...) {
  
  # gates are created for each unique group in groupBy
  # fs contains samples for current group to gate
  # we obtain group index by searching for group containing all samples in fs
  # group is identified by name to ensure correct matching
  # names used in gates must include group variables separated by " "
  
  # prepare groupBy - split out grouping variables
  if(is.character(groupBy)) {
    if(nchar(groupBy) == 0) {
      groupBy <- NA
    } else {
      groupBy <- unlist(strsplit(groupBy, ":"))
    }
  }
  
  # experiment details
  pd <- pData(gs)
  
  # grouping
  if(!all(is.na(groupBy))) {
    # variable names supplied
    if(all(grepl("^[A-Za-z]+$", groupBy, ignore.case = TRUE))) {
      # check for valid variable names
      if(!all(groupBy %in% colnames(pd))) {
        stop(
          "'groupBy' references variables that don't exist in the GatingSet!"
        )
      }
      # split experiment details by grouping variables
      groups <- split(
        pd, pd[, groupBy],
        sep = " ",
        lex.order = TRUE,
        drop = TRUE
      )
      # locate group containing all samples in fs -> group name
      ind <- names(groups)[
        unlist(
          lapply(
            groups,
            function(group) {
              all(sampleNames(fs) %in% c(rownames(group), group$name))
            }
          )
        )
      ]
    # numeric grouping is not supported - inconsistent gate order
    } else {
      if(groupBy == length(gs)) {
        ind <- 1
      } else {
        message(
          "Numeric groupBy is not supported. Grouping all samples."
        )
        ind <- 1
      }
    }
    # no grouping - use first set of gates
  } else {
    ind <- 1
  }
  
  # return group index or name (preferred)
  return(ind)
  
}

#' Gating plugin for manual/templated gating
#'
#' @param fr a \code{flowFrame} object containing the collapsed data for a
#'   sample group to be gated.
#' @param pp_res indicates the name of the current sample group as provided by
#'   \code{pp_gate_template()}.
#' @param channels names of the channels or markers in which the supplied
#'   \code{gate(s)} were constructed.
#' @param gate list of filters per sample group each named based on the levels
#'   of the grouping variables separated by \code{" "}.
#'   \code{generate_gate_template()} should be used to create the
#'   \code{gatingTemplate} so that the gates can be appropriately formatted and
#'   passed to \code{gate_template()} through the gatingTemplate
#'   \code{gating_args}.
#' @param ... not in use.
#'
#' @return a list of filters for the current sample group to be used for gating.
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@ozette.com}
#'
#' @noRd
.gate_template <- function(fr,
                          pp_res,
                          channels,
                          gate,
                          ...) {
  
  
  # gate = list of filters per sample group (all groups)
  # pp_res = name of sample group to extract gate(s) from gate
  
  # no sample groups - apply first set of gates to all samples
  if(is.null(pp_res)) {
    ind <- 1
  }
  
  # pp_res is index of gate(s) to extract for sample fr
  if(is.numeric(pp_res)) {
    ind <- pp_res
  # pp_res is the name of sample group
  } else if(is.character(pp_res)) {
    # attempt to locate gate for sample group
    ind <- match(
      pp_res,
      names(gate)
    )
    # sample group doesn't exist in the gatingTemplate
    if(all(is.na(ind))) {
      stop(
        "The gatingTemplate doesn't contain gates for the ",
        pp_res,
        " sample group!"
      )
    }
  }
  
  # extract and apply matching set of gates
  return(gate[[ind]])
  
}

#' Store gates from a GatingSet in a gatingTemplate for templated gating
#'
#' \code{generate_gate_template()} iterates through samples groups in a
#' GatingSet to extract and store existing gates into an openCyto
#' \code{gatingTemplate} that can be easily applied to new GatingSets using
#' \code{gt_gating()}.
#'
#' @param x a \code{GatingHierarchy} or \code{GatingSet} object.
#' @param groupBy a vector of variable names that exist in
#'   \code{colnames(pData(x))} to split samples into groups before extracting
#'   the gates specified by \code{nodes}, set to NA by default to produce
#'   consensus gates across all samples.
#' @param nodes a vector of full node paths that exist in \code{x} for which
#'   templated gating entries are required, set to all nodes in \code{x} by
#'   default.
#' @param gatingTemplate logical to indicate whether a \code{gatingTemplate}
#'   object should be returned, set to FALSE by default. \code{gatingTemplate}
#'   objects can only be created when all nodes in \code{x} are specified.
#' @param save_as name of the CSV file to which the gatingTemplate entries
#'   should be written, set to NULL to prevent the gatingTemplate from being
#'   written to a CSV file.
#' @param ... additional arguments passed to \code{\link{gatingTemplate-class}}.
#'
#' @return either a \code{data.table} or \code{gatingTemplate} object containing
#'   the requested gatingTemplate entries for downstream templated gating using
#'   \code{gt_gating()}.
#'
#' @importFrom flowCore polygonGate rectangleGate quadGate filters
#' @importFrom flowWorkspace pData gh_pop_get_gate gh_get_pop_paths
#'   gh_pop_get_data gh_pop_is_negated
#' @importFrom grDevices chull
#' @importFrom methods as
#' @importFrom data.table fwrite
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@ozette.com}
#'
#' @examples
#' \dontrun{
#' # gs contains variable BioSampleType
#' pData(gs)$BioSampleType <- rep(c("cells", "beads"), each = 4)
#' # create gatingTemplate object
#' gt <- generate_gate_template(
#'   gs,
#'   groupBy = "BioSampleType",
#'   gatingTemplate = TRUE,
#'   save_as = "gatingTemplate.csv"
#' )
#' # gt is a gatingTemplate object
#' gt
#' # apply gt to a new GatingSet
#' gt_gating(
#'   gt,
#'   gs_new
#' )
#' }
#'
#' @export
generate_gate_template <- function(x,
                                   groupBy = NA,
                                   nodes = NULL,
                                   gatingTemplate = FALSE,
                                   save_as = NULL,
                                   ...) {
  
  # NOTE: All samples in x must contain the same set of gates! We cannot
  # generate gatingTemplate entries for clusters (i.e. labelled events) as
  # these cannot be transferred to new GatingSets - only valid gate objects
  # (e.g. booleanFilter, rectangleGate, ellipsoidGate, quadrantGate etc.) are
  # supported.
  
  # experiment details
  pd <- pData(x)
  
  # extract full node paths
  nodes_full <- gh_get_pop_paths(
    x[[1]],
    path = "full"
  )
  
  # default nodes
  if(is.null(nodes)) {
    nodes <- nodes_full
  # manual nodes
  } else {
    if(!all(nodes %in% nodes_full)) {
      stop(
        "'nodes' must contain full paths to nodes existing in the ",
        class(x), "!"
      )
    }
  }
  
  # prepare groupBy - split samples into groups
  if(!all(is.na(groupBy))) {
    # split grouping variables
    if(length(groupBy) == 1 & grepl("\\:", groupBy)) {
      groupBy <- strsplit(
        groupBy,
        "\\:"
      )[[1]]
    }
    # check groupBy for valid variables
    if(!all(groupBy %in% colnames(pd))) {
      stop(
        "'groupBy' references variables that don't exist in pData(x)!"
      )
    }
  }
  
  # split details into groups by groupBy
  if(all(is.na(groupBy))) {
    pd_groups <- list(
      pd
    )
  # group order doesn't matter
  } else {
    pd_groups <- split(
      pd,
      pd[, groupBy],
      drop = TRUE,
      sep = " "
    )
  }
  
  # compute instrument limits - replace infinite co-ordinates
  rng <- apply(
    do.call(
      "rbind",
      lapply(
        seq_along(x),
        function(id) {
          cf <- gh_pop_get_data(
            x[[id]]
          )
          # determine range using instrument and data limits - data often lower
          rbind(
            range(cf, type = "instrument"),
            range(cf, type = "data")
          )
        }
      )
    ),
    2,
    "range"
  )
  
  # rename groupBy
  group_by <- groupBy
  
  # initialize gatingTemplate variables - CRAN R CMD CHECK NOTE
  parent <- NULL
  alias <- NULL
  pop <- NULL
  gating_method <- NULL
  gating_args <- NULL
  collapseDataForGating <- NULL
  groupBy <- NULL
  preprocessing_method <- NULL
  
  # generate a bare bones gatingTemplate 
  gt <- gh_generate_template(x[[1]])
  
  # split bare bones gatingTemplate by parent and dims
  gt_chunks <- split(
    gt,
    gt[, c("parent", "dims")],
    sep = "|",
    drop = TRUE
  )
  
  # order gt_chunks by nodes - ensures node order correct in gatingTemplate
  idx <- sapply(
    gt_chunks,
    function(gt_chunk) {
      parent <- gt_chunk$parent
      parent[parent == "root"] <- ""
      alias <- unlist(strsplit(gt_chunk$alias, ","))
      min(
        match(
          paste0(
            parent,
            "/",
            alias
          ),
          nodes
        )
      )
    }
  )
  
  # drop unwanted nodes
  gt_chunks <- gt_chunks[!is.na(idx)]
  idx <- idx[!is.na(idx)]
  
  # reorder gatingTemplate chunks
  gt_chunks <- gt_chunks[order(idx)]
  
  # construct gatingTemplate
  gt_new <- do.call(
    "rbind",
    structure(
      lapply(
        gt_chunks,
        function(gt_chunk) {
          # parent
          parent <- unique(gt_chunk$parent)
          # channels
          channels <- strsplit(unique(gt_chunk$dims), ",")[[1]]
          # alias
          alias <- gt_chunk$alias
          # REQUIRE: list of filters per population in each sample
          group_gate_list <- structure(
            lapply(
              pd_groups,
              function(pd_group) {
                # groupBy = NA -> all samples
                if(all(is.na(groupBy))) {
                  idx <- seq_along(x)
                # groupBy subset of samples
                } else {
                  idx <- match(pd_group$name, pd$name)
                }
                # extract gates per sample in group
                sample_gates <- lapply(
                  idx,
                  function(id) {
                    # extract gate objects
                    gate_list <- structure(
                      lapply(
                        alias,
                        function(pop) {
                          gate <- gh_pop_get_gate(
                            x[[id]],
                            paste0(
                              if(parent == "root") {
                                ""
                              } else {
                                parent
                              },
                              "/",
                              pop
                            )
                          )
                          # inherit negated flag
                          attributes(gate)$negated <- gh_pop_is_negated(
                            x[[id]],
                            paste0(
                              if(parent == "root") {
                                ""
                              } else {
                                parent
                              },
                              "/",
                              pop
                            )
                          )
                          return(gate)
                        }
                      ),
                      names = alias
                    )
                    # format quadrant gates
                    quad_idx <- which(
                      sapply(
                        gate_list,
                        function(gate){
                          if(inherits(gate, "rectangleGate")) {
                            any(grepl("quad", names(attributes(gate))))
                          } else {
                            FALSE
                          }
                        }
                      )
                    )
                    # quadrants located
                    if(length(quad_idx) > 0) {
                      # pull out quadrant rectangleGates
                      quads <- gate_list[quad_idx]
                      # drop quadrants from gate_list
                      gate_list <- gate_list[-quad_idx]
                      # quadrant co-ordinates
                      coords <- do.call(
                        "rbind",
                        lapply(
                          quads,
                          function(quad) {
                            structure(
                              c(quad@min, quad@max),
                              names = parameters(quad)
                            )
                          }
                        )
                      )
                      # quadrant center
                      coords <- unlist(
                        sapply(
                          channels,
                          function(z) {
                            unique(coords[, z][is.finite(coords[, z])])
                          }
                        )
                      )
                      names(coords) <- channels
                      # create quadrantGate
                      quads <- list(
                        quadGate(
                          filterId = paste0(
                            names(attributes(quads[[1]])[["quadrants"]]),
                            collapse = "|"
                          ),
                          .gate = coords
                        )
                      )
                      attributes(quads[[1]])$negated <- FALSE
                      gate_list <- c(quads, gate_list)
                    }
                    names(gate_list) <- sapply(
                      gate_list,
                      function(gate) {
                        gate@filterId
                      }
                    )
                    return(gate_list)
                  }
                )
                # construct gatingTemplate entries
                structure(
                  lapply(
                    names(sample_gates[[1]]),
                    function(gate) {
                      # extract unique gates across samples in group
                      gates <- unique(
                        lapply(
                          sample_gates,
                          `[[`,
                          gate
                        )
                      )
                      # samples in group have different gates
                      if(length(gates) > 1) {
                        # multiRangeGate or QuadGate - use first gate
                        if(inherits(gates[[1]], "multiRangeGate") |
                           inherits(gates[[1]], "quadGate")) {
                          return(gates[[1]])
                        # 1-D rectangleGate
                        } else if(inherits(gates[[1]], "rectangleGate") &
                                  length(parameters(gates[[1]])) == 1) {
                          # combine coords from each gate
                          coords <- do.call(
                            "rbind",
                            lapply(
                              gates,
                              function(pop) {
                                rbind(
                                  pop@min,
                                  pop@max
                                )
                              }
                            )
                          )
                        # coerce gates to polygon and chull
                        } else {
                          # combine coords from each gate
                          coords <- do.call(
                            "rbind",
                            lapply(
                              gates,
                              function(pop) {
                                as(pop, "polygonGate")@boundaries
                              }
                            )
                          )
                        }
                        # replace infinite co-ordinates with instrument limits
                        lapply(
                          colnames(coords),
                          function(chan) {
                            coords[, chan][
                              is.infinite(coords[, chan]) & coords[, chan] < 0
                            ] <<- min(rng[, chan])
                            coords[, chan][
                              is.infinite(coords[, chan]) & coords[, chan] > 0
                            ] <<- max(rng[, chan])
                          }
                        )
                        # 1-D gate - range
                        if(ncol(coords) == 1) {
                          # use range across all gates
                          chull <- rectangleGate(
                            .gate = apply(
                              coords,
                              2,
                              "range"
                            ),
                            filterId = gate
                          )
                        # 2-D gate - chull
                        } else {
                          # get convex hull
                          chull <- coords[chull(coords), ]
                          # chull is a rectangle
                          if(nrow(chull) < 3) {
                            # construct rectangle
                            chull <- rectangleGate(
                              chull,
                              filterId = gate
                            )
                            # convert to valid polygon
                            chull <- as(
                              chull, 
                              "polygonGate"
                            )
                          # chull is a polygon
                          } else {
                            # construct polygonGate
                            chull <- polygonGate(
                              .gate = chull,
                              filterId = gate
                            )
                          }
                        }
                        # inherit negated flag
                        attributes(chull)$negated <- 
                          attributes(gates[[1]])$negated
                        # return consensus gate
                        return(chull)
                      # boolean gates captured here - shared logic
                      } else {
                        return(gates[[1]])
                      }
                    }
                  ),
                  names = names(sample_gates[[1]])
                )
              }
            )
          )
          # create gatingTemplate entries
          do.call(
            "rbind",
            lapply(
              names(group_gate_list[[1]]),
              function(alias) {
                # gate flags
                bool <- FALSE
                quad <- FALSE
                negated <- FALSE
                # get gates from each group
                gates <- structure(
                  lapply(
                    group_gate_list,
                    function(gate_list) {
                      gate <- gate_list[[alias]]
                      # negated flag
                      negated <<- attributes(gate)$negated
                      # booleanFilter flag
                      if(inherits(gate, "booleanFilter")) {
                        bool <<- TRUE
                        return(
                          gate
                        )
                      # quadGates aren't wrapped in filters()
                      } else if(inherits(gate, "quadGate")) {
                        quad <<- TRUE
                        return(
                          gate_list[[alias]]
                        )
                      } else {
                        return(
                          filters(
                            gate_list[alias]
                          )
                        )
                      }
                    }
                  ),
                  names = names(group_gate_list)
                )
                # create gatingTemplate entry
                gt_entry <- as.data.table(gt_chunk[1, , drop = FALSE])
                # alias
                gate <- gsub("\\|", ",", alias)
                # pop
                if(negated) {
                  gt_entry[, pop := "-"]
                } else {
                  gt_entry[, pop := "+"]
                }
                # alias
                gt_entry[, alias := gate]
                # preprocessing method
                gt_entry[, preprocessing_method := "pp_gate_template"]
                # gating method
                gt_entry[, gating_method := "gate_template"]
                # collapse data - required or need gate for every sample
                gt_entry[, collapseDataForGating := TRUE]
                # groupBy
                if(all(is.na(group_by))) {
                  group_by <- NA
                } else {
                  group_by <- paste0(
                    group_by,
                    collapse = ":"
                  )
                }
                gt_entry[, groupBy := group_by]
                # booleanFilter
                if(bool) {
                  gt_entry[, groupBy := NA]
                  gt_entry[, preprocessing_method := NA]
                  gt_entry[, gating_method := "boolGate"]
                  gt_entry[, gating_args := gates[[1]]@deparse]
                  gt_entry[, collapseDataForGating := NA]
                # other gates handled above
                } else {
                  if(quad) {
                    gt_entry[, pop := "*"]
                  }
                  gt_entry[, gating_args := .argDeparser(
                    list(
                      gate = gates,
                      openCyto.minEvents  = -1
                    )
                  )]
                }
                return(gt_entry)
              }
            )
          )
        }
      ),
      names = names(gt_chunks)
    )
  )

  # write gatingTemplate to CSV file
  if(!is.null(save_as)) {
    fwrite(
      gt_new,
      file = save_as
    )
  }

  # convert to gatingTemplate object - only if its complete
  if(gatingTemplate & all(nodes %in% nodes_full)) {
    gt_new <- gatingTemplate(
      gt_new,
      ...
    )
  }

  # return constructed gatingTemplate
  return(gt_new)
  
}
