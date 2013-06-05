.update_ref_node <- function(ref_node, new_df, this_parent) {
  ind <- match(ref_node, new_df[, "alias"])
  if (!is.na(ind)) {
    if (this_parent == "root") {
      stop("Not able to to find unique reference for ", ref_node)
    } else {
      # assuming this_parent is uniquely identifiable
      ref_node <- paste(this_parent, ref_node, sep = "/")
    }
  }
  ref_node
}

# make sure the alias is unique under one parent
.check_alias <- function(new_df, alias, this_parent) {
  siblings <- subset(new_df, parent == this_parent)[, "alias"]
  matched_sibs <- match(alias, siblings)
  
  if (length(matched_sibs) >= 2) {
    stop(alias, " is not unique within ", this_parent)
  }
}

.preprocess_csv <- function(x) {
  df <- read.csv(x, as.is = T)
  new_df <- df[0, ]
  for (i in 1:nrow(df)) {
    this_row <- df[i, , drop = FALSE]
    
    popName <- this_row[1, "pop"]
    dims <- this_row[1, "dims"]
    gm <- this_row[1, "method"]
    
    if (!grepl("[+-]", popName)) {
      popName <- paste0(popName, "+")
    }
    
    dim_count <- length(strsplit(split = ",", dims)[[1]])
    if (!dim_count %in% c(1, 2)) {
      if (!(dim_count == 0 && gm %in% c("polyFunctions", "boolGate", "refGate"))) {
        stop(popName, " has invalid number of dimensions: ", dim_count)
      }
    }

    ######################################################################
    # validity check for pop and determine if expanded to multiple pops
    # three different types of expansion
    # 1.expand quadGate defined by refGate
    # 2.expand quadGate defined by 1d gating method
    # 3.expand two pops defined by +/-
    ######################################################################

    one_pop_token <- "[\\+-]"
    pop_name_pat <- "[^\\+-]+"
    one_pop_pat <- paste(pop_name_pat, one_pop_token, sep = "")
    
    two_pop_token <- "(\\+/-)|(-/\\+)"
    two_pop_pat <- paste(pop_name_pat, "(", two_pop_token, ")", sep = "")
    
    if (grepl(paste0("^", one_pop_pat, "$"), popName)) {
      # A+ no expansion(simply update flowClust gm)
      if (gm == "flowClust") {
        if (dim_count == 1) {
          this_row[1, "method"] <- "flowClust.1d"
        } else {
          this_row[1, "method"] <- "flowClust.2d"
        }
      }
      
      res <- this_row
      
    } else if (grepl(paste("^", two_pop_pat, "$", sep = ""), popName)) {
      # A+/-
      
      if (gm == "flowClust") {
        if (dim_count == 1) {
          this_row[1, "method"] <- "flowClust.1d"
        } else {
          this_row[1, "method"] <- "flowClust.2d"
        }
      }
      # expand to two rows
      cat("expanding pop: ", popName, "\n")
      cur_dim <- sub(two_pop_token, "", popName)
      new_pops <- paste(cur_dim, c("+", "-"), sep = "")
      
      # create 1d gate
      res_1d <- c(alias = new_pops[1], pop = new_pops[1], parent = this_row[1, "parent"], 
                    dims, this_row[1, "method"], this_row["args"])
      # create ref gate
      res_ref <- c(alias = new_pops[2], pop = new_pops[2], parent = this_row[1, "parent"], 
                        dims, "refGate", file.path(this_row[1, "parent"],new_pops[1]))
      res <- rbind(res_1d, res_ref)

    } else if (grepl(paste("^(", one_pop_pat, "){2}$", sep = ""), popName)) {
      # A+B+
      
      if (gm == "refGate") {
        # no expansion
        res <- this_row
        
      } else {
        
        if (gm == "flowClust") {
          if (dim_count == 2) {
          cat("expanding pop: ", popName, "\n")
          
          this_row[1, "method"] <- "flowClust.1d"
          } else {
          stop("dimensions '", dims, "' is not consistent with pop name '", 
            popName, "'")
          }
        }
        # needs to be split into two 1d gates and one refgate
        split_terms <- .splitTerms(pop_pat = one_pop_pat, two_pop_token, 
          popName)
        
        # create 1d gate for each dim
        res_1d <- .gen_1dgate(split_terms$terms, this_row, one_pop_token, 
          two_pop_token, new_df)
        
        res_ref <- .gen_refGate(split_terms$splitted_terms, this_row, ref_nodes = res_1d[, 
          "alias"], alias = this_row["alias"], new_df = new_df)
        res <- rbind(res_1d, res_ref)
        
      }
    } else if (grepl(paste0("^(", two_pop_pat, "){2}$"), popName) ||
               grepl(paste0("^", two_pop_pat, one_pop_pat, "$"), popName) ||
               grepl(paste0("^", one_pop_pat, two_pop_pat, "$"), popName)) {
      # A+/-B+/-
      cat("expanding pop: ", popName, "\n")
      pop_stat <- paste0("(", two_pop_pat, ")|(", one_pop_pat, ")")
      split_terms <- .splitTerms(pop_pat = pop_stat, two_pop_token, popName)
      if (gm == "refGate") {
        res <- .gen_refGate(split_terms$splitted_terms, this_row = this_row, 
          new_df = new_df)
      } else {
        
        if (gm == "flowClust") {
          if (dim_count == 2) {
          this_row[1, "method"] <- "flowClust.1d"
          
          } else {
          stop("dimensions '", dims, "' is not consistent with pop name '", 
            popName, "'")
          }
        }
        
        # create 1d gate for each dim
        res_1d <- .gen_1dgate(split_terms$terms, this_row, one_pop_token, 
          two_pop_token, new_df)
        res_ref <- .gen_refGate(split_terms$splitted_terms, this_row, ref_nodes = res_1d[, 
          "alias"], new_df = new_df)
        res <- rbind(res_1d, res_ref)
      }
      
    } else {
      stop("invalid population pattern '", popName, "'")
    }
    
    if (is.matrix(res)) {
      colnames(res) <- names(this_row)
    } else {
      names(res) <- names(this_row)
    }
    
    new_df <- rbind(new_df, res)
  }
  new_df
  
}

.splitTerms <- function(pop_pat, two_pop_token, popName) {
  term_pos <- gregexpr(pop_pat, popName)[[1]]
  x_term <- substr(popName, 1, term_pos[2] - 1)
  y_term <- substr(popName, term_pos[2], nchar(popName))
  terms <- c(x_term, y_term)
  
  splitted_terms <- lapply(terms, function(cur_term) {
    token_pos <- gregexpr(two_pop_token, cur_term)[[1]]
    if (token_pos > 0) {
      dim_name <- substr(cur_term, 1, token_pos - 1)
      splitted_term <- paste(dim_name, c("+", "-"), sep = "")
      
    } else {
      splitted_term <- cur_term
    }
    splitted_term
  })
  list(terms = terms, splitted_terms = splitted_terms)
}

.gen_1dgate <- function(terms, this_row, one_pop_token, two_pop_token, new_df) {
  res <- do.call(rbind, lapply(terms, function(cur_term) {
    toReplace <- paste("(", two_pop_token, ")|(", one_pop_token, ")", sep = "")
    cur_dim <- sub(toReplace, "", cur_term)
    new_pop_name <- paste(cur_dim, "+", sep = "")
    this_parent <- this_row[1, "parent"]
    .check_alias(new_df, new_pop_name, this_parent)
    
    c(alias = new_pop_name, pop = new_pop_name, parent = this_parent, dims = cur_dim, 
      this_row["method"], this_row["args"])
  }))
  rownames(res) <- NULL
  res
}
.gen_refGate <- function(splitted_terms, this_row, ref_nodes = NULL, alias = NULL, 
  new_df) {
  this_parent <- this_row[1, "parent"]
  
  if (is.null(ref_nodes)) {
    # simply copy ref args from the row
    ref_args <- this_row[1, "args"]
  } else {
    # use the new generated 1d pops to construct ref args
    # prepend the path to ref_nodes if needed
    ref_nodes <- lapply(ref_nodes, .update_ref_node, new_df = new_df, this_parent = this_parent)
    ref_args <- paste(ref_nodes, collapse = ":")
  }
  
  new_pops <- as.character(outer(splitted_terms[[1]], splitted_terms[[2]], paste, 
    sep = ""))
  
  if (is.null(alias)) {
    alias <- new_pops
  }
  # create ref gate for each new_pop )
  do.call(rbind, mapply(new_pops, alias, FUN = function(new_pop, cur_alias) {
    .check_alias(new_df, cur_alias, this_parent)
    c(alias = cur_alias, pop = new_pop, parent = this_parent, this_row["dims"], 
      method = "refGate", args = ref_args)
  }, SIMPLIFY = FALSE))
}

read.FCS.csv <- function(file, stains = NA) {
  mat <- as.matrix(read.csv(file, check.names = FALSE))
  
  fr <- new("flowFrame", exprs = mat)
  
  pd <- pData(parameters(fr))
  pd$desc <- as.character(pd$desc)
  pd$name <- as.character(pd$name)

  ## update the desc with marker name
  if (!is.na(stains)) {
    ind <- match(names(stains), pd$name)
    pd[ind, ]$desc <- as.character(stains)

    ## update SSC and FSC description with NA
    ind <- grepl("[F|S]SC", pd$desc)
    pd[ind, ]$desc <- NA
  } else pd$desc <- NA
  
  ## update minRange with -111 for proper display of the data
  pd$minRange[pd$minRange < (-111)] <- -111
  pData(parameters(fr)) <- pd
  fr
}

read.flowSet.csv <- function(files, ...) {
  fs <- flowSet(lapply(files, read.FCS.csv, ...))
  sampleNames(fs) <- basename(files)
  fs
}

#' Creates a matrix of points on an ellipse from a fitted flowClust model.
#'
#' The ellipse is constructed from a contour from the fitted flowClust model.
#'
#' By default, the contour level is extracted from the \code{ruleOutliers} slot
#' in the \code{filter}. The user can override this level by specifying value
#' between 0 and 1 in \code{quantile}.
#'
#' @param filter object containing the fitted \code{flowClust} model.
#' @param include the mixture component in the fitted \code{flowClust} model for
#' which the contour (ellipse) is returned
#' @param ecol TODO
#' @param elty TODO
#' @param quantile the contour level of the ellipse. See details.
#' @param npoints the number of points on the ellipse
#' @param subset the dimensions of the mixture component to return
#' @return matrix containing the points of the ellipse from the flowClust contour
.getEllipse <- function(filter = NULL, include = seq_len(filter@K), ecol = 1, elty = 1, 
  quantile = NULL, npoints = 501, subset = c(1, 2)) {
  
  # Sets the quantile of the ellipse.
  if (is.null(quantile)) {
    quantile <- filter@ruleOutliers[2]
  } else {
    if (!is.numeric(quantile) || quantile < 0 || quantile > 1) {
      stop("The 'quantile' must be a numeric value between 0 and 1.")
    }
  }
  
  # py is the degrees of freedom?
  py <- 2
  ecol <- matrix(ecol, length(include))
  elty <- matrix(elty, length(include))
  
  if (all(filter@nu != Inf)) {
    if (filter@ruleOutliers[1] == 0) {
      # 0 means quantile
      cc <- py * qf(p = quantile, py, filter@nu)
    } else {
      # 1 means u.cutoff
      cc <- ((filter@nu + py)/quantile - filter@nu)
    }
  } else {
    cc <- qchisq(p = quantile, py)
  }
  
  j <- 0
  if (length(filter@lambda) > 0) {
    lambda <- rep(filter@lambda, length.out = filter@K)
  } else {
    lambda <- numeric(0)
  }
  cc <- rep(cc, length.out = filter@K)
  for (i in include) {
    eigenPair <- eigen(filter@sigma[i, subset, subset])
    l1 <- sqrt(eigenPair$values[1]) * sqrt(cc)
    l2 <- sqrt(eigenPair$values[2]) * sqrt(cc)
    angle <- atan(eigenPair$vectors[2, 1]/eigenPair$vectors[1, 1]) * 180/pi
    
    if (length(lambda) > 0) {
      res <- rbox(flowClust:::.ellipsePoints(a = l1[i], b = l2[i], alpha = angle, 
        loc = filter@mu[i, subset], n = npoints), lambda[i])
    } else {
      res <- flowClust:::.ellipsePoints(a = l1[i], b = l2[i], alpha = angle, 
        loc = filter@mu[i, subset], n = npoints)
    }
  }
  res
}

.flowParamMatch <- function(pd, name, fix = FALSE, partial = FALSE) {
  # try to compelete word match by following with a space or the end of string
  if (partial) 
    pname <- name else pname <- paste0(name, "([ ]|$)")
  
  if (fix) {
    ind <- which(toupper(pd$name) %in% toupper(name))
  } else {
    ind <- which(grepl(pname, pd$name, ignore.case = T))
  }
  
  if (length(ind) == 0) {
    # try marker name
    ind <- which(unlist(lapply(pd$des, function(x) {
      # split by white space and then match each individual string
      if (fix) {
        any(unlist(lapply(strsplit(x, " "), function(y) toupper(y) %in% toupper(name))))
      } else {
        grepl(pattern = pname, x, ignore.case = T)
      }
    })))
  }
  ind
}

getChannelMarker <- function(frm, name, ...) {
  # try stain name
  pd <- pData(parameters(frm))

  # try complete match first
  ind <- .flowParamMatch(pd, name, ...)

  if (length(ind) > 1) {
    stop("multiple markers matched: ", name)
  }
  
  if (length(ind) == 0) {
    # if no match then give a second try to patial match
    ind <- .flowParamMatch(pd, name, partial = TRUE, ...)
    if (length(ind) == 0) 
      stop("can't find ", name) else if (length(ind) > 1) 
      stop("multiple markers matched: ", name) else warning(name, " is partially matched with ", pd[ind, c("name", "desc")])
  }
  
  pd[ind, c("name", "desc")]
}

#' For the given workflow, we look up the given markers and return the
#' corresponding channels.
#'
#' @param flow_frame object of type \code{flowFrame}
#' @param markers the markers from which we obtain the corresponding channel names
#' @return vector of channel names
markers2channels <- function(flow_frame, markers) {
  # First, we build a lookup table for the channels and markers.
  channel_markers <- lapply(colnames(flow_frame), function(channel) {
    marker_name_desc <- getChannelMarker(flow_frame, channel)
    marker <- with(marker_name_desc, ifelse(is.na(desc), name, desc))
    cbind(channel, marker = unname(marker))
  })
  channel_markers <- data.frame(do.call(rbind, channel_markers), stringsAsFactors = FALSE)
  
  # Now, we query the channels for the specified markers.
  channels <- sapply(markers, function(marker) {
    channel_markers$channel[grepl(marker, channel_markers$marker)]
  })
  
  as.vector(channels)
}

#' For the given flow frame, we look up the given markers and return the
#' corresponding channels.
#'
#' @param flow_frame object of type \code{flowFrame}
#' @param channels the channels from which we obtain the corresponding markers
#' @return vector of markers
channels2markers <- function(flow_frame, channels) {
  markers <- sapply(channels, function(channel) {
    marker <- getChannelMarker(flow_frame, channel)
    with(marker, ifelse(is.na(desc), name, desc))
  })
  unname(markers)
}

#' Removes any observation from the given flow frame that has values less than
#' (greater than) the minimum (maxim) value.
#'
#' The minimum/maximum values are ignored if \code{NULL}.
#'
#' @param flow_frame a \code{flowFrame} object
#' @param min a numeric value that sets the lower boundary for data filtering
#' @param max a numeric value that sets the upper boundary for data filtering
#' @return a \code{flowFrame} object
truncate_flowframe <- function(flow_frame, channel, min = NULL, max = NULL) {

  channel <- as.character(channel)
  if (length(channel) != 1) {
    stop("Only one 'channel' may be truncated.")
  }
  
  x_channel <- exprs(flow_frame)[, channel]

  # For comparison purposes, we update the min and max values to -Inf and Inf,
  # respectively, if one is NULL.
  min <- ifelse(is.null(min), -Inf, min)
  max <- ifelse(is.null(max), Inf, max)
  
  # Removes any observation that has an observation outside of the min and max
  # values specified.
  exprs_truncated <- exprs(flow_frame)[min < x_channel & x_channel < max, ]
  if (is.vector(exprs_truncated)) {
    exprs_truncated <- t(as.matrix(exprs_truncated))
  }
  exprs(flow_frame) <- exprs_truncated
  
  flow_frame
}

#' Computes the quantile from flowClust for a given vector of probabilties
#'
#' We estimate the quantile from a \code{flowClust} fit with a combination of
#' numerical integration and a root-finding method. We are effectively
#' estimating the cumulative distribution function (CDF) of the mixture density
#' estimated by \code{flowClust}.
#'
#' Because we are using numerical methods, we also need an \code{interval} of
#' values in which we will attempt to find the specified quantile.
#'
#' @param p vector of probabilities
#' @param object an object containing the \code{flowClust} fit
#' @param interval a vector of length 2 containing the end-points of the interval
#' of values to find the quantile
#' @param ... Additional arguments that are passed to \code{uniroot} to find the
#' quantile.
#' @return the quantile corresponding to the specified probabilities
quantile_flowClust <- function(p, object, interval, ...) {
  cdf_target <- function(x, p, object) {
    cdf_values <- sapply(seq_len(object@K), function(k) {
      nu <- ifelse(length(object@nu) == 1, object@nu, object@nu[k])
      lambda <- ifelse(length(object@lambda) == 1, object@lambda, object@lambda[k])
      
      # TODO: Incorporate the Box-Cox transformation (i.e., box(qt(...), lambda =
      # lambda)) into quantile The case of 'lambda = 1' may be not be trivial -- this
      # case is largely ignored in flowClust.
      pt((x - object@mu[k])/sqrt(object@sigma[k]), df = nu)
    })
    weighted.mean(cdf_values, w = object@w) - p
  }
  
  uniroot(cdf_target, interval = interval, p = p, object = object, ...)$root
}

#' Extracts the quadrants of a quadGate as a list of rectangleGates
#'
#' The quadrants are numbered in a clockwise manner with the top-left quadrant
#' numbered 1, the top-right quadrant numbered 2, and so on.
#'
#' @param quad_gate a \code{quadGate} object
#' @param markers character vector of the marker names for the x- and y-axes
#' @param channels character vector of the channel names for the x- and y-axes
#' @param quadrants a vector indicating the quadrants to extract
#' @return a \code{filters} object containing a list of the rectangle gates
quadGate2rectangleGates <- function(quad_gate, markers, channels, quadrants = 1:4) {
  x_gate <- quad_gate@boundary[1]
  y_gate <- quad_gate@boundary[2]
  
  gates_list <- list()
  
  # Top-left quadrant
  gates <- list(c(-Inf, x_gate), c(y_gate, Inf))
  names(gates) <- channels
  gates_list[[paste0(markers[1], "-", markers[2], "+")]] <- rectangleGate(gates)
  
  # Top-right quadrant
  gates <- list(c(x_gate, Inf), c(y_gate, Inf))
  names(gates) <- channels
  gates_list[[paste0(markers[1], "+", markers[2], "+")]] <- rectangleGate(gates)
  
  # Lower-right quadrant
  gates <- list(c(x_gate, Inf), c(-Inf, y_gate))
  names(gates) <- channels
  gates_list[[paste0(markers[1], "+", markers[2], "-")]] <- rectangleGate(gates)
  
  # Lower-left quadrant
  gates <- list(c(-Inf, x_gate), c(-Inf, y_gate))
  names(gates) <- channels
  gates_list[[paste0(markers[1], "-", markers[2], "-")]] <- rectangleGate(gates)
  
  filters(gates_list[quadrants])
}

#' Constructs a vector of all the combinations of A & B & C
#'
#' The \code{permutations} function is from the \code{gregmisc} package on CRAN.
#' @param markers character vector of marker names
#' @return vector containing all combinations of the markers
#' @examples
#' polyfunction_nodes(c('IFNg', 'IL2', 'TNFa', 'GzB', 'CD57'))
polyfunction_nodes <- function(markers) {
  
  markers <- paste0(markers, "+")
  num_markers <- length(markers)
  and_list <- as.vector(permutations(n = 1, r = num_markers - 1, c("&"), repeats = TRUE))
  isnot_list <- permutations(n = 2, r = num_markers, c("!", ""), repeats = TRUE)
  apply(isnot_list, 1, function(isnot_row) {
    isnot_row[-1] <- paste0(and_list, isnot_row[-1])
    paste(paste0(isnot_row, markers), collapse = "")
  })
} 
