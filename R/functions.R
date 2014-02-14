#' generate a partially complete csv template from the existing gating hierarchy 
#' 
#' To ease the process of replicating the existing (usually a manual one) gating schemes, 
#' this function populate an empty gating template with the 'alias', 'pop', 'parent' and 'dims' 
#' columns that exacted from an \code{GatingHierarchy}, and leave the other columns (e.g. `gating_method`) blank.
#' So users can make changes to that template instead of writing from scratch.
#' 
#' @param gh a \code{GatingHierarchy} likely parsed from a xml workspace
#' @return a gating template in \code{data.frame} format that requires further edition after output to csv 
#' @export 
templateGen <- function(gh){
  nodes <- getNodes(gh, order = "tsort")
  ldply(nodes[-1], function(thisNode){
        thisGate <- getGate(gh, thisNode)
        dims <- paste(as.vector(parameters(thisGate)), collapse = ",")
        parent <- getParent(gh, thisNode)
        alias <- basename(thisNode)
        pop <- alias
        c(alias = alias
          , pop = pop
          , parent = parent
          , dims = dims
          , gating_method = NA
          , gating_args = NA
          , collapseDataForGating = NA
          , groupBy = NA
          , preprocessing_method = NA
          , preprocessing_args = NA
          )
      })
  
  
}

.getFullPath <- function(ref_node, df){

  ref_node <- flowWorkspace:::trimWhiteSpace(ref_node)
  #prepend root if start with /
  if(substr(ref_node, 1, 1) == "/")
    ref_node <- paste0("root", ref_node)
  
  tokens <- strsplit(ref_node, "/")[[1]]
  res_path <- NULL
  
#browser()
  df_toSearch <- df
  # start to match the tokens in the path
  while (length(tokens) > 0) {
    #pop the current one
    curToken <- tokens[1]
    tokens <- tokens[-1]
    if(curToken == "root")
      res_path <- c(res_path, "root")
    else{
      toMatch <- gsub("\\+", "\\\\\\+", curToken)
      toMatch <- paste0("^",toMatch,"$")
      ind <- grep(toMatch, df_toSearch[, "alias"])
      if(length(ind) == 0)
        stop("Not able to to find reference to: ", curToken)
      else if(length(ind) > 1)
        stop("Non-unique reference to: ", curToken)
      else{
          # prepend its parent to make it full path
          thisParent <- df_toSearch[ind, "parent"]
          if(thisParent == "root")
            curToken <- paste0("/", curToken)
          else
            curToken <- paste(thisParent, curToken, sep = "/")
          
          #only save the full path of the first token 
          if(is.null(res_path))
            toSave <- curToken
          else
            toSave <- basename(curToken)
          
          res_path <- c(res_path, toSave)
      }
    }
#    browser()
    #subset the data frame by parent
    df_toSearch <- subset(df, parent == curToken)
          
  }
        
  
  res_path <- paste(res_path, collapse = "/")
#  browser()    
  if(res_path == "root")
    res_path <- "root"
  else
    res_path <- sub("root", "", res_path)
  res_path

}

# make sure the alias is unique under one parent
.check_alias <- function(new_df, alias, this_parent) {
  if(grepl("[\\|\\&|\\:|\\/]", alias))
    stop(alias , "contains illegal character: |,&,:,/")
  
  siblings <- subset(new_df, parent == this_parent)[, "alias"]
  toMatch <- gsub("\\+", "\\\\\\+", alias)
  toMatch <- paste0("^",toMatch,"$")
  
  matched_sibs <- grep(toMatch, siblings)
  
  if (length(matched_sibs) >= 2) {
    stop(alias, " is not unique within ", this_parent)
  }
}
#' preprocess the csv template
#' 
#' It expands the definition of gates or construct reference gates when necessary
#' @importFrom data.table fread
.preprocess_csv <- function(x) {
  df <- as.data.frame(fread(x))
  df <- df[, 1:10] #only parse first 10 columns and ignore the rest 
  new_df <- df[0, ]
  
  for (i in 1:nrow(df)) {
    this_row <- df[i, , drop = FALSE]
    
    .check_alias(new_df, this_row[1, "alias"], this_row[1, "parent"])
    
    popName <- this_row[1, "pop"]
    dims <- this_row[1, "dims"]
    gm <- this_row[1, "gating_method"]
    
#    browser()
    #update parent with full path
    this_row[1, "parent"] <- .getFullPath(this_row[1, "parent"], new_df)
    
    
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
#    browser()
    if (grepl(paste0("^", one_pop_pat, "$"), popName)) {
      # A+ no expansion(simply update flowClust gm)
      if (gm == "flowClust") {
        if (dim_count == 1) {
          this_row[1, "gating_method"] <- "flowClust.1d"
        } else {
          this_row[1, "gating_method"] <- "flowClust.2d"
        }
      }
      
      res <- this_row
      new_df <- .addToDf(res, this_row, new_df)
    } else if (grepl(paste("^", two_pop_pat, "$", sep = ""), popName)) {
      # A+/-
      
      if (gm == "flowClust") {
        if (dim_count == 1) {
          this_row[1, "gating_method"] <- "flowClust.1d"
        } else {
          this_row[1, "gating_method"] <- "flowClust.2d"
        }
      }
      # expand to two rows
      message("expanding pop: ", popName, "\n")
      cur_dim <- sub(two_pop_token, "", popName)
      
      new_pops <- paste(cur_dim, c("+", "-"), sep = "")
  
      # create 1d gate
    .check_alias(new_df, new_pops[1], this_row[1, "parent"])
      res_1d <- c(alias = new_pops[1], pop = new_pops[1], parent = this_row[1, "parent"], 
                    dims, this_row[1, "gating_method"], this_row["gating_args"]
                  , this_row["collapseDataForGating"], this_row["groupBy"], this_row["preprocessing_method"], this_row["preprocessing_args"])
      new_df <- .addToDf(res_1d, this_row, new_df)
      # create ref gate
  
      refNode <- file.path(this_row[1, "parent"], new_pops[1])
      .check_alias(new_df, new_pops[2], this_row[1, "parent"])
      res_ref <- c(alias = new_pops[2], pop = new_pops[2], parent = this_row[1, "parent"], 
                        dims, "refGate", refNode, this_row["collapseDataForGating"], this_row["groupBy"], NA,NA)
      new_df <- .addToDf(res_ref, this_row, new_df)
#      browser()
    } else if (grepl(paste("^(", one_pop_pat, "){2}$", sep = ""), popName)) {
      # A+B+
      
      if (gm == "refGate") {
        # no expansion
        res <- this_row
        new_df <- .addToDf(res, this_row, new_df)
      } else {
        
        if (gm == "flowClust") {
          if (dim_count == 2) {
          message("expanding pop: ", popName, "\n")
          
          this_row[1, "gating_method"] <- "flowClust.1d"
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
        new_df <- .addToDf(res_1d, this_row, new_df)
        
        res_ref <- .gen_refGate(split_terms$splitted_terms, this_row, ref_nodes = res_1d[, "alias"], alias = this_row["alias"], new_df = new_df)
        new_df <- .addToDf(res_ref, this_row, new_df)
        
      }
    } else if (grepl(paste0("^(", two_pop_pat, "){2}$"), popName) ||
               grepl(paste0("^", two_pop_pat, one_pop_pat, "$"), popName) ||
               grepl(paste0("^", one_pop_pat, two_pop_pat, "$"), popName)) {
      # A+/-B+/-
      message("expanding pop: ", popName, "\n")
      pop_stat <- paste0("(", two_pop_pat, ")|(", one_pop_pat, ")")
      split_terms <- .splitTerms(pop_pat = pop_stat, two_pop_token, popName)
      if (gm == "refGate") {
        res <- .gen_refGate(split_terms$splitted_terms, this_row = this_row, 
          new_df = new_df)
        new_df <- .addToDf(res, this_row, new_df)
      } else {
        
        if (gm == "flowClust") {
          if (dim_count == 2) {
          this_row[1, "gating_method"] <- "flowClust.1d"
          
          } else {
          stop("dimensions '", dims, "' is not consistent with pop name '", 
            popName, "'")
          }
        }
        
        # create 1d gate for each dim
        res_1d <- .gen_1dgate(split_terms$terms, this_row, one_pop_token, 
          two_pop_token, new_df)
#        browser()
        new_df <- .addToDf(res_1d, this_row, new_df)
        
        res_ref <- .gen_refGate(split_terms$splitted_terms, this_row, ref_nodes = res_1d[, "alias"], new_df = new_df)
        new_df <- .addToDf(res_ref, this_row, new_df)
      }
      
    } else {
      stop("invalid population pattern '", popName, "'")
    }

  }
#  browser()

   new_df[is.na(new_df)] <- ""
   new_df
  
}
.addToDf <- function(res,this_row, new_df){
#  browser()
  if (is.matrix(res)) {
    colnames(res) <- names(this_row)
  } else {
    names(res) <- names(this_row)
  }
  
  res <- as.data.frame(res)
  
  rbind(new_df, res)
}
#' split the population pattern into multiple population names 
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
#' convert to 1d gating based on the population pattern
.gen_1dgate <- function(terms, this_row, one_pop_token, two_pop_token, new_df) {
    
  res <- ldply(terms, function(cur_term) {
    toReplace <- paste("(", two_pop_token, ")|(", one_pop_token, ")", sep = "")
    cur_dim <- sub(toReplace, "", cur_term)
    new_pop_name <- paste(cur_dim, "+", sep = "")
    this_parent <- this_row[1, "parent"]
    .check_alias(new_df, new_pop_name, this_parent)
    
    data.frame(alias = new_pop_name, pop = new_pop_name, parent = this_parent, dims = cur_dim, 
      this_row["gating_method"], this_row["gating_args"], this_row["collapseDataForGating"], this_row["groupBy"], this_row["preprocessing_method"], this_row["preprocessing_args"])
  })
  rownames(res) <- NULL
  res
}
#' generate reference gate based on the splitted population patterns 
.gen_refGate <- function(splitted_terms, this_row, ref_nodes = NULL, alias = NULL, 
  new_df) {
  this_parent <- this_row[1, "parent"]
  
  if (is.null(ref_nodes)) {
    # simply copy ref args from the row
    ref_args <- this_row[1, "gating_args"]
  } else {
    # use the new generated 1d pops to construct ref args
    # prepend the path to ref_nodes 
    ref_nodes <- file.path(this_parent, ref_nodes)
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
    data.frame(alias = cur_alias, pop = new_pop, parent = this_parent, this_row["dims"], 
      method = "refGate", gating_args = ref_args, NA, NA, NA, NA)
  }, SIMPLIFY = FALSE))
}

#' reading flowFrame from csv spreadsheet.It is for internal usage.
#' @param file \code{character} csv file name that contains the fluorescence intensities
#' @param stains \code{character} a vector specifying the stains for channels. Default is \code{NA}
#' @return a \code{flowFrame} object
.read.FCS.csv <- function(file, stains = NA) {
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
#' reading flowSet from multiple csv spreadsheets.It is for internal usage.
#' @param files \code{character} csv file names
#' @param ... arguments passed to \link{read.FCS.csv}
#' @return a \code{flowSet} object
.read.flowSet.csv <- function(files, ...) {
  fs <- flowSet(lapply(files, .read.FCS.csv, ...))
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
#' @param ecol \code{numeric} to be documented
#' @param elty \code{numeric} to be documented
#' @param quantile the contour level of the ellipse. See details.
#' @param npoints the number of points on the ellipse
#' @param subset the dimensions of the mixture component to return
#' @param \code{...} additional parameters
#' @return matrix containing the points of the ellipse from the flowClust contour
#' @importFrom flowClust rbox
.getEllipse <- function(filter = NULL, include = seq_len(filter@K), ecol = 1, elty = 1, 
  quantile = NULL, npoints = 501, subset = c(1, 2),...) {
  
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
  
  #Does trans exist in the extra parameter list?
  #If not, set it to true by default
  ellipsis<-as.environment(list(...))
  if(exists("trans",envir=ellipsis)){
    trans<-get("trans",ellipsis)
  }else{
    trans<-1
  }

  #Test for trans==0 when lambda is defined to get around the off 
  #by one bug due to the reverse box-cox transformation
  if ((length(filter@lambda) > 0)&&trans==0) {
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
    
    if ((length(lambda) > 0)&trans==1) {
      res <- rbox(flowClust:::.ellipsePoints(a = l1[i], b = l2[i], alpha = angle, 
        loc = filter@mu[i, subset], n = npoints), lambda[i])
    } else {
      res <- flowClust:::.ellipsePoints(a = l1[i], b = l2[i], alpha = angle, 
        loc = filter@mu[i, subset], n = npoints)
    }
  }
  res
}



#' Removes any observation from the given flowFrame object that has values
#' outside the given range for the specified channels
#'
#' The minimum/maximum values are ignored if \code{NULL}.
#'
#' @param flow_frame a \code{flowFrame} object
#' @param min a numeric vector that sets the lower bounds for data filtering
#' @param max a numeric vector that sets the upper bounds for data filtering
#' @param channels \code{character} specifying which channel to operate on
#' @return a \code{flowFrame} object
#' @examples
#' \dontrun{
#'  library(flowClust)
#'  data(rituximab)
#'  # Consider the range of values for FSC.H and SSC.H
#'  summary(rituximab)
#' 
#'  # Truncates any observations with FSC.H outside [100, 950]
#'  rituximab2 <- .truncate_flowframe(rituximab, channels = "FSC.H", min = 100, max = 950)
#'  summary(rituximab2)
#'  # Next, truncates any observations with SSC.H outside [50, 1000]
#'  rituximab3 <- .truncate_flowframe(rituximab2, channels = "SSC.H", min = 50, max = 1000)
#'  summary(rituximab3)
#'
#'  # Instead, truncates both channels at the same time
#'  rituximab4 <- .truncate_flowframe(rituximab, channels = c("FSC.H", "SSC.H"),
#'  min = c(100, 50), max = c(950, 1000))
#'  summary(rituximab4)
#' }
#' 
.truncate_flowframe <- function(flow_frame, channels, min = NULL, max = NULL) {
  channels <- as.character(channels)
  num_channels <- length(channels)

  # For comparison purposes, we update the min and max values to -Inf and Inf,
  # respectively, if NULL.
  if (is.null(min)) {
    min <- rep(-Inf, num_channels)
  }
  if (is.null(max)) {
    max <- rep(Inf, num_channels)
  }

  if (!(num_channels == length(min) && num_channels == length(max))) {
    stop("The lengths of 'min' and 'max' must match the number of 'channels' given.")
  }
 
  gate_coordinates <- lapply(seq_len(num_channels), function(i) {
    c(min[i], max[i])
  })
  names(gate_coordinates) <- channels

  truncate_filter <- rectangleGate(gate_coordinates)
  Subset(flow_frame, truncate_filter)
}

#' Removes any observation from each flowFrame object within the flowSet that has
#' values outside the given range for the specified channels
#'
#' The minimum/maximum values are ignored if \code{NULL}.
#'
#' @param flow_set a \code{flowSet} object
#' @inheritParams .truncate_flowframe
#' @return a \code{flowSet} object
.truncate_flowset <- function(flow_set, channels, min = NULL, max = NULL) {
  channels <- as.character(channels)
  num_channels <- length(channels)

  # For comparison purposes, we update the min and max values to -Inf and Inf,
  # respectively, if NULL.
  if (is.null(min)) {
    min <- rep(-Inf, num_channels)
  }
  if (is.null(max)) {
    max <- rep(Inf, num_channels)
  }

  if (!(num_channels == length(min) && num_channels == length(max))) {
    stop("The lengths of 'min' and 'max' must match the number of 'channels' given.")
  }
 
  gate_coordinates <- lapply(seq_len(num_channels), function(i) {
    c(min[i], max[i])
  })
  names(gate_coordinates) <- channels

  truncate_filter <- rectangleGate(gate_coordinates)
  Subset(flow_set, truncate_filter)
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
.quantile_flowClust <- function(p, object, interval, ...) {
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
.quadGate2rectangleGates <- function(quad_gate, markers, channels, quadrants = 1:4) {
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
#' .polyfunction_nodes(c('IFNg', 'IL2', 'TNFa', 'GzB', 'CD57'))
.polyfunction_nodes <- function(markers) {
  
  markers <- paste0(markers, "+")
  num_markers <- length(markers)
  and_list <- as.vector(permutations(n = 1, r = num_markers - 1, c("&"), repeats.allowed = TRUE))
  isnot_list <- permutations(n = 2, r = num_markers, c("!", ""), repeats.allowed = TRUE)
  apply(isnot_list, 1, function(isnot_row) {
    isnot_row[-1] <- paste0(and_list, isnot_row[-1])
    paste(paste0(isnot_row, markers), collapse = "")
  })
} 

#' Finds values of a vector in an interval
#'
#' @param x numeric vector
#' @param interval numeric vector of length 2
#' @return numeric vector containing the values of \code{x} between
#' \code{interval}. If no values are found, \code{NA} is returned.
#' @examples
#' z <- seq.int(1, 9, by = 2)
#' .between_interval(z, interval = c(2, 8))
#' @export 
.between_interval <- function(x, interval) {
  x <- x[findInterval(x, interval) == 1]
  if (length(x) == 0) {
    x <- NA
  }
  x
}

#' Constructs a 2x2 rotation matrix for a given angle
#' @param theta \code{numeric} the degree of rotation that ensures the angle between the x-axis and the eigenvector is between 0 and pi
.rotation_matrix <- function(theta) {
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
}
