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

#' prepend all ancester nodes to construct the full path for the given node 
#' 
#' @param ref_node \code{character} the node to work with
#' @param dt \code{data.table} that stores the current preprocessed template
.getFullPath <- function(ref_node, dt){
  
  ref_node <- flowWorkspace:::trimWhiteSpace(ref_node)
  #prepend root if start with /
  if(substr(ref_node, 1, 1) == "/")
    ref_node <- paste0("root", ref_node)
  
  tokens <- strsplit(ref_node, "/")[[1]]
  res_path <- NULL
  
#browser()
  dt_toSearch <- dt
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
      ind <- grep(toMatch, dt_toSearch[, alias])
      if(length(ind) == 0)
        stop("Not able to to find reference to: ", curToken)
      else if(length(ind) > 1)
        stop("Non-unique reference to: ", curToken)
      else{
        # prepend its parent to make it full path
        thisParent <- dt_toSearch[ind, parent]
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
    dt_toSearch <- dt[parent == curToken, ]
    
  }
  
  
  res_path <- paste(res_path, collapse = "/")
#  browser()    
  if(res_path == "root")
    res_path <- "root"
  else
    res_path <- sub("root", "", res_path)
  res_path
  
}

#' validity check of alias column
.validity_check_alias <- function(alias){
  if(grepl("[\\|\\&|\\:|\\/]", alias))
    stop(alias , "contains illegal character: |,&,:,/")
  
}
#' check if alias is unique within template
#' 
#' @param new_dt \code{data.tab;e} contains the preprocessed rows
#' @param alias \code{character} the alias to be checked
#' @param this_parent \code{character} the full gating path of the parent node  
#' @return NULL when pass the check, otherwise, throw the error
.unique_check_alias <- function(new_dt, alias, this_parent) {
  
  siblings <- new_dt[parent == this_parent, alias]
  toMatch <- gsub("\\+", "\\\\\\+", alias)
  toMatch <- paste0("^",toMatch,"$")
  
  matched_sibs <- grep(toMatch, siblings)
  
  if (length(matched_sibs) > 0) {
    stop(alias, " is not unique within ", this_parent)
  }
}
#' preprocess the csv template
#' 
#' It parses the data table sequentially and does the valdidity checking and expansion row by row.
#'  
#' It expands the definition of gates or construct reference gates when necessary
#' @param dt \code{data.table} loaded directly from csv gating template
#' @return a preprocessed(expanded when applicable) \code{data.frame}
#' @import data.table
.preprocess_csv <- function(dt) {
  
  #only parse these columns(other columns may be used by user for other purpose e.g. comments)
  dt <- dt[, list(alias
          , pop
          , parent
          , dims
          , gating_method
          , gating_args
          , collapseDataForGating
          , groupBy
          , preprocessing_method
          , preprocessing_args
      )
  ]  
  
  new_dt <- dt[0, ]
  
  for (i in 1:nrow(dt)) {
    this_row <- dt[i, , drop = FALSE]
    .validity_check_alias(this_row[, alias])
    #update parent with full path
    this_row[, parent := .getFullPath(this_row[, parent], new_dt)]
    #preprocess the current row
    res <- .preprocess_row(this_row)
    
    #check the uniqueness of alias within res
    aliasVec <- res[,alias]
    dupVec <- duplicated(aliasVec)
    if(any(dupVec))
      stop(aliasVec[dupVec], "is duplicated!")
    
    #validity check of alias within the context of new_df
    apply(res, 1, function(row){
          thisAlias <- row[["alias"]]
          .validity_check_alias(thisAlias)
          .unique_check_alias(new_dt, thisAlias, row[["parent"]])
        })
    
    new_dt <- rbindlist(list(new_dt, res)) 
  }
#  browser()
  
  new_dt[is.na(new_dt)] <- ""
  new_dt
  
}    

#' The actual preprocessing logic for each entry of template
#' 
#' Here are the major preprocessing tasks:
#' 1. validity check for 'alias' (special character and uniqueness check)
#' 2. validity check for the numer of parameters('dim' column)
#' 3. dispatch 'flowClust' method to either 'flowClust.1d' or 'flowClust.2d' based on the 'pop' and 'dims' columns
#' 4. expand the single row to multiple rows when applicable. There are basically two types of expansion:
#'      4.1. expand to two 1d gates and one rectangelGate when 'pop' name is defined as quadrant pattern (e.g. "A+B+")  and 'gating_method' is not "refGate"
#'      4.2. expand to multiple gates when 'pop' is defined with '+/-' (e.g. "A+/-" or "A+/-B+/-") 
#' 
#' @param this_row a single-row \code{data.table}
#' @return  \code{data.table}
.preprocess_row <- function(this_row){
  #make sure it doesn't tamper the input
  this_row <- copy(this_row)
  
  popName <- this_row[1, pop]
  dims <- this_row[1, dims]
  gm <- this_row[1, gating_method]
  
  if (!grepl("[+-]$", popName)) {
    popName <- paste0(popName, "+")
    this_row[, pop := popName]
  }
  
  dim_count <- length(strsplit(split = ",", dims)[[1]])
  if (!dim_count %in% c(1, 2)) {
    if (!(dim_count == 0 && gm %in% c("polyFunctions", "boolGate", "refGate"))) {
      stop(popName, " has invalid number of dimensions: ", dim_count)
    }
  }
  
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
        this_row[1, gating_method := "flowClust.1d"] 
      } else {
        this_row[1, gating_method := "flowClust.2d"] 
      }
    }
    
    res <- this_row
    
  } else if (grepl(paste("^", two_pop_pat, "$", sep = ""), popName)) {
    # A+/-
    
    if (gm == "flowClust") {
      if (dim_count == 1) {
        this_row[1, gating_method := "flowClust.1d"]
      } else {
        this_row[1, gating_method := "flowClust.2d"]
      }
    }
    # expand to two rows
    message("expanding pop: ", popName, "\n")
    cur_dim <- sub(two_pop_token, "", popName)
    
    new_pops <- paste(cur_dim, c("+", "-"), sep = "")
    
    
    res <- rbindlist(list(this_row, this_row))
    
    # update 1d gate
    res[1, alias := new_pops[1]]
    res[1, pop := new_pops[1]]
    
    
    # update ref gate
    refNode <- file.path(this_row[, parent], new_pops[1])
    res[2, alias := new_pops[2]]
    res[2, pop := new_pops[2]]
    res[2, gating_method := "refGate"]
    res[2, gating_args := refNode]
    res[2, preprocessing_method := ""]
    res[2, preprocessing_args := ""]
    
  } else if (grepl(paste("^(", one_pop_pat, "){2}$", sep = ""), popName)) {
    # A+B+
    
    if (gm == "refGate") {
      # no expansion
      res <- this_row
      
    } else {
      
      if (gm == "flowClust") {
        if (dim_count == 2) {
          message("expanding pop: ", popName, "\n")
          
          this_row[1, gating_method := "flowClust.1d"]
        } else {
          stop("dimensions '", dims, "' is not consistent with pop name '", 
              popName, "'")
        }
      }
      # needs to be split into two 1d gates and one refgate
      split_terms <- .splitTerms(pop_pat = one_pop_pat, two_pop_token,popName)
      # create 1d gate for each dim
      res_1d <- .gen_1dgate(split_terms$terms, this_row, one_pop_token, two_pop_token)
      res_ref <- .gen_refGate(split_terms$splitted_terms, this_row, ref_nodes = res_1d[, alias], alias = this_row[, alias])
      res <- rbindlist(list(res_1d, res_ref))
    }
  } else if (grepl(paste0("^(", two_pop_pat, "){2}$"), popName) ||
      grepl(paste0("^", two_pop_pat, one_pop_pat, "$"), popName) ||
      grepl(paste0("^", one_pop_pat, two_pop_pat, "$"), popName)) {
    # A+/-B+/-
    message("expanding pop: ", popName, "\n")
    two_or_one_pop_pat <- paste0("(", two_pop_pat, ")|(", one_pop_pat, ")")
    split_terms <- .splitTerms(pop_pat = two_or_one_pop_pat, two_pop_token, popName)
    if (gm == "refGate") {
      res <- .gen_refGate(split_terms$splitted_terms, this_row = this_row)
      
    } else {
      
      if (gm == "flowClust") {
        if (dim_count == 2) {
          this_row[1, gating_method := "flowClust.1d"]
          
        } else {
          stop("dimensions '", dims, "' is not consistent with pop name '", 
              popName, "'")
        }
      }
      
      # create 1d gate for each dim
      res_1d <- .gen_1dgate(split_terms$terms, this_row, one_pop_token, two_pop_token)
      res_ref <- .gen_refGate(split_terms$splitted_terms, this_row, ref_nodes = res_1d[, alias])
      res <- rbindlist(list(res_1d, res_ref))
    }
    
  } else 
    stop("invalid population pattern '", popName, "'")
  res
}


#' split the population pattern into multiple population names
#' 
#' currently only works for A+/-B+/- .
#' TODO:  support A+/-  
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
#' convert to 1d gating based on the population pattern (A+/-B+/-)
.gen_1dgate <- function(terms, this_row, one_pop_token, two_pop_token) {
  
  res <- ldply(terms, function(cur_term) {
        toReplace <- paste("(", two_pop_token, ")|(", one_pop_token, ")", sep = "")
        cur_dim <- sub(toReplace, "", cur_term)
        new_pop_name <- paste(cur_dim, "+", sep = "")
        this_parent <- this_row[, parent]
        .validity_check_alias(new_pop_name)
        
        data.table(alias = new_pop_name
            , pop = new_pop_name
            , parent = this_parent
            , dims = cur_dim 
            , gating_method = this_row[, gating_method]
            , gating_args = this_row[, gating_args]
            , collapseDataForGating = this_row[, collapseDataForGating]
            , groupBy = this_row[, groupBy]
            , preprocessing_method = this_row[, preprocessing_method]
            , preprocessing_args = this_row[, preprocessing_args]
        )
      })
#  rownames(res) <- NULL
  as.data.table(res)
}
#' generate reference gate based on the splitted population patterns(A+/-B+/-) 
.gen_refGate <- function(splitted_terms, this_row, ref_nodes = NULL, alias = NULL) {
  this_parent <- this_row[, parent]
  
  if (is.null(ref_nodes)) {
    # simply copy ref args from the row
    ref_args <- this_row[, gating_args]
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
            .validity_check_alias(cur_alias)
            data.table(alias = cur_alias
                , pop = new_pop
                , parent = this_parent
                , dims = this_row[, dims] 
                , gating_method = "refGate"
                , gating_args = ref_args
                , collapseDataForGating = ""
                , groupBy = ""
                , preprocessing_method = ""
                , preprocessing_args = ""
            )
          }, SIMPLIFY = FALSE))
}

