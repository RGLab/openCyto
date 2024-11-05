#' @templateVar old templateGen
#' @templateVar new gh_generate_template
#' @template template-depr_pkg
NULL
#' generate a partially complete csv template from the existing gating hierarchy 
#' 
#' To ease the process of replicating the existing (usually a manual one) gating schemes, 
#' this function populate an empty gating template with the 'alias', 'pop', 'parent' and 'dims' 
#' columns that exacted from an \code{GatingHierarchy}, and leave the other columns (e.g. `gating_method`) blank.
#' So users can make changes to that template instead of writing from scratch.
#' 
#' @name gh_generate_template
#' @aliases templateGen 
#' @param gh a \code{GatingHierarchy} likely parsed from a xml workspace
#' @return a gating template in \code{data.frame} format that requires further edition after output to csv
#' @examples 
#' library(flowWorkspace)
#' dataDir <- system.file("extdata",package="flowWorkspaceData")
#' gs <- load_gs(list.files(dataDir, pattern = "gs_manual",full = TRUE))
#' gh_generate_template(gs[[1]])
#' @export 
gh_generate_template <- function(gh){
  nodes <- gs_get_pop_paths(gh, order = "tsort")
  dt = do.call(rbind, lapply(nodes[-1], function(thisNode){
        thisGate <- gh_pop_get_gate(gh, thisNode)
        dims <- paste(as.vector(parameters(thisGate)), collapse = ",")
        parent <- gs_pop_get_parent(gh, thisNode)
        alias <- basename(thisNode)
        pop <- alias
        data.frame(alias = alias
            , pop = "+"
            , parent = parent
            , dims = dims
            , gating_method = NA_character_
            , gating_args = NA_character_
            , collapseDataForGating = NA
            , groupBy = NA_character_
            , preprocessing_method = NA_character_
            , preprocessing_args = NA_character_
        )
      })
  )
  if(is.null(dt)||ncol(dt)==0){
    cn = c("alias","pop","parent","dims","gating_method","gating_args","collapseDataForGating","groupBy","preprocessing_method","preprocessing_args")
    dt = vector('list',10)
    names(dt) = cn
    dt = as.data.table(dt)
  }
  return(dt)
}

#' @export
templateGen <- function(gh){
  .Deprecated("gh_generate_template")
}

#' prepend all ancester nodes to construct the full path for the given node 
#' 
#' @param ref_node \code{character} the node to work with
#' @param dt \code{data.table} that stores the current preprocessed template
#' @noRd 
.getFullPath <- function(ref_node, dt){
  
  ref_node <- trimws(ref_node)
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
      # toMatch <- gsub("\\+", "\\\\\\+", curToken)
      toMatch <- paste0("\\Q",curToken,"\\E")
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
#' @noRd 
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
#' @noRd 
.unique_check_alias <- function(new_dt, alias, this_parent) {
  #skip the unexpanded multi-pop pattern 
  #since it is not supposed to be refered 
  if(alias == "*")
    return(NULL)
  siblings <- new_dt[parent == this_parent, alias]
  toMatch <- gsub("\\+", "\\\\\\+", alias)
  toMatch <- paste0("^",toMatch,"$")
  
  matched_sibs <- grep(toMatch, siblings)
  
  if (length(matched_sibs) > 0) {
    stop(alias, " is not unique within ", this_parent)
  }
}

#' split the rows that has multiple parents
#' @noRd 
.split_multi_parents <- function(dt)
{
  #can't use apply since it will coerce dt to df and prepend extra space to int as it convert the column to character
  res<- lapply(1:nrow(dt), function(i){
    
    row <- dt[i,]
    parent <- row[, parent]
    if(!grepl(",", parent))
      row
    else
    {
      message("splitting the row that has multiple parents: '", parent, "'")
      rbindlist(
        lapply(strsplit(parent, split = ",")[[1]], function(p){
          r <- copy(row)
          r[, parent := p]
          r
        })
      ) 
    }
     
  })
 
  rbindlist(res)
}

#' preprocess the csv template
#' 
#' It parses the data table sequentially and does the valdidity checking and expansion row by row.
#'  
#' It expands the definition of gates or construct reference gates when necessary
#' @param dt \code{data.table} loaded directly from csv gating template
#' @param strict whether validity check for pop alias. It is turned on in the normal template parsing. when used with add_pop it is turned off to bypass the check on some existing boolean gates (that has ! : symbols).  
#' @return a preprocessed(expanded when applicable) \code{data.frame}
#' @import data.table
#' @noRd 
.preprocess_csv <- function(dt, strict = TRUE) {
  dt <- .split_multi_parents(dt)
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
  
  #coerce most columns to character 
  #because data.table tends to read the empty columns as logical NAs which causes trouble in parsing later on
  
  for(col in colnames(dt))
  {
    if(col != "collapseDataForGating")
    {
      newCol <- dt[[col]]
      dt[, (col):= as.character(newCol)]
    }
  }
  
  new_dt <- dt[0, ]
  
  for (i in 1:nrow(dt)) {
    this_row <- dt[i, , drop = FALSE]
    if(strict)
      .validity_check_alias(this_row[, alias])
    #update parent with full path
    this_row[, parent := .getFullPath(this_row[, parent], new_dt)]
    #preprocess the current row
    res <- .preprocess_row(this_row, strict = strict)
    
    #check the uniqueness of alias within res
    aliasVec <- res[,alias]
    dupVec <- duplicated(aliasVec)
    if(any(dupVec))
      stop(aliasVec[dupVec], "is duplicated!")
    
    #validity check of alias within the context of new_df
    apply(res, 1, function(row){
          thisAlias <- row[["alias"]]
          if(strict)
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
#' 3. dispatch 'flowClust' method to either 'gate_flowclust_1d' or 'gate_flowclust_2d' based on the 'pop' and 'dims' columns
#' 4. expand the single row to multiple rows when applicable. There are basically two types of expansion:
#'      4.1. expand to two 1d gates and one rectangelGate when 'pop' name is defined as quadrant pattern (e.g. "A+B+")  and 'gating_method' is not "refGate"
#'      4.2. expand to multiple gates when 'pop' is defined with '+/-' (e.g. "A+/-" or "A+/-B+/-") 
#' 
#' @param this_row a single-row \code{data.table}
#' @return  \code{data.table}
#' @noRd 
.preprocess_row <- function(this_row, strict = TRUE){
  #make sure it doesn't tamper the input
  this_row <- copy(this_row)
  
  alias <- this_row[1, alias]
  popName <- this_row[1, pop]
  dims <- this_row[1, dims]
  gm <- this_row[1, gating_method]
  
  if(gm == "dummy") #skip dummy since it doesn't need any processing
    return(this_row)
  
  if(popName == "*" ){
    if(alias == "*"|| dims %in% c("", NA))
      return(this_row)
    else
      return(.gen_dummy_ref_gate(this_row))
    }
    
  if(isTRUE(getOption("openCyto")[["check.pop"]])&&!grepl("^(([\\+-])|(\\+/-)|(-/\\+)){1,2}$",popName))
    stop("'pop': ", popName, " should only contain the + and - symbols and must conform to the valid +/- or +/-+/- patterns!
			\n Please remove any other letters and correct the pop pattern!
			\n Or type ?openCyto.options to see how you can turn off the 'check.pop' flag in options('openCyto') to bypass this validiy check (Not recommended)")
  
  dim_count <- length(strsplit(split = ",", dims)[[1]])
  if (!dim_count %in% c(1, 2)) {
    if (!(dim_count == 0 && gm %in% c("polyFunctions", "boolGate", "refGate"))) {
      stop(popName, " has invalid number of dimensions: ", dim_count)
    }
  }
  
  one_pop_token <- "[\\+-]"
  pop_name_pat <- "[^\\+-]*" #now we allow empty pop name with pure +/- signs  "[^\\+-]+"
  one_pop_pat <- paste(pop_name_pat, one_pop_token, sep = "")
  
  two_pop_token <- "(\\+/-)|(-/\\+)"
  two_pop_pat <- paste(pop_name_pat, "(", two_pop_token, ")", sep = "")
#    browser()
  if (grepl(paste0("^", one_pop_pat, "$"), popName)) {
    # A+ no expansion(simply update flowClust gm)
    if (gm %in% c("flowClust", "gate_flowclust")) {
      if (dim_count == 1) {
        this_row[1, gating_method := "gate_flowclust_1d"] 
      } else {
        this_row[1, gating_method := "gate_flowclust_2d"] 
      }
    }
    res <- this_row
    #strip the other letters from legacy popName
    res[, pop := ifelse(grepl("-$", popName), "-", "+")]
    
  } else if (grepl(paste("^", two_pop_pat, "$", sep = ""), popName)) {
    # A+/-
    
    if (gm %in% c("flowClust", "gate_flowclust")) {
      if (dim_count == 1) {
        this_row[1, gating_method := "gate_flowclust_1d"]
      } else {
        this_row[1, gating_method := "gate_flowclust_2d"]
      }
    }
    # expand to two rows
    new_pops <- c("+", "-")
    
    if(dim_count == 1){
      new_alias <- paste(dims, new_pops, sep = "")  
      if(!trimws(alias) %in% c("", "*"))
        message("alias '", alias, "' is ignored since pop names are auto-generated from 'dims' when pop = '+/-' ")
    }else if(dim_count == 2)
    {
      
      if(alias == "*")
        stop("Please provide a proper alias name because dims(", dims, ") can't be used to derive the expanded pops!")
      new_alias <- paste(alias, new_pops, sep = "")  
      
    }else
      stop("Don't know how to handle ", popName, " for multi-dimensional gating! ")
    
    message("expanding pop: ", popName, " to ", paste(new_alias, collapse = "/"), "\n")
    
    
    
    res <- rbindlist(list(this_row, this_row))
    
    # update 1d gate
    res[1, alias := new_alias[1]]
    res[1, pop := new_pops[1]]
    
    
    # update ref gate
    refNode <- file.path(this_row[, parent], new_alias[1])
    res[2, alias := new_alias[2]]
    res[2, pop := new_pops[2]]
    res[2, gating_method := "refGate"]
    res[2, gating_args := refNode]
    res[2, preprocessing_method := ""]
    res[2, preprocessing_args := ""]
    
  } else if (grepl(paste("^(", one_pop_pat, "){2}$", sep = ""), popName)) {
    # A+B+
    split_terms <- .splitTerms(pop_pat = one_pop_pat, two_pop_token, popName, dims)
    new_pops <- as.character(outer(split_terms[[1]][["pop"]], split_terms[[2]][["pop"]], paste,sep = ""))
    new_alias <- as.character(outer(split_terms[[1]][["alias"]], split_terms[[2]][["alias"]], paste,sep = ""))
    
    if (gm == "refGate") {
      # no expansion
      res <- this_row
      res[, pop := new_pops]
    } else {
      if(gm %in% c("quadGate.seq", "gate_quad_sequential", "quadGate.tmix", "gate_quad_tmix")){
        if(!trimws(alias) %in% c("", "*"))
          message("alias '", alias, "' is ignored since pop names are auto-generated from 'dims' for gate_quad_ methods")
        this_row[1, alias := "*"]
        res <- this_row
      }else{
        if (gm %in% c("flowClust", "gate_flowclust")) {
          if (dim_count == 2) {
            message("expanding pop: ", popName, "\n")
            
            this_row[1, gating_method := "gate_flowclust_1d"]
          } else {
            stop("dimensions '", dims, "' is not consistent with pop name '", 
                 popName, "'")
          }
        }
        # needs to be split into two 1d gates and one refgate
        
        # create 1d gate for each dim
        res_1d <- .gen_1dgate(this_row, strict = strict)
        res_ref <- .gen_refGate(new_pops, this_row, ref_nodes = res_1d[, alias], alias = this_row[, alias], strict = strict)
        res <- rbindlist(list(res_1d, res_ref)) 
      }
    }
  } else if (grepl(paste0("^(", two_pop_pat, "){2}$"), popName) ||
      grepl(paste0("^", two_pop_pat, one_pop_pat, "$"), popName) ||
      grepl(paste0("^", one_pop_pat, two_pop_pat, "$"), popName)) {
    # A+/-B+/-
    message("expanding pop: ", popName, "\n")
    two_or_one_pop_pat <- paste0("(", two_pop_pat, ")|(", one_pop_pat, ")")
    split_terms <- .splitTerms(pop_pat = two_or_one_pop_pat, two_pop_token, popName, dims)
    new_pops <- as.character(outer(split_terms[[1]][["pop"]], split_terms[[2]][["pop"]], paste,sep = ""))
    new_alias <- as.character(outer(split_terms[[1]][["alias"]], split_terms[[2]][["alias"]], paste,sep = ""))
    if (gm == "refGate") {
      res <- .gen_refGate(new_pops, this_row = this_row, alias = new_alias, strict = strict)
      
    } else {
      if(gm %in% c("quadGate.seq", "gate_quad_sequential", "quadGate.tmix", "gate_quad_tmix")){
        if(!trimws(alias) %in% c("", "*"))
          message("alias '", alias, "' is ignored since pop names are auto-generated from 'dims' gate_quad_ methods")
        this_row[1, alias := "*"]
        res <- this_row
      }else{
        if (gm %in% c("flowClust", "gate_flowclust")) {
          if (dim_count == 2) {
            this_row[1, gating_method := "gate_flowclust_1d"]
            
          } else {
            stop("dimensions '", dims, "' is not consistent with pop name '", 
                 popName, "'")
          }
        }
        
        # create 1d gate for each dim
        res_1d <- .gen_1dgate(this_row, strict = strict)
        res_ref <- .gen_refGate(new_pops, this_row, ref_nodes = res_1d[, alias], alias = new_alias, strict = strict)
        res <- rbindlist(list(res_1d, res_ref))
      }
    }
    
  } else 
    stop("invalid population pattern '", popName, "'")
  res
}


#' split the population pattern into multiple population names
#' 
#' currently only works for A+/-B+/- .
#' TODO:  support A+/-  
#' @noRd 
.splitTerms <- function(pop_pat, two_pop_token, popName, dims) {
  dims_vec <- trimws(strsplit(split = ",", dims)[[1]])
  if(length(dims_vec)!=2)
    stop("Can't split the population pattern '", popName, "' into multiple population names due to invalid 'dims' column: '", dims, "'" )
  term_pos <- gregexpr(pop_pat, popName)[[1]]
  x_term <- substr(popName, 1, term_pos[2] - 1)
  y_term <- substr(popName, term_pos[2], nchar(popName))
  terms <- c(x_term, y_term)
  
  splitted_terms <- mapply(terms, dims_vec, FUN = function(cur_term, dim_name) {
        token_pos <- gregexpr(two_pop_token, cur_term)[[1]]
        if (token_pos > 0) {
          # dim_name <- substr(cur_term, 1, token_pos - 1)
          splitted_term <- list(alias = paste(dim_name, c("+", "-"), sep = "")
                                , pop = c("+", "-")
                            )
          
        } else {
          nlen <- nchar(cur_term)
          cur_token <- substr(cur_term, nlen, nlen)
          splitted_term <- list(alias = paste0(dim_name, cur_token)
                                , pop = cur_token
                            )
        }
        splitted_term
      }, SIMPLIFY = FALSE)
  splitted_terms
}

#' generate some dummy rows that just serves as reference without any gating method associated
#' @noRd 
.gen_dummy_ref_gate <- function(this_row){
  
  alias <- this_row[, alias]
  
  refNode <- file.path(this_row[, parent], alias)
  
  pops <- trimws(unlist(strsplit(split = ",", alias)))
  new_rows <- lapply(pops, function(thisPop){
                          dummy_row <- copy(this_row)
                          dummy_row[, alias := thisPop]
                          dummy_row[, pop := ""]
                          dummy_row[, gating_method := "dummy_gate"]
                          dummy_row[, gating_args := refNode]
                          dummy_row[, collapseDataForGating := ""]
                          dummy_row[, groupBy := ""]
                          dummy_row[, preprocessing_method := ""]
                          dummy_row[, preprocessing_args := ""]
                          dummy_row
                        })
  
  rbindlist(list(this_row,rbindlist(new_rows)))    
}
#' convert to 1d gating based on the population pattern (A+/-B+/-)
#' @noRd 
.gen_1dgate <- function(this_row, strict = TRUE) {
  
  dims.dims <- trimws(strsplit(split = ",", this_row[, dims])[[1]])
  
  res <- do.call(rbind, lapply(dims.dims, function(cur_dim) {
        new_pop_name <- paste(cur_dim, "+", sep = "")
        this_parent <- this_row[, parent]
        if(strict)
          .validity_check_alias(new_pop_name)
        
        data.table(alias = new_pop_name
            , pop = "+"
            , parent = this_parent
            , dims = cur_dim 
            , gating_method = this_row[, gating_method]
            , gating_args = this_row[, gating_args]
            , collapseDataForGating = this_row[, collapseDataForGating]
            , groupBy = this_row[, groupBy]
            , preprocessing_method = this_row[, preprocessing_method]
            , preprocessing_args = this_row[, preprocessing_args]
        )
      }))
#  rownames(res) <- NULL
  as.data.table(res)
}
#' generate reference gate based on the splitted population patterns(A+/-B+/-) 
#' @noRd 
.gen_refGate <- function(new_pops, this_row, ref_nodes = NULL, alias, strict = TRUE) {
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
  
  
  
  # create ref gate for each new_pop )
  do.call(rbind, mapply(new_pops, alias, FUN = function(new_pop, cur_alias) {
            if(strict)
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

