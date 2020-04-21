#' convert a gatingTemplate object to a data.table
#' 
#' It is the inverse function of gatingTemplate constructor.
#' 
#' @name as.data.table
#' @param x gatingTemplate object
#' @param keep.rownames not used
#' @return a data.table
#' @export
as.data.table.gatingTemplate <- function(x, keep.rownames = FALSE){
  
  
  # gate each node 
  gt_nodes <- gt_get_nodes(x, order = "tsort")[-1]
  
  
  res <- lapply(gt_nodes, function(node){
    
    # get parent node to gate
    nodePath <- node@id
    parent <- gt_get_parent(x, nodePath)
    # extract gate method from one edge(since multiple edge to the same node is
    # redudant)
    this_gate <- gt_get_gate(x, parent, nodePath)
    gating_args <- parameters(this_gate)
    #collapse into string when neccessary
    split <- !extends(class(this_gate), "refGate")
    gating_args <- .argDeparser(gating_args, split)
    
    #get preprocessing method
    this_ppm <- openCyto:::ppMethod(x, parent, nodePath)
    
    #preprocessing
    if(class(this_ppm) == "ppMethod")
    {
      preprocessing_method <- names(this_ppm)
      preprocessing_args <- parameters(this_ppm)
      preprocessing_args <- .argDeparser(preprocessing_args)
    }else{
      preprocessing_args <- preprocessing_method <- ""
    }
    
    data.table(alias = alias(node)
                          , pop = names(node)
                          , parent = parent
                          , dims = this_gate@dims
                          , gating_method = names(this_gate)
                          , gating_args =  gating_args
                          , collapseDataForGating = openCyto:::isCollapse(this_gate)
                          , groupBy = openCyto:::groupBy(this_gate)
                          , preprocessing_method = preprocessing_method
                          , preprocessing_args = preprocessing_args
                      )    
    
  })
  
  rbindlist(res)
  
}

#' deparse a list(named) of expression into a string 
#' inverse function of argParser
#' @noRd 
.argDeparser <- function(args, split = TRUE){
  if(split){
    args <- sapply(names(args), function(argn){
          
          argv <- deparse(args[[argn]])
          argv <- gsub("\"", "'", argv) #restore dquote to squote
          argv <- paste(argv, collapse = "")
          paste(argn, argv, sep = " = ")
        })
    
    
    paste(args, collapse = ", ")  
  }else
    as.character(args[[1]])
    
}
