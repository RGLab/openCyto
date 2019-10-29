#' @templateVar old remove_pop
#' @templateVar new gs_remove_gating_method
#' @template template-depr_pkg
NULL
#' Reverse the action of gating methods applied via \code{gs_add_gating_method}
#' 
#' This function provides an easy way to remove the gates and nodes created by the most
#' recent call to \code{\link{gs_add_gating_method}} on the specified \code{GatingSet} or \code{GatingSetList},
#' with a separate history being maintained for each such object. \code{gs_remove_gating_method} allows 
#' for repeated use, effectively serving as a multi-level undo function for \code{gs_add_gating_method}.
#' 
#' @name gs_remove_gating_method
#' @aliases remove_pop
#' @param gs The \code{GatingSet} or \code{GatingSetList} for which the most recent
#' \code{gs_add_gating_method} call should be reversed.
#' 
#' @usage 
#' gs_remove_gating_method(gs)
#' 
#' @seealso \code{\link{gs_add_gating_method}} \code{\link{gs_add_gating_method_init}}
#' @examples
#' \dontrun{
#'  # add quad gates 
#'  gs_add_gating_method(gs, gating_method = "mindensity", dims = "CCR7,CD45RA", parent = "cd4-cd8+", pop = "CCR7+/-CD45RA+/-")
#'  # Remove the gates and nodes resulting from that gs_add_gating_method call
#'  gs_remove_gating_method(gs)
#' } 
#' @importFrom flowWorkspace identifier
#' @export
gs_remove_gating_method <- function(gs){
  if(!(identifier(gs) %in% names(add_pop_history$records))){
    stop(paste("No gs_add_gating_method calls have been made for this GatingSet or GatingSetList"))
  }
  this_record <- add_pop_history$records[[identifier(gs)]]
  old_ids <- this_record[[length(this_record)-1]]
  new_ids <- this_record[[length(this_record)]]
  new_ids <- new_ids[!(new_ids %in% old_ids)] # setdiff won't work here
  lapply(new_ids, function(x) gs_pop_remove(x, gs = gs))
  
  ## If that succeeded, remove the records
  # If it's a GatingSetList, remove the records from all of its GatingSets
  if(is(gs, "GatingSetList")){
    lapply(gs, function(x){
      add_pop_history$records[[identifier(x)]][[length(add_pop_history$records[[identifier(x)]])]] <- NULL
      # If back to initial state, clear the add_pop history
      if(length(add_pop_history$records[[identifier(x)]]) <= 1){
        add_pop_history$records[[identifier(x)]] <- NULL
      }
    })
  }
  # Pop off the record
  add_pop_history$records[[identifier(gs)]][[length(this_record)]] <- NULL
  # If back to initial state, clear the add_pop history
  if(length(add_pop_history$records[[identifier(gs)]]) <= 1){
    add_pop_history$records[[identifier(gs)]] <- NULL
  }
  
  invisible(this_record)
} 

#' @export
remove_pop <- function(gs){
  .Deprecated("gs_remove_gating_method")
  gs_remove_gating_method(gs)
}