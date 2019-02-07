#' Reverse the action of gating methods applied via \code{add_pop}
#' 
#' This function provides an easy way to remove the gates and nodes created by the most
#' recent call to \code{\link{add_pop}} on the specified \code{GatingSet} or \code{GatingSetList},
#' with a separate history being maintained for each such object. \code{remove_pop} allows 
#' for repeated use, effectively serving as a multi-level undo function for \code{add_pop}.
#' 
#' @param gs The \code{GatingSet} or \code{GatingSetList} for which the most recent
#' \code{add_pop} call should be reversed.
#' 
#' @usage 
#' remove_pop(GatingSet)
#' remove_pop(GatingSetList)
#' 
#' @seealso \code{\link{add_pop}} \code{\link{add_pop_init}}
#' @examples
#' \dontrun{
#'  # add quad gates 
#'  add_pop(gs, gating_method = "mindensity", dims = "CCR7,CD45RA", parent = "cd4-cd8+", pop = "CCR7+/-CD45RA+/-")
#'  # Remove the gates and nodes resulting from that add_pop call
#'  remove_pop(gs)
#' } 
#' @export
remove_pop <- function(gs){
  if(!(gs@guid %in% names(add_pop_history$records))){
    stop(paste("No add_pop calls have been made for this GatingSet or GatingSetList"))
  }
  this_record <- add_pop_history$records[[gs@guid]]
  old_ids <- this_record[[length(this_record)-1]]
  new_ids <- this_record[[length(this_record)]]
  new_ids <- new_ids[!(new_ids %in% old_ids)] # setdiff won't work here
  lapply(new_ids, function(x) Rm(x, gs))
  
  ## If that succeeded, remove the records
  # If it's a GatingSetList, remove the records from all of its GatingSets
  if(is(gs, "GatingSetList")){
    lapply(gs, function(x){
      add_pop_history$records[[x@guid]][[length(add_pop_history$records[[x@guid]])]] <- NULL
      # If back to initial state, clear the add_pop history
      if(length(add_pop_history$records[[x@guid]]) <= 1){
        add_pop_history$records[[x@guid]] <- NULL
      }
    })
  }
  # Pop off the record
  add_pop_history$records[[gs@guid]][[length(this_record)]] <- NULL
  # If back to initial state, clear the add_pop history
  if(length(add_pop_history$records[[gs@guid]]) <= 1){
    add_pop_history$records[[gs@guid]] <- NULL
  }
  
  invisible(this_record)
}