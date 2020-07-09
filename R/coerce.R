# coerce method for ncdfFlowList
# 
# overwrite the version defined in ncdfFlow package
# It is at least 2 times (for 2-channel-9-sample ncfs collapsing) faster than origina version.
# @param from \code{ncdfFlowList}
#' @importClassesFrom ncdfFlow ncdfFlowList
#' @importFrom methods coerce
setAs(from = "ncdfFlowList", to = "flowFrame", def = function(from){
      selectMethod("coerce", signature = c("ncdfFlowSet", "flowFrame"))(from)      
      
    })

# fast collapsing
# 
# It is revised version of flowSet coerce method,  
# 1. without the overhead of validity checking for parameters 
# 2. without the extra column 'Original'
# 3. concatenate matrices within c++ to avoid copying
fast_coerce <- function(from){
  if(length(from) == 1)
  	from[[1, returnType = "flowFrame"]]
  else {
  	thisFr <- from[[1, use.exprs = FALSE]]
  	
  	params <- parameters(thisFr)
  	colnames <- colnames(thisFr)
  	
  	mat_list <- fsApply(from, function(fr)exprs(fr), simplify = FALSE)
  	
  	
  	desc  <- list(description="Synthetic Frame",sampleNames=sampleNames(from))
  	new("flowFrame",exprs = collapseData(mat_list, colnames)
  			, parameters = params
  			, description = desc
  	)
  }
}  
#' @importFrom flowCore colnames
#' @importClassesFrom ncdfFlow ncdfFlowSet
setAs(from="ncdfFlowSet", to="flowFrame", def=function(from)
			fast_coerce(from)
)
setAs(from="cytoset", to="flowFrame", def=function(from){
			fast_coerce(from)
		})