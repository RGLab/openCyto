#' overwrite the coerce method for ncdfFlowList defined in ncdfFlow package
#' 
#' It is at least 2 times (for 2-channel-9-sample ncfs collapsing) faster than origina version. 
setAs(from = "ncdfFlowList", to = "flowFrame", def = function(from){
      selectMethod("coerce", signature = c("ncdfFlowSet", "flowFrame"))(from)      
      
    })

#' fast collapsing
#' 
#' It is revised version of flowSet coerce method,  
#' 1. without the overhead of validity checking for parameters 
#' 2. without the extra column 'Original'
#' 3. concatenate matrices within c++ to avoid copying
#'  
#' @importFrom flowCore colnames
setAs(from="ncdfFlowSet", to="flowFrame", def=function(from)
    {
      if(length(from) == 1)
        from[[1]]
      else {
        thisFr <- from[[1, use.exprs = FALSE]]
        
        params <- parameters(thisFr)
        colnames <- colnames(thisFr)
        
        mat_list <- fsApply(from, function(fr)exprs(fr), simplify = FALSE)
        
        
        desc  <- list(description="Synthetic Frame",sampleNames=sampleNames(from))
        new("flowFrame",exprs = .Call('openCyto_collapseData', PACKAGE = 'openCyto', mat_list, colnames)
            , parameters = params
            , description = desc
        )
      }
    }
)
