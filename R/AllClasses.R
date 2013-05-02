###############################################################################
## a class representing the gating method and population info in a graph object 
###############################################################################
setClass("gatingTemplate"
		,contains="graphNEL"
		,representation(name="character")
	)
	
###############################################################################	
#extend flowCore gate classes to have extra slot to store posteriors from flowClust gating routine 	
###############################################################################
setClass("fcFilter",representation("VIRTUAL",priors="list",posteriors="list"))

#setClassUnion("fcRectangleGate",c("rectangleGate","fcFilter"))
setClass("fcRectangleGate",contains=c("fcFilter","rectangleGate"))


fcRectangleGate<-function(x,priors,posts)
{
	res<-as(x,"fcRectangleGate")
    res@priors<-priors
    res@posteriors<-posts
	res
}

setClass("fcPolygonGate",contains=c("fcFilter","polygonGate"))

fcPolygonGate<-function(x,priors,posts)
{
  res<-as(x,"fcPolygonGate")
  res@priors<-priors
  res@posteriors<-posts
  res
}

setClass("fcFilterList"
		,contains="filterList"
		)
fcFilterList<-function(x)
{
#  browser()
  
  if(!all(unlist(lapply(x,function(i)extends(class(i),"fcFilter")))))
    stop("not all filters are fcFilter!")
  if(class(x) == "list")
    x <- filterList(x)
  
  sname<-names(x)
  x <- as(x,"fcFilterList")

  attr(x,"names")<-sname
  x

  	
}
#for the purpose of method dispatching
#setClass("fcFilterList2d"
#		,contains="fcFilterList"
#)
###############################################################################	
##a flowClust tree is a container to hold priors and posteriors that can be visualized
## for the purpose of fine-tunning parameters for flowClust algorithm
###############################################################################
setClass("fcTree"
		,contains="gatingTemplate"
		)
fcTree<-function(gt){
		
		res<-as(gt,"fcTree")
		nodeDataDefaults(res,"fList")<-new("filterList")
		res
		}	
setClass("gtMethod"
		,representation(name="character"
						,dims="character"
						,args="list"
						)	
)

setClass("refGate"
    ,contains="gtMethod"
    ,representation(refNodes="character")
)

setClass("boolMethod"
		,contains="refGate"
		)
setClass("polyFunctions"
		,contains="boolMethod"
		)
        
        
setClass("gtPopulation"
		,representation(id="numeric" #consistent with node label in gating template graph
						,name="character"
						,alias="character"
                        ,parentID="numeric"
		)	
)
setClass("gtSubsets"
		,contains="gtPopulation"
		)	


isPolyfunctional<-function(gm){
#	browser()
#	grepl("^\\[\\:.+\\:\\]$",x)
	class(gm)=="polyFunctions"
}

#gating arguments parser
.argParser<-function(txt,split=TRUE){
	
	#trim whitespaces at beginning and the end of string
	txt<-gsub("^\\s*|\\s*$","",txt)
    if(split){
    	paired_args<-paste("c(",txt,")")
    	paired_args<-try(parse(text=paired_args),silent = TRUE)
    	if(class(paired_args)=="try-error")
    	{
    		errmsg<-attr(paired_args,"condition")
    		msg <- conditionMessage(errmsg)
    		stop("invalid gating argument:\n",msg)
    	}
        
        paired_args<-as.list(as.list(paired_args)[[1]])[-1]
        names(paired_args)<-tolower(names(paired_args))  
    }else
    {
      
#      browser()
      paired_args<-as.symbol(txt)
      
      paired_args<-list(paired_args)
    }
	
	paired_args
}

#match node ID by name with the subset
.getNodeID <- function(g, subset, this_node){
  node_id_ind <- which(unlist(lapply(subset,function(this_node_id){
                alias(getNodes(g,this_node_id)) == this_node
              })))
  
  
  if(length(node_id_ind) == 0){
    stop("Population '", this_node, "' not found under the current parent:", alias(getNodes(g,getParent(g,subset[2]))))
  }else if(length(node_id_ind) > 1){
    stop("Population '", this_node, "' not unique under the current parent:", alias(getNodes(g,getParent(g,subset[2]))))
  }else{
    node_id <- subset[node_id_ind]
  }
  node_id
}
#search for node ID by path with gating template tree and add the edge
.searchNode<-function(g, node_name){
#browser()
  if(node_name == "root"){
    node_id <- "1"
  }else{
    #travese the graph to get node list for matching later on
    discovered<-dfs(g)[["discovered"]]
    #split by "/" for each reference node
      
    node_name <- flowWorkspace:::trimWhiteSpace(node_name)
    tokens<-strsplit(node_name,"/")[[1]]
    
    ##locate the first token in traversed node list 
    firstToken<-tokens[1]
    tokens<-tokens[-1]

#    browser()
    node_id <- .getNodeID(g, subset = discovered, this_node =firstToken)
      
    #start from matchedID to match the rest of tokens in the path
    while(length(tokens)>0)
    {
      curToken<-tokens[1]
      tokens<-tokens[-1]
      #find the id within the edges sourced from current ancester
      dests<-edges(g,node_id)[[1]]
      node_id <- .getNodeID(g, subset = dests, this_node =curToken)        
    }
    
  }
      
  node_id
    
  
}
setGeneric("gatingTemplate", function(x, ...) standardGeneric("gatingTemplate"))
#constructor from csv
setMethod("gatingTemplate",signature(x="character"),function(x,name){
	
	df <- .preprocess_csv(x)
#    browser()
	#create graph with root node
	g<-graphNEL(nodes="1",edgemode="directed")
	g<-as(g,"gatingTemplate")
	nodeDataDefaults(g,"pop")<-""
	edgeDataDefaults(g,"gtMethod")<-""
	#add default root
	nodeData(g,"1","pop")<-new("gtPopulation",id=1,name="root",alias="root")
	#parse each row
	nEdges<-nrow(df)
	edgs<-vector("list",nEdges)
	for(i in 1:nEdges){
		
		curNodeID<-as.character(i+1)
        #extract info from dataframe
		parent<-as.character(df[i,"parent"])
        #get parent ID
        parentID <- .searchNode(g, parent)
#        browser()
		curPop<-as.character(df[i,"alias"])
		curPopName<-as.character(df[i,"pop"])
        #create pop object
		curNode<-new("gtPopulation" 
                      , id = as.numeric(curNodeID)
                      , name = curPopName
                      , alias = curPop
                      , parentID=as.numeric(parentID)
                      )
					
        #create gating method object
        cur_method <- as.character(df[i, "method"])
        cur_args <- as.character(df[i, "args"])
        cur_dims <- as.character(df[i, "dims"])
        
        #do not parse args for refGate-like gate since they might break 
        #the current parse due to the +/- | &,! symbols
        if(cur_method%in%c("boolGate","polyfunctions","refGate")){
          split_args <- FALSE
        }else{
          split_args <- TRUE
        }
        cur_args <- .argParser(cur_args,split_args)
#          browser()
          
		gm<-new("gtMethod"
				, name = cur_method
				, dims = cur_dims
				, args = cur_args
				)
		#specialize gtMethod as needed		
		if(names(gm) == "boolGate"){
          gm <- as(gm, "boolMethod")
        }else if(names(gm) == "polyFunctions"){
          gm <- as(gm, "polyFunctions")
        }else if(names(gm) == "refGate"){
          gm <- as(gm, "refGate")
        }
           
          
		cat("Adding population:", curPop, "\n")
        #add current node to graph
        g_updated<-graph::addNode(curNodeID,g)
        
        if(!extends(class(gm), "refGate")){
          #add edge from parent
          g_updated<-addEdge(parentID, curNodeID,g_updated)
          #add the gm object to the edge
          edgeData(g_updated,parentID, curNodeID, "gtMethod") <- gm
        }else{
          ##########################################
          #refGate-like methods need extra parsing
          ##########################################
          
          #get argument 
          args <- gm@args[[1]]
          args <- deparse(args)
          
          #parsing reference nodes
          if(class(gm) == "boolMethod"){
            #strip ! symbols
            args <- gsub("!", "", args) 
            #split by logical operator when regular boolean gates
            refNodes <- strsplit(args, "&|\\|")[[1]]
          }else{
            #split by colon for refGate or polyfunctional boolean gates
            refNodes <- strsplit(args, "\\:")[[1]] 
          }
            
           #update refNodes slot of gm object 
           gm@refNodes <- refNodes                
            
          #specialize the node type for polyfunctions
          if(class(gm)=="polyFunctions"){
            curNode<-as(curNode,"gtSubsets")
          }
            
          
          #add edges from reference nodes
          for(ref_node in refNodes)
          {
            #get node id for reference node
            #(using the old graph object since the new graph has unconnected new node
            ref_id <- .searchNode(g, ref_node)
             #add the edge from it 
            g_updated<-addEdge(ref_id, curNodeID, g_updated)
            #append the gm object to the edge
            edgeData(g_updated,ref_id, curNodeID, "gtMethod") <- gm
          }
                  
        }
        
        #add the current population object to the current node
        nodeData(g_updated, curNodeID, "pop") <- curNode 
        #update graph
        g<-g_updated
   }

		
	g@name<-name
	g
})

