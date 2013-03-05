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
setClass("boolMethod"
		,contains="gtMethod"
		)
setClass("polyFunctions"
		,contains="boolMethod"
		)
setClass("gtPopulation"
		,representation(id="numeric" #consistent with node label in gating template graph
						,name="character"
						,alias="character"
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
.argParser<-function(txt){
	
	#trim whitespaces at beginning and the end of string
	txt<-gsub("^\\s*|\\s*$","",txt)
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
	paired_args
}
setGeneric("gatingTemplate", function(x, ...) standardGeneric("gatingTemplate"))
#constructor from csv
setMethod("gatingTemplate",signature(x="character"),function(x,name){
	
	#read csv
	df<-read.csv(x,as.is=T)
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
		parent<-df[i,"parent"]
		curPop<-df[i,"alias"]
		curPopName<-df[i,"pop"]
		curNode<-new("gtPopulation",id=as.numeric(curNodeID),name=curPopName,alias=curPop)
					
#		browser()
		gm<-new("gtMethod"
				,name=df[i,"method"]
				,dims=df[i,"dims"]
				,args=.argParser(df[i,"args"])
				)
				
		if(names(gm)=="boolGate")
			gm<-as(gm,"boolMethod")
		else if(names(gm)=="polyFunctions")
			gm<-as(gm,"polyFunctions")
					
		cat("Adding population:",curPop,"\n")
		if(parent=="root")
		{
			src<-"1"
			g<-graph::addNode(curNodeID,g)
			
			nodeData(g,curNodeID,"pop")<-curNode
			g<-addEdge(src,curNodeID,g)

			edgeData(g,src,curNodeID,"gtMethod")<-gm
		}else
		{
			#travese the graph to get node list for matching later on
			discovered<-dfs(g)[[1]]
			g<-graph::addNode(curNodeID,g)
			
#			TODO:deal with not symbol!
#			isPoly<-isPolyfunctional(gm)
			if(class(gm)=="polyFunctions")
				refNodes<-strsplit(parent,"\\:")[[1]] #split by colon when polyfunctional boolean gates
			else if(class(gm)=="boolMethod")
				refNodes<-strsplit(parent,"&|\\|")[[1]]#split by logical operator when regular boolean gates
			else
				refNodes<-parent
			if(class(gm)=="polyFunctions")
				curNode<-as(curNode,"gtSubsets")
			nodeData(g,curNodeID,"pop")<-curNode
			
			for(refNode in refNodes)
			{
						#split by "/" for each reference node
						
				
				tokens<-strsplit(refNode,"/")[[1]]
				##locate the first token in traversed node list 
				firstToken<-tokens[1]
				tokens<-tokens[-1]
				curAncesterID<-NULL
				for(ancesterID in discovered)
				{
#										browser()
					if(alias(getNodes(g,ancesterID))==firstToken)
					{
						curAncesterID<-ancesterID
						
					}
					
				}
				if(is.null(curAncesterID))
					stop("Population '",firstToken,"' not found")
				#start from matchedID to match the rest of tokens in the path
				while(length(tokens)>0)
				{
					curToken<-tokens[1]
					tokens<-tokens[-1]
					#find the id within the edges sourced from current ancester
					dests<-edges(g,curAncesterID)[[1]]
					curAncesterID<-NULL
					for(ancesterID in dests)
					{
#										browser()
						if(alias(getNodes(g,ancesterID))==curToken)
						{
							curAncesterID<-ancesterID
							
						}
						
					}
					if(is.null(curAncesterID))
						stop("Population '",curToken,"' not found")
				}
#				browser()
				g<-addEdge(curAncesterID,curNodeID,g)
				edgeData(g,curAncesterID,curNodeID,"gtMethod")<-gm
			}
			
		}
	
	}
	
	
	
#	browser()
	
	g@name<-name
	g
})

