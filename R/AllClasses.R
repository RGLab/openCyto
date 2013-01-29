# TODO: Add comment
# 
# Author: wjiang2
###############################################################################
## exsits for the purpose of dispatching 
setClass("gatingTemplate",
		contains="graphNEL"
		,representation(name="character")
#		,validity=function(object){
#			#check if the parent is ambiguous or not valid			
#			parent_matches<-sapply(levels(object[,"parent"]),function(x){
#						matches<-length(grep(x,object[,"pop"],fix=T))
#						if(matches==0&&x=="root")
#							matches<-1
#						matches
#					})
#			
#			invalid_parent<-which(parent_matches==0)
#			ambiguous_parent<-which(parent_matches>1)
#			if(length(invalid_parent)>0)
#				return(paste("The following parent population names are not valid:"
#							,paste(names(invalid_parent)
#								,collapse=",")
#							)
#						)
#			else if(length(ambiguous_parent)>0)
#				return(paste("The following parent population names are ambiguous:"
#								,paste(names(ambiguous_parent)
#										,collapse=",")
#						)
#				)
#			else 
#				TRUE
#		}
	)

setClass("gatingMethod"
		,representation(name="character"
						,dims="character"
						,args="character"
						)	
)

gatingTemplate<-function(file,name){
	#read csv
	df<-read.csv(file,as.is=T)
	#create graph with root node
	g<-graphNEL(nodes="0",edgemode="directed")
	nodeDataDefaults(g,"label")<-"root"
	edgeDataDefaults(g,"gatingMethod")<-""
	#parse each row
	nEdges<-nrow(df)
	edgs<-vector("list",nEdges)
	for(i in 1:nEdges){
		
		curNodeID<-as.character(i)
		parent<-df[i,"parent"]
		curPop<-df[i,"alias"]
#		browser()
		gm<-new("gatingMethod"
				,name=df[i,"method"]
				,dims=df[i,"dims"]
				,args=df[i,"args"]
				)
		cat("Adding population:",curPop,"\n")
		if(parent=="root")
		{
			src<-"0"
			g<-graph::addNode(curNodeID,g)
			nodeData(g,curNodeID,"label")<-curPop
			g<-addEdge(src,curNodeID,g)
#			browser()
			edgeData(g,src,curNodeID,"gatingMethod")<-gm
		}else
		{
			#travese the graph to get node list for matching later on
			discovered<-dfs(g)[[1]]
			g<-graph::addNode(curNodeID,g)
			nodeData(g,curNodeID,"label")<-curPop
			#split by logical operator in case of boolean gates
			refNodes<-strsplit(parent,"&|\\|")[[1]]
			for(refNode in refNodes)
			{
						#split by "/" for each reference node
						
#				browser()
				tokens<-strsplit(refNode,"/")[[1]]
				##locate the first token in traversed node list 
				firstToken<-tokens[1]
				tokens<-tokens[-1]
				curAncesterID<-NULL
				for(ancesterID in discovered)
				{
#										browser()
					if(nodeData(g,ancesterID,"label")==firstToken)
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
						if(nodeData(g,ancesterID,"label")==curToken)
						{
							curAncesterID<-ancesterID
							
						}
						
					}
					if(is.null(curAncesterID))
						stop("Population '",curToken,"' not found")
				}
#				browser()
				g<-addEdge(curAncesterID,curNodeID,g)
				edgeData(g,curAncesterID,curNodeID,"gatingMethod")<-gm
			}
			
		}
	
	}
	
	
	
#	browser()
	gt<-as(g,"gatingTemplate")
	gt@name<-name
	gt
}

