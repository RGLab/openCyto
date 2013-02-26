# TODO: Add comment
# 
# Author: wjiang2
###############################################################################

setMethod("getNodes",signature=c("fcTree"),definition=function(x,y)
		{

			if(missing(y))
				nodeData(x)
			else
				nodeData(x,y)
		})
setMethod("plot",sig=c("fcTree","character"),definition=function(x,y,channel=NULL,...){
			#get fcobject
#			browser()
			nodes<-getNodes(x)
			
				
			allAlias<-lapply(nodes,function(curNode)alias(curNode$pop))
			ind<-which(allAlias%in%y)
			if(length(ind)>1)
				stop("Population '",y,"' is ambiguous!")
			else if(is.na(ind))
				stop("Population '",y,"' is not found!")
			else
			{
				matchedNode<-nodes[[ind]]
			}
		
			
			fcObj<-matchedNode$fcObj
			
			plot(x=fcObj,y=channel,main=matchedNode$pop@name,...)
		})

setMethod("plot",sig=c("fcTree","numeric"),definition=function(x,y,channel=NULL,...){
#			browser()
			y<-as.character(y)
			nodes<-getNodes(x)
			matchedNode<-nodes[[y]]
			fcObj<-matchedNode$fcObj
			
			plot(x=fcObj,y=channel,main=matchedNode$pop@name,...)
		})