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
setMethod("plot",sig=c("fcTree","character"),definition=function(x,y,isAlias=TRUE,...){
			#get fcobject
#			browser()
			nodes<-getNodes(x)
			if(isAlias)#The node index is a pop alias
			{
				
				allAlias<-lapply(nodes,function(curNode)alias(curNode$pop))
				ind<-match(y,allAlias)
				if(length(ind)>1)
					stop("Population '",y,"' is ambiguous!")
				else if(length(ind)==0)
					stop("Population '",y,"' is not found!")
				else
				{
					matchedNode<-nodes[[ind]]
				}
			}else
				matchedNode<-nodes[[y]]
			
			fcObj<-matchedNode$fcObj
			plot(x=fcObj,main=matchedNode$pop@name,...)
		})

