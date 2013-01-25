# TODO: Add comment
# 
# Author: wjiang2
###############################################################################
## exsits for the purpose of dispatching 
setClass("gatingTemplate",
		contains="data.frame"
		,representation(name="character")
		,validity=function(object){
			#check if the parent is ambiguous or not valid			
			parent_matches<-sapply(levels(object[,"parent"]),function(x){
						matches<-length(grep(x,object[,"pop"],fix=T))
						if(matches==0&&x=="root")
							matches<-1
						matches
					})
			
			invalid_parent<-which(parent_matches==0)
			ambiguous_parent<-which(parent_matches>1)
			if(length(invalid_parent)>0)
				return(paste("The following parent population names are not valid:"
							,paste(names(invalid_parent)
								,collapse=",")
							)
						)
			else if(length(ambiguous_parent)>0)
				return(paste("The following parent population names are ambiguous:"
								,paste(names(ambiguous_parent)
										,collapse=",")
						)
				)
			else 
				TRUE
		}
	)



gatingTemplate<-function(file,name){
	gt<-read.csv(file)
	gt<-as(gt,"gatingTemplate")
	gt@name<-name
	gt
}