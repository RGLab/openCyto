# 
# ICS
# 
###############################################################################
unloadNamespace("openCyto")
library(flowWorkspace)
#library(openCyto)
source("/home/wjiang2/rglab/workspace/openCyto/R/AllClasses.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gatingTemplate-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gating-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gating-functions.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gtMethod-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gtPopulation-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/fcTree-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/fcFilterList-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/fcFilter-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/bayes-flowClust.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/median-logicle-transform.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/functions.R")




path<-"/home/wjiang2/rglab/workspace/openCyto"
load(file.path(path,"data/065_fs.rda"))

gt<-gatingTemplate(file.path(path,"data/ICS_GatingTemplate.csv"),"ICS")
gt
#getNodes(gt,"14")
getChildren(obj=gt,y="2")
getGate(gt,"14","17")
#getNodes(gt)
#png("openCyto/gatingTemplate.png")
plot(gt)
class(gt)
#dev.off()

##transform the ICS data
paramters<-colnames(fs[[1]])
trans <- estimateLogicle(fs[[1]], channels = paramters[!grepl("[F|S]SC|[T|t]ime",paramters)])
fs_trans<-transform(fs,trans)
gs<-GatingSet(fs_trans)
env1<-new.env(parent=emptyenv())
gating(gt,gs,env1, prior_group='VISITNO')
plot(gs[[1]])
plot(gs[[1]],bool=T)
getNodes(gs[[1]])
plotGate(gs[[1]],xbin=128,margin=T)


#plot priors
plotGate(x=gs,3)
X11()
plot(env1$fct,"v",posteriors=T)
plotGate(x=gs,4)
plot(env1$fct,"nonDebris",posteriors=T)
plot(env1$fct,"cd4",posteriors=T,channel="PE Cy55-A")
plot(env1$fct,"cd4",posteriors=T,channel="FITC-A")
plot(env1$fct,"TNFa",posteriors=T)
#Tcell is already transformed
load(file.path(path,"data/fs_tcell.rda"))

gt1<-gatingTemplate(file.path(path,"data/Cytotrol_Tcell_GatingTemplate.csv"),"Tcell")
plot(gt1)

fs_tcell[[1]]
gs1<-GatingSet(fs_tcell)
env1<-new.env(parent=emptyenv())
gating(gt1,gs1,env1)
plotGate(gs1[[1]],xbin=128)
plot(env1$fct,"nonDebris",post=T)
plot(env1$fct,"cd3",post=T)
plot(env1$fct,"cd4",post=T,channel="<B710-A>")
plot(env1$fct,"activated cd4",post=T,channel="<R660-A>")
plot(env1$fct,2,post=T)
plot(env1$fct,"lymph",post=T)

#xyplot(`<B710-A>`~`<R660-A>`,fs_tcell)
#densityplot(~.,fs_tcell[[1]])

#Bcell
load(file.path(path,"data/fs_bcell.rda"))

gt2<-gatingTemplate(file.path(path,"data/Cytotrol_Bcell_GatingTemplate.csv"),"Bcell")
plot(gt2)

fs_bcell[[1]]
gs2<-GatingSet(fs_bcell)
env1<-new.env(parent=emptyenv())
gating(gt2,gs2,env1)
plot(gs2[[1]],bool=T)
plotGate(gs2[[1]],bool=T,xbin=64)


plot(env1$fct,"nonDebris",post=T)
plot(env1$fct,"cd19",post=T)
plot(env1$fct,"IgD-cd27+",channel="<G780-A>",post=T)
plot(env1$fct,"IgD-cd27+",post=T,channel="<V545-A>")


