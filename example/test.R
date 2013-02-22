# 
# ICS
# 
###############################################################################
unloadNamespace("openCyto")
library(flowWorkspace)
library(openCyto)
source("/home/wjiang2/rglab/workspace/openCyto/R/AllGenerics.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/AllClasses.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gatingTemplate-methods.R")
lapply(list.files("/home/wjiang2/rglab/workspace/openCyto/R",full=T),source)

path<-"openCyto"
load(file.path(path,"data/065_fs.rda"))

gt<-gatingTemplate(file.path(path,"data/ICS_GatingTemplate.csv"),"ICS")
gt
#getNodes(gt,"14")
getChildren(obj=gt,y="2")
getGate(gt,"14","17")
#getNodes(gt)
#png("openCyto/gatingTemplate.png")
plot(gt)
#dev.off()

##transform the ICS data
paramters<-colnames(fs[[1]])
trans <- estimateLogicle(fs[[1]], channels = paramters[!grepl("[F|S]SC|[T|t]ime",paramters)])
fs_trans<-transform(fs,trans)
gs<-GatingSet(fs_trans)
#gatingTemplate(gs)<-gt
#gatingTemplate(gs)
env1<-new.env(parent=emptyenv())
gating(gt,gs,env1)
plot(gs[[1]])
getNodes(gs[[1]])
plotGate(gs[[1]],xbin=128)
getData(gs1)[[1]]

#plot priors
plot(env1$fct,"v",posteriors=T)
plot(env1$fct,"nonDebris",posteriors=T)

#Tcell is already transformed
load(file.path(path,"data/fs_tcell.rda"))

gt1<-gatingTemplate(file.path(path,"data/Cytotrol_Tcell_GatingTemplate.csv"),"Tcell")
plot(gt1)

fs_tcell[[1]]
gs1<-GatingSet(fs_tcell)
gating(gt1,gs1)
plotGate(gs1[[1]],xbin=128)
#xyplot(`<B710-A>`~`<R660-A>`,fs_tcell)
#densityplot(~.,fs_tcell[[1]])

#Bcell
load(file.path(path,"data/fs_bcell.rda"))

gt2<-gatingTemplate(file.path(path,"data/Cytotrol_Bcell_GatingTemplate.csv"),"Bcell")
plot(gt2)

fs_bcell[[1]]
gs2<-GatingSet(fs_bcell)
gating(gt2,gs2)
plotGate(gs2[[1]],bool=T,xbin=128)

#
library(openCyto)
archive_path <- '/loc/no-backup/ramey'

gt <- gatingTemplate(file.path(archive_path, "HVTN065-GatingTemplate.csv"), "HVTN065")
gs <- unarchive(file = file.path(archive_path, "test-HVTN065.tar"), archive_path)
gating(gt, gs)
plotGate(gs[[1]],bool=T,xbin=128,margin=T)

