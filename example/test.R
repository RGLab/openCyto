# 
# ICS
# 
###############################################################################
unloadNamespace("flowIncubator")
unloadNamespace("openCyto")
unloadNamespace("QUALIFIER")
unloadNamespace("flowStats")
unloadNamespace("flowWorkspace")


library(openCyto)
library(flowIncubator)
library(flowWorkspace)
library(flowClust)
library(flowStats)
library(MASS)
library(plyr)
library(clue)
library(gtools)

source("/home/wjiang2/rglab/workspace/openCyto/R/AllClasses.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gatingTemplate-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gating-cytokines.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gating-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gating-functions.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gtMethod-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gtPopulation-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/fcTree-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/fcFilterList-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/fcFilter-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/bayes-flowClust.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/functions.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/bayes-flowClust.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/wrapper-functions.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/preprocessing-method.R")

path<-"/home/wjiang2/rglab/workspace/openCyto"

###############
##ICS
###############

gt<-gatingTemplate(file.path(path,"data/ICS_expanded.csv"),"ICS")
gt<-gatingTemplate(file.path(path,"data/ICS.csv"),"ICS")
gt
#getNodes(gt,"14")
#getChildren(obj=gt,y="2")
#getGate(gt,"14","17")
#getNodes(gt)
#png("openCyto/gatingTemplate.png")
plot(gt)

#dev.off()

##transform the ICS data
load(file.path(path,"data/065_fs.rda"))
paramters<-colnames(fs[[1]])
trans <- estimateLogicle(fs[[1]], channels = paramters[!grepl("[F|S]SC|[T|t]ime",paramters)])
fs_trans<-transform(fs,trans)

gs<-GatingSet(fs_trans)
env1<-new.env(parent=emptyenv())
gating(gt,gs,env1, prior_group='VISITNO')
plot(gs[[1]])
plot(gs[[1]],bool=T)
getNodes(gs[[1]])
plotGate(gs[[1]],xbin=64,margin=T)
plotGate(gs[[1]],c(8,22),xbin=128,margin=T,digits=3)
getPopStats(gs[[1]])[22,]

#plot priors
plotGate(x=gs,3)
X11()
plot(env1$fct,"v",posteriors=T)
plotGate(x=gs,4)
plot(env1$fct,"nonDebris",posteriors=T)
plot(env1$fct,"cd4",posteriors=T,channel="PE Cy55-A")
plot(env1$fct,"cd4",posteriors=T,channel="FITC-A")
plot(env1$fct,"TNFa",posteriors=T)

###############
#Tcell is already transformed
###############


gt1<-gatingTemplate(file.path(path,"data/Cytotrol_Tcell_expanded.csv"),"Tcell")
gt1<-gatingTemplate(file.path(path,"data/Cytotrol_Tcell.csv"),"Tcell")
plot(gt1)

load(file.path(path,"data/fs_tcell.rda"))
gs1<-GatingSet(fs_tcell)
env1<-new.env(parent=emptyenv())
gating(gt1,gs1,env1, num_cores=4,parallel_type="MPI")
plotGate(gs1[[1]],xbin=64)
plot(env1$fct,"nonDebris",post=T)
plot(env1$fct,"cd3",post=T)
plot(env1$fct,"cd4+",post=T,channel="<B710-A>")
plot(env1$fct,"activated cd4",post=T,channel="<R660-A>")
plot(env1$fct,2,post=T)
plot(env1$fct,"lymph",post=T)

#xyplot(`<B710-A>`~`<R660-A>`,fs_tcell)
#densityplot(~.,fs_tcell[[1]])

###############
#Bcell is already transformed
###############


gt2<-gatingTemplate(file.path(path,"data/Cytotrol_Bcell_expanded.csv"),"Bcell")
gt2<-gatingTemplate(file.path(path,"data/Cytotrol_Bcell.csv"),"Bcell")
plot(gt2)

getNodes(gs2[[1]])
load(file.path(path,"data/fs_bcell.rda"))
gs2<-GatingSet(fs_bcell)
env1<-new.env(parent=emptyenv())
gating(gt2,gs2,env1)
getGate(gs2,4)
plot(gs2[[1]],bool=T)
plotGate(gs2[[1]],bool=T,xbin=64)
getData(gs2[[1]],"cd19&cd20")
getNodes(gs2[[1]])
densityplot(~`<G780-A>`,getData(gs2,"cd19&!cd20"))

getGate(env1$fct,"nonDebris")
plot(env1$fct,"nonDebris",post=T)
plot(env1$fct,"cd19",post=T)
plot(env1$fct,"IgD-cd27+",channel="<G780-A>",post=T)
plot(env1$fct,"IgD-cd27+",post=T,channel="<V545-A>")

###debug
library(openCyto)
gs_HVTN065 <- load_gs("/loc/no-backup/ramey/HVTN/065/gating-results")
gating_template <- gatingTemplate("/home/jramey/rglab/papers/paper-opencyto/gt-HVTN065.csv", "HVTN065")
Rm("cd4", gs_HVTN065)
Rm("cd8", gs_HVTN065)
gs <- clone(gs_HVTN065[1:12])
getData(gs)
Rm("cd4", gs)
Rm("cd8", gs)

gating(gating_template, gs, prior_group = 'Stim'
    , num_cores = 6, parallel_type = "MPI"
)

