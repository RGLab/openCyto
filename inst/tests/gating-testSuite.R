flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
gs <- load_gs(file.path(flowDataPath,"gs_manual"))

fcsFiles <- list.files(pattern = "CytoTrol", flowDataPath, full = TRUE)
ncfs  <- read.ncdfFlowSet(fcsFiles)
compMat <- getCompensationMatrices(gs[[1]])
ncfs_comp <- compensate(ncfs, compMat)


chnls <- parameters(compMat)
transFuncts <- estimateLogicle(ncfs[[1]], channels = chnls)
ncfs_trans <- transform(ncfs_comp, transFuncts)

gs <- GatingSet(ncfs_trans)

gt_tcell <- gatingTemplate(gtFile)
gating(gt_tcell, gs)
