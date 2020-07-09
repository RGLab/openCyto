library(data.table)
library(utils)
library(tools)
library(flowCore)
library(ncdfFlow)
gtFile <- system.file("extdata/gating_template/tcell.csv", package = "openCyto")
expectResults <- readRDS("expectResults.rds")
dataDir <- system.file("extdata",package="flowWorkspaceData")
gs_dir <- list.files(dataDir, pattern = "gs_manual",full = TRUE)
# set_default_backend("tile")
if(get_default_backend()=="tile")
{
  tmp <- tempfile()
  convert_backend(gs_dir, tmp)
  gs_dir <- tmp
}
gs <- load_gs(gs_dir)
data("GvHD")

