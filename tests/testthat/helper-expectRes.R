library(data.table)
library(utils)
library(tools)
library(flowCore)
library(ncdfFlow)
gtFile <- system.file("extdata/gating_template/tcell.csv", package = "openCyto")
expectResults <- readRDS("expectResults.rds")
dataDir <- system.file("extdata",package="flowWorkspaceData")
gs_dir <- list.files(dataDir, pattern = "gs_manual",full = TRUE)
load_gs_local <- function(gs_dir, ...)
{
  
  
  # set_default_backend("tile")
  if(get_default_backend()=="tile")
  {
    tmp <- tempfile()
    convert_backend(gs_dir, tmp)
    if(use_on_disk_idx())
    {
      gs <- load_gs(tmp)
      gs_convert_idx_to_ondisk(gs)
      expect_equal(basename(gh_idx_get_uri(gs[[1]])), paste0(sampleNames(gs[[1]]), ".idx"))
      tmp1 <- tempfile()
      save_gs(gs, tmp1)
      unlink(tmp, recursive = TRUE)
      tmp <- tmp1
    }
    gs_dir <- tmp
  }
  gs <- load_gs(gs_dir, ...)
  gs
}
gs <- load_gs_local(gs_dir)
data("GvHD")

