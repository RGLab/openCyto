context("gating...")

gatingResults <- readRDS(system.file("tests/gatingResults.rds", package = "openCyto"))

localPath <- "~/rglab/workspace/openCyto"

test_that("tcell", {
      
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
      
      thisRes <- getPopStats(gs)
      expectRes <- gatingResults[["gating_tcell"]]
      expect_equal(thisRes, expectRes, tol = 0.006)
      
    })

test_that("ICS", {
      
      
      gtfile <- system.file("extdata/gating_template/ICS.csv", package = "openCyto")
      gt <- gatingTemplate(gtfile)
      
      load(file.path(localPath,"misc/testSuite/fs_080.rda"))
      gs <- GatingSet(fs[1:2])
      gating(gt, gs, mc.core = 2, parallel_type = "multicore")
      
      thisRes <- getPopStats(gs)
      expectRes <- gatingResults[["gating_ICS"]]
      expect_equal(thisRes, expectRes, tol = 0.006)
      
    })

test_that("treg", {
      
      
      gtfile <- system.file("extdata/gating_template/treg.csv", package = "openCyto")
      gt <- gatingTemplate(gtfile)
      
      gs <- load_gs(file.path(localPath,"misc/testSuite/gs-treg"))
      Rm("boundary", gs)
      gating(gt, gs, mc.core = 3, parallel_type = "multicore")
      
      thisRes <- getPopStats(gs)
      expectRes <- gatingResults[["gating_treg"]]
      expect_equal(thisRes, expectRes, tol = 0.006)
      
    })

test_that("bcell", {
      
      
      gtfile <- system.file("extdata/gating_template/bcell.csv", package = "openCyto")
      gt <- gatingTemplate(gtfile)
      
      gs <- load_gs(path = file.path(localPath,"misc/testSuite/gs-bcell"))
      Rm("boundary", gs)
      gating(gt, gs, mc.core = 3, parallel_type = "multicore")
      
      thisRes <- getPopStats(gs)
      expectRes <- gatingResults[["gating_bcell"]]
      expect_equal(thisRes, expectRes, tol = 0.006)
      
    })

