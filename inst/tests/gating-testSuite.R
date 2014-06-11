context("gating...")

gatingResults <- readRDS(system.file("tests/gatingResults.rds", package = "openCyto"))

localPath <- "~/rglab/workspace/openCyto"

test_that("tcell", {
      
      gt_tcell <- gatingTemplate(gtFile, autostart = 1L)
      
      gs <- load_gs(file.path(localPath,"misc/testSuite/gs-tcell"))
      
      gating(gt_tcell, gs, mc.core = 2, parallel_type = "multicore")
      
      thisRes <- getPopStats(gs, path = "full")
      expectRes <- gatingResults[["gating_tcell"]]
      expect_equal(thisRes, expectRes, tol = 0.04)
       
    })

test_that("ICS", {
      
      
      gtfile <- system.file("extdata/gating_template/ICS.csv", package = "openCyto")
      gt <- gatingTemplate(gtfile)
      
      
      gs <- load_gs(file.path(localPath,"misc/testSuite/gs-ICS"))
      Rm("s", gs)
      gating(gt, gs, mc.core = 2, parallel_type = "multicore")
      
      thisRes <- getPopStats(gs, path = "full")
      expectRes <- gatingResults[["gating_ICS"]]
      expect_equal(thisRes, expectRes, tol = 0.05)
      
    })

test_that("treg", {
      
      
      gtfile <- system.file("extdata/gating_template/treg.csv", package = "openCyto")
      gt <- gatingTemplate(gtfile)
      
      gs <- load_gs(file.path(localPath,"misc/testSuite/gs-treg"))
      Rm("boundary", gs)
      gating(gt, gs, mc.core = 3, parallel_type = "multicore")
      
      thisRes <- getPopStats(gs, path = "full")
      expectRes <- gatingResults[["gating_treg"]]
      expect_equal(thisRes, expectRes, tol = 0.25)
      
    })

test_that("bcell", {
      
      
      gtfile <- system.file("extdata/gating_template/bcell.csv", package = "openCyto")
      gt <- gatingTemplate(gtfile, autostart = 1L)
      
      gs <- load_gs(path = file.path(localPath,"misc/testSuite/gs-bcell"))
      Rm("boundary", gs)
      gating(gt, gs, mc.core = 3, parallel_type = "multicore")
      
      thisRes <- getPopStats(gs, path = "full")
      expectRes <- gatingResults[["gating_bcell"]]
      expect_equal(thisRes, expectRes, tol = 0.08)
      
    })

