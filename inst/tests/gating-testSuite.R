context("gating...")

gatingResults <- readRDS(system.file("tests/gatingResults.rds", package = "openCyto"))

localPath <- "~/rglab/workspace/openCyto"

test_that("tcell", {
      
      gt_tcell <- gatingTemplate(gtFile, autostart = 1L)
      
      gs <- load_gs(file.path(localPath,"misc/testSuite/gs-tcell"))
      
      gating(gt_tcell, gs, mc.core = 2, parallel_type = "multicore")
      
      thisRes <- getPopStats(gs, path = "full", format = "wide")
      expectRes <- gatingResults[["gating_tcell"]]
      expect_equal(thisRes, expectRes, tol = 0.04)
      
      #test the interactive gating API
      nodes <- getChildren(gs[[1]], "cd4-cd8+")
      for(node in nodes)
        Rm(node, gs)
      add_pop(gs, gating_method = "tailgate", dims = "CD38,HLA", parent = "cd4-cd8+", pop = "CD38+HLA+", alias = "activated cd8", preprocessing_method = "standardize_flowset")
      add_pop(gs, gating_method = "mindensity", dims = "CCR7,CD45RA", parent = "cd4-cd8+", pop = "CCR7+/-CD45RA+/-")
      thisRes <- getPopStats(gs, path = "full", format = "wide")
      expect_equal(thisRes, expectRes, tol = 0.04)
      
    })

test_that("ICS", {
      
      
      gtfile <- system.file("extdata/gating_template/ICS.csv", package = "openCyto")
      gt <- gatingTemplate(gtfile)
      
      
      gs <- load_gs(file.path(localPath,"misc/testSuite/gs-ICS"))
      Rm("s", gs)
      gating(gt, gs, mc.core = 2, parallel_type = "multicore")
      
      thisRes <- getPopStats(gs, path = "full", format = "wide")
      expectRes <- gatingResults[["gating_ICS"]]
      expect_equal(thisRes, expectRes, tol = 0.05)
      
      #test add_pop
      nodes <- getChildren(gs[[1]], "cd8")[-(1:4)]
      for(node in nodes)
        Rm(node, gs)
      add_pop(gs, gating_method = "polyFunctions", parent = "cd8", gating_args = "cd8/IFNg:cd8/IL2:cd8/TNFa")
      
      thisRes <- getPopStats(gs, path = "full", format = "wide")
      expect_equal(thisRes, expectRes, tol = 0.04)
      
    })

test_that("treg", {
      
      
      gtfile <- system.file("extdata/gating_template/treg.csv", package = "openCyto")
      gt <- gatingTemplate(gtfile)
      
      gs <- load_gs(file.path(localPath,"misc/testSuite/gs-treg"))
      Rm("boundary", gs)
      gating(gt, gs, mc.core = 3, parallel_type = "multicore")
      
      thisRes <- getPopStats(gs, path = "full", format = "wide")
      expectRes <- gatingResults[["gating_treg"]]
      expect_equal(thisRes, expectRes, tol = 0.25)
      
    })

test_that("bcell", {
      
      
      gtfile <- system.file("extdata/gating_template/bcell.csv", package = "openCyto")
      gt <- gatingTemplate(gtfile, autostart = 1L)
      
      gs <- load_gs(path = file.path(localPath,"misc/testSuite/gs-bcell"))
      Rm("boundary", gs)
      gating(gt, gs, mc.core = 3, parallel_type = "multicore")
      
      thisRes <- getPopStats(gs, path = "full", format = "wide")
      expectRes <- gatingResults[["gating_bcell"]]
      expect_equal(thisRes, expectRes, tol = 0.08)
      
    })

