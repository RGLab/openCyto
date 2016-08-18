context("gating...")

gatingResults <- readRDS("gatingResults.rds")

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
      add_pop(gs, gating_method = "tailgate", dims = "CD38,HLA"
              , parent = "cd4-cd8+", pop = "CD38+HLA+"
              , alias = "activated cd8", preprocessing_method = "standardize_flowset")
      #test new .mindensity2 wrapper
      add_pop(gs, gating_method = "mindensity2", dims = "CCR7,CD45RA", parent = "cd4-cd8+", pop = "CCR7+/-CD45RA+/-")
      thisRes <- getPopStats(gs, path = "full", format = "wide")
      expect_equal(thisRes, expectRes, tol = 0.04)
      
      #use channel in pop and stains in dims
      for(node in nodes[7:9])
        Rm(node, gs)
      add_pop(gs, gating_method = "tailgate"
              , dims = "CD38,HLA"
              , parent = "cd4-cd8+"
              , pop = "R660-HLA+"
              , alias = "activated cd8"
              , preprocessing_method = "standardize_flowset"
      )      
      thisRes <- getPopStats(gs, path = "full", format = "wide")
      expect_equal(thisRes, expectRes, tol = 0.04)
      
      #pure +-
      for(node in nodes[7:9])
        Rm(node, gs)
      add_pop(gs, gating_method = "tailgate"
              , dims = "CD38,HLA"
              , parent = "cd4-cd8+"
              , pop = "-+"
              , alias = "activated cd8"
              , preprocessing_method = "standardize_flowset"
      )      
      thisRes <- getPopStats(gs, path = "full", format = "wide")
      expect_equal(thisRes, expectRes, tol = 0.04)
                    
    })

test_that("tcell--asinhtGml2", {
  
  fcs <- list.files(pattern = "CytoTrol*", system.file("extdata", package = "flowWorkspaceData"), full.names = TRUE)
  fs <- read.ncdfFlowSet(fcs)  
  gs <- GatingSet(fs)
  #compensate
  comp <- compensation(spillover(fs[[1]])[["SPILL"]])
  gs <- compensate(gs, comp)
  #transform
  trans <- asinhtGml2_trans()
  chnls <- as.vector(parameters(comp))
  trans.list <- transformerList(chnls, trans)
  gs <- transform(gs, trans.list)
  
  # autoplot(getData(gs), chnls[1])
  
  #load gating template
  localPath <- "~/rglab/workspace/openCyto"
  gtFile <- system.file("extdata/gating_template/tcell.csv", package = "openCyto")
  #modify scale-specific parameters
  dt <- fread(gtFile)
  dt[c(5, 8), gating_args := NA]
  tmp <- tempfile()
  write.csv(dt, file = tmp, row.names = FALSE)
  gt_tcell <- gatingTemplate(tmp, autostart = 1L)
  gating(gt_tcell, gs, mc.core = 2, parallel_type = "multicore")
  # autoplot(gs[[1]])  
  
  thisRes <- getPopStats(gs, path = "full")
  expectRes <- gatingResults[["gating_tcell_asinhtGml2"]]
  expect_equal(thisRes, expectRes, tol = 0.005)

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
      
      #test add_pop with polyFunctions
      nodes <- getChildren(gs[[1]], "cd8")[-(1:4)]
      for(node in nodes)
        Rm(node, gs)
      add_pop(gs, gating_method = "polyFunctions", parent = "cd8", gating_args = "cd8/IFNg:cd8/IL2:cd8/TNFa")
      
      thisRes <- getPopStats(gs, path = "full", format = "wide")
      expect_equal(thisRes, expectRes, tol = 0.04)
      
      #test add_pop with boolean method
      pop <- "IL2orIFNg"
      Rm(pop, gs)
      add_pop(gs, alias = pop, gating_method = "boolGate", parent = "cd4", gating_args = "cd4/IL2|cd4/IFNg")
      
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

