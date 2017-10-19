context("gating...")

gatingResults <- readRDS("gatingResults.rds")
# gatingResults <- readRDS("tests/testthat/gatingResults.rds")

localPath <- "~/rglab/workspace/openCyto"

test_that("tcell", {
      
      gt_tcell <- gatingTemplate(gtFile, autostart = 1L)
      
      gs <- load_gs(file.path(localPath,"misc/testSuite/gs-tcell"))
      
      expect_warning(gating(gt_tcell, gs, mc.core = 2, parallel_type = "multicore"), regexp = "HLA is partially matched")
      
      #test toggle helperGates
      
      expect_equal(length(getNodes(gs)), 29)
      helperGates <- get.helperGates(gt_tcell, gs)
      expect_true(setequal(helperGates, c('/nonDebris/singlets/lymph/cd3/cd4+',
                                          '/nonDebris/singlets/lymph/cd3/cd8+',
                                          '/nonDebris/singlets/lymph/cd3/cd4+cd8-/CD45_neg',
                                          '/nonDebris/singlets/lymph/cd3/cd4+cd8-/CD45_neg/CCR7_gate',
                                          '/nonDebris/singlets/lymph/cd3/cd4+cd8-/CD38+',
                                          '/nonDebris/singlets/lymph/cd3/cd4+cd8-/HLA+',
                                          '/nonDebris/singlets/lymph/cd3/cd4-cd8+/CCR7+',
                                          '/nonDebris/singlets/lymph/cd3/cd4-cd8+/CD45RA+',
                                          '/nonDebris/singlets/lymph/cd3/cd4-cd8+/CD38+',
                                          '/nonDebris/singlets/lymph/cd3/cd4-cd8+/HLA+')))
      #hide
      toggle.helperGates(gt_tcell, gs)
      expect_equal(length(getNodes(gs)), 19)
      expect_equal(length(getNodes(gs, showHidden = TRUE)), 29)
      #recover
      toggle.helperGates(gt_tcell, gs)
      expect_equal(length(getNodes(gs)), 29)
    
      #rm helper gates
      gs1 <- clone(gs,isNew = FALSE)
      delete.helperGates(gt_tcell, gs1)
      expect_equal(length(getNodes(gs1, showHidden = TRUE)), 19)
      
      
      thisRes <- getPopStats(gs, path = "full", format = "wide")
      expectRes <- gatingResults[["gating_tcell"]]
      expect_equal(thisRes, expectRes, tol = 0.04)
      
      #test the interactive gating API
      nodes <- getChildren(gs[[1]], "cd4-cd8+")
      for(node in nodes)
        Rm(node, gs)
      opt <- getOption("openCyto")
      opt[["check.pop"]] <- FALSE
      options(openCyto = opt)
      expect_warning(add_pop(gs, gating_method = "tailgate", dims = "CD38,HLA"
              , parent = "cd4-cd8+", pop = "CD38+HLA+"
              , alias = "activated cd8", preprocessing_method = "standardize_flowset")
        , regexp = "HLA is partially matched")
      opt[["check.pop"]] <- TRUE
      options(openCyto = opt)
      #test new .mindensity2 wrapper
      add_pop(gs, gating_method = "mindensity2", dims = "CCR7,CD45RA", parent = "cd4-cd8+", pop = "+/-+/-")
      thisRes <- getPopStats(gs, path = "full", format = "wide")
      expect_equal(thisRes, expectRes, tol = 0.04)
      
      expect_true(all(c("cd4-cd8+/CD38+", "cd4-cd8+/HLA+") %in% getNodes(gs, path = "auto")))
      
     
      #use channel in pop and stains in dims
      for(node in nodes[7:9])
        Rm(node, gs)
      expect_warning(add_pop(gs, gating_method = "tailgate"
                      , dims = "CD38,HLA"
                      , parent = "cd4-cd8+"
                      , pop = "-+"
                      , alias = "activated cd8"
                      , preprocessing_method = "standardize_flowset"
                        )      
      , regexp = "HLA is partially matched")
      thisRes <- getPopStats(gs, path = "full", format = "wide")
      expect_equal(thisRes, expectRes, tol = 0.04)
      
      #pure +-
      for(node in nodes[7:9])
        Rm(node, gs)
      expect_warning(
        add_pop(gs, gating_method = "tailgate"
              , dims = "CD38,HLA"
              , parent = "cd4-cd8+"
              , pop = "-+"
              , alias = "activated cd8"
              , preprocessing_method = "standardize_flowset"
            )      
      , regexp = "HLA is partially matched")
      thisRes <- getPopStats(gs, path = "full", format = "wide")
      expect_equal(thisRes, expectRes, tol = 0.04)
      
      
      #test keep.helperGates = FALSE
      nodes <- getChildren(gs[[1]], "cd4-cd8+")
      for(node in nodes)
        Rm(node, gs)
      expect_warning(
        add_pop(gs, gating_method = "tailgate", dims = "CD38,HLA"
                , parent = "cd4-cd8+", pop = "++"
                , alias = "activated cd8", keep.helperGates = FALSE)
        , regexp = "HLA is partially matched")
      expect_false(any(c("cd4-cd8+/CD38+", "cd4-cd8+/HLA+") %in% getNodes(gs, path = "auto")))
      expect_equal(getChildren(gs[[1]], "cd4-cd8+", path = "auto"), "activated cd8")
      
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
  expect_warning(gating(gt_tcell, gs, mc.core = 2, parallel_type = "multicore"), regexp = "HLA is partially matched")
  # autoplot(gs[[1]])  
  
  thisRes <- getPopStats(gs, path = "full")
  expectRes <- gatingResults[["gating_tcell_asinhtGml2"]]
  expect_equal(thisRes, expectRes, tol = 0.006)

})

test_that("ICS", {
      
      
      gtfile <- system.file("extdata/gating_template/ICS.csv", package = "openCyto")
      gt <- gatingTemplate(gtfile)
      
      
      gs <- load_gs(file.path(localPath,"misc/testSuite/gs-ICS"))
      Rm("s", gs)
      expect_warning(gating(gt, gs, mc.core = 2, parallel_type = "multicore"), regexp = "Pacific Blue-A is partially matched")
      
      thisRes <- getPopStats(gs, path = "full", format = "wide")
      expectRes <- gatingResults[["gating_ICS"]]
      expect_equal(thisRes, expectRes, tol = 0.05)
      
      #test add_pop with polyFunctions
      nodes <- getChildren(gs[[1]], "cd8")[-(1:4)]
      for(node in nodes)
        Rm(node, gs)
      expect_warning(
        add_pop(gs, gating_method = "polyFunctions", parent = "cd8", gating_args = "cd8/IFNg:cd8/IL2:cd8/TNFa")
      , regexp = "is replaced with")
      
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

