context("gating...")

gatingResults <- readRDS("gatingResults.rds")
# gatingResults <- readRDS("tests/testthat/gatingResults.rds")

localPath <- "~/rglab/workspace/openCyto"

skip_if_not(dir.exists(file.path(localPath, "misc")))
library(parallel)
register_plugins(flowStats:::.tailgate, "tailgate")
test_that("tcell", {
      
      gt_tcell <- gatingTemplate(gtFile)
      
      gs <- load_gs(file.path(localPath,"misc/testSuite/gs-tcell"))
      
      expect_warning(gating(gt_tcell, gs, mc.core = 1, parallel_type = "multicore"), regexp = "HLA is partially matched")
      
      #test toggle helperGates
      
      expect_equal(length(gs_get_pop_paths(gs)), 29)
      helperGates <- gt_get_helpergates(gt_tcell, gs)
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
      gt_toggle_helpergates(gt_tcell, gs)
      expect_equal(length(gs_get_pop_paths(gs)), 19)
      expect_equal(length(gs_get_pop_paths(gs, showHidden = TRUE)), 29)
      #recover
      gt_toggle_helpergates(gt_tcell, gs)
      expect_equal(length(gs_get_pop_paths(gs)), 29)
    
      #rm helper gates
      gs1 <- gs_copy_tree_only(gs)
      gt_delete_helpergates(gt_tcell, gs1)
      expect_equal(length(gs_get_pop_paths(gs1, showHidden = TRUE)), 19)
      
      
      thisRes <- gs_pop_get_stats(gs[1], path = "full", type = "percent")
      expectRes <- gatingResults[["gating_tcell"]]
      expect_equivalent(thisRes[,percent], expectRes[thisRes[,pop],1], tol = 0.04)
      expect_equivalent(gs_pop_get_stats(gs[2], path = "full", type = "percent")[,percent], expectRes[thisRes[,pop],2], tol = 0.04)
      #test the interactive gating API
      nodes <- gs_pop_get_children(gs[[1]], "cd4-cd8+")
      for(node in nodes)
        gs_pop_remove(gs, node)
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
      gs_add_gating_method(gs, gating_method = "gate_mindensity2", dims = "CCR7,CD45RA", parent = "cd4-cd8+", pop = "+/-+/-")
      thisRes <- gs_pop_get_stats(gs[1], path = "full", type = "percent")
      expect_equivalent(thisRes[,percent], expectRes[thisRes[,pop],1], tol = 0.04)
      
      expect_true(all(c("cd4-cd8+/CD38+", "cd4-cd8+/HLA+") %in% gs_get_pop_paths(gs, path = "auto")))
      
     
      #use channel in pop and stains in dims
      for(node in nodes[7:9])
        gs_pop_remove(gs, node)
      expect_warning(add_pop(gs, gating_method = "tailgate"
                      , dims = "CD38,HLA"
                      , parent = "cd4-cd8+"
                      , pop = "-+"
                      , alias = "activated cd8"
                      , preprocessing_method = "standardize_flowset"
                        )      
      , regexp = "HLA is partially matched")
      thisRes <- gs_pop_get_stats(gs[1], path = "full", type = "percent")
      expect_equivalent(thisRes[,percent], expectRes[thisRes[,pop],1], tol = 0.04)
      
      
      #pure +-
      for(node in nodes[7:9])
        gs_pop_remove(gs, node)
      expect_warning(
        add_pop(gs, gating_method = "tailgate"
              , dims = "CD38,HLA"
              , parent = "cd4-cd8+"
              , pop = "-+"
              , alias = "activated cd8"
              , preprocessing_method = "standardize_flowset"
            )      
      , regexp = "HLA is partially matched")
      thisRes <- gs_pop_get_stats(gs[1], path = "full", type = "percent")
      expect_equivalent(thisRes[,percent], expectRes[thisRes[,pop],1], tol = 0.04)
      
      
      #test keep.helperGates = FALSE
      nodes <- gs_pop_get_children(gs[[1]], "cd4-cd8+")
      for(node in nodes)
        gs_pop_remove(gs, node)
      expect_warning(
        add_pop(gs, gating_method = "tailgate", dims = "CD38,HLA"
                , parent = "cd4-cd8+", pop = "++"
                , alias = "activated cd8", keep.helperGates = FALSE)
        , regexp = "HLA is partially matched")
      expect_false(any(c("cd4-cd8+/CD38+", "cd4-cd8+/HLA+") %in% gs_get_pop_paths(gs, path = "auto")))
      expect_equal(gs_pop_get_children(gs[[1]], "cd4-cd8+", path = "auto"), "activated cd8")
      
    })

test_that("tcell--asinhtGml2", {
  
  fcs <- list.files(pattern = "CytoTrol*", system.file("extdata", package = "flowWorkspaceData"), full.names = TRUE)
  fs <- load_cytoset_from_fcs(fcs)  
  gs <- GatingSet(fs)
  #compensate
  comp <- compensation(spillover(fs[[1]])[["SPILL"]])
  gs <- compensate(gs, comp)
  #transform
  trans <- asinhtGml2_trans()
  chnls <- as.vector(parameters(comp))
  trans.list <- transformerList(chnls, trans)
  gs <- transform(gs, trans.list)
  
  # autoplot(gs_pop_get_data(gs), chnls[1])
  
  #load gating template
  localPath <- "~/rglab/workspace/openCyto"
  gtFile <- system.file("extdata/gating_template/tcell.csv", package = "openCyto")
  #modify scale-specific parameters
  dt <- fread(gtFile)
  dt[c(5, 8), gating_args := NA]
  tmp <- tempfile()
  write.csv(dt, file = tmp, row.names = FALSE)
  gt_tcell <- gatingTemplate(tmp)
  library(parallel)
  expect_warning(gating(gt_tcell, gs, mc.core = 1, parallel_type = "multicore"), regexp = "HLA is partially matched")
  # autoplot(gs[[1]])  
  
  thisRes <- gs_pop_get_count_fast(gs, path = "full")
  expectRes <- gatingResults[["gating_tcell_asinhtGml2"]]
  expect_equal(thisRes, expectRes, tol = 0.018)

})

test_that("ICS", {
      
      
      gtfile <- system.file("extdata/gating_template/ICS.csv", package = "openCyto")
      gt <- gatingTemplate(gtfile)
      
      
      gs <- load_gs(file.path(localPath,"misc/testSuite/gs-ICS"))[1]
	  gs_pop_remove(gs, "s")
	  #TODO:investigate segfault associated with multicore and L#125 asinhtGml2_trans
      
      expect_warning(gating(gt, gs
                            # , mc.core = 2
                            # , parallel_type = "multicore"
                            ), regexp = "Pacific Blue-A is partially matched")
      
      expectRes <- gatingResults[["gating_ICS"]]
      thisRes <- gs_pop_get_stats(gs[1], path = "full", type = "percent")
      expect_equivalent(thisRes[,percent], expectRes[thisRes[,pop],1], tol = 0.04)
      
      #test add_pop with polyFunctions
      
      #test add_pop with boolean method
      pop <- "IL2orIFNg"
      gs_pop_remove(gs, pop)
      gs_add_gating_method(gs, alias = pop, gating_method = "boolGate", parent = "cd4", gating_args = "cd4/IL2|cd4/IFNg")
      
      thisRes <- gs_pop_get_stats(gs[1], path = "full", type = "percent")
      expect_equivalent(thisRes[,percent], expectRes[thisRes[,pop],1], tol = 0.04)
      
    })

test_that("treg", {
      
      skip("TO investigate the failure")
      gtfile <- system.file("extdata/gating_template/treg.csv", package = "openCyto")
      gt <- gatingTemplate(gtfile)
      
      gs <- load_gs(file.path(localPath,"misc/testSuite/gs-treg"))
      gs_pop_remove(gs, "boundary")
      #TODO: fix mc.core=3 error
      expect_warning(gating(gt, gs, mc.core = 1, parallel_type = "multicore"), "did not converge")
      
      expectRes <- gatingResults[["gating_treg"]]
      thisRes <- gs_pop_get_stats(gs[1], path = "full", type = "percent")
      expect_equivalent(thisRes[,percent], expectRes[thisRes[,pop],1], tol = 0.025)
      
    })

test_that("bcell", {
      
      
      gtfile <- system.file("extdata/gating_template/bcell.csv", package = "openCyto")
      gt <- gatingTemplate(gtfile)
      
      gs <- load_gs(path = file.path(localPath,"misc/testSuite/gs-bcell"))
      gs_pop_remove(gs, "boundary")
      expect_warning(gating(gt, gs, mc.core = 1, parallel_type = "multicore"), "did not converge")
      
      expectRes <- gatingResults[["gating_bcell"]]
      thisRes <- gs_pop_get_stats(gs[1], path = "full", type = "percent")
      expect_equivalent(thisRes[,percent], expectRes[thisRes[,pop],1], tol = 0.08)
      
    })

