context("refgate")

test_that("refGate", {
  dataDir <- system.file("extdata",package="flowWorkspaceData")
  gs <- load_gs(list.files(dataDir, pattern = "gs_manual",full = TRUE))
  #clear existing gates
  Rm("not debris", gs)
  
  #add 1d gates
  add_pop(gs, alias = "*", pop = "+/-", parent = "root", dims = "CD4", gating_method = "mindensity")
  add_pop(gs, alias = "CD8", pop = "+", parent = "root", dims = "CD8", gating_method = "mindensity")
  #add quadgates through refGates
  add_pop(gs, pop = "+/-+/-", parent = "root", dims = "CD4,CD8", gating_method = "refGate", gating_args = "CD4+:CD8")
  nodes <- getNodes(gs)[5:8]
  stats1 <- getPopStats(gs[[1]], subpopulations = nodes)
  for(node in nodes)
    Rm(node, gs)
  #refer to cd4-
  add_pop(gs, pop = "+/-+/-", parent = "root", dims = "CD4,CD8", gating_method = "refGate", gating_args = "CD4-:CD8")
  stats2 <- getPopStats(gs[[1]], subpopulations = nodes)
  
  expect_equal(stats1, stats2)
  
})
