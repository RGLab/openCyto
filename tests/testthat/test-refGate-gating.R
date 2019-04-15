context("refgate")

test_that("refGate", {
  dataDir <- system.file("extdata",package="flowWorkspaceData")
  gs <- load_gs(list.files(dataDir, pattern = "gs_manual",full = TRUE))
  #clear existing gates
  gs_remove_gate("not debris", gs)
  
  #add 1d gates
  gs_add_pop(gs, alias = "*", pop = "+/-", parent = "root", dims = "CD4", gating_method = "mindensity")
  gs_add_pop(gs, alias = "CD8", pop = "+", parent = "root", dims = "CD8", gating_method = "mindensity")
  #add quadgates through refGates
  gs_add_pop(gs, pop = "+/-+/-", parent = "root", dims = "CD4,CD8", gating_method = "refGate", gating_args = "CD4+:CD8")
  nodes <- gs_get_pop_paths(gs)[5:8]
  stats1 <- gh_get_pop_stats(gs[[1]], subpopulations = nodes)
  for(node in nodes)
    gs_remove_gate(node, gs)
  #refer to cd4-
  gs_add_pop(gs, pop = "+/-+/-", parent = "root", dims = "CD4,CD8", gating_method = "refGate", gating_args = "CD4-:CD8")
  stats2 <- gh_get_pop_stats(gs[[1]], subpopulations = nodes)
  
  expect_equal(stats1, stats2)
  
})
