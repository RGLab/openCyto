context("refgate")

test_that("refGate", {
  gs <- gs_copy_tree_only(gs)
  #clear existing gates
  gs_pop_remove("not debris", gs = gs)
  
  #add 1d gates
  gs_add_gating_method(gs, alias = "*", pop = "+/-", parent = "root", dims = "CD4", gating_method = "mindensity")
  gs_add_gating_method(gs, alias = "CD8", pop = "+", parent = "root", dims = "CD8", gating_method = "mindensity")
  #add quadgates through refGates
  gs_add_gating_method(gs, pop = "+/-+/-", parent = "root", dims = "CD4,CD8", gating_method = "refGate", gating_args = "CD4+:CD8")
  nodes <- gs_get_pop_paths(gs)[5:8]
  stats1 <- gh_pop_compare_stats(gs[[1]], subpopulations = nodes)
  for(node in nodes)
    gs_pop_remove(node, gs = gs)
  #refer to cd4-
  gs_add_gating_method(gs, pop = "+/-+/-", parent = "root", dims = "CD4,CD8", gating_method = "refGate", gating_args = "CD4-:CD8")
  stats2 <- gh_pop_compare_stats(gs[[1]], subpopulations = nodes)
  
  expect_equal(stats1, stats2)
  
})
