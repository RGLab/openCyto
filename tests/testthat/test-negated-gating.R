context("negated gate")

test_that("negated ellipse gate", {
  gs <- gs_copy_tree_only(gs)
  #grab an existing ellipse gate for the demo purpose
  ellipse_gate <- gh_pop_get_gate(gs[[1]], "CD4")
  #create a dummy gating function that simply return this static gate
  dummy_func <- function(fr, pp_res, channels, ...)
  {
    return (ellipse_gate)
  }
  #register it as the plugin gating method
  register_plugins(dummy_func, "cd4gate")
  #apply it to the gs
  gs_add_gating_method(gs, alias = "cd4_negated", parent = "CD3+", pop = "-", dims = "CD4,CD8", gating_method = "cd4gate")
  
  expect_true(flowWorkspace::gh_pop_is_negated(gs[[1]],"cd4_negated"))
  expect_equal(gh_pop_get_proportion(gs[[1]], "cd4_negated"), 0.3753648)
  
  #add 2d ref gate  
  gs_add_gating_method(gs, alias = "cd4_ref", parent = "CD3+", pop = "+", dims = "CD4,CD8", gating_method = "refGate",gating_args = "CD4")
  expect_false(flowWorkspace::gh_pop_is_negated(gs[[1]],"cd4_ref"))
  expect_equal(gh_pop_get_proportion(gs[[1]], "CD4"), gh_pop_get_proportion(gs[[1]], "cd4_ref"))
    
  gs_add_gating_method(gs, alias = "cd4_ref_negated", parent = "CD3+", pop = "-", dims = "CD4,CD8", gating_method = "refGate",gating_args = "CD4")
  expect_true(flowWorkspace::gh_pop_is_negated(gs[[1]],"cd4_ref_negated"))
  expect_equal(gh_pop_get_proportion(gs[[1]], "cd4_negated"), gh_pop_get_proportion(gs[[1]], "cd4_ref_negated"))
})
