context("negated gate")

test_that("negated ellipse gate", {
  dataDir <- system.file("extdata",package="flowWorkspaceData")
  gs <- load_gs(list.files(dataDir, pattern = "gs_manual",full = TRUE))
  #grab an existing ellipse gate for the demo purpose
  ellipse_gate <- gh_get_gate(gs[[1]], "CD4")
  #create a dummy gating function that simply return this static gate
  dummy_func <- function(fr, pp_res, channels, ...)
  {
    return (ellipse_gate)
  }
  #register it as the plugin gating method
  registerPlugins(dummy_func, "cd4gate")
  #apply it to the gs
  add_pop(gs, alias = "cd4_negated", parent = "CD3+", pop = "-", dims = "CD4,CD8", gating_method = "cd4gate")
  
  expect_true(flowWorkspace:::gh_is_negated(gs[[1]],"cd4_negated"))
  expect_equal(gh_get_proportion(gs[[1]], "cd4_negated"), 0.3753648)
  
})
