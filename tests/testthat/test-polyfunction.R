context("polyfunction")

test_that("polyfunction testing", {
  skip("no longer supported")
  gs <- GatingSet(GvHD[1])
  gs <- transform(gs, estimateLogicle(gs[[1]], colnames(gs)[3:7]))
  gs_add_gating_method(gs, alias = "FL1+", parent = "root", dims = "FL1", gating_method = "tailgate")
  gs_add_gating_method(gs, alias = "FL3+", parent = "root", dims = "FL3", gating_method = "tailgate")
  gs_add_gating_method(gs, alias = "FL4+", parent = "root", dims = "FL4", gating_method = "tailgate")
  gs_add_gating_method(gs, alias = "*", parent = "root", dims = "*", gating_method = "polyfunctions", gating_args = "FL1+:FL3+:FL4+")
  expect_equal(gs_get_pop_paths(gs)[-(1:4)], c("/FL1+&FL3+&FL4+", "/FL1+&FL3+&!FL4+", "/FL1+&!FL3+&FL4+", "/FL1+&!FL3+&!FL4+", "/!FL1+&FL3+&FL4+", "/!FL1+&FL3+&!FL4+", "/!FL1+&!FL3+&FL4+", "/!FL1+&!FL3+&!FL4+"))
  node <- "/FL1+&!FL3+&!FL4+"
  expect_equal(gh_pop_get_count(gs[[1]], node), sum(gh_pop_get_indices(gs[[1]], "FL1+")&!gh_pop_get_indices(gs[[1]], "FL3+")&!gh_pop_get_indices(gs[[1]], "FL4+")))
  })