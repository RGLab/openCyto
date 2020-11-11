context("gating functions")
data(GvHD)
fr <- GvHD[[1]]
fr <- transform(fr, transformList(c("FL1-H"),  log))
test_that("gate_mindensity2", {
      
      expect_equal(gate_mindensity2(fr, "SSC-H"), rectangleGate(`SSC-H` = c(901.120, Inf), filterId = ""), tol = 5e-07)
      expect_equal(gate_mindensity2(fr, "FL1-H"), rectangleGate(`FL1-H` = c(5.216064, Inf), filterId = ""), tol = 5e-07)
      #the latest pull from phu seems to change the cutpoint for this call
      #      expect_equal(mindensity2(fr, "FL2-H"), rectangleGate(`FL2-H` = c(511.0989417, Inf), filterId = ""), tol = 5e-07)
      expect_equal(gate_mindensity2(fr, "FL2-H"), rectangleGate(`FL2-H` = c(4777.69089, Inf), filterId = ""), tol = 5e-07)
      
      
    })

test_that("gate_tail", {
  
  # Right tail
  gs1 <- load_gs(system.file("extdata", "gs_manual", package = "flowWorkspaceData"))
  gs2 <- gs_clone(gs1)
  
  # flowWorkspace::gs_pop_add approach
  right_tail_gate <- gate_tail(gh_pop_get_data(gs1[[1]], "CD4"), channel = "<G560-A>", side = "right")
  right_tail_gate@filterId <- "CCR7_right_tail"
  right_tail_gate
  gs_pop_add(gs1, right_tail_gate, parent = "CD4", name = "CCR7_right_tail_pos")
  gs_pop_add(gs1, right_tail_gate, parent = "CD4", name = "CCR7_right_tail_neg", negated = TRUE)
  recompute(gs1)
  
  # openCyto::gs_add_gating_method approach
  openCyto::gs_add_gating_method(gs2, alias = "CCR7_right_tail_pos", pop = "+", parent = "CD4",
                                 dims = "<G560-A>", gating_method = "gate_tail", 
                                 gating_args = "side='right'")
  openCyto::gs_add_gating_method(gs2, alias = "CCR7_right_tail_neg", pop = "-", parent = "CD4",
                                 dims = "<G560-A>", gating_method = "gate_tail", 
                                 gating_args = "side='right'")
  
  g1_pos <- gh_pop_get_gate(gs1[[1]], "CCR7_right_tail_pos")
  g2_pos <- gh_pop_get_gate(gs2[[1]], "CCR7_right_tail_pos")
  g1_neg <- gh_pop_get_gate(gs1[[1]], "CCR7_right_tail_pos")
  g2_neg <- gh_pop_get_gate(gs2[[1]], "CCR7_right_tail_pos")
  
  # For each, check that the cut point is on the right side of the peak
  test_bound <- c("<G560-A>" = 2661)
  expect_equal(g1_pos@min, test_bound, tolerance = 1)
  expect_equal(g2_pos@min, test_bound, tolerance = 1)
  expect_equal(g1_neg@min, test_bound, tolerance = 1)
  expect_equal(g2_neg@min, test_bound, tolerance = 1)
  
  # And check that their infinite bounds are in the appropriate direction
  # Negation is done at the logical level, so the bounds should be untouched
  test_bound <- c("<G560-A>" = Inf)
  expect_equal(g1_pos@max, test_bound)
  expect_equal(g2_pos@max, test_bound)
  expect_equal(g1_pos@max, test_bound)
  expect_equal(g2_pos@max, test_bound)
  
  # Check that the resulting stats of both approaches agree
  expect_equal(gh_pop_get_stats(gs1[[1]], "CCR7_right_tail_pos"), gh_pop_get_stats(gs2[[1]], "CCR7_right_tail_pos"))
  expect_equal(gh_pop_get_stats(gs1[[1]], "CCR7_right_tail_neg"), gh_pop_get_stats(gs2[[1]], "CCR7_right_tail_neg"))
  
  
  # Left tail
  
  gs1 <- load_gs(system.file("extdata", "gs_manual", package = "flowWorkspaceData"))
  gs2 <- gs_clone(gs1)
  
  # flowWorkspace::gs_pop_add approach
  left_tail_gate <- gate_tail(gh_pop_get_data(gs1[[1]], "CD4"), channel = "<G560-A>", side = "left")
  left_tail_gate@filterId <- "CCR7_left_tail"
  left_tail_gate
  gs_pop_add(gs1, left_tail_gate, parent = "CD4", name = "CCR7_left_tail_pos")
  gs_pop_add(gs1, left_tail_gate, parent = "CD4", name = "CCR7_left_tail_neg", negated = TRUE)
  recompute(gs1)
  
  # openCyto::gs_add_gating_method approach
  openCyto::gs_add_gating_method(gs2, alias = "CCR7_left_tail_pos", pop = "+", parent = "CD4",
                                 dims = "<G560-A>", gating_method = "gate_tail", 
                                 gating_args = "side='left'")
  openCyto::gs_add_gating_method(gs2, alias = "CCR7_left_tail_neg", pop = "-", parent = "CD4",
                                 dims = "<G560-A>", gating_method = "gate_tail", 
                                 gating_args = "side='left'")
  
  g1_pos <- gh_pop_get_gate(gs1[[1]], "CCR7_left_tail_pos")
  g2_pos <- gh_pop_get_gate(gs2[[1]], "CCR7_left_tail_pos")
  g1_neg <- gh_pop_get_gate(gs1[[1]], "CCR7_left_tail_pos")
  g2_neg <- gh_pop_get_gate(gs2[[1]], "CCR7_left_tail_pos")
  
  # For each, check that the cut point is on the left side of the peak
  test_bound <- c("<G560-A>" = 1991)
  expect_equal(g1_pos@min, test_bound, tolerance = 1)
  expect_equal(g2_pos@min, test_bound, tolerance = 1)
  expect_equal(g1_neg@min, test_bound, tolerance = 1)
  expect_equal(g2_neg@min, test_bound, tolerance = 1)
  
  # And check that their infinite bounds are in the appropriate direction
  # Negation is done at the logical level, so the bounds should be untouched
  test_bound <- c("<G560-A>" = Inf)
  expect_equal(g1_pos@max, test_bound)
  expect_equal(g2_pos@max, test_bound)
  expect_equal(g1_pos@max, test_bound)
  expect_equal(g2_pos@max, test_bound)
  
  # Check that the resulting stats of both approaches agree
  expect_equal(gh_pop_get_stats(gs1[[1]], "CCR7_left_tail_pos"), gh_pop_get_stats(gs2[[1]], "CCR7_left_tail_pos"))
  expect_equal(gh_pop_get_stats(gs1[[1]], "CCR7_left_tail_neg"), gh_pop_get_stats(gs2[[1]], "CCR7_left_tail_neg"))
})
