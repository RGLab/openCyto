library(testthat)

require(flowWorkspace)
require(openCyto)
require(flowWorkspaceData)

file <- system.file("extdata/gs_auto", package="flowWorkspaceData")
gs <- load_gs(file)
ff <- getData(gs[[1]])
x <- exprs(ff)[, "FSC-A"]
test_that( "mindensity finds the valley it should", {
  gate <- mindensity(ff, "FSC-A")
  expect_that( unname(gate@min), equals(57500, tol=100) )
})
