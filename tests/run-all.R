library(testthat)
library(openCyto)
library(data.table)
library(utils)
library(tools)

gtFile <- system.file("extdata/gating_template/tcell.csv", package = "openCyto")
test_package("openCyto")

#taking quite some time , thus only for internal testing
#test_file(system.file("tests/gating-testSuite.R", package = "openCyto"))
