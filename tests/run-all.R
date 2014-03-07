library(testthat)
library(openCyto)
library(data.table)
library(utils)
library(tools)

resultDir <- system.file("tests/expect_result",package="openCyto")
gtFile <- system.file("extdata/tcell.csv", package = "openCyto")
test_package("openCyto")


#test_file("/home/wjiang2/rglab/workspace/flowWorkspace/inst/tests/test-archive.R")
