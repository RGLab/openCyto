library(data.table)
library(utils)
library(tools)
library(flowCore)
library(ncdfFlow)
gtFile <- system.file("extdata/gating_template/tcell.csv", package = "openCyto")
expectResults <- readRDS("expectResults.rds")

