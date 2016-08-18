library(data.table)
library(utils)
library(tools)

gtFile <- system.file("extdata/gating_template/tcell.csv", package = "openCyto")
expectResults <- readRDS("expectResults.rds")