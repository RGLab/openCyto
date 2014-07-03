context("misc functions")

test_that("fast collapsing for ncdfFlowSet", {
      library(openCyto)
      data(GvHD)
      fs <- GvHD[pData(GvHD)$Patient %in% 6:7][1:4]
      suppressMessages(ncfs <- ncdfFlowSet(fs))
      
      fr <- as(fs, "flowFrame")
      fr1 <- as(ncfs, "flowFrame")
      expect_is(fr1, "flowFrame")
      expect_equal(fr@exprs[, -9], fr1@exprs, tol = 2e-7)
      
    })
