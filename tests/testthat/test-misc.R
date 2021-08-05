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

test_that("robust mean", {
  data(GvHD)
  fr <- GvHD[[1]]
  x <- exprs(fr)[,1]
  sd <- mad(x)
  res <- robust_m_estimator(x, sd)
  
  expect_equal(sd, 148.26, tol = 2e-7)
  expect_equal(res, 229.625, tol = 2e-7)
  
  set.seed(1)
  x = rnorm(100)
  sd <- mad(x)
  res <- robust_m_estimator(x, sd)
  expect_equal(sd, 0.8700003, tol = 2e-7)
  expect_equal(res, 0.1153234, tol = 2e-7)
  
})
