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

test_that("solve_LSAP", {
  
  x <- matrix(c(5, 1, 4, 3, 5, 2, 2, 4, 4), nrow = 3)
  x
  y <- as.vector(solve_LSAP(x))
  expect_equal(y, c(3,1,2))  
  
  x = structure(c(120.19673111448, 109.212815628447
                  , 693.431895464961, 611.636808127566)
                , .Dim = c(2L, 2L)
                , .Dimnames = list(c("1", "2"), c("1", "2")))
  y <- as.vector(solve_LSAP(x))
  expect_equal(y, c(1,2))  
  
  x = structure(c(49.9764377352581, 150.456415722488, 998.26447559313, 882.252993034892)
                , .Dim = c(2L, 2L)
                , .Dimnames = list(c("1", "2"), c("1", "2")))
  y <- as.vector(solve_LSAP(x))
  expect_equal(y, c(1,2))  
  
})
