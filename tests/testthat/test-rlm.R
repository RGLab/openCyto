test_that("fast_rlm", {

  
  set.seed(1)
  n <- 1e3
  x <- seq_len(n)
  y <- x * 2.5 - 1.3 + rnorm(n, sd = 30)
  df <- data.frame(x, y)
  #convert to rlm data format
  names(y) <- x
  x <- cbind(1, x)
  # r1 <- MASS::rlm(y~x, df)
  system.time(r1 <- MASS::rlm(x, y))
  system.time(r2 <- fast_rlm(x, y))
  # r2$coefficients is unnamed
  expect_equivalent(r1$coefficients, r2$coefficients)
  sel = c("coefficients"
          , "residuals"
          , "fitted.values"
          , "weights"
          , "rank"
          , "converged"
          , "wresid"
          , "w" ,"x", "s"
  )
  expect_equivalent(r1[sel], r2[sel])
  
  
  })