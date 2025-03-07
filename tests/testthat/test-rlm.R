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
  expect_equal(r1$coefficients, r2$coefficients, check.attributes = FALSE)
  
  #del unrelevant parameters before comparison
  r1 = r1[names(r2)]
  expect_equivalent(r1, r2)
  # library(ggplot2)
  # ggplot(df, aes(x,y)) +geom_point() + geom_abline(slope = r1$coefficients)
  # plot(y~x, data = data.frame(x, y)) + abline(r1)
  
  
  })