context("gating functions")
data(GvHD)
fr <- GvHD[[1]]
fr <- transform(fr, transformList(c("FL1-H"),  log))
test_that("mindensity2", {
      
      expect_equal(mindensity2(fr, "SSC-H"), rectangleGate(`SSC-H` = c(901.120, Inf), filterId = ""), tol = 5e-07)
      expect_equal(mindensity2(fr, "FL1-H"), rectangleGate(`FL1-H` = c(5.216064, Inf), filterId = ""), tol = 5e-07)
      #the latest pull from phu seems to change the cutpoint for this call
      #      expect_equal(mindensity2(fr, "FL2-H"), rectangleGate(`FL2-H` = c(511.0989417, Inf), filterId = ""), tol = 5e-07)
      expect_equal(mindensity2(fr, "FL2-H"), rectangleGate(`FL2-H` = c(4777.69089, Inf), filterId = ""), tol = 5e-07)
      
      
    })
