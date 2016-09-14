context("gating functions")
data(GvHD)
test_that("mindensity2", {
      fr <- GvHD[[1]]
      fr <- transform(fr, transformList(c("FL1-H"),  log))
  
      expect_equal(mindensity2(fr, "SSC-H"), rectangleGate(`SSC-H` = c(901.120, Inf), filterId = ""), tol = 5e-07)
      expect_equal(mindensity2(fr, "FL1-H"), rectangleGate(`FL1-H` = c(5.216064, Inf), filterId = ""), tol = 5e-07)
      #the latest pull from phu seems to change the cutpoint for this call
      #      expect_equal(mindensity2(fr, "FL2-H"), rectangleGate(`FL2-H` = c(511.0989417, Inf), filterId = ""), tol = 5e-07)
      expect_equal(mindensity2(fr, "FL2-H"), rectangleGate(`FL2-H` = c(4777.69089, Inf), filterId = ""), tol = 5e-07)
      
      
    })

fr <- transform(GvHD[[1]], transformList(c("FL1-H"),  log))
trans <- estimateLogicle(fr, c("FL1-H", "FL2-H", "FL3-H", "FL4-H", "FL2-A"))
fr <- transform(fr, trans)

test_that("mindensity2", {
    
  expect_equal(mindensity(fr, "SSC-H")@min, 901.120, tol = 5e-07, check.attributes = FALSE)
  expect_equal(mindensity2(fr, "FL1-H")@min, 0.856, tol = 4e-04, check.attributes = FALSE)
  expect_equal(mindensity(fr, "FL2-H")@min, 3.0569, tol = 5e-06, check.attributes = FALSE)
  expect_equal(mindensity(fr, "FL2-H", gate_range=c(1,2))@min, 1.2839, tol = 5e-06, check.attributes = FALSE)
  expect_equal(mindensity(fr, "FL3-H")@min, 0.538, tol = 5e-04, check.attributes = FALSE)
  expect_equal(mindensity(fr, "FL4-H")@min, 0.538, tol = 5e-04, check.attributes = FALSE)
  expect_equal(mindensity(fr, "FL2-A")@min, 1.047, tol = 2e-04, check.attributes = FALSE)
  expect_equal(mindensity(fr, "FSC-H")@min, 892.65, tol = 5e-06, check.attributes = FALSE)
  
})

