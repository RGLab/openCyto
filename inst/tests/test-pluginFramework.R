context("plugin Framework")

test_that("listgtMethods", {
      thisRes <- capture.output(listgtMethods())
      expect_equal(thisRes, .DEFAULTS)
      
    })

test_that(".isRegistered", {
      
      expect_true(.isRegistered("mindensity"))
      expect_false(.isRegistered("myFunc"))
      expect_true(.isRegistered(".mindensity"))
      
    })
test_that("registerGatingFunction", {
      
      
      #depdency not installed
      expect_false(.isRegistered("flowDensity"))
      expect_message(thisRes <- registerGatingFunction(fun=.flowDensity,methodName="flowDensity",dep="flowDensity"), "Can't register flowDensity with dependency on flowDensity, because dependency is not installed.")
      expect_false(thisRes)
      expect_false(.isRegistered("flowDensity"))
      
      #invalid formal
      myFunc <- function(fr1, pp_res, yChannel, filterId, ...){
        cat("dummy gating function")
      }
      expect_message(thisRes <- registerGatingFunction(fun=myFunc,methodName="myFunc"), "Formals of function don't match expected template")
      expect_false(thisRes)
      expect_false(.isRegistered("myFunc"))

      #successfully registered without rep
      myFunc <- function(fr, pp_res, yChannel, filterId, ...){
        cat("dummy gating function")
      }
      
      
      expect_message(thisRes <- registerGatingFunction(fun = myFunc, methodName = "myFunc"), "Registered myFunc")
      expect_true(thisRes)
      expect_true(.isRegistered("myFunc"))
      expect_true(.unregister("myFunc"))
      expect_false(.isRegistered("myFunc"))
      
    })


