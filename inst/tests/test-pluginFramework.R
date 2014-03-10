context("plugin Framework")

test_that(".isRegistered", {
      
      expect_true(.isRegistered("mindensity"))
      expect_true(.isRegistered("warpSet"))
      
      expect_true(.isRegistered(".mindensity"))
      expect_true(.isRegistered(".warpSet"))
      
      expect_false(.isRegistered("myFunc"))
    })
test_that("registerPlugins", {
      
      
      #depdency not installed
      expect_false(.isRegistered("flowDensity"))
      expect_message(thisRes <- registerPlugins(fun=.flowDensity,methodName="flowDensity",dep="flowDensity"), "Can't register flowDensity with dependency on flowDensity, because dependency is not installed.")
      expect_false(thisRes)
      expect_false(.isRegistered("flowDensity"))
      
      #invalid formal
      myGatingFunc <- function(fr1, pp_res, yChannel, ...){
        cat("dummy gating function")
      }
      expect_message(thisRes <- registerPlugins(fun = myGatingFunc,methodName="myGatingFunc"), "Formals of function don't match expected template")
      expect_false(thisRes)
      expect_false(.isRegistered("myGatingFunc"))

      #successfully registered without rep
      myGatingFunc <- function(fr, pp_res, yChannel, ...){
        cat("dummy gating function")
      }
      
      expect_message(thisRes <- registerPlugins(fun = myGatingFunc, methodName = "myGatingFunc"), "Registered myGatingFunc")
      expect_true(thisRes)
      expect_true(.isRegistered("myGatingFunc"))
      expect_true(.unregister("myGatingFunc"))
      expect_false(.isRegistered("myGatingFunc"))

      
      #pp method
      myPreprocessingFunc <- function(fr, pp_res, yChannel, ...){
        cat("dummy preprocessing function")
      }
      
      expect_error(thisRes <- registerPlugins(type = "pp", fun = myPreprocessingFunc, methodName = "myPreprocessingFunc"), "'arg' should be one of ")
      expect_message(thisRes <- registerPlugins(type = "preprocess", fun = myPreprocessingFunc, methodName = "myPreprocessingFunc"), "Formals of function don't match expected template")
      expect_false(.isRegistered("myPreprocessingFunc"))
      
      myPreprocessingFunc <- function(fs, gs, gm, xChannel, yChannel, ...){
        cat("dummy preprocessing function")
      }
      expect_message(thisRes <- registerPlugins(type = "preprocess", fun = myPreprocessingFunc, methodName = "myPreprocessingFunc"), "Registered myPreprocessingFunc")
      expect_true(.isRegistered("myPreprocessingFunc"))
      expect_true(.unregister("myPreprocessingFunc", type = "pre"))
      expect_false(.isRegistered("myPreprocessingFunc"))
      
    })


