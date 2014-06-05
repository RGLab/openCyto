context("plugin Framework")

test_that(".isRegistered", {
      
      expect_true(.isRegistered("mindensity"))
      expect_true(.isRegistered("warpSet"))
      
      expect_true(.isRegistered(".mindensity"))
      expect_true(.isRegistered(".warpSet"))
      
      expect_false(.isRegistered("myFunc"))
    })
test_that("registerPlugins", {
      
      dummyMethod <- function(){}
      #depdency not installed
      expect_false(openCyto:::.isRegistered("dummyMethod"))
      expect_message(thisRes <- registerPlugins(fun=dummyMethod,methodName="dummyMethod",dep="dummyPkg")
                    , "Can't register dummyMethod with dependency on dummyPkg, because dependency is not installed.")
      expect_false(thisRes)
      expect_false(openCyto:::.isRegistered("dummyMethod"))
      
      #invalid formal
      myGatingFunc <- function(fr1, pp_res, yChannel, ...){
        cat("dummy gating function")
      }
      expect_message(thisRes <- registerPlugins(fun = myGatingFunc,methodName="myGatingFunc"), "Formals of function don't match expected template")
      expect_false(thisRes)
      expect_false(openCyto:::.isRegistered("myGatingFunc"))

      #successfully registered without rep
      myGatingFunc <- function(fr, pp_res, yChannel, ...){
        cat("dummy gating function")
      }
      
      expect_message(thisRes <- registerPlugins(fun = myGatingFunc, methodName = "myGatingFunc"), "Registered myGatingFunc")
      expect_true(thisRes)
      expect_true(openCyto:::.isRegistered("myGatingFunc"))
      expect_true(openCyto:::.unregister("myGatingFunc"))
      expect_false(openCyto:::.isRegistered("myGatingFunc"))

      
      #pp method
      myPreprocessingFunc <- function(fr, pp_res, yChannel, ...){
        cat("dummy preprocessing function")
      }
      
      expect_error(thisRes <- registerPlugins(type = "pp", fun = myPreprocessingFunc, methodName = "myPreprocessingFunc"), "'arg' should be one of ")
      expect_message(thisRes <- registerPlugins(type = "preprocess", fun = myPreprocessingFunc, methodName = "myPreprocessingFunc"), "Formals of function don't match expected template")
      expect_false(openCyto:::.isRegistered("myPreprocessingFunc"))
      
      
      myPreprocessingFunc <- function(fs, gs, gm, xChannel, yChannel, ...){
        cat("dummy preprocessing function")
      }
      expect_message(thisRes <- registerPlugins(type = "preprocess", fun = myPreprocessingFunc, methodName = "myPreprocessingFunc")
                    , "Formals of function don't match")      
      expect_false(openCyto:::.isRegistered("myPreprocessingFunc"))
                
                
      myPreprocessingFunc <- function(fs, gs, gm, xChannel, yChannel, groupBy, isCollapse, ...){
        cat("dummy preprocessing function")
      }                
      expect_message(thisRes <- registerPlugins(type = "preprocess", fun = myPreprocessingFunc, methodName = "myPreprocessingFunc")
                    , "Registered myPreprocessingFunc")
      expect_true(openCyto:::.isRegistered("myPreprocessingFunc"))
      expect_true(openCyto:::.unregister("myPreprocessingFunc", type = "pre"))
      expect_false(openCyto:::.isRegistered("myPreprocessingFunc"))
      
    })


