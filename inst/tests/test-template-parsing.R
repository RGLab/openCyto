context("Gating Template")

expectResults <- readRDS(system.file("tests/expectResults.rds", package = "openCyto"))

one_pop_token <- "[\\+-]"
pop_name_pat <- "[^\\+-]+"
one_pop_pat <- paste(pop_name_pat, one_pop_token, sep = "")

two_pop_token <- "(\\+/-)|(-/\\+)"
two_pop_pat <- paste(pop_name_pat, "(", two_pop_token, ")", sep = "")
two_or_one_pop_pat <- paste0("(", two_pop_pat, ")|(", one_pop_pat, ")")
#test_that("templateGen", {
#      
#      templateGen(gh)
#    })

dt <- fread(gtFile, autostart = 1L)


test_that(".preprocess_csv", {
      
      suppressMessages(
          preprocessed_dt <<- .preprocess_csv(dt)
      )
      expect_equivalent(preprocessed_dt, expectResults[["preprocess_csv"]])
      
    })


test_that(".gatingTemplate", {
      suppressMessages(thisRes <- .gatingTemplate(preprocessed_dt, name="default"))
      expect_equal(thisRes, expectResults[["tcell"]])
    })

test_that(".getFullPath", {
      ref_node <- "root"
      expect_equal(.getFullPath(ref_node, preprocessed_dt), "root")
      
      ref_node <- "cd3"
      expect_equal(.getFullPath(ref_node, preprocessed_dt), "/nonDebris/singlets/lymph/cd3")
      
      ref_node <- "lymph/cd3"
      expect_equal(.getFullPath(ref_node, preprocessed_dt), "/nonDebris/singlets/lymph/cd3")
      
      ref_node <- "/lymph/cd3"
      expect_error(.getFullPath(ref_node, preprocessed_dt), "Not able to to find reference")
      
      ref_node <- "cd31"
      expect_error(.getFullPath(ref_node, preprocessed_dt), "Not able to to find reference")
      
    })

test_that(".validity_check_alias", {
      
      #illegal character |,&,:,/
      errMsg <- "contains illegal character:"
      expect_error(.validity_check_alias("cd4+:"), errMsg)
      expect_error(.validity_check_alias("cd3/cd4+"), errMsg)
      expect_error(.validity_check_alias("cd8&cd4+"), errMsg)
      expect_error(.validity_check_alias("cd4+|cd8"), errMsg)
      
    })

test_that(".unique_check_alias", {
      
      
      parentPath <- "/nonDebris/singlets/lymph/cd3"
      
      expect_null(.unique_check_alias(preprocessed_dt, alias = "cd4", this_parent = parentPath))
      expect_error(.unique_check_alias(preprocessed_dt, alias = "cd4+", this_parent = parentPath), , "not unique")
    })

test_that(".splitTerms", {
      
      # A+B+
      expectRes <- list(terms = c("A+", "B+")
                        , splitted_terms = list("A+", "B+")
                        )
      expect_equal(.splitTerms(one_pop_pat, two_pop_token, "A+B+"), expectRes)
    
      #A+/-B+
      
      expect_equal(.splitTerms(two_or_one_pop_pat, two_pop_token, "A+B+"), expectRes)
    
      expectRes <- list(terms = c("A+/-", "B+")
                      , splitted_terms = list(c("A+", "A-")
                                              , "B+"
                                              )
                      ) 
      expect_equal(.splitTerms(two_or_one_pop_pat, two_pop_token, "A+/-B+"), expectRes)                      
      
    })

#used by multiple test cases
template_row <- data.table(alias = "TH", pop = "cd4+cd8-"
                          , parent = "cd3", dims = "cd4,cd8"
                          , gating_method = "mindensity", gating_args = ""
                          , collapseDataForGating = "TRUE", groupBy = "4"
                          , preprocessing_method = "", preprocessing_args = ""
                      )
                      
test_that(".gen_1dgate", {
      
      #A+B+                            
      this_row <- copy(template_row)
      expectRes <- rbindlist(list(this_row, this_row))
      expectRes[1, alias := "cd4+"]
      expectRes[1, pop := "cd4+"]
      expectRes[1, dims := "cd4"]
      expectRes[2, alias := "cd8+"]
      expectRes[2, pop := "cd8+"]
      expectRes[2, dims := "cd8"]
      
      thisRes <- .gen_1dgate(c("cd4+", "cd8+"), this_row, one_pop_token, two_pop_token)
      expect_equal(thisRes, expectRes)

      #A-/B+/-
      this_row <- copy(template_row)
      this_row[, pop := "cd4-cd8+/-"]
      thisRes <- .gen_1dgate(c("cd4-", "cd8+/-"), this_row, one_pop_token, two_pop_token)
      expect_equal(thisRes, expectRes)
      
    })

test_that(".gen_refGate", {
      
      
      #A+B+                            
      this_row <- copy(template_row)
      expectRes <- copy(this_row)
      expectRes[, gating_method := "refGate"]
      expectRes[, gating_args := "cd3/cd4+:cd3/cd8+"]
      expectRes[, collapseDataForGating := ""]
      expectRes[, groupBy := ""]
      
      thisRes <- .gen_refGate(c("cd4+", "cd8-"), this_row, c("cd4+", "cd8+"), this_row[, alias])
      expect_equal(thisRes, expectRes)
      
      #A-/B+/-
      this_row <- copy(template_row)
      this_row[, pop := "cd4-cd8+/-"]
      expectRes <- rbindlist(list(this_row, this_row))
      expectRes[1, alias := "cd4-cd8+"]
      expectRes[1, pop := "cd4-cd8+"]
      expectRes[1, gating_method := "refGate"]
      expectRes[1, gating_args := "cd3/cd4+:cd3/cd8+"]
      expectRes[1, collapseDataForGating := ""]
      expectRes[1, groupBy := ""]
      expectRes[2, alias := "cd4-cd8-"]
      expectRes[2, pop := "cd4-cd8-"]
      expectRes[2, gating_method := "refGate"]
      expectRes[2, gating_args := "cd3/cd4+:cd3/cd8+"]
      expectRes[2, collapseDataForGating := ""]
      expectRes[2, groupBy := ""]
      thisRes <- .gen_refGate(list("cd4-", c("cd8+" ,"cd8-")), this_row, c("cd4+", "cd8+"))
      expect_equivalent(thisRes, expectRes)
      
      
    })

test_that(".preprocess_row", {
      
      #return as it is
      this_row <- copy(template_row)
      this_row[, pop := "cd4+"]
      expectRes <- this_row
      expect_equal(.preprocess_row(this_row), expectRes)
      
      #invalid pop name
      this_row <- copy(template_row)
      this_row[, pop := "cd3++"]
      expect_error(.preprocess_row(this_row), "invalid population pattern ")
      
      #A+/-
      this_row <- copy(template_row)
      this_row[, pop := "cd4+/-"]
      
      expectRes <- rbindlist(list(this_row, this_row))
      expectRes[1, alias := "cd4+"]
      expectRes[1, pop := "cd4+"]
      expectRes[2, alias := "cd4-"]
      expectRes[2, pop := "cd4-"] 
      expectRes[2, gating_method := "refGate"]
      expectRes[2, gating_args := "cd3/cd4+"]
      
      expect_equal(suppressMessages(.preprocess_row(this_row)), expectRes)
      
      #A+B+ 
      this_row <- copy(template_row)
      this_row[, pop := "cd4+cd8-"]
      this_row[, dims := "cd4,cd8"]
      
      expectRes <- rbindlist(list(this_row, this_row, this_row))
      expectRes[1, alias := "cd4+"]
      expectRes[1, pop := "cd4+"]
      expectRes[1, dims := "cd4"]
      expectRes[2, alias := "cd8+"]
      expectRes[2, pop := "cd8+"]
      expectRes[2, dims := "cd8"]
      expectRes[3, gating_method := "refGate"]
      expectRes[3, gating_args := "cd3/cd4+:cd3/cd8+"]
      expectRes[3, collapseDataForGating := ""]
      expectRes[3, groupBy := ""]
      expect_equal(suppressMessages(.preprocess_row(this_row)), expectRes)
      
      #A+B+ (refGate)
      
      this_row <- copy(template_row)
      this_row[, pop := "cd4+cd8-"]
      this_row[, dims := "cd4,cd8"]
      this_row[, gating_method := "refGate"]
      expect_equal(.preprocess_row(this_row), this_row)
      
      #A+/-B+/-
      this_row <- copy(template_row)
      this_row[, pop := "cd4+/-cd8-"]
      this_row[, dims := "cd4,cd8"]
      
      expectRes <- rbindlist(list(this_row, this_row, this_row, this_row))
      expectRes[1, alias := "cd4+"]
      expectRes[1, pop := "cd4+"]
      expectRes[1, dims := "cd4"]
      expectRes[2, alias := "cd8+"]
      expectRes[2, pop := "cd8+"]
      expectRes[2, dims := "cd8"]
      expectRes[3, alias := "cd4+cd8-"]
      expectRes[3, pop := "cd4+cd8-"]
      expectRes[3, gating_method := "refGate"]
      expectRes[3, gating_args := "cd3/cd4+:cd3/cd8+"]
      expectRes[3, collapseDataForGating := ""]
      expectRes[3, groupBy := ""]
      expectRes[4, alias := "cd4-cd8-"]
      expectRes[4, pop := "cd4-cd8-"]
      expectRes[4, gating_method := "refGate"]
      expectRes[4, gating_args := "cd3/cd4+:cd3/cd8+"]
      expectRes[4, collapseDataForGating := ""]
      expectRes[4, groupBy := ""]
      expect_equal(suppressMessages(.preprocess_row(this_row)), expectRes)
      
      #TODO: invalid parent
#      this_row[["parent"]] <- "lymph1"
      
      
      #flowClust.1d
      this_row <- copy(template_row)
      this_row[, pop := "cd4"]
      this_row[, dims := "cd4"]
      this_row[, gating_method := "flowClust"]
      
      expectRes <- copy(this_row)
      expectRes[, gating_method := "flowClust.1d"]
      expectRes[, pop := "cd4+"]
      expect_equal(.preprocess_row(this_row), expectRes)
      
      #flowClust.2d
      this_row <- copy(template_row)
      this_row[, pop := "cd4+"]
      this_row[, gating_method := "flowClust"]
      expectRes <- copy(this_row)
      expectRes[, gating_method := "flowClust.2d"]
      expect_equal(.preprocess_row(this_row), expectRes)
      
    })


test_that("gatingTemplate constructor", {
  thisPath <- system.file(package = "openCyto")
#  thisPath <- file.path(thisPath, "inst")
  gtfiles <- list.files(file.path(thisPath, "extdata/gating_template"), full = TRUE)
  for(thisFile in gtfiles){
    
    templateName <- file_path_sans_ext(basename(thisFile))
    message(templateName)
    suppressWarnings(suppressMessages(
        thisRes <- gatingTemplate(thisFile, autostart = 1L)
    ))

    expect_equal(thisRes, expectResults[[templateName]])
  }
  
}) 


