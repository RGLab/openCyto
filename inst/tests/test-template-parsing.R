context("Gating Template")

expectResults <- readRDS(system.file("tests/expectResults.rds", package = "openCyto"))

one_pop_token <- "[\\+-]"
pop_name_pat <- "[^\\+-]*"
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
      preprocessed_dt[, isMultiPops := FALSE]
      
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

#used by multiple test cases
template_row <- data.table(alias = "TH", pop = "cd4+cd8-"
                          , parent = "cd3", dims = "cd4,cd8"
                          , gating_method = "mindensity", gating_args = ""
                          , collapseDataForGating = "TRUE", groupBy = "4"
                          , preprocessing_method = "", preprocessing_args = ""
                      )

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
      expect_error(.unique_check_alias(preprocessed_dt, alias = "cd4+", this_parent = parentPath),  "not unique")

      #multi-pops
      #with pop names specified
      this_row <- copy(template_row)
      this_row[, alias := "cd4,cd8"]
      this_row[, pop := "*"]
      
      expect_error(.unique_check_alias(this_row, alias = "cd4,cd8", this_parent = "cd3"),  "not unique")
      expect_null(.unique_check_alias(this_row, alias = "cd4", this_parent = "cd3"))
      
      #no pop names specified
      this_row[, alias := "*"]
      expect_null(.unique_check_alias(this_row, alias = "*", this_parent = "cd3"))
    })

test_that(".splitTerms", {
      
      dims <- "A,B"
      # A+B+
      expectRes <- c("A+", "B+")
      expect_equivalent(.splitTerms(one_pop_pat, two_pop_token, "A+B+", dims), expectRes)
      expect_equivalent(.splitTerms(one_pop_pat, two_pop_token, "++", dims), expectRes)
    
      #A+/-B+
      
      expect_equivalent(.splitTerms(two_or_one_pop_pat, two_pop_token, "A+B+", dims), expectRes)
    
      expectRes <- list(c("A+", "A-") , "B+") 
      expect_equivalent(.splitTerms(two_or_one_pop_pat, two_pop_token, "A+/-B+", dims), expectRes)                      
    
      #+/-+  
      expect_equivalent(.splitTerms(two_or_one_pop_pat, two_pop_token, "+/-+", dims), expectRes)                      
    })

test_that(".gen_dummy_ref_gate", {
      
      this_row <- copy(template_row)
      this_row[, alias := "cd4,cd8"]
      this_row[, pop := "*"]
      
      dummy_row1 <- copy(this_row)
      dummy_row1[, alias := "cd4"]
      dummy_row1[, pop := ""]
      dummy_row1[, gating_method := "dummy_gate"]
      dummy_row1[, gating_args := file.path("cd3/cd4,cd8")]
      dummy_row1[, collapseDataForGating := ""]
      dummy_row1[, groupBy := ""]
      dummy_row1[, preprocessing_method := ""]
      dummy_row1[, preprocessing_args := ""]
      
      dummy_row2 <- copy(dummy_row1)
      dummy_row2[, alias := "cd8"]
      
      expectRes <- rbindlist(list(this_row, dummy_row1, dummy_row2)) 
      expect_equal(.gen_dummy_ref_gate(this_row), expectRes)
    })

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
      
      thisRes <- .gen_1dgate(this_row)
      expect_equal(thisRes, expectRes)
      
      #A-/B+/-
      this_row <- copy(template_row)
      this_row[, pop := "cd4-cd8+/-"]
      thisRes <- .gen_1dgate(this_row)
      expect_equal(thisRes, expectRes)
      
      #inconsistency between dims and pops
      this_row <- copy(template_row)
      this_row[, pop := "A-cd8+"]
      expect_equal(.gen_1dgate(this_row), expectRes)
    })

test_that(".gen_refGate", {
      
      
      #A+B+                            
      this_row <- copy(template_row)
      expectRes <- copy(this_row)
      expectRes[, gating_method := "refGate"]
      expectRes[, gating_args := "cd3/cd4+:cd3/cd8+"]
      expectRes[, collapseDataForGating := ""]
      expectRes[, groupBy := ""]
      
      thisRes <- .gen_refGate("cd4+cd8-", this_row, c("cd4+", "cd8+"), this_row[, alias])
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
      thisRes <- .gen_refGate(c("cd4-cd8+" ,"cd4-cd8-"), this_row, c("cd4+", "cd8+"))
      expect_equivalent(thisRes, expectRes)
      
      
    })

test_that(".preprocess_row", {
      
      #return as it is
      this_row <- copy(template_row)
      this_row[, pop := "cd4+"]
      expectRes <- this_row
      expect_equal(.preprocess_row(this_row), expectRes)
      this_row[, pop := "+"]
      expect_equal(.preprocess_row(this_row), expectRes)
      
      #invalid pop name
      this_row <- copy(template_row)
      this_row[, pop := "cd3++-"]
      expect_error(.preprocess_row(this_row), "invalid population pattern ")
      
      
      #A+/- 1d
      ########################
      this_row <- copy(template_row)
      this_row[, pop := "cd4+/-"]
      this_row[, dims := "cd4"]
      
      expectRes <- rbindlist(list(this_row, this_row))
      expectRes[1, alias := "cd4+"]
      expectRes[1, pop := "cd4+"]
      expectRes[2, alias := "cd4-"]
      expectRes[2, pop := "cd4-"] 
      expectRes[2, gating_method := "refGate"]
      expectRes[2, gating_args := "cd3/cd4+"]
      
      expect_equal(suppressMessages(.preprocess_row(this_row)), expectRes)
      #pure +/-
      this_row[, pop := "+/-"]
      expect_equal(suppressMessages(.preprocess_row(this_row)), expectRes)
      #generic name
      this_row[, pop := "A+/-"]
      expect_equal(suppressMessages(.preprocess_row(this_row)), expectRes)
      
      #A+/- 2d
      ########################
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
      
      #+/- 2d
      this_row[, pop := "+/-"]
      expect_error(suppressMessages(.preprocess_row(this_row)), regexp = "Please provide a proper pop name")
      
      #3d
      this_row[, dims := "cd3,cd4,cd8"]
      expect_error(suppressMessages(.preprocess_row(this_row)), regexp = "invalid number of dim")
      
      #A+B+
      ########################
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
      
      #generic
      this_row[, pop := "A+B-"]
      expect_equal(suppressMessages(.preprocess_row(this_row)), expectRes)
      #pure +-
      this_row[, pop := "+-"]
      expect_equal(suppressMessages(.preprocess_row(this_row)), expectRes)
      
      #A+B+ (refGate)
      ########################
      this_row <- copy(template_row)
      this_row[, pop := "cd4+cd8-"]
      this_row[, dims := "cd4,cd8"]
      this_row[, gating_method := "refGate"]
      expect_equal(.preprocess_row(this_row), this_row)
      #generic
      expectRes <- copy(this_row)
      this_row[, pop := "A+B-"]
      expect_equal(suppressMessages(.preprocess_row(this_row)), expectRes)
      #pure +-
      this_row[, pop := "+-"]
      expect_equal(suppressMessages(.preprocess_row(this_row)), expectRes)
      
      #A+/-B+/-
      ########################
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
      
      #generic
      this_row[, pop := "A+/-B-"]
      expect_equal(suppressMessages(.preprocess_row(this_row)), expectRes)
      #pure +-
      this_row[, pop := "+/--"]
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
      
      #multi-pops
      #no expansion    
      this_row <- copy(template_row)
      this_row[, alias := "*"]
      this_row[, pop := "*"]
      expect_equal(.preprocess_row(this_row), this_row)
      
      #with expansion
      this_row[, alias := "cd4,cd8"]
      expect_equal(.preprocess_row(this_row), .gen_dummy_ref_gate(this_row))
    })


test_that("gatingTemplate constructor", {
  thisPath <- system.file(package = "openCyto")
#  thisPath <- file.path(thisPath, "inst")
  gtfiles <- list.files(file.path(thisPath, "extdata/gating_template"), full = TRUE)
  options(scipen=999)
  for(thisFile in gtfiles){
    
    templateName <- file_path_sans_ext(basename(thisFile))
    message(templateName)
    suppressWarnings(suppressMessages(
        thisRes <- gatingTemplate(thisFile, autostart = 1L)
    ))

    expect_equal(thisRes, expectResults[[templateName]])
    
    
    #test as.data.table method
    dt.orig <- fread(thisFile, autostart = 1L)
    dt.orig <- openCyto:::.preprocess_csv(dt.orig)
    
    dt.new <- as.data.table(thisRes)
    
    #sort by alias
    keys <- c("alias", "parent")
    setkeyv(dt.orig, keys)
    setkeyv(dt.new, keys)
    #standardize collapseDataForGating col
    dt.orig[, collapseDataForGating:=as.logical(ifelse(collapseDataForGating == "", FALSE, collapseDataForGating))]
    
    
    #skip gating_args (for tcell panel) since there is format changes due to the scientifc notation (1e3 -->1000)
    #and we don't want to change this csv since it is used elsewehere
    if(templateName == "tcell")        
      expect_equal(dt.orig[, -6, with = FALSE], dt.new[, -6, with = FALSE])
    else
      expect_equal(dt.orig, dt.new)

  }
  
}) 


        
