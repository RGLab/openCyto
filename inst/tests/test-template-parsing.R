context("Gating Template")

expectResults <- readRDS(system.file("tests/expectResults.rds", package = "openCyto"))
#test_that("templateGen", {
#      
#      templateGen(gh)
#    })

df <- as.data.frame(fread(gtFile))

test_that(".preprocess_csv", {
      
      suppressMessages(
          pre_csv <<- .preprocess_csv(df)
      )
      expect_equivalent(pre_csv, expectResults[["preprocess_csv"]])
      
    })

test_that(".gatingTemplate", {
      suppressMessages(thisRes <- .gatingTemplate(pre_csv, name="default"))
      expect_equal(thisRes, expectResults[["tcell"]])
    })

test_that(".getFullPath", {
      
      .getFullPath(ref_node, df)
    })

test_that(".check_alias", {
      
      this_df <- pre_csv
      parentPath <- "/nonDebris/singlets/lymph/cd3"
      #unique
      expect_null(.check_alias(this_df, alias = "cd4+", this_parent = parentPath))
      #not found
      expect_null(.check_alias(this_df, alias = "cd4", this_parent = parentPath))
      #illegal character |,&,:,/
      errMsg <- "contains illegal character:"
      expect_error(.check_alias(this_df, alias = "cd4+:", this_parent = parentPath), errMsg)
      expect_error(.check_alias(this_df, alias = "cd3/cd4+", this_parent = parentPath), errMsg)
      expect_error(.check_alias(this_df, alias = "cd8&cd4+", this_parent = parentPath), errMsg)
      expect_error(.check_alias(this_df, alias = "cd4+|cd8", this_parent = parentPath), errMsg)
      
      #duplicated entry
      tmp <- subset(this_df, alias == "cd8+") 
      tmp[["alias"]] <- "cd4+"
      new_df <- rbind(this_df, tmp)
      
      expect_error(.check_alias(new_df, alias = "cd4+", this_parent = parentPath), "not unique") 
      
    })

test_that(".splitTerms", {
      one_pop_token <- "[\\+-]"
      pop_name_pat <- "[^\\+-]+"
      one_pop_pat <- paste(pop_name_pat, one_pop_token, sep = "")
      
      two_pop_token <- "(\\+/-)|(-/\\+)"
      two_pop_pat <- paste(pop_name_pat, "(", two_pop_token, ")", sep = "")
      two_or_one_pop_pat <- paste0("(", two_pop_pat, ")|(", one_pop_pat, ")")
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
      
      #A+/-
      .splitTerms(two_or_one_pop_pat, two_pop_token, popName = "A+/-B-")
    })

test_that(".gen_1dgate", {
      
      .preprocess_csv(df)
    })

test_that(".gen_refGate", {
      
      .preprocess_csv(df)
    })

test_that(".gen_refGate", {
      
      .preprocess_csv(df)
    })

test_that(".addToDf", {
      
      .addToDf(res,this_row, new_df)
    })

expectResults[["gt1"]] <- thisRes

test_that("gatingTemplate constructor", {
    
  for(tn in c("bcell", "treg", "ICS", "tcell")){
    thisPath <- system.file(package = "openCyto")
#    thisPath <- file.path(thisPath, "inst")
    thisFile <- paste0("extdata/", tn, ".csv")
    suppressMessages(thisRes <- gatingTemplate(file.path(thisPath, thisFile)))
    expect_equal(thisRes, expectResults[[tn]])
  }

}) 

