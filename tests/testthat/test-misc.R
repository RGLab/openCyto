context("misc functions")

test_that("fast collapsing for ncdfFlowSet", {
      library(openCyto)
      data(GvHD)
      fs <- GvHD[pData(GvHD)$Patient %in% 6:7][1:4]
      suppressMessages(ncfs <- ncdfFlowSet(fs))
      
      fr <- as(fs, "flowFrame")
      fr1 <- as(ncfs, "flowFrame")
      expect_is(fr1, "flowFrame")
      expect_equal(fr@exprs[, -9], fr1@exprs, tol = 2e-7)
      
    })

test_that("start/stop.at argument in gt_gating function", {
  library(openCyto)
  
  dataDir <- system.file("extdata",package="flowWorkspaceData")
  files <- list.files(dataDir, pattern = "a2004_O1T2pb05i_A1_A01.fcs$",full = TRUE)
  fs <- read.flowSet(files[1])
  
  compMat <- spillover(fs[[1]]) %>% compensation
  fs_comp <- compensate(fs, compMat)
  chnls <- parameters(compMat)
  transFuncts <- estimateLogicle(fs[[1]], channels = chnls)
  fs_trans <- transform(fs_comp, transFuncts)
  
  gs <- GatingSet(fs_trans)
  gt <- gatingTemplate(system.file("extdata/gating_template/testGatingTemplate.csv", package = "openCyto"))
  gt_gating(gt, gs, stop.at = "/boundary/singlet/Alexa700_gate")
  expect_setequal(gs_get_pop_paths(gs),c("root", "/boundary", "/boundary/singlet", "/boundary/singlet/PerCPCY5_5_gate", "/boundary/singlet/APC_gate", "/boundary/singlet/AmCyan_gate"))  
  gt_gating(gt, gs, start = "/boundary/singlet/Alexa700_gate")
  expect_setequal(gs_get_pop_paths(gs),c("root", "/boundary", "/boundary/singlet"
                                        , "/boundary/singlet/PerCPCY5_5_gate"
                                        , "/boundary/singlet/APC_gate"
                                        , "/boundary/singlet/AmCyan_gate"
                                        , "/boundary/singlet/Lymph"
                                        , "/boundary/singlet/Alexa700_gate"
                                        , "/boundary/singlet/Lymph/Alexa700pPerCPCY5_5p", "/boundary/singletRefGate", "/boundary/singletRefGate/singletRefGate_LymphRefGate"))  
})

test_that("robust mean", {
  data(GvHD)
  fr <- GvHD[[1]]
  x <- exprs(fr)[,1]
  sd <- mad(x)
  res <- robust_m_estimator(x, sd)
  
  expect_equal(sd, 148.26, tol = 2e-7)
  expect_equal(res, 229.625, tol = 2e-7)
  
  set.seed(1)
  x = rnorm(100)
  sd <- mad(x)
  res <- robust_m_estimator(x, sd)
  expect_equal(sd, 0.8700003, tol = 2e-7)
  expect_equal(res, 0.1153234, tol = 2e-7)
  
})

test_that("solve_LSAP", {
  
  x <- matrix(c(5, 1, 4, 3, 5, 2, 2, 4, 4), nrow = 3)
  x
  y <- as.vector(solve_LSAP(x))
  expect_equal(y, c(3,1,2))  
  
  x = structure(c(120.19673111448, 109.212815628447
                  , 693.431895464961, 611.636808127566)
                , .Dim = c(2L, 2L)
                , .Dimnames = list(c("1", "2"), c("1", "2")))
  y <- as.vector(solve_LSAP(x))
  expect_equal(y, c(1,2))  
  
  x = structure(c(49.9764377352581, 150.456415722488, 998.26447559313, 882.252993034892)
                , .Dim = c(2L, 2L)
                , .Dimnames = list(c("1", "2"), c("1", "2")))
  y <- as.vector(solve_LSAP(x))
  expect_equal(y, c(1,2))  
  
})
