context("gating functions")
data(GvHD)
fr <- GvHD[[1]]
fr <- transform(fr, transformList(c("FL1-H"),  log))

test_that("gate_singlet", {
  fr <- read.FCS(system.file("extdata/CytoTrol_CytoTrol_1.fcs", package = "flowWorkspaceData"))
  g <- gate_singlet(fr, "FSC-A", "FSC-H"
                    # , robust = F
                    # , prediction_level =0.9
                    )
  expect_equal(g, polygonGate(structure(c(24601, 24601
                                          , 262143, 262143
                                          , 11476.9417996388
                                          , 48163.8684594367
                                          , 236614.412689126, 199926.505500422)
                                          , .Dim = c(4L, 2L)
                                          , .Dimnames = list(NULL, c("FSC-A", "FSC-H")))
                                , filterId = "singlet"), tol = 5e-07)
  g <- gate_singlet(fr, "FSC-A", "FSC-H"
                    , sidescatter = "SSC-A"
                    # , prediction_level =0.9
  )
  expect_equal(g, polygonGate(structure(c(24601, 24601
                                          , 262143, 262143
                                          , 16988.0208566085, 41393.3256789415
                                          , 227799.565338899, 203392.940551338)
                                        , .Dim = c(4L, 2L)
                                        , .Dimnames = list(NULL, c("FSC-A", "FSC-H")))
                              , filterId = "singlet"), tol = 5e-06)
  # library(ggcyto)
  # autoplot(fr, "FSC-A", "FSC-H") + geom_gate(g) + geom_stats()
})

test_that("gate_mindensity2", {
      
      expect_equal(gate_mindensity2(fr, "SSC-H"), rectangleGate(`SSC-H` = c(901.120, Inf), filterId = ""), tol = 5e-07)
      expect_equal(gate_mindensity2(fr, "FL1-H"), rectangleGate(`FL1-H` = c(5.216064, Inf), filterId = ""), tol = 5e-07)
      #the latest pull from phu seems to change the cutpoint for this call
      #      expect_equal(mindensity2(fr, "FL2-H"), rectangleGate(`FL2-H` = c(511.0989417, Inf), filterId = ""), tol = 5e-07)
      expect_equal(gate_mindensity2(fr, "FL2-H"), rectangleGate(`FL2-H` = c(4777.69089, Inf), filterId = ""), tol = 5e-07)
      
      
    })
