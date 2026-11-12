# Regression tests for gl.run.assignpop (review of PRs #80/#81).
#
# As merged, split and two-object modes wrote AnalysisInfo.txt,
# AssignmentResult.txt and UsedLoci.txt to getwd() by default.

test_that("single-object mode returns the assignPOP data structure", {
  dat <- gl.run.assignpop(possums.gl, verbose = 0)
  expect_named(dat, c("DataMatrix", "SampleID", "LocusName"))
  expect_equal(ncol(dat$DataMatrix), 2 * nLoc(possums.gl) + 1)
})

test_that("split mode writes to tempdir(), not the working directory", {
  skip_if_not_installed("assignPOP")
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  out <- gl.run.assignpop(possums.gl, unknown.id = "1", verbose = 0)
  expect_length(list.files(wd), 0)
  expect_equal(out$results$Ind.ID, "1")
})
