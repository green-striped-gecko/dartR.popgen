# Regression tests for gl.assign.pca (review of PRs #80/#81).
#
# As merged, the function returned the imputed genlight, labelled the unknown
# "Unknown" (its siblings use "unknown"), and stopped at verbose >= 3 when a
# population could not be tested ("missing value where TRUE/FALSE needed").

skip_if_not_installed("SIBER")

test_that("returns observed genotypes and labels the unknown 'unknown'", {
  r <- gl.assign.pca(testset.gl, unknown = "UC_00146", plot.out = FALSE,
                     verbose = 0)
  expect_true("unknown" %in% popNames(r))
  expect_false("Unknown" %in% popNames(r))
  expect_gt(sum(is.na(as.matrix(r))), 0)
  expect_true("EmmacMaclGeor" %in% kept_pops(r))
})

test_that("an untestable population is retained without stopping", {
  gl <- sim_pops(n.pop = 2, n.ind = 15, n.loc = 100)
  gl <- gl[c(1:15, 16:17), ]      # B keeps 2 individuals: no ellipse
  out <- utils::capture.output(
    r <- gl.assign.pca(gl, unknown = "A1", nmin = 1, plot.out = FALSE,
                       verbose = 3)
  )
  expect_true("B" %in% kept_pops(r))
  expect_true(any(grepl("not tested", out)))
})
