# Regression tests for gl.assign.mahalanobis (review of PRs #80/#81).
#
# As merged: singular covariance matrices were inverted with tol = 1e-20,
# giving negative D2 scored p = 1 (6 populations on the documented example,
# best assignment EmmacBrisWive instead of EmmacMaclGeor); with one dimension
# the "unknown" label fell on the first reference individual; p-values used
# chi-square instead of the small-sample F test; unknown was the 6th
# argument; and the imputed genlight was returned.

test_that("unknown is the second argument", {
  r <- quiet(gl.assign.mahalanobis(testset.gl, "UC_00146", verbose = 0))
  expect_s4_class(r, "genlight")
})

test_that("singular populations are not given negative distances", {
  out <- utils::capture.output(
    r <- gl.assign.mahalanobis(testset.gl, unknown = "UC_00146", verbose = 3)
  )
  rows <- grep("^ *[0-9]+ +Emmac", out, value = TRUE)
  d2 <- suppressWarnings(as.numeric(vapply(strsplit(trimws(rows), " +"),
                                           `[`, "", 3)))
  expect_true(all(d2[!is.na(d2)] >= 0))
  expect_true(any(grepl("not tested", out)))
  expect_true("EmmacMaclGeor" %in% kept_pops(r))
})

test_that("with one dimension the unknown's own distance is tested", {
  gl <- sim_pops(n.pop = 2, n.ind = 20, n.loc = 200, fst = 0.2, seed = 6)
  r <- gl.assign.mahalanobis(gl, unknown = "A1", verbose = 0)
  # dim.limit = nPop - 1 = 1; A1 lies far from B
  expect_equal(kept_pops(r), "A")
})

test_that("members of a population are not rejected from it", {
  # Hotelling's F test holds its size for a new individual; the chi-square
  # approximation rejected true sources at n = 10
  gl <- sim_pops(n.pop = 3, n.ind = 12, n.loc = 300, fst = 0.1, seed = 7)
  hits <- vapply(c("A1", "B13", "C25"), function(u) {
    r <- gl.assign.mahalanobis(gl, unknown = u, verbose = 0)
    substr(u, 1, 1) %in% kept_pops(r)
  }, logical(1))
  expect_true(all(hits))
})

test_that("returns observed, not imputed, genotypes", {
  r <- quiet(gl.assign.mahalanobis(testset.gl, unknown = "UC_00146",
                                   verbose = 0))
  expect_gt(sum(is.na(as.matrix(r))), 0)
  expect_true("unknown" %in% popNames(r))
})
