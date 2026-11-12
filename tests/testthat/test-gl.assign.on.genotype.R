# Regression tests for gl.assign.on.genotype (review of PRs #80/#81).
#
# As merged, a locus with no genotype in a population gave a NaN term that
# sum(na.rm = TRUE) dropped for that population only, so a population with
# more missing loci summed fewer negative terms and won the AIC weighting;
# absent alleles were scored log(1e-10).

test_that("missing data in one population does not make it more likely", {
  # A and B drawn from the same frequencies; B has no data at half the loci
  set.seed(4)
  q <- stats::runif(400, 0.1, 0.9)
  m <- t(replicate(41, stats::rbinom(400, 2, q)))
  m[21:40, 1:200] <- NA
  gl <- make_gl(m, c(rep("A", 20), rep("B", 20), "A"),
                ind = c(paste0("A", 1:20), paste0("B", 1:20), "U"))
  r <- gl.assign.on.genotype(gl, unknown = "U", verbose = 0)
  # As merged: B AIC.wt 1.0, A 6e-89 and excluded
  expect_true("A" %in% kept_pops(r))
})

test_that("an allele absent from a population sample does not exclude it", {
  # U carries the alternate allele at one locus where A's sample has none
  set.seed(5)
  q <- stats::runif(100, 0.2, 0.8)
  m <- t(replicate(41, stats::rbinom(100, 2, q)))
  m[1:20, 1] <- 0L
  m[41, 1] <- 1L
  gl <- make_gl(m, c(rep("A", 20), rep("B", 20), "A"),
                ind = c(paste0("A", 1:20), paste0("B", 1:20), "U"))
  r <- gl.assign.on.genotype(gl, unknown = "U", verbose = 0)
  expect_true("A" %in% kept_pops(r))
})

test_that("SilicoDArT data are refused", {
  expect_error(
    gl.assign.on.genotype(testset.gs, unknown = indNames(testset.gs)[1],
                          verbose = 0)
  )
})

test_that("n.best larger than the number of populations keeps them all", {
  r <- gl.assign.on.genotype(testset.gl, unknown = "UC_00146", n.best = 50,
                             verbose = 0)
  expect_gt(nPop(r), 10)
})

test_that("verbose = 0 prints nothing and the true source is kept", {
  out <- utils::capture.output(
    r <- gl.assign.on.genotype(testset.gl, unknown = "UC_00146", verbose = 0)
  )
  expect_length(out, 0)
  expect_true("EmmacMaclGeor" %in% kept_pops(r))
})
