# Regression tests for gl.assign.pa (review of PRs #80/#81).
#
# As merged, count.pa() counted only alternate-allele private alleles, so an
# unknown carrying a reference allele absent from a population scored zero
# against it; it counted loci with no data in the population as private; it
# stopped on an individual with no non-zero genotype ("invalid 'type' (list)")
# and on a single-member population; and n.best larger than the number of
# populations stopped in gl.keep.pop().

# Population A fixed for the alternate allele, B for the reference allele,
# at 20 loci; 30 further loci vary at random in both
toy_fixed <- function(unknown.geno, unknown.pop) {
  set.seed(1)
  m <- rbind(matrix(2L, 11, 20), matrix(0L, 11, 20), rep(unknown.geno, 20))
  m <- cbind(m, matrix(sample(0:2, 23 * 30, replace = TRUE), 23))
  make_gl(m, c(rep("A", 11), rep("B", 11), unknown.pop),
          ind = c(paste0("A", 1:11), paste0("B", 1:11), "U"))
}

test_that("private reference alleles exclude a population", {
  # U is 0/0 at the 20 loci where A is fixed 2/2
  r <- gl.assign.pa(toy_fixed(0L, "B"), unknown = "U", verbose = 0)
  expect_equal(kept_pops(r), "B")
})

test_that("private alternate alleles exclude a population", {
  r <- gl.assign.pa(toy_fixed(2L, "A"), unknown = "U", verbose = 0)
  expect_equal(kept_pops(r), "A")
})

test_that("loci with no data in a population are not private alleles", {
  set.seed(2)
  m <- matrix(sample(0:2, 23 * 40, replace = TRUE), 23)
  m[23, ] <- m[1, ]          # U is a copy of A1
  m[1:11, 1:10] <- NA        # A not genotyped at 10 loci
  m[23, 1:10] <- 2L
  gl <- make_gl(m, c(rep("A", 11), rep("B", 11), "A"),
                ind = c(paste0("A", 1:11), paste0("B", 1:11), "U"))
  r <- gl.assign.pa(gl, unknown = "U", verbose = 0)
  expect_true("A" %in% kept_pops(r))
})

test_that("an individual with no non-zero genotype does not stop the run", {
  set.seed(3)
  m <- rbind(matrix(sample(0:2, 22 * 10, replace = TRUE), 22), rep(0L, 10))
  gl <- make_gl(m, c(rep("A", 11), rep("B", 11), "A"),
                ind = c(paste0("A", 1:11), paste0("B", 1:11), "U"))
  expect_s4_class(gl.assign.pa(gl, unknown = "U", verbose = 0), "genlight")
})

test_that("a single-member population does not stop the run", {
  r <- quiet(gl.assign.pa(testset.gl, unknown = "UC_00146", nmin = 1,
                          verbose = 0))
  expect_s4_class(r, "genlight")
})

test_that("n.best larger than the number of populations keeps them all", {
  r <- gl.assign.pa(testset.gl, unknown = "UC_00146", n.best = 50,
                    verbose = 0)
  expect_gt(nPop(r), 10)
})

test_that("verbose = 0 prints nothing", {
  out <- utils::capture.output(
    r <- gl.assign.pa(testset.gl, unknown = "UC_00146", verbose = 0)
  )
  expect_length(out, 0)
})

test_that("the documented example keeps the true source", {
  r <- gl.assign.pa(testset.gl, unknown = "UC_00146", verbose = 0)
  expect_true("EmmacMaclGeor" %in% kept_pops(r))
  expect_true("unknown" %in% popNames(r))
})
