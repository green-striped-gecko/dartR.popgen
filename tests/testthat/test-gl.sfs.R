# Characterization tests for gl.sfs -- Phase A baseline of the dartR function
# review (dartR.popgen at 461af97), updated in Phase C: assertions tagged
# [approved n] were changed with the matching approved change; see
# function-review/reports/dartR.popgen/gl.sfs.md (dartR.base). Reference values are computed independently from the
# genotype matrix.

sfs_2pop <- function() {
  x <- gl.keep.pop(possums.gl, pop.list = popNames(possums.gl)[1:2],
                   verbose = 0)
  x <- x[c(1:4, 31:33), ]
  x <- gl.filter.callrate(x, threshold = 1, verbose = 0)
  gl.filter.monomorphs(x, verbose = 0)
}

sfs_counts <- function(x) {
  m <- as.matrix(x)
  p <- pop(x)
  list(n = nInd(x), cs = colSums(m, na.rm = TRUE),
       n1 = sum(p == levels(p)[1]), n2 = sum(p == levels(p)[2]),
       c1 = colSums(m[p == levels(p)[1], , drop = FALSE]),
       c2 = colSums(m[p == levels(p)[2], , drop = FALSE]))
}

test_that("single-population SFS matches an independent count [baseline]", {
  x <- sfs_2pop()
  k <- sfs_counts(x)
  s <- gl.sfs(x, singlepop = TRUE, plot.out = FALSE, verbose = 0)
  expect_equal(names(s), paste0("d", 0:k$n))
  expect_equal(unname(s), tabulate(pmin(k$cs, 2 * k$n - k$cs) + 1, k$n + 1))
  su <- gl.sfs(x, singlepop = TRUE, folded = FALSE, plot.out = FALSE,
               verbose = 0)
  expect_equal(unname(su), tabulate(k$cs + 1, 2 * k$n + 1))
  s2 <- gl.sfs(x, singlepop = TRUE, minbinsize = 2, plot.out = FALSE,
               verbose = 0)
  expect_equal(names(s2)[1], "d2")
})

test_that("two-population unfolded SFS matches an independent count [baseline]", {
  x <- sfs_2pop()
  k <- sfs_counts(x)
  s <- gl.sfs(x, folded = FALSE, plot.out = FALSE, verbose = 0)
  ref <- table(factor(k$c1, 0:(2 * k$n1)), factor(k$c2, 0:(2 * k$n2)))
  expect_equal(as.vector(dim(s)), c(2 * k$n1 + 1, 2 * k$n2 + 1))
  expect_true(all(unclass(s) == unclass(ref)))
})

test_that("two-population folded SFS folds on the combined minor allele [approved 2]", {
  x <- sfs_2pop()
  k <- sfs_counts(x)
  s <- gl.sfs(x, folded = TRUE, plot.out = FALSE, verbose = 0)
  # folding on the minor allele over both populations (fastsimcoal2 MAF)
  flip <- (k$c1 + k$c2) > (k$n1 + k$n2)
  g1 <- ifelse(flip, 2 * k$n1 - k$c1, k$c1)
  g2 <- ifelse(flip, 2 * k$n2 - k$c2, k$c2)
  ref <- table(factor(g1, 0:(2 * k$n1)), factor(g2, 0:(2 * k$n2)))
  expect_true(all(unclass(s) == unclass(ref)))
  expect_equal(sum(s), nLoc(x))
})

test_that("two-population minbinsize keeps private polymorphisms [approved 3]", {
  x <- sfs_2pop()
  k <- sfs_counts(x)
  s0 <- gl.sfs(x, folded = FALSE, minbinsize = 0, plot.out = FALSE,
               verbose = 0)
  s <- gl.sfs(x, folded = FALSE, minbinsize = 1, plot.out = FALSE,
              verbose = 0)
  expect_equal(sum(k$c1 == 0 | k$c2 == 0), 10)
  expect_equal(dim(s), dim(s0))
  expect_equal(sum(s), nLoc(x)) # no monomorphic loci in x
  s2 <- gl.sfs(x, folded = FALSE, minbinsize = 2, plot.out = FALSE,
               verbose = 0)
  tot <- outer(0:(2 * k$n1), 0:(2 * k$n2), `+`)
  expect_equal(sum(s2), sum(s0[tot >= 2]))
  expect_true(all(s2[tot < 2] == 0))
})

test_that("loci with missing calls are excluded and reported [approved 1]", {
  x <- gl.keep.pop(platypus.gl, pop.list = popNames(platypus.gl)[1],
                   verbose = 0)
  m <- as.matrix(x)
  full <- colSums(is.na(m)) == 0
  expect_equal(sum(!full), 237)
  out <- capture.output(
    s <- gl.sfs(x, singlepop = TRUE, plot.out = FALSE, verbose = 1)
  )
  expect_match(out, "237 of 1000 loci have missing calls", all = FALSE)
  cs <- colSums(m[, full])
  n <- nInd(x)
  expect_equal(sum(s), sum(full))
  expect_equal(unname(s), tabulate(pmin(cs, 2 * n - cs) + 1, n + 1))
})

test_that("plot.file saves without displaying the plot [approved 4]", {
  x <- sfs_2pop()
  f <- file.path(tempdir(), "sfs_test.RDS")
  unlink(f)
  s <- gl.sfs(x, singlepop = TRUE, minbinsize = 2, plot.out = FALSE,
              plot.file = "sfs_test", plot.dir = tempdir(), verbose = 0)
  expect_true(file.exists(f))
  gp <- readRDS(f)
  # bars at the true classes
  expect_equal(gp$data$names[1], 2)
  # a multidimensional sfs is returned, the save is skipped
  expect_true(is.array(
    gl.sfs(x, plot.out = FALSE, plot.file = "sfs_test2",
           plot.dir = tempdir(), verbose = 0)
  ))
  unlink(f)
})

test_that("SilicoDArT input is rejected [approved 5]", {
  expect_error(
    gl.sfs(testset.gs[1:20, 1:50], singlepop = TRUE, plot.out = FALSE,
           verbose = 0),
    "SilicoDArT"
  )
})

test_that("no population message is silent at verbose 0 [approved 6]", {
  x <- sfs_2pop()
  pop(x) <- NULL
  out <- capture.output(s <- gl.sfs(x, plot.out = FALSE, verbose = 0))
  expect_length(out, 0)
  expect_false(is.array(s))
})

test_that("too many dimensions stops with a message [approved 6]", {
  expect_error(
    gl.sfs(possums.gl, plot.out = FALSE, verbose = 0),
    "Cannot create a multidimensional sfs"
  )
})
