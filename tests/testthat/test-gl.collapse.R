# Regression tests for gl.collapse (function-review campaign).
#
# Before the review (dartR.popgen 1.2.2) the grouping made one pass, so a
# chain of similar populations could leave a population in two groups and
# gl.merge.pop() stopped ("not present in the dataset"); a single group made
# gl.fixed.diff() stop; and the returned matrices were recomputed with
# gl.collapse's tloc (default 0) whatever tloc built 'fd'.
#
# Grouping is tested on hand-made matrices through the internal helper
# utils.collapse.groups(); the rest end to end on testset.gl.

groups <- function(m, tpop) {
  vapply(utils.collapse.groups(m, tpop), paste, "", collapse = ",")
}

sym <- function(v, nm) {
  m <- matrix(0, length(nm), length(nm), dimnames = list(nm, nm))
  m[lower.tri(m)] <- v
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}

quiet <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressWarnings(expr))
  out
}

test_that("groups: A-B and C-D at fd 0, E apart", {
  m <- sym(c(0, 5, 5, 5,
             5, 5, 5,
             0, 5,
             5), LETTERS[1:5])
  expect_equal(groups(m, 0), c("A,B", "C,D", "E"))
  expect_equal(groups(m, 5), "A,B,C,D,E")
})

test_that("groups follow chains that the old single pass missed", {
  nm <- c("L", "S", "M", "Q", "A", "Z", "I", "R", "B")
  m <- matrix(c(0,5,3,0,8,8,0,3,8,
                5,0,0,5,5,3,8,3,0,
                3,0,0,0,5,0,3,3,5,
                0,5,0,0,5,3,8,8,0,
                8,5,5,5,0,8,8,8,3,
                8,3,0,3,8,0,5,5,8,
                0,8,3,8,8,5,0,5,8,
                3,3,3,8,8,5,5,0,8,
                8,0,5,0,3,8,8,8,0), 9, byrow = TRUE,
              dimnames = list(nm, nm))
  expect_equal(groups(m, 0), c("A", "B,I,L,M,Q,S,Z", "R"))
})

test_that("groups equal connected components on random matrices", {
  comp <- function(m, t) {
    n <- nrow(m)
    g <- seq_len(n)
    repeat {
      ch <- FALSE
      for (i in 1:n) for (j in 1:n) {
        if (m[i, j] <= t && g[i] != g[j]) {
          g[g == max(g[i], g[j])] <- min(g[i], g[j])
          ch <- TRUE
        }
      }
      if (!ch) break
    }
    sort(vapply(split(rownames(m), g), function(v) {
      paste(sort(v), collapse = ",")
    }, ""))
  }
  withr::with_seed(2, {
    for (r in 1:100) {
      n <- sample(4:9, 1)
      m <- matrix(sample(c(0, 3, 5, 8, 8), n * n, TRUE), n)
      m[lower.tri(m)] <- t(m)[lower.tri(m)]
      diag(m) <- 0
      rownames(m) <- colnames(m) <- sample(LETTERS, n)
      expect_equal(sort(groups(m, 0)), unname(comp(m, 0)))
    }
  })
})

test_that("NA distances do not join populations", {
  m <- sym(c(NA, 5, 12), c("A", "B", "C"))  # A-B NA, A-C 5, B-C 12
  expect_equal(groups(m, 10), c("A,C", "B"))
})

test_that("testset.gl example: tpop = 1 collapses to three groups", {
  fd <- quiet(gl.fixed.diff(dartR.data::testset.gl, tloc = 0.05, verbose = 0))
  fd2 <- quiet(gl.collapse(fd, tpop = 1, tloc = 0.05, verbose = 0))
  expect_s3_class(fd2, "fd")
  expect_equal(popNames(fd2$gl),
               c("EmmacBrisWive+", "EmsubRopeMata", "EmvicVictJasp"))
  expect_equal(nInd(fd2$gl), nInd(dartR.data::testset.gl))
  expect_equal(names(fd2), c("gl", "fd", "pcfd", "nobs", "nloc", "expfpos",
                             "sdfpos", "pval"))
  # returned matrices use the same tloc as the input
  at05 <- quiet(gl.fixed.diff(fd2$gl, tloc = 0.05, verbose = 0))
  expect_equal(as.matrix(fd2$fd), as.matrix(at05$fd))
})

test_that("a tloc that does not match fd stops with a clear error", {
  fd <- quiet(gl.fixed.diff(dartR.data::testset.gl, tloc = 0.05, verbose = 0))
  expect_error(quiet(gl.collapse(fd, tpop = 1, verbose = 0)),
               "do not match those computed with tloc = 0")
})

test_that("all populations in one group returns a one-population fd", {
  fd <- quiet(gl.fixed.diff(dartR.data::testset.gl, tloc = 0, verbose = 0))
  o <- quiet(gl.collapse(fd, tpop = 100, verbose = 0))
  expect_s3_class(o, "fd")
  expect_equal(nPop(o$gl), 1)
  expect_equal(popNames(o$gl), "EmmacBrisWive+")
  expect_equal(nInd(o$gl), nInd(dartR.data::testset.gl))
  expect_equal(dim(as.matrix(o$fd)), c(1L, 1L))
})

test_that("no amalgamation returns the input fd unchanged", {
  x <- dartR.data::testset.gl
  x <- dartR.base::gl.keep.pop(x, pop.list = c("EmsubRopeMata",
                                                "EmvicVictJasp",
                                                "EmmacCoopAvin"),
                               verbose = 0)
  fd <- quiet(gl.fixed.diff(x, tloc = 0, verbose = 0))
  expect_true(all(as.matrix(fd$fd)[lower.tri(as.matrix(fd$fd))] > 0))
  o <- quiet(gl.collapse(fd, tpop = 0, verbose = 0))
  expect_identical(o, fd)
})

test_that("invalid arguments stop with clear errors", {
  fd <- quiet(gl.fixed.diff(dartR.data::testset.gl, tloc = 0, verbose = 0))
  expect_error(quiet(gl.collapse(fd, tpop = "1", verbose = 0)), "tpop")
  expect_error(quiet(gl.collapse(fd, tpop = -1, verbose = 0)), "tpop")
  expect_error(quiet(gl.collapse(fd, tloc = 0.7, verbose = 0)), "tloc")
  expect_error(quiet(gl.collapse(fd, tloc = c(0, 0.1), verbose = 0)), "tloc")
})

test_that("verbose = 0 is silent", {
  fd <- quiet(gl.fixed.diff(dartR.data::testset.gl, tloc = 0, verbose = 0))
  expect_silent(gl.collapse(fd, tpop = 1, verbose = 0))
})
