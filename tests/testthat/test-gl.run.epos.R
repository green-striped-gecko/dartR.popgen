# Characterization tests for gl.run.epos -- Phase A baseline of the dartR
# function review (dartR.popgen at 1e90a1e), updated in Phase C: assertions
# tagged [approved n] were changed with the matching approved change; see
# function-review/reports/dartR.popgen/gl.run.epos.md (dartR.base).
#
# The real-run tests need the EPOS binaries: set EPOS_DIR to a folder with
# epos and epos2plot (and bootSfs for the bootstrap test), e.g. the folder
# made by gl.download.binary("epos").

ep_genlight <- function() {
  gl.filter.monomorphs(possums.gl[1:30, 1:200], verbose = 0)
}

ep_dir <- function(boot = FALSE) {
  d <- Sys.getenv("EPOS_DIR")
  progs <- c("epos", "epos2plot", if (boot) "bootSfs")
  if (.Platform$OS.type == "windows") progs <- paste0(progs, ".exe")
  if (!nzchar(d) || !all(file.exists(file.path(d, progs)))) return(NULL)
  normalizePath(d)
}

ep_run <- function(x, ...) {
  out <- NULL
  capture.output(
    out <- suppressWarnings(
      gl.run.epos(x, epos.path = ep_dir(), plot.display = FALSE,
                  outpath = tempdir(), verbose = 0, ...)
    )
  )
  out
}

test_that("missing binaries name the files and the download [approved 3]", {
  expect_error(
    capture.output(
      gl.run.epos(ep_genlight(), epos.path = tempfile("nope"), L = 1e5,
                  u = 1e-8, verbose = 0)
    ),
    "Cannot find epos(\\.exe)?, epos2plot(\\.exe)?.*gl.download.binary"
  )
})

test_that("unknown method stops [baseline]", {
  expect_error(
    gl.run.epos(ep_genlight(), epos.path = tempdir(), method = "foo",
                verbose = 0),
    "method must be one of"
  )
})

test_that("default run returns history, sfs and diagnostics [baseline]", {
  skip_if(is.null(ep_dir()), "EPOS_DIR not set")
  x <- ep_genlight()
  r <- ep_run(x, L = 1e5, u = 1e-8, seed = 1) # [approved 6]
  expect_named(r, c("history", "plot", "sfs", "diagnostics"))
  expect_equal(names(r$history), c("generation", "low", "median", "high"))
  expect_equal(r$history$median, c(1, 1, 29700, 29700))
  expect_equal(r$history$generation, c(0, 1.27, 1.27, 79200))
  expect_equal(r$sfs$r, 1:30)
  expect_equal(sum(r$sfs$fr), nLoc(x))
  # EPOS doubles the middle bin of a folded SFS on input: 176 + 4
  expect_equal(r$diagnostics$polymorphic_sites, 180)
  expect_equal(r$diagnostics$monomorphic_sites, 1e5 - 180)
  expect_equal(r$diagnostics$likelihood, -80, tolerance = 1e-3)
})

test_that("verbose = 0 is silent [approved 7]", {
  skip_if(is.null(ep_dir()), "EPOS_DIR not set")
  out <- capture.output(
    r <- gl.run.epos(ep_genlight(), epos.path = ep_dir(), L = 1e5, u = 1e-8,
                     plot.display = FALSE, outpath = tempdir(), verbose = 0)
  )
  expect_length(out, 0)
})

# expected values below come from epos run directly on correctly labelled
# SFS files (class 0 row; -x 1; -U with classes 1..59), seed 1
test_that("minbinsize = 0 uses the zero class instead of L [approved 1]", {
  skip_if(is.null(ep_dir()), "EPOS_DIR not set")
  r <- ep_run(possums.gl[1:30, 1:200], u = 1e-8, minbinsize = 0, seed = 1)
  expect_equal(r$diagnostics$monomorphic_sites, 24)
  expect_equal(r$history$median[3], 1.5e7)
  expect_equal(r$sfs$r[1], 0)
})

test_that("minbinsize = 2 excludes singletons with their true classes [approved 1]", {
  skip_if(is.null(ep_dir()), "EPOS_DIR not set")
  r <- ep_run(ep_genlight(), L = 1e5, u = 1e-8, minbinsize = 2, seed = 1)
  expect_equal(r$sfs$r[1], 2)
  expect_equal(r$diagnostics$sfs[[1]]$r[1], 2)
  expect_equal(max(r$diagnostics$sfs[[1]]$r), 30)
  expect_equal(r$history$median[3], 31100)
  expect_equal(max(r$history$generation), 83100)
})

test_that("folded = FALSE runs unfolded on classes 1..2n-1 [approved 1]", {
  skip_if(is.null(ep_dir()), "EPOS_DIR not set")
  r <- ep_run(ep_genlight(), L = 1e5, u = 1e-8, folded = FALSE, seed = 1)
  expect_equal(max(r$diagnostics$sfs[[1]]$r), 59)
  expect_equal(r$history$median, c(9440, 9440))
  expect_equal(r$history$generation, c(0, 37100))
})

test_that("L required, u optional, options need no leading space [approved 3]", {
  skip_if(is.null(ep_dir()), "EPOS_DIR not set")
  x <- ep_genlight()
  expect_error(ep_run(x, u = 1e-8), "L \\(sequence length\\) is required")
  expect_equal(nrow(ep_run(x, L = 1e5)$history), 4)
  expect_equal(nrow(ep_run(x, L = 1e5, u = 1e-8, other.options = "-m 10")$history), 4)
  # an epos failure stops with its exit status
  expect_error(ep_run(x, L = 1e5, u = 1e-8, other.options = "-L abc"),
               "epos failed \\(exit status")
})

test_that("SilicoDArT input is rejected [approved 4]", {
  expect_error(
    capture.output(
      gl.run.epos(testset.gs[1:20, 1:60], epos.path = tempdir(), L = 1e5,
                  verbose = 0)
    ),
    "SilicoDArT"
  )
})

test_that("an unnamed user sfs gives the same result [approved 5]", {
  skip_if(is.null(ep_dir()), "EPOS_DIR not set")
  x <- ep_genlight()
  s <- as.numeric(gl.sfs(x, minbinsize = 1, singlepop = TRUE,
                         plot.out = FALSE, verbose = 0))
  a <- ep_run(x, sfs = s, L = 1e5, u = 1e-8, seed = 1)
  b <- ep_run(x, L = 1e5, u = 1e-8, seed = 1)
  expect_identical(a$history, b$history)
})

test_that("bootstrap run returns one sfs per replicate [baseline]", {
  skip_if(is.null(ep_dir(boot = TRUE)), "EPOS_DIR without bootSfs")
  r <- ep_run(ep_genlight(), L = 1e5, u = 1e-8, boot = 5, seed = 3)
  expect_length(r$diagnostics$sfs, 5)
  expect_true(all(sapply(r$diagnostics$sfs, nrow) == 30))
  expect_s3_class(r$diagnostics$sfs_plot, "ggplot")
  expect_equal(names(r$history), c("generation", "low", "median", "high"))
  # [approved 6] same seed, same result
  r2 <- ep_run(ep_genlight(), L = 1e5, u = 1e-8, boot = 5, seed = 3)
  expect_identical(r$history, r2$history)
  # [approved 2] narrower quantiles give narrower limits
  r3 <- ep_run(ep_genlight(), L = 1e5, u = 1e-8, boot = 5, seed = 3,
               upper = 0.75, lower = 0.25)
  expect_false(identical(r$history, r3$history))
  expect_true(all(r3$history$low >= r$history$low))
  expect_true(all(r3$history$high <= r$history$high))
})

# ---- Follow-up review (2026-11-11), baseline at dartR.popgen 8262192 --------
# Baseline captured before any change, then updated in Phase C: assertions
# tagged [approved n] were changed with the matching approved change of the
# addendum in function-review/reports/dartR.popgen/gl.run.epos.md (dartR.base).

# one platypus population: 41 individuals, 484 polymorphic loci with call
# rate >= 0.95, 104 of them with at least one missing call
ep_missing <- function(polymorphic = TRUE) {
  x <- platypus.gl[pop(platypus.gl) == "TENTERFIELD", ]
  m <- as.matrix(x)
  poly <- apply(m, 2, function(g) length(unique(g[!is.na(g)])) > 1)
  x[, colMeans(!is.na(m)) >= 0.95 & (poly | !polymorphic)]
}

# folded SFS of 30 sequences, the EPOS example data/testF.dat
# (github.com/EvolBioInf/epos)
ep_testF <- c(9526, 3998, 2315, 1487, 1075, 873, 689, 570, 567, 627, 547,
              605, 460, 573, 236)

ep_out <- function(x, ...) {
  r <- NULL
  out <- capture.output(
    r <- suppressWarnings(
      gl.run.epos(x, epos.path = ep_dir(), plot.display = FALSE,
                  outpath = tempdir(), seed = 1, ...)
    )
  )
  list(r = r, out = out)
}

test_that("loci with missing calls are reported and L is scaled to the loci used [approved 9]", {
  skip_if(is.null(ep_dir()), "EPOS_DIR not set")
  x <- ep_missing()
  m <- as.matrix(x)
  complete <- colSums(is.na(m)) == 0
  expect_equal(c(nInd(x), nLoc(x), sum(!complete)), c(41, 484, 104))
  o <- ep_out(x, L = nLoc(x) * 69, u = 1e-8, verbose = 2)
  # [approved 9] was: no message
  expect_true(any(grepl("104 of 484 loci have missing calls", o$out)))
  expect_true(any(grepl("33396 x 380 / 484 = 26220", o$out)))
  expect_equal(sum(o$r$sfs$fr), 380)
  # [approved 9] was: 484 * 69 sites and Ne 8600, 8600, 76700, 76700
  expect_equal(o$r$diagnostics$monomorphic_sites +
                 o$r$diagnostics$polymorphic_sites, 380 * 69)
  expect_equal(o$r$history$median, c(11000, 11000, 97700, 97700))
  # the same history as the complete loci run with their own L
  ref <- ep_out(x[, complete], L = 380 * 69, u = 1e-8, verbose = 0)
  expect_identical(o$r$history, ref$r$history)
  # the message is shown at verbose = 1 and not at verbose = 0
  expect_true(any(grepl("loci have missing calls",
                        ep_out(x, L = nLoc(x) * 69, u = 1e-8, verbose = 1)$out)))
  expect_length(ep_out(x, L = nLoc(x) * 69, u = 1e-8, verbose = 0)$out, 0)
  # with the zero class instead of L there is nothing to scale (the
  # monomorphic loci of this population are kept for the zero class)
  z <- ep_out(ep_missing(polymorphic = FALSE), u = 1e-8, minbinsize = 0,
              verbose = 1)
  expect_true(any(grepl("loci have missing calls", z$out)))
  expect_false(any(grepl("L scaled", z$out)))
  expect_equal(z$r$sfs$r[1], 0)
  # a user-supplied sfs leaves L as given
  s <- gl.sfs(x, singlepop = TRUE, plot.out = FALSE, verbose = 0)
  u <- ep_out(x, sfs = s, L = nLoc(x) * 69, u = 1e-8, verbose = 2)
  expect_false(any(grepl("missing", u$out)))
  expect_equal(u$r$diagnostics$monomorphic_sites +
                 u$r$diagnostics$polymorphic_sites, 484 * 69)
})

test_that("complete data: no message and L as given [approved 9]", {
  skip_if(is.null(ep_dir()), "EPOS_DIR not set")
  o <- ep_out(ep_genlight(), L = 1e5, u = 1e-8, verbose = 2)
  expect_false(any(grepl("missing", o$out)))
  expect_equal(o$r$diagnostics$monomorphic_sites +
                 o$r$diagnostics$polymorphic_sites, 1e5)
})

test_that("an Ne = 1 step is returned with a warning [approved 10]", {
  skip_if(is.null(ep_dir()), "EPOS_DIR not set")
  o <- ep_out(ep_genlight(), L = 1e5, u = 1e-8, verbose = 2)
  # the history is unchanged
  expect_equal(o$r$history$median, c(1, 1, 29700, 29700))
  expect_equal(o$r$history$generation, c(0, 1.27, 1.27, 79200))
  # [approved 10] was: no message
  expect_true(any(grepl("Ne = 1 for generations 0 to 1.27", o$out)))
  expect_true(any(grepl("Ne = 1", ep_out(ep_genlight(), L = 1e5, u = 1e-8,
                                         verbose = 1)$out)))
  # no warning when no interval is 1
  n <- ep_out(ep_genlight(), sfs = ep_testF, L = 1e7, u = 1.2e-8, verbose = 2)
  expect_true(all(n$r$history$median > 1))
  expect_false(any(grepl("Ne = 1|lower limit", n$out)))
})

test_that("a bootstrap run warns when the median or the lower limit is 1 [approved 10]", {
  skip_if(is.null(ep_dir(boot = TRUE)), "EPOS_DIR without bootSfs")
  o <- ep_out(ep_genlight(), L = 1e5, u = 1e-8, boot = 20, verbose = 1)
  expect_true(any(o$r$history$low == 1))
  expect_true(any(grepl("Ne = 1 for generations|lower limit of Ne is 1",
                        o$out)))
})

test_that("exhaustive with a depth below 3 is the greedy search and says so [approved 11]", {
  skip_if(is.null(ep_dir()), "EPOS_DIR not set")
  x <- ep_genlight()
  g <- ep_out(x, sfs = ep_testF, L = 1e7, u = 1.2e-8, verbose = 2)
  e <- ep_out(x, sfs = ep_testF, L = 1e7, u = 1.2e-8, verbose = 2,
              method = "exhaustive")
  e3 <- ep_out(x, sfs = ep_testF, L = 1e7, u = 1.2e-8, verbose = 2,
               method = "exhaustive", depth = 3)
  # results are unchanged: the default depth is still 2
  expect_identical(g$r$history, e$r$history)
  expect_equal(g$r$diagnostics$likelihood, -106.137405)
  expect_equal(e3$r$diagnostics$likelihood, -99.505349)
  # [approved 11] was: no message
  expect_true(any(grepl("depth = 2 is the greedy search", e$out)))
  expect_false(any(grepl("greedy search", g$out)))
  expect_false(any(grepl("greedy search", e3$out)))
  expect_length(ep_out(x, sfs = ep_testF, L = 1e7, u = 1.2e-8, verbose = 0,
                       method = "exhaustive")$out, 0)
})

test_that("an unfolded sfs from gl.sfs gives the computed result [approved 13]", {
  skip_if(is.null(ep_dir()), "EPOS_DIR not set")
  x <- ep_genlight()
  su <- gl.sfs(x, folded = FALSE, singlepop = TRUE, plot.out = FALSE,
               verbose = 0)
  a <- ep_out(x, L = 1e5, u = 1e-8, folded = FALSE, verbose = 0)
  b <- ep_out(x, sfs = su, L = 1e5, u = 1e-8, folded = FALSE, verbose = 0)
  expect_equal(max(a$r$diagnostics$sfs[[1]]$r), 59)
  expect_equal(a$r$history$median, c(9440, 9440))
  # [approved 13] was: classes 1..60 sent, Ne 9400
  expect_equal(max(b$r$diagnostics$sfs[[1]]$r), 59)
  expect_identical(a$r$history, b$r$history)
  # the same for an unnamed vector of the classes 1..60
  d <- ep_out(x, sfs = as.numeric(su)[-1], L = 1e5, u = 1e-8, folded = FALSE,
              verbose = 0)
  expect_identical(a$r$history, d$r$history)
  # an unfolded sfs of another sample size is sent as given
  f <- ep_out(x, sfs = c(ep_testF, rev(ep_testF)[-1]), L = 1e7, u = 1.2e-8,
              folded = FALSE, verbose = 0)
  expect_equal(max(f$r$diagnostics$sfs[[1]]$r), 29)
})

test_that("several populations are pooled with a message [approved 12]", {
  skip_if(is.null(ep_dir()), "EPOS_DIR not set")
  x <- gl.filter.monomorphs(possums.gl[c(1:15, 31:45), ], verbose = 0)
  expect_equal(nPop(x), 2)
  # [approved 12] was: no message
  o <- ep_out(x, L = 1e5, u = 1e-8, verbose = 2)
  expect_true(any(grepl("x has 2 populations", o$out)))
  expect_false(any(grepl("populations", ep_out(x, L = 1e5, u = 1e-8,
                                               verbose = 1)$out)))
  # one population: no message
  expect_false(any(grepl("populations", ep_out(ep_genlight(), L = 1e5,
                                               u = 1e-8, verbose = 2)$out)))
})
