# Tests for gl.run.bayescan.
#
# The real-run tests need the BayeScan executable: set BAYESCAN_DIR to the
# folder holding it (e.g. the binaries folder of the BayeScan 2.1
# distribution, or a folder with an executable named "bayescan").

bs_dir <- function() {
  d <- Sys.getenv("BAYESCAN_DIR")
  if (!nzchar(d) || !dir.exists(d)) return(NULL)
  normalizePath(d)
}

# Two populations, 80 loci, with three edited loci at known positions:
# locus 5 monomorphic, locus 10 a fixed difference between the populations,
# and locus 30 without calls in population B
bs_genlight <- function() {
  x <- gl.keep.pop(possums.gl[, 1:80], pop.list = c("A", "B"), verbose = 0)
  m <- as.matrix(x)
  p <- as.character(pop(x))
  m[, 5] <- 0
  m[, 10] <- ifelse(p == "A", 0, 2)
  m[p == "B", 30] <- NA
  x@gen <- methods::new("genlight", m, ploidy = 2)@gen
  x
}

bs_run <- function(x, ...) {
  out <- NULL
  capture.output(
    out <- gl.run.bayescan(x, bayescan.path = bs_dir(), n.iter = 1000,
                           thin = 5, n.pilot = 5, pilot.length = 1000,
                           burnin = 2000, n.cores = 2, outpath = tempdir(),
                           plot.display = FALSE, verbose = 0, ...)
  )
  out
}

test_that("a missing executable names the download page", {
  expect_error(
    gl.run.bayescan(bs_genlight(), bayescan.path = tempfile("nope"),
                    verbose = 0),
    "does not exist.*cmpg.unibe.ch"
  )
  empty <- tempfile("empty")
  dir.create(empty)
  expect_error(
    gl.run.bayescan(bs_genlight(), bayescan.path = empty, verbose = 0),
    "Cannot find a BayeScan executable"
  )
})

test_that("one population, SilicoDArT data and bad settings stop", {
  x <- bs_genlight()
  expect_error(
    gl.run.bayescan(gl.keep.pop(x, pop.list = "A", verbose = 0),
                    bayescan.path = tempdir(), verbose = 0),
    "at least two populations"
  )
  expect_error(
    gl.run.bayescan(testset.gs, bayescan.path = tempdir(), verbose = 0)
  )
  expect_error(
    gl.run.bayescan(x, bayescan.path = tempdir(), burnin = -1, verbose = 0),
    "burnin must be whole numbers"
  )
  expect_error(
    gl.run.bayescan(x, bayescan.path = tempdir(), q.threshold = 1,
                    verbose = 0),
    "q.threshold must be"
  )
})

test_that("results have one row per locus, mapped to the right loci", {
  skip_if(is.null(bs_dir()), "BAYESCAN_DIR not set")
  x <- bs_genlight()
  r <- bs_run(x)
  expect_named(r, c("results", "fst.pop", "trace", "plot", "trace.plot"))
  expect_equal(nrow(r$results), nLoc(x))
  expect_equal(r$results$locus, locNames(x))
  # the monomorphic locus is not tested
  expect_true(all(is.na(r$results[5, -1])))
  expect_equal(sum(is.na(r$results$prob)), 1)
  # the fixed difference is the strongest outlier
  expect_equal(r$results$selection[10], factor("diversifying",
               levels = c("neutral", "diversifying", "balancing")))
  expect_equal(which.max(r$results$fst), 10)
  # the classification follows q.threshold and the sign of alpha
  t <- r$results[!is.na(r$results$q.value), ]
  expect_true(all(t$selection[t$q.value > 0.05] == "neutral"))
  expect_true(all(t$selection[t$q.value <= 0.05 & t$alpha > 0] ==
                    "diversifying"))
  expect_equal(r$fst.pop$population, c("A", "B"))
  expect_equal(colnames(r$trace), c("logL", "A", "B"))
  expect_s3_class(r$plot, "ggplot")
})

test_that("the trace plot has one panel per traced parameter", {
  skip_if(is.null(bs_dir()), "BAYESCAN_DIR not set")
  r <- bs_run(bs_genlight(), plot.file = "bs_trace_test",
              plot.dir = tempdir())
  expect_s3_class(r$trace.plot, "ggplot")
  d <- r$trace.plot$data
  expect_equal(nrow(d), 3 * nrow(r$trace))
  expect_equal(nlevels(d$parameter), 3)
  # the values of each panel are the columns of the trace, in order
  expect_equal(d$value[d$parameter == levels(d$parameter)[2]], r$trace$A)
  # iterations run from the end of the burn-in in steps of thin
  it <- unique(d$iteration)
  expect_true(all(diff(it) == 5) && it[1] > 2000)
  expect_true(file.exists(file.path(tempdir(), "bs_trace_test_trace.RDS")))
  expect_true("ess" %in% names(r$fst.pop))
})

test_that("the effective sample sizes match coda", {
  skip_if(is.null(bs_dir()), "BAYESCAN_DIR not set")
  skip_if_not_installed("coda")
  r <- bs_run(bs_genlight())
  ess <- coda::effectiveSize(coda::mcmc(as.matrix(r$trace)))
  expect_equal(r$fst.pop$ess, unname(ess[-1]))
  expect_match(levels(r$trace.plot$data$parameter)[1],
               paste0("^Log-likelihood \\(ESS ", round(ess[1]), "\\)$"))
})

test_that("the results are the _fst.txt file BayeScan wrote", {
  skip_if(is.null(bs_dir()), "BAYESCAN_DIR not set")
  r <- bs_run(bs_genlight(), outfile = "bs_test")
  raw <- read.table(file.path(tempdir(), "bs_test_fst.txt"), skip = 1)
  tested <- r$results[!is.na(r$results$prob), ]
  expect_equal(tested$q.value, raw[[4]])
  expect_equal(tested$alpha, raw[[5]])
  expect_true(file.exists(file.path(tempdir(), "bs_test_input.txt")))
})

test_that("a stricter q.threshold never adds outliers", {
  skip_if(is.null(bs_dir()), "BAYESCAN_DIR not set")
  r <- bs_run(bs_genlight(), q.threshold = 0.001)
  t <- r$results[!is.na(r$results$q.value), ]
  expect_true(all(t$selection[t$q.value > 0.001] == "neutral"))
})
