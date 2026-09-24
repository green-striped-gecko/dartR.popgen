# Regression tests for gl.outflank (function-review campaign).
#
# Before the review (dartR.popgen 1.2.2) every SNP entered the analysis
# twice (one genind column per allele), which doubled the outlier counts and
# moved df and q-values away from the OutFLANK package; all-NA loci were
# dropped, shifting the index; names with a dot failed; there was no verbose.
#
# Data: Balding-Nichols simulation, 5 populations x 20 diploids, 600
# neutral loci (F = 0.05) and 10 selected loci (F = 0.5), 5% missing calls.

sim_outflank_gl <- function(seed = 1) {
  withr::with_seed(seed, {
    npop <- 5; nper <- 20; L0 <- 600; L1 <- 10
    L <- L0 + L1
    Fv <- c(rep(0.05, L0), rep(0.5, L1))
    p0 <- stats::runif(L, 0.1, 0.9)
    G <- matrix(NA_integer_, npop * nper, L)
    for (l in seq_len(L)) {
      for (k in seq_len(npop)) {
        pk <- stats::rbeta(1, p0[l] * (1 - Fv[l]) / Fv[l],
                           (1 - p0[l]) * (1 - Fv[l]) / Fv[l])
        G[((k - 1) * nper + 1):(k * nper), l] <- stats::rbinom(nper, 2, pk)
      }
    }
    G[sample(length(G), 0.05 * length(G))] <- NA
  })
  nm <- c(sprintf("neutral-%d-A/G", seq_len(L0)), sprintf("sel-%d-A/G", seq_len(L1)))
  x <- new("genlight", as.data.frame(G), ploidy = 2,
           ind.names = sprintf("ind%d", seq_len(npop * nper)), loc.names = nm,
           pop = factor(rep(paste0("P", seq_len(npop)), each = nper)))
  suppressWarnings(dartR.base::gl.compliance.check(x, verbose = 0))
}

run_quiet <- function(...) {
  out <- NULL
  utils::capture.output(out <- suppressWarnings(gl.outflank(...)))
  out
}

test_that("returns index and outflank list; selected loci are flagged", {
  skip_if_not_installed("qvalue")
  x <- sim_outflank_gl()
  res <- run_quiet(x, plot = FALSE)
  expect_named(res, c("index", "outflank"))
  expect_length(res$index, nLoc(x))
  o <- res$outflank
  expect_equal(nrow(o$results), nLoc(x))
  expect_equal(o$results$LocusName, locNames(x))
  flagged <- o$results$LocusName[which(o$results$OutlierFlag)]
  expect_true(all(grepl("^sel-", flagged)))
  expect_gte(length(flagged), 5)
})

test_that("index is TRUE for non-outliers, FALSE for outliers, NA untested", {
  skip_if_not_installed("qvalue")
  x <- sim_outflank_gl()
  res <- run_quiet(x, plot = FALSE)
  expect_identical(res$index, !res$outflank$results$OutlierFlag)
})

test_that("results match the OutFLANK package exactly", {
  skip_if_not_installed("qvalue")
  skip_if_not_installed("OutFLANK")
  x <- sim_outflank_gl()
  o <- run_quiet(x, plot = FALSE)$outflank
  g <- as.matrix(x)
  g[is.na(g)] <- 9
  utils::capture.output(
    fm <- OutFLANK::MakeDiploidFSTMat(g, locNames(x), as.character(pop(x)))
  )
  suppressMessages(library(qvalue))
  up <- OutFLANK::OutFLANK(fm, LeftTrimFraction = 0.05,
                           RightTrimFraction = 0.05, Hmin = 0.1,
                           NumberOfSamples = nPop(x), qthreshold = 0.05)
  expect_equal(o$results$FSTNoCorr, up$results$FSTNoCorr)
  expect_equal(o$FSTbar, up$FSTbar)
  expect_equal(o$dfInferred, up$dfInferred)
  expect_equal(o$results$qvalues, up$results$qvalues)
  expect_equal(o$results$OutlierFlag, up$results$OutlierFlag)
  expect_equal(o$numberHighFstOutliers, up$numberHighFstOutliers)
  expect_equal(o$numberHighFstOutliers, sum(o$results$OutlierFlag, na.rm = TRUE))
})

test_that("meanAlleleFreq is the reference-allele frequency", {
  skip_if_not_installed("qvalue")
  x <- sim_outflank_gl()
  o <- run_quiet(x, plot = FALSE)$outflank
  ref_freq <- 1 - colMeans(as.matrix(x), na.rm = TRUE) / 2
  expect_equal(unname(o$results$meanAlleleFreq), unname(ref_freq))
})

test_that("all-NA loci stay in the output as NA, in input order", {
  skip_if_not_installed("qvalue")
  x <- sim_outflank_gl()
  m <- as.matrix(x)
  m[, 5] <- NA
  y <- new("genlight", as.data.frame(m), ploidy = 2, ind.names = indNames(x),
           loc.names = locNames(x), pop = pop(x))
  y <- suppressWarnings(dartR.base::gl.compliance.check(y, verbose = 0))
  res <- run_quiet(y, plot = FALSE)
  expect_length(res$index, nLoc(y))
  expect_equal(res$outflank$results$LocusName, locNames(y))
  expect_true(is.na(res$index[5]))
})

test_that("locus names with dots work", {
  skip_if_not_installed("qvalue")
  x <- sim_outflank_gl()
  locNames(x) <- sub("neutral-", "scaf.", locNames(x))
  res <- run_quiet(x, plot = FALSE)
  expect_equal(res$outflank$results$LocusName, locNames(x))
})

test_that("genind input gives the same result as the genlight", {
  skip_if_not_installed("qvalue")
  x <- sim_outflank_gl()
  gi <- suppressWarnings(dartR.base::gl2gi(x, verbose = 0))
  a <- run_quiet(gi, plot = FALSE)
  b <- run_quiet(x, plot = FALSE)
  expect_equal(a$index, b$index)
})

test_that("invalid input stops with clear errors", {
  skip_if_not_installed("qvalue")
  x <- sim_outflank_gl()
  one <- x
  pop(one) <- factor(rep("A", nInd(one)))
  expect_error(run_quiet(one, plot = FALSE), "at least two populations")
  nopop <- x
  pop(nopop) <- NULL
  expect_error(run_quiet(nopop, plot = FALSE), "assigned to a population")
  expect_error(run_quiet(dartR.data::testset.gs, plot = FALSE), "SilicoDArT")
  expect_error(gl.outflank(x, plot = FALSE, NumberOfSamples = 2),
               "unused argument")
})

test_that("verbose = 0 is silent; verbose = 3 reports a summary", {
  skip_if_not_installed("qvalue")
  x <- sim_outflank_gl()
  expect_silent(gl.outflank(x, plot = FALSE, verbose = 0))
  out <- utils::capture.output(gl.outflank(x, plot = FALSE, verbose = 3))
  expect_true(any(grepl("Loci tested:", out)))
  expect_true(any(grepl("Inferred df:", out)))
  expect_true(any(grepl("Completed:", out)))
})

test_that("plot = TRUE draws without changing the result", {
  skip_if_not_installed("qvalue")
  x <- sim_outflank_gl()
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  a <- run_quiet(x, plot = TRUE, Hmin = 0.2)
  b <- run_quiet(x, plot = FALSE, Hmin = 0.2)
  expect_identical(a$index, b$index)
})
