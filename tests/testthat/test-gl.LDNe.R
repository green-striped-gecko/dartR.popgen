# Characterization baseline for gl.LDNe (dartr-function-review, 2026-09-18).
# Snapshots what the function does today on possums.gl[1:60, 1:100] (the
# roxygen example: two populations of 30, 100 SNPs) with NeEstimator V2.
# The pre-review baseline pinned the defects of the reviewed state; those
# expectations were flipped to the approved behaviour recorded in
# function-review/reports/dartR.popgen/gl.LDNe.md (changes 1 to 10). The
# estimates themselves are unchanged from the baseline.
#
# The whole file is skipped unless NEEST_DIR names a directory holding the
# NeEstimator executable for this OS (Ne2-1.exe / Ne2-1M / Ne2-1L).

neest_dir <- function() {
  d <- Sys.getenv("NEEST_DIR", "")
  prog <- switch(Sys.info()[["sysname"]],
                 Windows = "Ne2-1.exe", Darwin = "Ne2-1M", Linux = "Ne2-1L",
                 "")
  if (d == "" || prog == "" || !file.exists(file.path(d, prog))) return("")
  d
}

quiet_ldne <- function(...) {
  # the binary prints to the console; keep the test log readable
  res <- NULL
  invisible(capture.output(res <- suppressMessages(gl.LDNe(...))))
  res
}

stat_rows <- c("Lowest Allele Frequency Used", "Harmonic Mean Sample Size",
               "Independent Comparisons", "OverAll r^2",
               "Expected r^2 Sample", "Estimated Ne^", "CI low Parametric",
               "CI high Parametric", "CI low JackKnife", "CI high JackKnife")

test_that("two populations, critical c(0, 0.05): one table per population", {
  skip_if(neest_dir() == "", "NeEstimator binary not available (set NEEST_DIR)")
  pops <- possums.gl[1:60, 1:100]
  res <- quiet_ldne(pops, neest.path = neest_dir(), critical = c(0, 0.05),
                    plot.out = FALSE, verbose = 0)
  expect_type(res, "list")
  expect_named(res, c("A", "B"))
  expect_equal(dim(res$A), c(10, 4))
  expect_equal(res$A$Statistic, stat_rows)
  expect_equal(colnames(res$A),
               c("Statistic", "Frequency 1", "Frequency 2", "Frequency 3"))
  expect_equal(unname(unlist(res$A[1, -1])), c("0.050", "0+", "No S*"))
  expect_equal(as.numeric(res$A[6, -1]), c(13.5, 15, 14.2))
  expect_equal(as.numeric(res$A[7, -1]), c(11.5, 12.9, 12.2))
  expect_equal(as.numeric(res$A[9, -1]), c(8.3, 9.9, 9.4))
  expect_equal(as.numeric(res$B[6, -1]), c(16.1, 17.4, 16.1))
})

test_that("mating accepts the documented 'monogamy' and the abbreviation 'mono'", {
  skip_if(neest_dir() == "", "NeEstimator binary not available (set NEEST_DIR)")
  pops <- possums.gl[1:60, 1:100]
  mono <- quiet_ldne(pops, neest.path = neest_dir(), mating = "monogamy",
                     plot.out = FALSE, verbose = 0)
  expect_equal(as.numeric(mono$A[6, -1]), c(31.6, 30.1))
  abbrev <- quiet_ldne(pops, neest.path = neest_dir(), mating = "mono",
                       plot.out = FALSE, verbose = 0)
  expect_equal(as.numeric(abbrev$A[6, -1]), as.numeric(mono$A[6, -1]))
  random <- quiet_ldne(pops, neest.path = neest_dir(), mating = "random",
                       plot.out = FALSE, verbose = 0)
  expect_equal(as.numeric(random$A[6, -1]), c(15, 14.2))
  expect_error(quiet_ldne(pops, neest.path = neest_dir(), mating = "clonal",
                          plot.out = FALSE, verbose = 0),
               "'arg' should be one of")
})

test_that("verbose gates every message: nothing at 0, tables at 2", {
  skip_if(neest_dir() == "", "NeEstimator binary not available (set NEEST_DIR)")
  pops <- possums.gl[1:60, 1:100]
  quiet <- capture.output(res <- suppressMessages(
    gl.LDNe(pops, neest.path = neest_dir(), plot.out = FALSE, verbose = 0)))
  expect_length(quiet, 0)
  expect_named(res, c("A", "B"))
  loud <- capture.output(suppressMessages(
    gl.LDNe(pops, neest.path = neest_dir(), plot.out = FALSE, verbose = 2)))
  expect_true(any(grepl("Estimated Ne", loud)))
  expect_false(any(grepl("Completed: gl2genepop", loud)))
})

test_that("Waples correction applies eq 1a and rejects invalid arguments", {
  skip_if(neest_dir() == "", "NeEstimator binary not available (set NEEST_DIR)")
  pops <- possums.gl[1:60, 1:100]
  res <- quiet_ldne(pops, neest.path = neest_dir(), critical = c(0, 0.05),
                    Waples.correction = "nChromosomes",
                    Waples.correction.value = 10, plot.out = FALSE,
                    verbose = 0)
  expect_equal(nrow(res$A), 15)
  expect_equal(res$A$Statistic[11], "Waples' corrected Ne")
  expect_equal(as.numeric(res$A[11, -1]),
               round(as.numeric(res$A[6, -1]) / (0.098 + 0.219 * log(10)), 1))
  gl <- quiet_ldne(pops, neest.path = neest_dir(),
                   Waples.correction = "genomeLength",
                   Waples.correction.value = 3000, plot.out = FALSE,
                   verbose = 0)
  expect_equal(as.numeric(gl$A[11, -1]),
               round(as.numeric(gl$A[6, -1]) / (-0.910 + 0.219 * log(3000)), 1))
  expect_error(quiet_ldne(pops, neest.path = neest_dir(),
                          Waples.correction = "foo",
                          Waples.correction.value = 10, plot.out = FALSE,
                          verbose = 0),
               "can only be either 'nChromosomes' or")
  expect_error(quiet_ldne(pops, neest.path = neest_dir(),
                          Waples.correction = "nChromosomes",
                          plot.out = FALSE, verbose = 0),
               "single positive number")
  expect_error(quiet_ldne(pops, neest.path = neest_dir(),
                          Waples.correction = "nChromosomes",
                          Waples.correction.value = c(10, 20),
                          plot.out = FALSE, verbose = 0),
               "single positive number")
})

test_that("naive = TRUE appends a naive Ne row and keeps the population names", {
  skip_if(neest_dir() == "", "NeEstimator binary not available (set NEEST_DIR)")
  pops <- possums.gl[1:60, 1:100]
  res <- quiet_ldne(pops, neest.path = neest_dir(), critical = c(0, 0.05),
                    naive = TRUE, plot.out = FALSE, verbose = 0)
  expect_named(res, c("A", "B"))
  expect_equal(nrow(res$A), 11)
  expect_equal(res$A$Statistic[11], "Naive Estimated Ne^")
  expect_equal(as.numeric(res$A[11, -1]), c(14.2, 15.3, 14.7))
})

test_that("pairing = 'separate' restricts comparisons and an unknown value is rejected", {
  skip_if(neest_dir() == "", "NeEstimator binary not available (set NEEST_DIR)")
  pops <- possums.gl[1:60, 1:100]
  set.seed(1)
  pops@chromosome <- as.factor(sample(1:10, size = nLoc(pops), replace = TRUE))
  res <- quiet_ldne(pops, neest.path = neest_dir(), critical = c(0, 0.05),
                    pairing = "separate", plot.out = FALSE, verbose = 0)
  expect_equal(as.numeric(res$A[3, -1]), c(2986, 3526, 3363))
  expect_error(quiet_ldne(pops, neest.path = neest_dir(), pairing = "foo",
                          plot.out = FALSE, verbose = 0),
               "'arg' should be one of")
})

test_that("a population of two individuals keeps its own missing jackknife CIs wherever it sorts", {
  skip_if(neest_dir() == "", "NeEstimator binary not available (set NEEST_DIR)")
  # small population first: the defect this replaces gave A the values of B
  small_first <- rbind(possums.gl[1:2, 1:100], possums.gl[31:60, 1:100])
  res <- quiet_ldne(small_first, neest.path = neest_dir(),
                    critical = c(0, 0.05), plot.out = FALSE, verbose = 0)
  expect_equal(as.numeric(res$A[2, -1]), c(2, 2, 2))
  expect_equal(as.numeric(res$A[6, -1]), c(Inf, Inf, Inf))
  expect_true(all(is.na(res$A[9, -1])))
  expect_true(all(is.na(res$A[10, -1])))
  expect_equal(as.numeric(res$B[9, -1]), c(11, 11.7, 11))
  expect_equal(as.numeric(res$B[10, -1]), c(24.7, 27.6, 24.8))

  # small population second: unchanged from the baseline
  small_second <- rbind(possums.gl[1:30, 1:100], possums.gl[31:32, 1:100])
  res2 <- quiet_ldne(small_second, neest.path = neest_dir(),
                     critical = c(0, 0.05), plot.out = FALSE, verbose = 0)
  expect_equal(as.numeric(res2$A[9, -1]), c(8.3, 9.9, 9.4))
  expect_true(all(is.na(res2$B[9, -1])))

  # small population in the middle of three
  small_mid <- rbind(possums.gl[1:30, 1:100], possums.gl[31:32, 1:100],
                     possums.gl[61:90, 1:100])
  res3 <- quiet_ldne(small_mid, neest.path = neest_dir(), critical = 0,
                     plot.out = FALSE, verbose = 0)
  expect_equal(as.numeric(res3$A[9, -1][1]), 9.9)
  expect_true(all(is.na(res3$B[9, -1])))
  expect_equal(as.numeric(res3$C[9, -1][1]), 9.1)
})

test_that("plot.file saves without plot.out, and colours are taken one per population", {
  skip_if(neest_dir() == "", "NeEstimator binary not available (set NEEST_DIR)")
  pops <- possums.gl[1:60, 1:100]
  dir <- withr::local_tempdir()
  res <- quiet_ldne(pops, neest.path = neest_dir(), plot.out = FALSE,
                    plot.file = "ldne", plot.dir = dir, verbose = 0)
  expect_named(res, c("A", "B"))
  expect_true(file.exists(file.path(dir, "ldne.RDS")))
  expect_true(file.exists(file.path(dir, "ldne_tab.RDS")))

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  # more colours than populations: the extra colours are ignored
  many <- quiet_ldne(pops, neest.path = neest_dir(), critical = c(0, 0.05),
                     plot_colors_pop = c("red", "blue", "green", "black"),
                     verbose = 0)
  expect_named(many, c("A", "B"))
  # a palette function is accepted
  pal <- quiet_ldne(pops, neest.path = neest_dir(),
                    plot_colors_pop = function(n) rep("red", n), verbose = 0)
  expect_named(pal, c("A", "B"))
  # fewer colours than populations is reported as such
  expect_error(quiet_ldne(pops, neest.path = neest_dir(),
                          plot_colors_pop = "red", verbose = 0),
               "at least one colour per population")
})

test_that("the saved output file is refreshed by a second run", {
  skip_if(neest_dir() == "", "NeEstimator binary not available (set NEEST_DIR)")
  pops <- possums.gl[1:60, 1:100]
  od <- withr::local_tempdir()
  quiet_ldne(pops, neest.path = neest_dir(), critical = 0, outpath = od,
             outfile = "ne.txt", plot.out = FALSE, verbose = 0)
  s1 <- file.size(file.path(od, "ne.txt"))
  quiet_ldne(pops, neest.path = neest_dir(), critical = c(0, 0.05, 0.1),
             outpath = od, outfile = "ne.txt", plot.out = FALSE, verbose = 0)
  expect_gt(file.size(file.path(od, "ne.txt")), s1)
})

test_that("each call runs in its own directory and cleans it up", {
  skip_if(neest_dir() == "", "NeEstimator binary not available (set NEEST_DIR)")
  pops <- possums.gl[1:60, 1:100]
  before <- list.files(tempdir(), pattern = "^LDNe_")
  quiet_ldne(pops, neest.path = neest_dir(), plot.out = FALSE, verbose = 0)
  expect_equal(list.files(tempdir(), pattern = "^LDNe_"), before)
})

test_that("missing executable and SilicoDArT input stop with a readable message", {
  skip_if(neest_dir() == "", "NeEstimator binary not available (set NEEST_DIR)")
  pops <- possums.gl[1:60, 1:100]
  expect_error(quiet_ldne(pops, neest.path = withr::local_tempdir(),
                          plot.out = FALSE, verbose = 0),
               "Cannot find")
  expect_error(quiet_ldne(testset.gs[1:30, 1:50], neest.path = neest_dir(),
                          plot.out = FALSE, verbose = 0),
               "Only SNP")
})
