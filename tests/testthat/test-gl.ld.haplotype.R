# Characterization tests for gl.ld.haplotype -- Phase A baseline of the
# dartR function review (dartR.popgen at 31fb40c), updated in Phase C:
# assertions tagged [approved Fn] were flipped with the matching approved
# finding. See function-review/reports/dartR.popgen/gl.ld.haplotype.md
# (dartR.base repo) for the findings.

skip_if_not_installed("snpStats")
skip_if_not_installed("raster")
skip_if_not_installed("sp")
skip_if_not_installed("scales")
skip_if_not_installed("viridis")

# The documented example: TENTERFIELD, first 15 individuals, chromosome slots
# filled from the platypus locus metrics.
ld_fixture <- function() {
  x <- gl.filter.callrate(platypus.gl, threshold = 1, verbose = 0)
  x <- gl.keep.pop(x, pop.list = "TENTERFIELD", verbose = 0)[1:15, ]
  x$chromosome <- as.factor(x$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1)
  x$position <- x$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
  x
}
chr1 <- "NC_041728.1_chromosome_1"
chr2 <- "NC_041729.1_chromosome_2"
haplo_cols <- c(
  "population", "chromosome", "haplotype", "start", "end",
  "start_ld_plot", "end_ld_plot", "midpoint", "midpoint_ld_plot", "labels"
)

# Quiet the deprecated-fortify warning (ggplot2 >= 3.4.4, see report F15) so
# the test log stays readable; the pins below are on values.
run_quiet <- function(expr) {
  suppressWarnings(suppressMessages(
    capture.output(res <- withVisible(expr))
  ))
  res
}

test_that("default call on the documented example returns an empty 10-column table, invisibly, and writes nothing", {
  withr::local_options(dartR_wd = NULL)
  pd <- withr::local_tempdir()
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- ld_fixture()
  expect_equal(nLoc(x), 569)
  expect_equal(nInd(x), 15)

  r <- run_quiet(gl.ld.haplotype(x, chrom_name = chr1, ld_max_pairwise = 1e7,
                                 plot.dir = pd, verbose = 0))
  expect_false(r$visible)
  res <- r$value
  expect_s3_class(res, "data.frame")
  expect_equal(dim(res), c(0L, 10L))
  expect_named(res, haplo_cols)

  # [approved F2] plot.save = FALSE (the default) writes no PDF
  expect_false(file.exists(file.path(pd, paste0("TENTERFIELD_", chr1, ".pdf"))))
  # [approved F5] the PLINK intermediates are removed after reading
  expect_length(list.files(tempdir(), pattern = "^gl_plink_"), 0)
})

test_that("plot.save controls the PDF in both branches", {
  withr::local_options(dartR_wd = NULL)
  pd <- withr::local_tempdir()
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- ld_fixture()
  f <- file.path(pd, paste0("TENTERFIELD_", chr1, ".pdf"))
  # [approved F2] no haplotype branch: nothing written unless asked
  run_quiet(gl.ld.haplotype(x, chrom_name = chr1, plot.out = FALSE,
                            plot.save = FALSE, plot.dir = pd, verbose = 0))
  expect_false(file.exists(f))
  run_quiet(gl.ld.haplotype(x, chrom_name = chr1, plot.out = FALSE,
                            plot.save = TRUE, plot.dir = pd, verbose = 0))
  expect_true(file.exists(f))
  unlink(f)
  # [approved F2] haplotype branch: written when asked
  run_quiet(gl.ld.haplotype(x, chrom_name = chr1, haplo_id = TRUE,
                            min_snps = 2, ld_threshold_haplo = 0.5,
                            plot.out = FALSE, plot.save = TRUE,
                            plot.dir = pd, verbose = 0))
  expect_true(file.exists(f))
})

test_that("haplotypes follow the adjacent-pair LD threshold and min_snps", {
  withr::local_options(dartR_wd = NULL)
  pd <- withr::local_tempdir()
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- ld_fixture()
  # [approved F1] On chromosome 1 (22 SNPs after the MAF filter) a single
  # adjacent pair has r2 >= 0.5 (0.72, SNPs 165892285 and 171211062), so
  # min_snps = 3 finds nothing; the pre-fix code reported one block from
  # 10300000 to 166000000.
  res3 <- run_quiet(gl.ld.haplotype(x, chrom_name = chr1, haplo_id = TRUE,
                                    min_snps = 3, ld_threshold_haplo = 0.5,
                                    plot.out = FALSE, plot.dir = pd,
                                    verbose = 0))$value
  expect_equal(nrow(res3), 0L)
  expect_named(res3, haplo_cols)

  # min_snps = 2 finds exactly that pair, with exact bp coordinates
  # [approved F1, F10]
  res2 <- run_quiet(gl.ld.haplotype(x, chrom_name = chr1, haplo_id = TRUE,
                                    min_snps = 2, ld_threshold_haplo = 0.5,
                                    plot.out = FALSE, plot.dir = pd,
                                    verbose = 0))$value
  expect_equal(nrow(res2), 1L)
  expect_equal(res2$population, "TENTERFIELD")
  expect_equal(res2$chromosome, chr1)
  expect_equal(res2$haplotype, 1L)
  expect_equal(res2$start, 165892285)
  expect_equal(res2$end, 171211062)
  expect_equal(res2$midpoint, (165892285 + 171211062) / 2)
  expect_equal(res2$labels, "166-171")
  expect_true(res2$start_ld_plot < res2$end_ld_plot)
  expect_equal(res2$midpoint_ld_plot,
               (res2$start_ld_plot + res2$end_ld_plot) / 2)

  # a threshold above the largest adjacent r2 finds nothing
  res73 <- run_quiet(gl.ld.haplotype(x, chrom_name = chr1, haplo_id = TRUE,
                                     min_snps = 2, ld_threshold_haplo = 0.73,
                                     plot.out = FALSE, plot.dir = pd,
                                     verbose = 0))$value
  expect_equal(nrow(res73), 0L)

  # threshold 0: every adjacent pair qualifies, one block over the whole
  # chromosome from the first to the last retained SNP
  res0 <- run_quiet(gl.ld.haplotype(x, chrom_name = chr1, haplo_id = TRUE,
                                    min_snps = 2, ld_threshold_haplo = 0,
                                    plot.out = FALSE, plot.dir = pd,
                                    verbose = 0))$value
  expect_equal(nrow(res0), 1L)
  expect_equal(res0$start, 10320334)
  expect_equal(res0$end, 179168634)
  expect_equal(res0$start_ld_plot, 0)
})

test_that("haplotype blocks agree with an independent adjacent-pair r2 computation", {
  withr::local_options(dartR_wd = NULL)
  pd <- withr::local_tempdir()
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- ld_fixture()
  # replicate the function's per-population preparation for chromosome 1
  y <- gl.keep.loc(x, loc.list = locNames(x)[x$chromosome == chr1], verbose = 0)
  y <- y[, order(y$position)]
  y <- gl.recalc.metrics(y, verbose = 0)
  y <- gl.filter.maf(y, threshold = 0.05, verbose = 0)
  y <- y[, !duplicated(y$position)]
  g <- as.matrix(y)
  sm <- methods::new("SnpMatrix",
    matrix(as.raw(ifelse(is.na(g), 0, g + 1)), nrow = nrow(g),
           dimnames = list(indNames(y), as.character(y$position))))
  adj <- as.matrix(snpStats::ld(sm, depth = 1, stats = "R.squared"))
  adj <- adj[cbind(1:(ncol(g) - 1), 2:ncol(g))]
  thr <- 0.3
  runs <- rle(!is.na(adj) & adj >= thr)
  run_end <- cumsum(runs$lengths)
  run_start <- run_end - runs$lengths + 1
  expected <- data.frame(
    start = y$position[run_start[runs$values]],
    end = y$position[run_end[runs$values] + 1]
  )
  expected <- expected[(run_end[runs$values] - run_start[runs$values] + 2) >= 2, ]
  expect_gt(nrow(expected), 1)

  res <- run_quiet(gl.ld.haplotype(x, chrom_name = chr1, haplo_id = TRUE,
                                   min_snps = 2, ld_threshold_haplo = thr,
                                   plot.out = FALSE, plot.dir = pd,
                                   verbose = 0))$value
  expect_equal(res$start, expected$start)
  expect_equal(res$end, expected$end)
  expect_equal(res$haplotype, seq_len(nrow(expected)))
})

test_that("verbose = 0 is silent", {
  withr::local_options(dartR_wd = NULL)
  pd <- withr::local_tempdir()
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- ld_fixture()
  # [approved F7, F8] no warnings, no table, no 'no haplotypes' note
  out <- suppressWarnings(suppressMessages(capture.output(
    gl.ld.haplotype(x, chrom_name = chr1, haplo_id = TRUE, min_snps = 3,
                    plot.out = FALSE, plot.dir = pd, verbose = 0)
  )))
  expect_length(out, 0)
})

test_that("a population with exactly ind.limit individuals is analysed", {
  withr::local_options(dartR_wd = NULL)
  pd <- withr::local_tempdir()
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- ld_fixture()[1:10, ]
  # [approved F11]
  out <- suppressWarnings(suppressMessages(capture.output(
    res <- gl.ld.haplotype(x, chrom_name = chr1, ind.limit = 10,
                           plot.out = FALSE, plot.dir = pd, verbose = 2)
  )))
  expect_false(any(grepl("Skipping population", out)))
  expect_true(any(grepl("Analysing chromosome", out)))
  expect_equal(nrow(res), 0L)
  # one fewer individual is skipped, at verbose 1 too; with no population
  # left the function stops instead of returning an empty table
  out1 <- character(0)
  expect_error(
    out1 <- capture.output(
      gl.ld.haplotype(x[1:9, ], chrom_name = chr1, ind.limit = 10,
                      plot.out = FALSE, plot.dir = pd, verbose = 1),
      type = "output"
    ),
    "No population was analysed.*TENTERFIELD \\(9 individuals"
  )
  expect_error(
    gl.ld.haplotype(x[1:9, ], chrom_name = chr1, ind.limit = 10,
                    plot.out = FALSE, plot.dir = pd, verbose = 0),
    "fewer than ind.limit = 10"
  )
})

test_that("skipping some populations warns and analyses the rest", {
  withr::local_options(dartR_wd = NULL)
  pd <- withr::local_tempdir()
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- ld_fixture()
  x2 <- x[1:5, ]
  pop(x2) <- factor(rep("SMALL", 5))
  x <- rbind(x, x2)
  x$chromosome <- as.factor(x$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1)
  x$position <- x$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
  w <- NULL
  out <- suppressMessages(capture.output(
    res <- withCallingHandlers(
      gl.ld.haplotype(x, chrom_name = chr1, ind.limit = 10,
                      plot.out = FALSE, plot.dir = pd, verbose = 2),
      warning = function(cnd) {
        if (grepl("Populations skipped", conditionMessage(cnd))) {
          w <<- conditionMessage(cnd)
        }
        invokeRestart("muffleWarning")
      }
    )
  ))
  expect_match(w, "SMALL \\(5 individuals")
  expect_false(grepl("TENTERFIELD", w))
  expect_true(any(grepl("Calculating pairwise LD in population TENTERFIELD", out)))
})

test_that("the all-chromosome default runs on the example and skips small chromosomes", {
  withr::local_options(dartR_wd = NULL)
  pd <- withr::local_tempdir()
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- ld_fixture()
  # [approved F3] chromosome X2 (2 SNPs after the MAF filter) aborted the
  # run before; it is now skipped with a warning
  out <- suppressWarnings(suppressMessages(capture.output(
    res <- gl.ld.haplotype(x, plot.out = FALSE, plot.dir = pd, verbose = 1)
  )))
  expect_equal(dim(res), c(0L, 10L))
  expect_true(any(grepl("Skipping chromosome NC_041750.1_chromosome_X2", out)))
})

test_that("input validation stops with a clear message", {
  withr::local_options(dartR_wd = NULL)
  x <- ld_fixture()
  # [approved F12]
  expect_error(
    suppressWarnings(suppressMessages(capture.output(
      gl.ld.haplotype(testset.gl, verbose = 0)
    ))),
    "x@chromosome"
  )
  expect_error(
    suppressWarnings(suppressMessages(capture.output(
      gl.ld.haplotype(testset.gs, verbose = 0)
    ))),
    "SilicoDArT"
  )
  # [approved F3]
  expect_error(
    suppressWarnings(suppressMessages(capture.output(
      gl.ld.haplotype(x, chrom_name = "nope", verbose = 0)
    ))),
    "not found in x@chromosome"
  )
  expect_error(
    suppressWarnings(suppressMessages(capture.output(
      gl.ld.haplotype(x, pop_name = "nope", verbose = 0)
    ))),
    "not found in x"
  )
})

test_that("the function runs after gl.set.wd() and leaves no PLINK files behind", {
  pd <- withr::local_tempdir()
  wd <- withr::local_tempdir()
  withr::local_options(dartR_wd = wd)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- ld_fixture()
  # [approved F5] the PLINK files used to be written to dartR_wd and read
  # from tempdir()
  res <- run_quiet(gl.ld.haplotype(x, chrom_name = c(chr1, chr2),
                                   plot.out = FALSE, plot.save = TRUE,
                                   plot.dir = pd, verbose = 0))$value
  expect_equal(nrow(res), 0L)
  expect_true(file.exists(file.path(pd, paste0("TENTERFIELD_", chr1, ".pdf"))))
  expect_true(file.exists(file.path(pd, paste0("TENTERFIELD_", chr2, ".pdf"))))
  expect_length(list.files(wd, pattern = "gl_plink"), 0)
  expect_length(list.files(tempdir(), pattern = "^gl_plink_"), 0)
})
