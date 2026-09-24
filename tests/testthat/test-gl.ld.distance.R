# Regression tests for gl.ld.distance (function-review campaign).
#
# Before the review (dartR.popgen 1.2.2) a resolution above the largest
# distance failed inside fields::stats.bin, verbose = 0 still printed the
# table, and palette functions for pop.colors failed.
#
# Input is a synthetic ld.report with the columns gl.report.ld.map returns
# that gl.ld.distance uses (pop, distance, ld.stat): two populations, pairs
# at known distances, so the binned means can be worked out by hand.

make_ld_report <- function() {
  data.frame(
    pop = rep(c("A", "B"), each = 6),
    distance = rep(c(100, 900, 1500, 1800, 2500, 3000), 2),
    ld.stat = c(0.8, 0.6, 0.4, 0.2, 0.1, 0.1,
                0.5, 0.3, 0.3, 0.1, 0.0, 0.2)
  )
}

run_quiet <- function(...) {
  out <- NULL
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  utils::capture.output(out <- suppressWarnings(gl.ld.distance(...)))
  out
}

test_that("bins are (1, res], ... labelled by their upper edge; means by pop", {
  res <- run_quiet(make_ld_report(), ld.resolution = 1000, verbose = 0)
  # breaks: 1, 1001, 2001, 3000
  expect_equal(res$distance, rep(c(1001, 2001, 3000), 2))
  expect_equal(as.character(res$pop), rep(c("A", "B"), each = 3))
  expect_equal(res$ld.stat, c(0.7, 0.3, 0.1, 0.4, 0.2, 0.1))
  expect_equal(res$n.pairs, rep(c(2L, 2L, 2L), 2))
})

test_that("finer bins give the hand-computed means and counts", {
  res <- run_quiet(make_ld_report(), ld.resolution = 500, verbose = 0)
  a <- res[res$pop == "A", ]
  # breaks 1, 501, 1001, 1501, 2001, 2501, 3000; one pair per bin, as the
  # old fields::stats.bin() binning gave on the same data
  expect_equal(a$distance, c(501, 1001, 1501, 2001, 2501, 3000))
  expect_equal(a$ld.stat, c(0.8, 0.6, 0.4, 0.2, 0.1, 0.1))
  expect_equal(a$n.pairs, rep(1L, 6))
})

test_that("empty bins give NA means and zero pairs", {
  res <- run_quiet(make_ld_report(), ld.resolution = 250, verbose = 0)
  a <- res[res$pop == "A", ]
  # (251, 501] and (501, 751] hold no pair
  expect_true(all(is.na(a$ld.stat[a$distance %in% c(501, 751)])))
  expect_equal(a$n.pairs[a$distance %in% c(501, 751)], c(0L, 0L))
  expect_true(all(a$n.pairs[is.na(a$ld.stat)] == 0))
  expect_equal(sum(a$n.pairs), 6)
})

test_that("a resolution above the largest distance gives one bin", {
  res <- run_quiet(make_ld_report(), ld.resolution = 1e6, verbose = 0)
  expect_equal(res$distance, c(3000, 3000))
  expect_equal(res$ld.stat, c(mean(c(0.8, 0.6, 0.4, 0.2, 0.1, 0.1)),
                              mean(c(0.5, 0.3, 0.3, 0.1, 0.0, 0.2))))
  expect_equal(res$n.pairs, c(6L, 6L))
})

test_that("no repeated last break when the sequence ends on the maximum", {
  res <- run_quiet(make_ld_report(), ld.resolution = 2999, verbose = 0)
  expect_equal(nrow(res), 2)
  expect_false(anyNA(res$ld.stat))
})

test_that("invalid input stops with clear errors", {
  expect_error(run_quiet(data.frame(a = 1), verbose = 0), "gl.report.ld.map")
  expect_error(run_quiet(make_ld_report(), ld.resolution = -5, verbose = 0),
               "positive")
  expect_error(run_quiet(make_ld_report(), ld.resolution = "1e5",
                         verbose = 0), "positive")
})

test_that("returns a data.table; verbose = 0 is silent, 3 prints the table", {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  expect_silent(res <- gl.ld.distance(make_ld_report(), ld.resolution = 1000,
                                      plot.out = FALSE, verbose = 0))
  expect_s3_class(res, "data.table")
  out <- utils::capture.output(
    gl.ld.distance(make_ld_report(), ld.resolution = 1000, plot.out = FALSE,
                   verbose = 3)
  )
  expect_true(any(grepl("n.pairs", out)))
})

test_that("pop.colors accepts a palette function; too few colours error", {
  expect_s3_class(
    run_quiet(make_ld_report(), ld.resolution = 1000, pop.colors = rainbow,
              verbose = 0),
    "data.table"
  )
  expect_error(
    run_quiet(make_ld_report(), ld.resolution = 1000, pop.colors = "blue",
              verbose = 0),
    "1 colors but there are 2 populations"
  )
})

test_that("plot.file saves a ggplot; the threshold line is in the legend", {
  dir <- withr::local_tempdir()
  run_quiet(make_ld_report(), ld.resolution = 1000, plot.file = "ldd",
            plot.dir = dir, verbose = 0)
  f <- list.files(dir, "ldd", full.names = TRUE)
  expect_length(f, 1)
  p <- readRDS(f)
  expect_s3_class(p, "ggplot")
  expect_setequal(ggplot2::get_guide_data(p, "colour")$.label, c("A", "B"))
  expect_equal(ggplot2::get_guide_data(p, "linetype")$.label,
               "R.squared = 0.2 (unlinked threshold)")
})
