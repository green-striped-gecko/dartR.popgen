# Characterization tests for gl.run.stairway2 -- Phase A baseline of the
# dartR function review (dartR.popgen at 7ddbca9), updated in Phase C:
# assertions tagged [approved n] were changed with the matching approved
# change; see function-review/reports/dartR.popgen/gl.run.stairway2.md
# (dartR.base).
#
# The run = FALSE tests need no external software. The real-run tests need
# Java and the Stairway Plot 2 folder `stairway_plot_es`: set STAIRWAY2_DIR to
# the directory that contains it (e.g. ~/programs).

sw_genlight <- function() {
  gl.filter.monomorphs(possums.gl[1:30, 1:200], verbose = 0)
}

sw_dir <- function() {
  d <- Sys.getenv("STAIRWAY2_DIR")
  if (!nzchar(d) || !dir.exists(file.path(d, "stairway_plot_es")) ||
      !nzchar(Sys.which("java"))) {
    return(NULL)
  }
  normalizePath(d)
}

# run = FALSE still checks for the binary folder, so a fake one is enough
sw_fake_dir <- function() {
  d <- file.path(tempfile("swfake"))
  dir.create(file.path(d, "stairway_plot_es"), recursive = TRUE)
  d
}

test_that("run = FALSE returns the run folder with the blueprint [approved 2]", {
  x <- sw_genlight()
  out <- capture.output(
    r <- gl.run.stairway2(x, L = 1e5, mu = 1e-8,
                          stairway2.path = sw_fake_dir(), run = FALSE,
                          seed = 1, verbose = 0)
  )
  expect_length(out, 0) # [approved 4] gl.sfs no longer prints at verbose 0
  expect_named(r, c("history", "plot", "run.dir"))
  expect_null(r$history)
  expect_null(r$plot)
  expect_true(dir.exists(r$run.dir))
  expect_true(dir.exists(file.path(r$run.dir, "stairway_plot_es")))
  # [approved 1] a subfolder of tempdir(), not tempdir() itself
  expect_equal(normalizePath(dirname(r$run.dir)), normalizePath(tempdir()))
  bp <- readLines(file.path(r$run.dir, "blueprint"))
  expect_true(any(grepl("^nseq: 60 ", bp)))
  expect_true(any(grepl("^L: 1e\\+05 ", bp)))
  expect_true(any(grepl("^mu: 1e-08 ", bp)))
  expect_true(any(grepl("^random_seed: 1$", bp)))
  expect_true(any(grepl("^nrand: 14 29 44 58 ", bp)))
  expect_true(any(grepl("^largest_size_of_SFS_bin_used_for_estimation: 30 ",
                        bp)))
  expect_true(any(grepl("^ninput: 200 ", bp)))
  # SFS line = gl.sfs(folded, minbinsize = 1): bins 1..nInd
  sfs_line <- grep("^SFS: ", bp, value = TRUE)
  sfs_vals <- as.numeric(strsplit(sub(" #.*", "", sub("^SFS: ", "", sfs_line)),
                                  " ")[[1]])
  expect_length(sfs_vals, nInd(x))
  expect_equal(sum(sfs_vals), nLoc(x))
})

test_that("missing mu stops with a message [baseline]", {
  expect_error(
    capture.output(
      gl.run.stairway2(sw_genlight(), L = 1e5, stairway2.path = sw_fake_dir(),
                       run = FALSE, verbose = 0)
    ),
    "Mutation rate per site per generation not specified"
  )
})

test_that("missing binary folder names the folder and the download [approved 7]", {
  expect_error(
    capture.output(
      gl.run.stairway2(sw_genlight(), L = 1e5, mu = 1e-8,
                       stairway2.path = tempfile("nope"), verbose = 0)
    ),
    "Cannot find the folder stairway_plot_es.*gl.download.binary"
  )
})

test_that("SilicoDArT is rejected at every verbosity [approved 3]", {
  gs <- testset.gs[1:20, 1:50]
  for (v in c(0, 2)) {
    expect_error(
      capture.output(
        gl.run.stairway2(gs, L = 1e5, mu = 1e-8,
                         stairway2.path = sw_fake_dir(), run = FALSE,
                         verbose = v)
      ),
      "SilicoDArT"
    )
  }
})

test_that("a failing Stairway run stops with a clear error [approved 6]", {
  skip_if(is.null(sw_dir()), "STAIRWAY2_DIR / java not available")
  # the fake folder has no Java classes, so Stairbuilder fails
  expect_error(
    capture.output(
      gl.run.stairway2(sw_genlight(), L = 1e5, mu = 1e-8,
                       stairway2.path = sw_fake_dir(), verbose = 0)
    ),
    "could not create the run script"
  )
})

test_that("real run returns the Stairway summary table [approved 8]", {
  d <- sw_dir()
  skip_if(is.null(d), "STAIRWAY2_DIR / java not available")
  x <- sw_genlight()
  capture.output(
    res <- gl.run.stairway2(x, L = 1e5, mu = 1e-8, stairway2.path = d,
                            nreps = 2, seed = 1, cleanup = FALSE,
                            plot.display = FALSE, verbose = 0)
  )
  expect_named(res, c("history", "plot", "run.dir")) # [approved 2]
  expect_s3_class(res$plot, "ggplot")
  expect_true(dir.exists(res$run.dir)) # cleanup = FALSE keeps it
  h <- res$history
  expect_equal(ncol(h), 11)
  expect_equal(names(h),
               c("mutation_per_site", "n_estimation", "theta_per_site_median",
                 "theta_per_site_2.5", "theta_per_site_97.5", "year",
                 "Ne_median", "low95", "high95", "low75", "high75"))
  expect_equal(nrow(h), 233)
  # year = mutation_per_site / mu * gentime
  expect_equal(h$year, h$mutation_per_site / 1e-8, tolerance = 1e-6)
  unlink(res$run.dir, recursive = TRUE)
})

test_that("cleanup, plot.file and parallel [approved 1, 4, 5]", {
  d <- sw_dir()
  skip_if(is.null(d), "STAIRWAY2_DIR / java not available")
  marker <- file.path(tempdir(), "sw_marker.txt")
  writeLines("x", marker)
  on.exit(unlink(c(marker, file.path(tempdir(), "swplot.RDS"))), add = TRUE)
  plan_before <- class(future::plan())
  # Stairway Plot 2's own output goes to the console directly, not through R,
  # so capture.output only sees R messages; at verbose = 0 there are none
  out <- capture.output(
    r <- gl.run.stairway2(sw_genlight(), L = 1e5, mu = 1e-8,
                          stairway2.path = d, nreps = 2, seed = 1,
                          plot.display = FALSE, plot.dir = tempdir(),
                          plot.file = "swplot", parallel = 2, verbose = 0)
  )
  expect_length(out, 0)                    # [approved 4]
  expect_true(file.exists(marker))         # [approved 1]
  expect_null(r$run.dir)                   # [approved 2]
  expect_true(file.exists(file.path(tempdir(), "swplot.RDS"))) # [approved 1]
  expect_equal(nrow(r$history), 233)
  expect_identical(class(future::plan()), plan_before) # [approved 5]
})
