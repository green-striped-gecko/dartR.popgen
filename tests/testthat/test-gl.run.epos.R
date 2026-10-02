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
