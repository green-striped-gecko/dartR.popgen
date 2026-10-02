# Characterization tests for gl.evanno and utils.structure.evanno -- Phase A
# baseline of the dartR function review (dartR.popgen at 95fde36), updated
# in Phase C: assertions tagged [approved n] were flipped with the matching
# approved change; see
# function-review/reports/dartR.popgen/gl.evanno.md (dartR.base repo).
#
# structure.result objects are built by hand with controlled LnP(K)
# values, so no STRUCTURE executable is needed.

# lnp: named list, K -> est.ln.prob of each replicate
ev_sr <- function(lnp) {
  runs <- list()
  for (k in names(lnp)) {
    for (r in seq_along(lnp[[k]])) {
      lab <- paste0("k", k, ".r", r)
      runs[[lab]] <- list(
        summary = c(k = as.numeric(k), est.ln.prob = lnp[[k]][r],
                    mean.lnL = 0, var.lnL = 0),
        q.mat = NULL, prior.anc = NULL, label = lab)
    }
  }
  class(runs) <- c("structure.result", "list")
  runs
}

ev_lnp <- list("1" = c(-1000, -1002, -998), "2" = c(-800, -805, -795),
               "3" = c(-780, -790, -770), "4" = c(-775, -785, -765))

ev_quiet <- function(expr) {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  suppressWarnings(suppressMessages(expr))
}

test_that("delta K matches a hand calculation on consecutive K", {
  r <- ev_quiet(gl.evanno(ev_sr(ev_lnp), plot.out = FALSE, verbose = 0))
  expect_named(r, c("df", "plots"))
  expect_named(r$df, c("k", "reps", "mean.ln.k", "sd.ln.k", "ln.pk",
                       "ln.ppk", "delta.k"))
  m <- vapply(ev_lnp, mean, 0)
  s <- vapply(ev_lnp, sd, 0)
  expect_equal(r$df$k, 1:4)
  expect_equal(r$df$reps, rep(3, 4))
  expect_equal(r$df$mean.ln.k, unname(m))
  expect_equal(r$df$sd.ln.k, unname(s))
  expect_equal(r$df$ln.pk, c(NA, unname(diff(m))))
  expect_equal(r$df$ln.ppk, c(NA, abs(diff(diff(unname(m)))), NA))
  expect_equal(r$df$delta.k,
               unname(c(NA, abs(m[3] - 2 * m[2] + m[1]) / s[2],
                        abs(m[4] - 2 * m[3] + m[2]) / s[3], NA)))
  # [approved 6] the combined figure is returned
  expect_named(r$plots, c("mean.ln.k", "ln.pk", "ln.ppk", "delta.k",
                          "combined"))
  expect_s3_class(r$plots$delta.k, "ggplot")
  expect_s3_class(r$plots$combined, "patchwork")
})

test_that("gl.evanno returns what utils.structure.evanno returns", {
  sr <- ev_sr(ev_lnp)
  a <- ev_quiet(gl.evanno(sr, plot.out = FALSE, verbose = 0))
  b <- ev_quiet(utils.structure.evanno(sr, plot = FALSE))
  expect_equal(a$df, b$df)
})

test_that("K values 8:11 keep numeric order", {
  r <- ev_quiet(gl.evanno(ev_sr(setNames(ev_lnp, 8:11)), plot.out = FALSE, verbose = 0))
  expect_equal(r$df$k, 8:11)
})

test_that("K values with gaps: statistics only from true neighbours [approved 2]", {
  sr <- ev_sr(setNames(ev_lnp, c(1, 3, 4, 5)))
  out <- capture.output(r <- ev_quiet(gl.evanno(sr, plot.out = FALSE,
                                                verbose = 1)))
  expect_equal(r$df$k, c(1, 3, 4, 5))
  m <- vapply(ev_lnp, mean, 0)
  s <- vapply(ev_lnp, sd, 0)
  # K = 3 has no K - 1: LnP'(3), LnP''(3) and delta K(3) are NA
  expect_equal(r$df$ln.pk, unname(c(NA, NA, m[3] - m[2], m[4] - m[3])))
  expect_equal(r$df$ln.ppk, unname(c(NA, NA, abs(m[4] - 2 * m[3] + m[2]),
                                     NA)))
  expect_equal(r$df$delta.k, unname(c(NA, NA,
                                      abs(m[4] - 2 * m[3] + m[2]) / s[3],
                                      NA)))
  expect_true(any(grepl("not consecutive.*K = 3", out)))
  r7 <- ev_quiet(gl.evanno(ev_sr(setNames(ev_lnp, c(1, 3, 5, 7))),
                           plot.out = FALSE, verbose = 0))
  expect_true(all(is.na(r7$df$delta.k)))
  expect_true(all(is.na(r7$df$ln.pk)))
})

test_that("identical replicates give delta K = NA with a warning [approved 1]", {
  lnp <- lapply(ev_lnp, function(v) rep(v[1], 3))
  out <- capture.output(r <- ev_quiet(gl.evanno(ev_sr(lnp),
                                                plot.out = FALSE,
                                                verbose = 1)))
  expect_equal(r$df$sd.ln.k, rep(0, 4))
  expect_true(all(is.na(r$df$delta.k)))
  expect_true(any(grepl("sd = 0\\) for K = 2, 3", out)))
})

test_that("one replicate per K: delta K NA with a warning [approved 3]", {
  lnp <- lapply(ev_lnp, `[`, 1)
  out <- capture.output(r <- ev_quiet(gl.evanno(ev_sr(lnp),
                                                plot.out = FALSE,
                                                verbose = 1)))
  expect_true(any(grepl("at least two replicates per K", out)))
  expect_true(all(is.na(r$df$delta.k)))
  expect_named(r$plots, c("mean.ln.k", "ln.pk", "ln.ppk", "combined"))
  out0 <- capture.output(r0 <- ev_quiet(gl.evanno(ev_sr(lnp),
                                                  plot.out = FALSE,
                                                  verbose = 0)))
  expect_length(out0, 0)
})

test_that("input errors [approved 4]", {
  expect_error(ev_quiet(gl.evanno(ev_sr(ev_lnp[1:2]), plot.out = FALSE,
                                  verbose = 0)),
               "needs at least three values of K; 2 found")
  expect_error(ev_quiet(gl.evanno(list(), plot.out = FALSE, verbose = 0)),
               "not a structure.result object")
  expect_error(ev_quiet(utils.structure.evanno(list(), plot = FALSE)),
               "not a structure.result object")
})

test_that("plot.out = TRUE prints the patchwork; plot.file saves it [approved 5, 6]", {
  sr <- ev_sr(ev_lnp)
  r <- ev_quiet(gl.evanno(sr, verbose = 0))
  expect_named(r, c("df", "plots"))
  d <- withr::local_tempdir()
  ev_quiet(gl.evanno(sr, plot.out = FALSE, plot.file = "ev", plot.dir = d,
                     verbose = 0))
  expect_equal(list.files(d), "ev.RDS")
  # [approved 5] dartR theme on every panel
  r2 <- ev_quiet(gl.evanno(sr, plot.out = FALSE, verbose = 0))
  expect_equal(r2$plots$mean.ln.k$theme$panel.background,
               theme_dartR()$panel.background)
  out <- capture.output(r3 <- ev_quiet(gl.evanno(sr, plot.out = FALSE,
                                                 verbose = 3)))
  expect_true(any(grepl("Largest delta K at K = 2", out)))
})

test_that("new arguments; no aes_string() [approved 5, 7]", {
  expect_named(formals(gl.evanno), c("sr", "plot.out", "plot.theme",
                                     "plot.dir", "plot.file", "verbose"))
  src <- deparse(body(utils.structure.evanno))
  expect_false(any(grepl("aes_string(", src, fixed = TRUE)))
  expect_false(any(grepl("gridExtra", src, fixed = TRUE)))
})
