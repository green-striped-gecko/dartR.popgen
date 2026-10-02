# Characterization tests for gl.run.snmf, gl.plot.snmf and gl.map.snmf --
# Phase A baseline of the dartR function review (dartR.popgen at a16da26),
# updated in Phase C: assertions tagged [approved n] were flipped with the
# matching approved change; see
# function-review/reports/dartR.popgen/gl.run.snmf.md (dartR.base repo).
#
# sNMF runs through the LEA package, so these tests need no external binary.

sn_genlight <- function() {
  x <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1:3],
                   verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  gl.filter.callrate(x, threshold = 0.9, verbose = 0)
}

# one sNMF run shared by the tests of this file
sn_cache <- new.env()
sn_run <- function() {
  if (is.null(sn_cache$r)) {
    x <- sn_genlight()
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    utils::capture.output(sn_cache$r <- suppressWarnings(
      gl.run.snmf(x, minK = 1, maxK = 3, rep = 2, plot.out = FALSE,
                  verbose = 0, seed = 42)))
    sn_cache$x <- x
  }
  sn_cache
}

sn_quiet <- function(expr) {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  suppressWarnings(suppressMessages(expr))
}

sn_calls <- function(m, method) {
  cl <- m$x$calls
  cl[vapply(cl, `[[`, "", "method") == method]
}

test_that("run: list with best runs, cross-entropy plot and Q per K", {
  s <- sn_run()
  r <- s$r
  expect_named(r, c("best_run", "cross_entropy", "matrix"))
  expect_named(r$matrix, c("K1", "K2", "K3"))
  expect_named(r$matrix$K3, c("Pop_1", "Pop_2", "Pop_3", "Cluster", "Pop",
                              "Label", "Order"))
  expect_setequal(r$matrix$K3$Label, indNames(s$x))
  expect_equal(unname(rowSums(r$matrix$K3[, 1:3])), rep(1, nInd(s$x)),
               tolerance = 1e-3)
  expect_s3_class(r$cross_entropy[[1]], "recordedplot")
  expect_length(r$best_run, 3)
  # [approved 2] files removed by default; best_run holds run names
  expect_match(r$best_run, "^K[1-3]/run[12]$")
})

test_that("run: gl2geno writes individuals in indNames(x) order, so Q rows match their labels", {
  x <- sn_genlight()
  d <- withr::local_tempdir()
  utils::capture.output(gl2geno(x, outpath = d, outfile = "t", verbose = 0))
  g <- readLines(file.path(d, "t.geno"))
  G <- do.call(rbind, lapply(strsplit(g, ""), as.integer))
  G[G == 9] <- NA
  expect_true(isTRUE(all.equal(G, t(as.matrix(x)), check.attributes = FALSE)))
})

test_that("run: cleanup removes the files; verbose 0 is silent [approved 2, 4]", {
  x <- sn_genlight()
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  before <- list.files(tempdir())
  out <- utils::capture.output(r <- suppressWarnings(
    gl.run.snmf(x, minK = 2, maxK = 2, plot.out = FALSE, verbose = 0,
                cleanup = TRUE)))
  expect_length(out, 0)
  expect_length(setdiff(list.files(tempdir()), before), 0)
  out2 <- utils::capture.output(r2 <- suppressWarnings(
    gl.run.snmf(x, minK = 2, maxK = 2, plot.out = FALSE, verbose = 0,
                cleanup = FALSE)))
  expect_true(dir.exists(r2$best_run))
})

test_that("run: argument checks [approved 5]", {
  x <- sn_genlight()
  expect_error(gl.run.snmf(x, minK = 3, maxK = 2, verbose = 0),
               "minK <= maxK")
  expect_error(gl.run.snmf(x, minK = 1, maxK = 2, rep = 0, verbose = 0),
               "rep must be")
  expect_false(any(grepl("build = ", deparse(body(gl.run.snmf)),
                         fixed = TRUE)))
})

test_that("run: plot.file goes to plot.dir, not the working directory [approved 3]", {
  x <- sn_genlight()
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  withr::local_options(dartR_wd = NULL)
  pd <- withr::local_tempdir()
  utils::capture.output(r <- suppressWarnings(
    gl.run.snmf(x, minK = 1, maxK = 2, plot.out = FALSE, plot.file = "ce",
                verbose = 0)))
  expect_false(file.exists(file.path(wd, "ce.RDS")))
  expect_true(file.exists(file.path(tempdir(), "ce.RDS")))
  unlink(file.path(tempdir(), "ce.RDS"))
  utils::capture.output(suppressWarnings(
    gl.run.snmf(x, minK = 1, maxK = 2, plot.out = FALSE, plot.file = "ce",
                plot.dir = pd, verbose = 0)))
  expect_true(file.exists(file.path(pd, "ce.RDS")))
  # plot.out = FALSE still returns the recorded plot
  expect_s3_class(r$cross_entropy[[1]], "recordedplot")
})

test_that("plot: bars and axis labels agree; plot.K must be one K [approved 6]", {
  s <- sn_run()
  q <- sn_quiet(gl.plot.snmf(s$r, plot.K = 3, verbose = 0))
  p <- ggplot2::last_plot()
  b <- ggplot2::ggplot_build(p)
  brks <- unlist(lapply(b$layout$panel_params, function(z) z$x$get_breaks()))
  labs <- unlist(lapply(b$layout$panel_params, function(z) {
    as.character(z$x$get_labels())
  }))
  bars <- unique(p$data[, c("Order", "Label")])
  expect_equal(as.character(bars$Label), labs[match(bars$Order, brks)])
  expect_error(sn_quiet(gl.plot.snmf(s$r, plot.K = 2:3, verbose = 0)),
               "plot.K must be one of the K values")
})

test_that("plot: verbose follows gl.set.verbosity; silent at 0 [approved 6]", {
  s <- sn_run()
  expect_null(formals(gl.plot.snmf)$verbose)
  out <- utils::capture.output(q <- sn_quiet(gl.plot.snmf(s$r, plot.K = 3,
                                                          verbose = 0)))
  expect_length(out, 0)
})

test_that("plot: clear input errors, palette function [approved 6]; hclust on distances, aes() [approved 7]", {
  s <- sn_run()
  expect_error(sn_quiet(gl.plot.snmf(s$r, plot.K = 5, verbose = 0)),
               "plot.K must be one of the K values in snmf.result: 1, 2, 3")
  q <- sn_quiet(gl.plot.snmf(s$r, plot.K = 3, color.clusters = rainbow,
                             verbose = 0))
  expect_equal(nrow(q), nInd(s$x))
  expect_error(sn_quiet(gl.plot.snmf(s$r, plot.K = 3, color.clusters = "red",
                                     verbose = 0)),
               "color.clusters has 1 colours but 3 are needed")
  src <- deparse(body(gl.plot.snmf))
  expect_false(any(grepl("aes_(", src, fixed = TRUE)))
  expect_false(any(grepl("stats::dist(res)", src, fixed = TRUE)))
  # den = TRUE: bars follow hclust of the Manhattan distances
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  suppressWarnings(gl.plot.snmf(s$r, plot.K = 3, den = TRUE, x = s$x,
                                verbose = 0))
  p <- ggplot2::last_plot()
  d <- gl.dist.ind(s$x, method = "Manhattan", plot.display = FALSE,
                   verbose = 0)
  dd <- stats::reorder(stats::as.dendrogram(stats::hclust(d)), TRUE,
                       agglo.FUN = mean)
  bars <- unique(p[[1]]$data[, c("Order", "Label")])
  bars <- bars[order(bars$Order), ]
  expect_equal(as.character(bars$Label), indNames(s$x)[order.dendrogram(dd)])
})

test_that("map: non-alphabetical pop levels, bars at their own centre [approved 1]", {
  skip_if_not_installed("leaflet")
  s <- sn_run()
  x <- s$x
  pop(x) <- factor(as.character(pop(x)),
                   levels = c("EmmacBurnBara", "EmmacBrisWive",
                              "EmmacBurdMist"))
  m <- sn_quiet(gl.map.snmf(x, s$r$matrix$K3))
  rc <- sn_calls(m$map, "addRectangles")
  lat <- vapply(rc, function(z) z$args[[1]], 0)
  grp <- rep(names(m$Q_name), vapply(m$Q_name, nrow, 0L) * 3)
  cen <- apply(x@other$latlon, 2, function(v) tapply(v, pop(x), mean))
  expect_equal(as.vector(tapply(lat, grp, min)[names(m$Q_name)]),
               unname(cen[names(m$Q_name), "lat"]), tolerance = 1e-8)
  expect_named(m, c("Q_name", "map"))
  expect_true(all(c("Pop_1", "Label", "Pop_name") %in% names(m$Q_name[[1]])))
  expect_equal(sum(vapply(m$Q_name, nrow, 0L)), nInd(x))
})

test_that("map: movepops by name; extra pops in x ignored [approved 1]", {
  skip_if_not_installed("leaflet")
  s <- sn_run()
  x <- s$x
  mp <- data.frame(lon = c(1, 0, 0), lat = c(0, 0, 0))
  m <- sn_quiet(gl.map.snmf(x, s$r$matrix$K3, movepops = mp))
  z <- sn_calls(m$map, "addMarkers")[[1]]
  cen <- apply(x@other$latlon, 2, function(v) tapply(v, pop(x), mean))
  expect_equal(z$args[[2]][1], unname(cen[1, "lon"] + 1))
  x4 <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1:4],
                    verbose = 0)
  m4 <- sn_quiet(gl.map.snmf(x4, s$r$matrix$K3, plot.out = FALSE,
                             verbose = 0))
  expect_named(m4$Q_name, popNames(x))
  out <- utils::capture.output(m0 <- sn_quiet(gl.map.snmf(x, s$r$matrix$K3,
                                                    plot.out = FALSE,
                                                    verbose = 0)))
  expect_length(out, 0)
})
