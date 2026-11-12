# Characterization tests for gl.run.popcluster, gl.plot.popcluster and
# gl.map.popcluster -- Phase A baseline of the dartR function review
# (dartR.popgen at 7ddbca9), updated in Phase C: assertions tagged
# [approved n] were flipped with the matching approved change; see
# function-review/reports/dartR.popgen/gl.run.popcluster.md (dartR.base).
#
# Run tests need the PopCluster binary: POPCLUSTER_DIR = folder holding
# PopClusterMac / PopClusterLnx. Plot and map tests use a hand-built run
# object in the layout gl.run.popcluster returns.

pc_genlight <- function() {
  x <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1:3],
                   verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  gl.filter.callrate(x, threshold = 0.9, verbose = 0)
}

pc_dir <- function() {
  d <- Sys.getenv("POPCLUSTER_DIR", unset = "")
  if (!nzchar(d)) skip("POPCLUSTER_DIR not set")
  path.expand(d)
}

# the layout of gl.run.popcluster()$matrix[[run]], K = 3
pc_q <- function(x, seed = 1) {
  set.seed(seed)
  m <- matrix(runif(nInd(x) * 3), nInd(x), 3)
  m <- m / rowSums(m)
  q <- data.frame(Index = seq_len(nInd(x)), Order = seq_len(nInd(x)),
                  Label = indNames(x), PercentMiss = "0.0",
                  Cluster = as.character(apply(m, 1, which.max)),
                  Pop_1 = m[, 1], Pop_2 = m[, 2], Pop_3 = m[, 3],
                  Pop = as.character(pop(x)), stringsAsFactors = FALSE)
  q <- q[order(q$Pop, as.numeric(q$Cluster)), ]
  q$Order <- seq_len(nrow(q))
  q
}

pc_result <- function(x) {
  list(output_path = tempdir(),
       best_run = data.frame(K = "3", BestRun = "run_K_3_R_1",
                             stringsAsFactors = FALSE),
       plots = list(),
       matrix = list(run_K_3_R_1 = pc_q(x)))
}

pc_quiet <- function(expr) {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  suppressWarnings(suppressMessages(expr))
}

pc_calls <- function(m, method) {
  cl <- m$x$calls
  cl[vapply(cl, `[[`, "", "method") == method]
}

test_that("run: result layout; labels and populations match (binary)", {
  pd <- pc_dir()
  x <- pc_genlight()
  out <- withr::local_tempdir()
  r <- pc_quiet(gl.run.popcluster(x, popcluster.path = pd, output.path = out,
                                  minK = 1, maxK = 3, rep = 2,
                                  plot.out = FALSE, verbose = 0))
  expect_named(r, c("output_path", "best_run", "plots", "matrix"))
  expect_length(r$matrix, 3)
  q <- r$matrix[[3]]
  expect_equal(q$Pop, as.character(pop(x))[match(q$Label, indNames(x))])
  expect_equal(unname(rowSums(q[, c("Pop_1", "Pop_2", "Pop_3")])),
               rep(1, nInd(x)), tolerance = 0.01)
  expect_named(r$plots, c("LogL_Mean", "DLK1", "DLK2", "FST.FIS"))
})

test_that("run: numeric likelihood table, plots in the right order [approved 1, 3]", {
  pd <- pc_dir()
  x <- pc_genlight()
  out <- withr::local_tempdir()
  r <- pc_quiet(gl.run.popcluster(x, popcluster.path = pd, output.path = out,
                                  minK = 1, maxK = 3, rep = 2,
                                  plot.out = FALSE, verbose = 0))
  expect_type(r$best_run$LogL_Mean, "double")
  expect_type(r$best_run$K, "double")
  expect_true(is.na(r$best_run$DLK1[1]))
  b <- ggplot2::ggplot_build(r$plots$LogL_Mean)
  # the plotted heights are the likelihoods themselves
  expect_equal(as.numeric(b$data[[2]]$y), r$best_run$LogL_Mean)
  expect_equal(r$plots$LogL_Mean$theme$panel.background,
               theme_dartR()$panel.background)
  r2 <- pc_quiet(gl.run.popcluster(x, popcluster.path = pd,
                                   output.path = out, minK = 1, maxK = 2,
                                   plot.out = FALSE, verbose = 0,
                                   plot_theme = ggplot2::theme_bw()))
  expect_equal(r2$plots$DLK1$theme$panel.background,
               ggplot2::theme_bw()$panel.background)
})

test_that("run: nothing in the working directory; run folder removed; silent at verbose 0 [approved 2, 4]", {
  pd <- pc_dir()
  x <- pc_genlight()
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  before <- list.files(tempdir())
  out <- utils::capture.output(r <- pc_quiet(gl.run.popcluster(
    x, popcluster.path = pd, minK = 2, maxK = 2, plot.out = FALSE,
    verbose = 0)))
  expect_length(list.files(wd), 0)
  expect_length(out, 0)
  # only the input files, in output.path = tempdir()
  expect_setequal(setdiff(list.files(tempdir()), before),
                  c("output.popcluster.dat", "output.popcluster.PcPjt"))
  expect_equal(r$output_path, tempdir())
})

test_that("run: clear argument and binary errors [approved 5]", {
  x <- pc_genlight()
  d <- withr::local_tempdir()
  expect_error(gl.run.popcluster(x, popcluster.path = d, output.path = d,
                                 minK = 1, maxK = 2, verbose = 0),
               "Cannot find PopCluster")
  pd <- pc_dir()
  expect_error(gl.run.popcluster(x, popcluster.path = pd, output.path = d,
                                 minK = 3, maxK = 2, verbose = 0),
               "minK <= maxK")
  expect_error(gl.run.popcluster(x, popcluster.path = pd, output.path = d,
                                 minK = 2, maxK = 3, model = 4,
                                 verbose = 0),
               "migration model")
  expect_length(list.files(d), 0)
})

test_that("plot: bars under their labels; verbose NULL; clear errors; palette function [approved 7]", {
  x <- pc_genlight()
  r <- pc_result(x)
  q <- pc_quiet(gl.plot.popcluster(r, plot.K = 3, verbose = 0))
  expect_equal(q, r$matrix[[1]])
  p <- ggplot2::last_plot()
  b <- ggplot2::ggplot_build(p)
  labs <- unlist(lapply(b$layout$panel_params, function(z) {
    as.character(z$x$get_labels())
  }))
  brks <- unlist(lapply(b$layout$panel_params, function(z) {
    as.character(z$x$get_breaks())
  }))
  bars <- unique(p$data[, c("Order", "Label")])
  expect_equal(as.character(bars$Label), labs[match(as.character(bars$Order),
                                                     brks)])
  expect_null(formals(gl.plot.popcluster)$verbose)
  expect_error(pc_quiet(gl.plot.popcluster(r, plot.K = 5, verbose = 0)),
               "plot.K must be one of the K values in pop_cluster_result: 3")
  q2 <- pc_quiet(gl.plot.popcluster(r, plot.K = 3, color_clusters = rainbow,
                                    verbose = 0))
  expect_equal(nrow(q2), nInd(x))
  expect_error(pc_quiet(gl.plot.popcluster(r, plot.K = 3,
                                           color_clusters = "red",
                                           verbose = 0)),
               "color_clusters has 1 colours but 3 are needed")
  out <- utils::capture.output(q3 <- pc_quiet(gl.plot.popcluster(
    r, plot.K = 3, plot.out = FALSE, verbose = 0)))
  expect_length(out, 0)
})

test_that("map: bars at their own centre; movepops by name [approved 8]", {
  skip_if_not_installed("leaflet")
  x <- pc_genlight()
  q <- pc_q(x)
  x2 <- x
  pop(x2) <- factor(as.character(pop(x)),
                    levels = c("EmmacBurnBara", "EmmacBrisWive",
                               "EmmacBurdMist"))
  m <- pc_quiet(gl.map.popcluster(x2, q))
  rc <- pc_calls(m$map, "addRectangles")
  lat <- vapply(rc, function(z) z$args[[1]], 0)
  grp <- rep(names(m$Q_name), vapply(m$Q_name, nrow, 0L) * 3)
  cen <- apply(x2@other$latlon, 2, function(v) tapply(v, pop(x2), mean))
  expect_equal(as.vector(tapply(lat, grp, min)[names(m$Q_name)]),
               unname(cen[names(m$Q_name), "lat"]), tolerance = 1e-8)
  expect_named(m, c("Q_name", "map"))
  mp <- data.frame(lon = c(1, 0, 0), lat = c(0, 0, 0))
  m2 <- pc_quiet(gl.map.popcluster(x, q, movepops = mp))
  z <- pc_calls(m2$map, "addMarkers")[[1]]
  cen1 <- apply(x@other$latlon, 2, function(v) tapply(v, pop(x), mean))
  expect_equal(z$args[[2]][1], unname(cen1[1, "lon"] + 1))
})
