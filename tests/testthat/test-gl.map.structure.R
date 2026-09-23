# Characterization tests for gl.map.structure -- Phase A baseline of the
# dartR function review (dartR.popgen at 95fde36), updated in Phase C:
# assertions tagged [approved n] were flipped with the matching approved
# change; see
# function-review/reports/dartR.popgen/gl.map.structure.md (dartR.base repo).
#
# The q-matrices are built by hand in the layout gl.plot.structure returns
# (Label, cluster1..K, K, orig.pop, ord), so no STRUCTURE executable is
# needed. The map is inspected through the leaflet widget's call list.

# 31 individuals from three populations of testset.gl, with lat/lon
ms_genlight <- function() {
  x <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1:3],
                   verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  gl.filter.callrate(x, threshold = 0.9, verbose = 0)
}

# one q-matrix table per entry of k (a vector of K labels)
ms_qmat <- function(x, k = "2", seed = 1) {
  set.seed(seed)
  out <- lapply(k, function(kl) {
    nk <- as.integer(sub("\\..*", "", kl))
    m <- matrix(runif(nInd(x) * nk), nInd(x), nk)
    m <- m / rowSums(m)
    d <- data.frame(Label = indNames(x))
    for (j in seq_len(nk)) d[[paste0("cluster", j)]] <- m[, j]
    d$K <- kl
    d$orig.pop <- pop(x)
    d$ord <- seq_len(nInd(x))
    d
  })
  names(out) <- as.character(seq_along(out))
  out
}

ms_calls <- function(r, method) {
  cl <- r$map$x$calls
  cl[vapply(cl, `[[`, "", "method") == method]
}

# mean longitude and minimum latitude of the bars drawn for each element
# of r$qmats (rectangles are added population by population)
ms_bar_positions <- function(r, K) {
  rc <- ms_calls(r, "addRectangles")
  lat <- vapply(rc, function(z) z$args[[1]], 0)
  lng <- vapply(rc, function(z) z$args[[2]], 0)
  grp <- factor(rep(names(r$qmats), vapply(r$qmats, nrow, 0L) * K),
                levels = names(r$qmats))
  data.frame(lng = as.vector(tapply(lng, grp, mean)),
             lat = as.vector(tapply(lat, grp, min)),
             row.names = names(r$qmats))
}

ms_labels <- function(r) {
  z <- ms_calls(r, "addMarkers")[[1]]
  data.frame(label = z$args[[11]], lat = z$args[[1]], lng = z$args[[2]])
}

ms_centres <- function(x) {
  apply(x@other$latlon, 2, function(v) tapply(v, pop(x), mean))
}

ms_run <- function(expr) {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  suppressWarnings(suppressMessages(expr))
}

test_that("default call: one table per population, bars at their own centre", {
  skip_if_not_installed("leaflet")
  x <- ms_genlight()
  q <- ms_qmat(x, c("2", "3"))
  r <- ms_run(gl.map.structure(q, x, K = 2, verbose = 0))
  expect_named(r, c("qmats", "map"))
  expect_s3_class(r$map, "leaflet")
  expect_named(r$qmats, levels(pop(x)))
  expect_equal(unname(vapply(r$qmats, nrow, 0L)), c(10L, 10L, 11L))
  # K = 2 table chosen, rows sorted by cluster1 within population
  expect_true(all(vapply(r$qmats, function(d) unique(d$K), "") == "2"))
  expect_false(is.unsorted(r$qmats[[1]]$cluster1))
  cen <- ms_centres(x)
  pos <- ms_bar_positions(r, 2)
  expect_equal(pos$lat, unname(cen[, "lat"]), tolerance = 1e-8)
  expect_lt(max(abs(pos$lng - cen[, "lon"])), 0.1)
  lb <- ms_labels(r)
  expect_equal(lb$label, rownames(cen))
  expect_equal(lb$lng, unname(cen[, "lon"]))
  # one rectangle per individual and cluster
  expect_length(ms_calls(r, "addRectangles"), nInd(x) * 2)
  # [approved 6] the gl.plot.structure default palette, not rainbow()
  cols <- unique(unlist(lapply(ms_calls(r, "addRectangles"), function(z) {
    grep("^#", unlist(z$args), value = TRUE)
  })))
  expect_setequal(cols, gl.select.colors(ncolors = 2, verbose = 0))
})

test_that("non-alphabetical pop levels: bars at their own centre [approved 1]", {
  skip_if_not_installed("leaflet")
  x <- ms_genlight()
  q <- ms_qmat(x)
  pop(x) <- factor(as.character(pop(x)),
                   levels = c("EmmacBurnBara", "EmmacBrisWive",
                              "EmmacBurdMist"))
  r <- ms_run(gl.map.structure(q, x, K = 2, verbose = 0))
  cen <- ms_centres(x)
  pos <- ms_bar_positions(r, 2)
  expect_equal(pos$lat, unname(cen[rownames(pos), "lat"]), tolerance = 1e-8)
  lb <- ms_labels(r)
  expect_equal(lb$label, rownames(pos))
  expect_equal(lb$lng, unname(cen[rownames(pos), "lon"]))
})

test_that("movepops is added by column name; named rows match populations [approved 3]", {
  skip_if_not_installed("leaflet")
  x <- ms_genlight()
  q <- ms_qmat(x)
  expect_equal(colnames(x@other$latlon), c("lat", "lon"))
  cen <- ms_centres(x)
  mp <- data.frame(lon = c(1, 0, 0), lat = c(0, 0, 0))
  lb <- ms_labels(ms_run(gl.map.structure(q, x, K = 2, movepops = mp,
                                          verbose = 0)))
  expect_equal(lb$lng[1], unname(cen[1, "lon"] + 1))
  expect_equal(lb$lat[1] - lb$lat[2],
               unname(cen[1, "lat"] - cen[2, "lat"]), tolerance = 1e-8)
  # one named row moves only that population
  mp2 <- data.frame(lon = 0, lat = 2, row.names = "EmmacBurnBara")
  lb2 <- ms_labels(ms_run(gl.map.structure(q, x, K = 2, movepops = mp2,
                                           verbose = 0)))
  expect_equal(lb2$lat[3] - lb2$lat[2],
               unname(cen[3, "lat"] + 2 - cen[2, "lat"]), tolerance = 1e-8)
  expect_error(ms_run(gl.map.structure(q, x, K = 2,
                                       movepops = mp[1:2, ], verbose = 0)),
               "one row per population")
})

test_that("extra populations in x are ignored; missing ones error [approved 2]", {
  skip_if_not_installed("leaflet")
  x <- ms_genlight()
  q <- ms_qmat(x)
  x4 <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1:4],
                    verbose = 0)
  out <- capture.output(r <- ms_run(gl.map.structure(q, x4, K = 2,
                                                     verbose = 2)))
  expect_named(r$qmats, levels(pop(x)))
  expect_true(any(grepl("not mapped: EmmacClarJack", out)))
  x2 <- gl.keep.pop(x, pop.list = popNames(x)[1:2], verbose = 0)
  expect_error(ms_run(gl.map.structure(q, x2, K = 2, verbose = 0)),
               "No coordinates in x for these populations of qmat: EmmacBurnBara")
})

test_that("one population and K = 1 are mapped [approved 5]", {
  skip_if_not_installed("leaflet")
  x <- ms_genlight()
  xs <- gl.keep.pop(x, pop.list = popNames(x)[1], verbose = 0)
  r <- ms_run(gl.map.structure(ms_qmat(xs), xs, K = 2, verbose = 0))
  rc <- ms_calls(r, "addRectangles")
  expect_length(rc, nInd(xs) * 2)
  expect_equal(rc[[1]]$args[[4]] - rc[[1]]$args[[2]], 0.01)
  r1 <- ms_run(gl.map.structure(ms_qmat(x, "1"), x, K = 1, verbose = 0))
  expect_length(ms_calls(r1, "addRectangles"), nInd(x))
})

test_that("inputs are checked up front [approved 4]", {
  skip_if_not_installed("leaflet")
  x <- ms_genlight()
  q <- ms_qmat(x)
  xn <- x
  xn@other$latlon[] <- NA
  expect_error(ms_run(gl.map.structure(q, xn, K = 2, verbose = 0)),
               "No coordinates in x")
  xl <- x
  colnames(xl@other$latlon) <- c("y", "x")
  expect_error(ms_run(gl.map.structure(q, xl, K = 2, verbose = 0)),
               "columns named lon and")
  expect_error(ms_run(gl.map.structure(list(1), x, K = 2, verbose = 0)),
               "list of q-matrices")
  expect_error(ms_run(gl.map.structure(q, x, K = 2:3, verbose = 0)),
               "single value")
  expect_error(ms_run(gl.map.structure(q, testset.gl@other$latlon, K = 2,
                                       verbose = 0)))
  expect_error(ms_run(gl.map.structure(q, x, K = 3, plot.colors = "red",
                                       verbose = 0)),
               "No entries for K = 3")
  expect_error(ms_run(gl.map.structure(q, x, K = 2, plot.colors = "red",
                                       verbose = 0)),
               "plot.colors has 1 colours but 2 are needed")
  r <- ms_run(gl.map.structure(q, x, K = 2, plot.colors = rainbow,
                               verbose = 0))
  expect_s3_class(r$map, "leaflet")
})

test_that("absent K errors; K picks a mode by label [approved 7]", {
  skip_if_not_installed("leaflet")
  x <- ms_genlight()
  expect_error(ms_run(gl.map.structure(ms_qmat(x), x, K = 5, verbose = 0)),
               "No entries for K = 5")
  q <- ms_qmat(x, c("2.1", "2.2"))
  out <- capture.output(r <- ms_run(gl.map.structure(q, x, K = 2,
                                                     verbose = 2)))
  expect_true(all(vapply(r$qmats, function(d) unique(d$K), "") == "2.1"))
  expect_true(any(grepl("matches 2 modes", out)))
  r2 <- ms_run(gl.map.structure(q, x, K = "2.2", verbose = 0))
  expect_true(all(vapply(r2$qmats, function(d) unique(d$K), "") == "2.2"))
})

test_that("verbose = 0 is silent and plot.out = FALSE still returns [approved 8]", {
  skip_if_not_installed("leaflet")
  x <- ms_genlight()
  q <- ms_qmat(x)
  out <- capture.output(r <- ms_run(gl.map.structure(q, x, K = 2,
                                                     plot.out = FALSE,
                                                     verbose = 0)))
  expect_length(out, 0)
  expect_named(r, c("qmats", "map"))
  out3 <- capture.output(ms_run(gl.map.structure(q, x, K = 2,
                                                 plot.out = FALSE,
                                                 verbose = 3)))
  expect_true(any(grepl("Completed: gl.map.structure", out3)))
})

test_that("individuals without a population: dropped with a warning, all NA is an error (follow-up)", {
  skip_if_not_installed("leaflet")
  x <- ms_genlight()
  q <- ms_qmat(x)
  q[[1]]$orig.pop <- as.character(q[[1]]$orig.pop)
  q[[1]]$orig.pop[1:3] <- NA
  out <- capture.output(r <- ms_run(gl.map.structure(q, x, K = 2,
                                                     plot.out = FALSE,
                                                     verbose = 1)))
  expect_true(any(grepl("3 individual\\(s\\) in qmat have no population",
                        out)))
  expect_equal(sum(vapply(r$qmats, nrow, 0L)), nInd(x) - 3)
  expect_false(anyNA(unlist(lapply(r$qmats, `[[`, "Label"))))
  expect_length(ms_calls(r, "addRectangles"), (nInd(x) - 3) * 2)
  q[[1]]$orig.pop <- NA
  expect_error(ms_run(gl.map.structure(q, x, K = 2, verbose = 0)),
               "No individual in qmat has a population")
})
