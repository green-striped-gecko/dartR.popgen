# Regression tests for gl.select.panel and gl.check.panel
# (function-review campaign).
#
# Before the review (dartR.popgen 1.2.2) dapc returned every locus on
# possums.gl, pahigh and monopop stopped with errors, individuals came back
# sorted by population, gl.check.panel failed on parameter = "Nall" and
# accepted a full data set holding different individuals.
#
# possums.gl (locus names X1..X200, no SNP metadata, 10 populations) and
# bandicoot.gl (DArT locus names and metadata, 5 populations).

sel <- function(...) {
  out <- NULL
  utils::capture.output(out <- suppressWarnings(gl.select.panel(...)))
  out
}
chk <- function(...) {
  out <- NULL
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  utils::capture.output(
    out <- suppressWarnings(suppressMessages(gl.check.panel(...)))
  )
  out
}
# Private-allele loci over all population pairs, computed independently
# of gl.report.pa() (which needs networkD3/tibble/tidyr)
all_pa <- function(x) {
  q <- lapply(seppop(x), function(p) colMeans(as.matrix(p), na.rm = TRUE) / 2)
  com <- utils::combn(length(q), 2)
  unique(unlist(lapply(seq_len(ncol(com)), function(i) {
    a <- q[[com[1, i]]]
    b <- q[[com[2, i]]]
    locNames(x)[which((b == 0 & a != 0) | (b == 1 & a != 1) |
                      (a == 0 & b != 0) | (a == 1 & b != 1))]
  })))
}
report_pa <- function(x) {
  pa <- gl.report.pa(x, loc.names = TRUE, verbose = 0)$names_loci
  unique(unlist(lapply(pa, function(z) c(z$pop1_pop2_pa, z$pop2_pop1_pa))))
}
maf <- function(g) {
  a <- gl.alf(g)[, 1]
  pmin(a, 1 - a)
}

test_that("ranking methods return nl loci with the top scores", {
  x <- dartR.data::possums.gl
  r <- sel(x, method = "hafall", nl = 20, verbose = 0)
  expect_equal(nLoc(r), 20)
  # hafall ranks by minor allele frequency (closest to 0.5)
  m <- maf(x)
  expect_gte(min(maf(r)), sort(m, decreasing = TRUE)[20] - 1e-12)
  r <- sel(x, method = "pic", nl = 20, verbose = 0)
  expect_equal(nLoc(r), 20)
  expect_equal(nrow(r@other$loc.metrics), 20)
})

test_that("random, stratified and hafpop return exactly nl loci", {
  x <- dartR.data::possums.gl
  for (m in c("random", "stratified", "hafpop", "picdart")) {
    set.seed(1)
    r <- sel(x, method = m, nl = 50, verbose = 0)
    expect_equal(nLoc(r), 50, info = m)
  }
})

test_that("bandicoot: dapc and pahigh selections (seeded snapshot)", {
  b <- dartR.data::bandicoot.gl
  set.seed(1)
  d <- sel(b, method = "dapc", nl = 30, verbose = 0)
  expect_equal(nLoc(d), 30)
  set.seed(1)
  p <- sel(b, method = "pahigh", nl = 30, verbose = 0, exact = FALSE)
  expect_true(all(locNames(p) %in% all_pa(b)))
  expect_equal(nLoc(p), 24)
})

test_that("dapc, pahigh and monopop return nl loci on possums.gl", {
  x <- dartR.data::possums.gl
  for (m in c("dapc", "pahigh", "monopop")) {
    set.seed(1)
    r <- sel(x, method = m, nl = 50, verbose = 0)
    expect_equal(nLoc(r), 50, info = m)
    expect_equal(nrow(r@other$loc.metrics), 50, info = m)
  }
})

test_that("dapc takes the top-contributing loci by position", {
  x <- dartR.data::possums.gl
  x <- x[pop(x) %in% levels(pop(x))[1:2], ]
  pop(x) <- droplevels(pop(x))
  dd <- adegenet::dapc(x, n.pca = 20, n.da = 5)
  top <- locNames(x)[order(dd$var.contr[, 1], decreasing = TRUE)[1:5]]
  r <- sel(x, method = "dapc", nl = 5, verbose = 0)
  expect_setequal(locNames(r), top)
})

test_that("pahigh on possums.gl selects private-allele loci only", {
  x <- dartR.data::possums.gl
  set.seed(1)
  r <- sel(x, method = "pahigh", nl = 30, verbose = 0, exact = FALSE)
  expect_gt(nLoc(r), 0)
  expect_true(all(locNames(r) %in% all_pa(x)))
})

test_that("monopop selects loci monomorphic in some population", {
  x <- dartR.data::possums.gl
  set.seed(1)
  r <- sel(x, method = "monopop", nl = 50, verbose = 0, exact = FALSE)
  mono_any <- sapply(seppop(r), function(p) {
    cm <- colMeans(as.matrix(p), na.rm = TRUE)
    cm == 0 | cm == 2
  })
  expect_true(all(apply(mono_any, 1, any)))
})

test_that("individuals keep the input order", {
  x <- dartR.data::possums.gl
  set.seed(3)
  xs <- x[sample(nInd(x)), ]
  for (m in c("random", "hafpop", "dapc")) {
    set.seed(1)
    r <- sel(xs, method = m, nl = 10, verbose = 0)
    expect_identical(indNames(r), indNames(xs), info = m)
  }
})

test_that("invalid arguments stop with clear errors", {
  x <- dartR.data::possums.gl
  expect_error(sel(x, method = "Random", nl = 10, verbose = 0),
               "'method' must be one of")
  expect_error(sel(x, method = "random", nl = 5000, verbose = 0),
               "'nl' must be a whole number between 1 and nLoc")
  expect_error(sel(x, method = "random", nl = 2.5, verbose = 0), "'nl'")
  y <- x
  pop(y) <- NULL
  expect_error(sel(y, method = "hafpop", nl = 10, verbose = 0),
               "assigned to a population")
  expect_equal(nLoc(sel(y, method = "random", nl = 10, verbose = 0)), 10)
  expect_error(sel(dartR.data::testset.gs, method = "hafall", nl = 10,
                   verbose = 0), "SilicoDArT")
  expect_error(gl.select.panel(x, nl = 10, plot.out = FALSE),
               "unused argument")
})

test_that("verbose = 0 is silent; top-up warns at verbose 1", {
  x <- dartR.data::possums.gl
  expect_silent(gl.select.panel(x, method = "hafpop", nl = 10, verbose = 0))
  set.seed(1)
  out <- utils::capture.output(
    gl.select.panel(x, method = "monopop", nl = 50, verbose = 1)
  )
  expect_true(any(grepl("adding .* random loci to reach nl = 50", out)))
})

test_that("check.panel returns orig vs panel values per population", {
  x <- dartR.data::possums.gl
  set.seed(1)
  p <- sel(x, method = "random", nl = 100, verbose = 0)
  r <- chk(p, x, parameter = "He", verbose = 0)
  expect_equal(names(r), c("het_orig", "het_panel"))
  expect_equal(nrow(r), nPop(x))
  expect_equal(r$het_orig, gl.report.heterozygosity(x, verbose = 0)$He)
  r <- chk(p, x, parameter = "Fst", verbose = 0)
  expect_equal(nrow(r), choose(nPop(x), 2))
})

test_that("check.panel accepts 'Na' and 'Nall'; rejects unknown values", {
  x <- dartR.data::possums.gl
  set.seed(1)
  p <- sel(x, method = "random", nl = 100, verbose = 0)
  a <- chk(p, x, parameter = "Na", verbose = 0)
  b <- chk(p, x, parameter = "Nall", verbose = 0)
  expect_equal(a, b)
  expect_equal(nrow(a), nPop(x))
  expect_error(chk(p, x, parameter = "fst", verbose = 0),
               "'parameter' must be one of")
})

test_that("check.panel matches individuals by name", {
  x <- dartR.data::possums.gl
  set.seed(1)
  p <- sel(x, method = "random", nl = 100, verbose = 0)
  x2 <- x
  indNames(x2) <- paste0("other", seq_len(nInd(x2)))
  expect_error(chk(p, x2, parameter = "He", verbose = 0),
               "same individuals")
  # a panel in another individual order gives the same result
  set.seed(4)
  ps <- p[sample(nInd(p)), ]
  expect_equal(chk(ps, x, parameter = "He", verbose = 0),
               chk(p, x, parameter = "He", verbose = 0))
})

test_that("check.panel: plot.out = FALSE draws nothing; plot.file saves", {
  x <- dartR.data::possums.gl
  set.seed(1)
  p <- sel(x, method = "random", nl = 100, verbose = 0)
  dir <- withr::local_tempdir()
  n_before <- length(grDevices::dev.list())
  utils::capture.output(suppressMessages(
    gl.check.panel(p, x, parameter = "Ho", plot.out = FALSE,
                   plot.file = "panel_ho", plot.dir = dir, verbose = 0)
  ))
  expect_equal(length(grDevices::dev.list()), n_before)
  f <- list.files(dir, "panel_ho", full.names = TRUE)
  expect_length(f, 1)
  expect_s3_class(readRDS(f), "ggplot")
  expect_silent(gl.check.panel(p, x, parameter = "Ho", plot.out = FALSE,
                               verbose = 0))
})

test_that("check.panel Ne (needs NeEstimator in NEEST_DIR)", {
  neest <- Sys.getenv("NEEST_DIR", "")
  skip_if(neest == "", "NEEST_DIR not set")
  x <- dartR.data::possums.gl
  x <- x[pop(x) %in% levels(pop(x))[1:3], ]
  set.seed(1)
  p <- sel(x, method = "random", nl = 100, verbose = 0)
  r <- chk(p, x, parameter = "Ne", neest.path = neest, verbose = 0)
  expect_equal(names(r), c("nes_orig", "nes_panel"))
  expect_equal(nrow(r), 3)
  expect_true(all(is.finite(r$nes_orig) | is.na(r$nes_orig)))
})

test_that("private-allele rule agrees with gl.report.pa()", {
  skip_if_not_installed("networkD3")
  skip_if_not_installed("tibble")
  skip_if_not_installed("tidyr")
  for (x in list(dartR.data::possums.gl, dartR.data::bandicoot.gl)) {
    expect_setequal(all_pa(x), report_pa(x))
  }
})

test_that("exact = FALSE stops when the method selects no loci", {
  # the 175 possums.gl loci polymorphic in all 3 populations hold no
  # monomorphic or private-allele loci; gl.keep.loc() given no loci
  # returned the whole object
  x <- dartR.data::possums.gl
  x <- x[pop(x) %in% levels(pop(x))[1:3], ]
  pop(x) <- droplevels(pop(x))
  poly <- sapply(seppop(x), function(p) {
    cm <- colMeans(as.matrix(p), na.rm = TRUE)
    !is.na(cm) & cm > 0 & cm < 2
  })
  x <- x[, apply(poly, 1, all)]
  for (m in c("pahigh", "monopop")) {
    expect_error(sel(x, method = m, nl = 5, exact = FALSE, verbose = 0),
                 "selected no loci", info = m)
    set.seed(1)
    expect_equal(nLoc(sel(x, method = m, nl = 5, verbose = 0)), 5, info = m)
  }
})

test_that("dapc copes with loci that have no calls in a population pair", {
  x <- dartR.data::platypus.gl
  set.seed(1)
  r <- sel(x, method = "dapc", nl = 10, verbose = 0)
  expect_equal(nLoc(r), 10)
  set.seed(1)
  r <- sel(x, method = "dapc", nl = 10, exact = FALSE, verbose = 0)
  expect_gte(nLoc(r), 10)
  expect_true(all(locNames(r) %in% locNames(x)))
})
