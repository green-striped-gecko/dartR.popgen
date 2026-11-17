# File-backed objects: gl.sfs, gl.select.panel and gl.amova read the
# genotypes in blocks of loci (dartR.base::utils.fbm.colwise,
# utils.fbm.reader, utils.fbm.pair.counts) instead of converting the object
# to a genlight. Results and messages must be identical to those of the
# gen-backed object for any block size, and no call may decode the whole FBM
# (gl.select.panel(method = "dapc") converts one population at a time, as
# adegenet needs a genlight).

fbm_pair <- function() {
  x <- gl.compliance.check(testset.gl[1:60, 1:120], verbose = 0)
  m <- as.matrix(x)
  m[, 5] <- NA                               # an all-NA locus
  m[pop(x) == levels(pop(x))[1], 7] <- NA    # all-NA in one population
  m[2, 20:40] <- NA                          # a poorly scored individual
  x@gen <- new("genlight", m, ploidy = ploidy(x))@gen
  x <- gl.recalc.metrics(x, verbose = 0)
  list(x = x, f = suppressWarnings(gl.gen2fbm(x, verbose = 0)))
}

# the multidimensional sfs needs few populations and individuals: two
# populations of four (missing loci are excluded by gl.sfs)
sfs_pair <- function(xs) {
  n <- table(pop(xs$x))
  p <- names(n)[n >= 4][1:2]
  rows <- unlist(lapply(p, function(k) which(pop(xs$x) == k)[1:4]))
  x <- gl.keep.pop(xs$x[rows, ], pop.list = p, verbose = 0)
  list(x = x, f = suppressWarnings(gl.gen2fbm(x, verbose = 0)))
}

panel <- function(method) {
  function(x) gl.select.panel(x, method = method, nl = 8, plot.out = FALSE,
                              verbose = 3)
}

calls <- list(
  sfs = function(x) gl.sfs(x, plot.out = FALSE, verbose = 3),
  sfs.single = function(x) gl.sfs(x, singlepop = TRUE, folded = FALSE,
                                  plot.out = FALSE, verbose = 3),
  panel.pahigh = panel("pahigh"),
  panel.monopop = panel("monopop"),
  panel.stratified = panel("stratified"),
  panel.hafall = panel("hafall"),
  panel.hafpop = panel("hafpop"),
  panel.pic = panel("pic"),
  amova = function(x) gl.amova(x, permutations = 9, verbose = 3)
)

strip <- function(v) {
  if (methods::is(v, "genlight")) {
    return(list(m = as.matrix(v), ind = indNames(v), loc = locNames(v),
                pop = pop(v), lm = v@other$loc.metrics))
  }
  if (is.list(v) && !is.data.frame(v)) {
    v[vapply(v, inherits, logical(1), "ggplot")] <- NULL
    v <- lapply(v, strip)
  }
  v
}

run <- function(FUN, x) {
  set.seed(1)
  out <- utils::capture.output(val <- suppressWarnings(FUN(x)))
  list(val = strip(val), out = out)
}

test_that("gen- and file-backed objects give identical results", {
  skip_if_not_installed("bigstatsr")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  old <- options(bigstatsr.block.sizeGB = getOption("bigstatsr.block.sizeGB"))
  on.exit(options(old), add = TRUE)
  set.seed(3)
  xs <- fbm_pair()
  ss <- sfs_pair(xs)
  for (bs in c(1, 1e-6)) {                   # 1e-6: one locus per block
    options(bigstatsr.block.sizeGB = bs)
    for (nm in names(calls)) {
      d <- if (grepl("^sfs", nm)) ss else xs
      g <- run(calls[[nm]], d$x)
      f <- run(calls[[nm]], d$f)
      expect_identical(f$val, g$val, label = paste(nm, bs))
      expect_identical(f$out, g$out, label = paste(nm, bs, "output"))
    }
  }
})

test_that("gl.select.panel(method = 'dapc') gives identical panels", {
  skip_if_not_installed("bigstatsr")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  # three populations of possums.gl (dapc needs variation within groups)
  x <- gl.keep.pop(possums.gl, pop.list = popNames(possums.gl)[1:3],
                   verbose = 0)
  f <- suppressWarnings(gl.gen2fbm(x, verbose = 0))
  g.res <- run(panel("dapc"), x)
  f.res <- run(panel("dapc"), f)
  expect_identical(f.res$val, g.res$val)
  expect_identical(f.res$out, g.res$out)
})

test_that("gl.amova distances equal StAMPP::stamppNeisD", {
  skip_if_not_installed("bigstatsr")
  set.seed(3)
  xs <- fbm_pair()
  y <- xs$x
  class(y) <- "genlight"
  expect_identical(.neis_d_ind_fbm(xs$f), StAMPP::stamppNeisD(y, FALSE))
})

test_that("file-backed objects are never decoded as a whole", {
  skip_if_not_installed("bigstatsr")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  set.seed(3)
  xs <- fbm_pair()
  ss <- sfs_pair(xs)
  local_mocked_bindings(.fbm_as_matrix = function(x) {
    stop("the whole FBM was decoded")
  }, .package = "dartR.base")
  for (nm in names(calls)) {
    d <- if (grepl("^sfs", nm)) ss else xs
    set.seed(1)
    expect_no_error(utils::capture.output(suppressWarnings(
      calls[[nm]](d$f))), message = nm)
  }
})

test_that("gl.check.panel and gl.check.future.panel: file-backed full data", {
  skip_if_not_installed("bigstatsr")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  old <- options(bigstatsr.block.sizeGB = getOption("bigstatsr.block.sizeGB"))
  on.exit(options(old), add = TRUE)
  set.seed(3)
  xs <- fbm_pair()
  # the panel: the first 20 loci of the full data, in both storage types
  p.g <- xs$x[, 1:20]
  p.f <- xs$f[, 1:20]
  chk <- function(p, x) {
    gl.check.panel(p, x, parameter = c("Ho_ind", "Fis_ind",
                                       "drift_resistance"),
                   plot.out = FALSE, verbose = 3)
  }
  fut <- function(p, x) {
    gl.check.future.panel(p, x, parameter = "He", n_gen = 2,
                          replicates = 2, plot.out = FALSE, verbose = 3)
  }
  for (bs in c(1, 1e-6)) {
    options(bigstatsr.block.sizeGB = bs)
    for (FUN in list(chk = chk, fut = fut)) {
      set.seed(1)
      g.out <- utils::capture.output(g <- suppressWarnings(FUN(p.g, xs$x)))
      set.seed(1)
      f.out <- utils::capture.output(f <- suppressWarnings(FUN(p.f, xs$f)))
      expect_identical(strip(f), strip(g))
      expect_identical(f.out, g.out)
    }
  }
  local_mocked_bindings(.fbm_as_matrix = function(x) {
    if (nLoc(x) > 20) stop("the full data were decoded")
    as.matrix(dartR.base::gl.fbm2gen(x, verbose = 0))
  }, .package = "dartR.base")
  set.seed(1)
  expect_no_error(utils::capture.output(suppressWarnings(chk(p.f, xs$f))))
  set.seed(1)
  expect_no_error(utils::capture.output(suppressWarnings(fut(p.f, xs$f))))
})
