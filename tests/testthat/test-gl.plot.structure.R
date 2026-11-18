# Characterization tests for gl.plot.structure -- Phase A baseline of the
# dartR function review (dartR.popgen at 95fde36), updated in Phase C:
# assertions tagged [approved n] were flipped with the matching approved
# change; see
# function-review/reports/dartR.popgen/gl.plot.structure.md (dartR.base repo).
#
# The structure.result objects are built by hand with the layout that
# gl.run.structure returns, so no STRUCTURE executable is needed.

# 31 individuals from three populations of testset.gl
ps_genlight <- function() {
  x <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1:3],
                   verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  gl.filter.callrate(x, threshold = 0.9, verbose = 0)
}

# one run in the layout of gl.run.structure: summary + q.mat
ps_run <- function(q, x, k, r) {
  q.mat <- data.frame(id = indNames(x), pct.miss = 0,
                      orig.pop = pop(x))
  for (j in seq_len(ncol(q))) q.mat[[paste0("Group.", j)]] <- q[, j]
  list(summary = c(k = k, est.ln.prob = -100 * k, mean.lnL = -90 * k,
                   var.lnL = 10),
       q.mat = q.mat, prior.anc = NULL, label = paste0("k", k, ".r", r))
}

ps_as_sr <- function(runs) {
  names(runs) <- vapply(runs, `[[`, "", "label")
  class(runs) <- c("structure.result", "list")
  runs
}

# normalised q matrix: population blocks plus a little noise
ps_q <- function(x, k, seed) {
  set.seed(seed)
  grp <- as.integer(pop(x))
  q <- t(vapply(grp, function(g) {
    v <- runif(k, 0, 0.1)
    v[(g - 1) %% k + 1] <- v[(g - 1) %% k + 1] + 1
    v
  }, numeric(k)))
  if (k == 1) q <- matrix(1, nrow = nInd(x), ncol = 1)
  q / rowSums(q)
}

# K = 1..3, three near-identical replicates each (a single mode per K)
ps_sr <- function(x) {
  runs <- list()
  for (k in 1:3) for (r in 1:3) {
    runs[[length(runs) + 1]] <- ps_run(ps_q(x, k, seed = 10 * k), x, k, r)
  }
  ps_as_sr(runs)
}

# K = 2, six replicates in two clearly different modes (3 + 3)
ps_sr_two_modes <- function(x) {
  set.seed(42)
  n <- nInd(x)
  mk <- function(split) {
    base <- cbind(c(rep(1, split), rep(0, n - split)),
                  c(rep(0, split), rep(1, n - split)))
    m <- base + matrix(runif(n * 2, 0, 0.1), n, 2)
    m / rowSums(m)
  }
  reps <- c(lapply(1:3, function(i) mk(11)), lapply(1:3, function(i) mk(21)))
  ps_as_sr(lapply(seq_along(reps), function(i) ps_run(reps[[i]], x, 2, i)))
}

quiet <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("K = 2:3 returns one data.table per K, rows sorted by Label", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- ps_genlight()
  sr <- ps_sr(x)
  set.seed(1)
  q <- quiet(gl.plot.structure(sr, K = 2:3, plot.out = FALSE, verbose = 0))
  expect_type(q, "list")
  expect_named(q, c("1", "2"))
  expect_s3_class(q[[1]], "data.table")
  expect_equal(unname(sapply(q, function(z) unique(z$K))), c("2", "3"))
  expect_named(q[[2]], c("Label", "cluster1", "cluster2", "cluster3",
                         "K", "orig.pop", "ord"))
  expect_equal(nrow(q[[1]]), nInd(x))
  expect_equal(q[[1]]$Label, sort(indNames(x)))
  # every K panel shares the individual order of the first panel
  expect_identical(q[[1]]$ord, q[[2]]$ord)
  # one mode per K: the returned q values are the mean of the replicates
  id <- sr[[4]]$q.mat$id
  got <- as.matrix(q[[1]][match(id, q[[1]]$Label), c("cluster1", "cluster2")])
  expect_equal(unname(rowSums(got)), rep(1, nInd(x)), tolerance = 1e-8)
})

test_that("K = NULL plots every K present in sr", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  sr <- ps_sr(ps_genlight())
  q <- quiet(gl.plot.structure(sr, plot.out = FALSE, verbose = 0))
  expect_equal(unname(sapply(q, function(z) unique(z$K))), c("1", "2", "3"))
})

test_that("two modes return the mean of each mode [approved 1]", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- ps_genlight()
  sr <- ps_sr_two_modes(x)
  set.seed(1)
  q <- quiet(gl.plot.structure(sr, K = 2, plot.out = FALSE, verbose = 0))
  expect_equal(unname(sapply(q, function(z) unique(z$K))), c("2.1", "2.2"))
  reps <- lapply(sr, function(z) as.matrix(z$q.mat[, 4:5]))
  means <- list(Reduce("+", reps[1:3]) / 3, Reduce("+", reps[4:6]) / 3)
  id <- sr[[1]]$q.mat$id
  near <- function(got, L) {
    min(sapply(L, function(r) min(max(abs(got - r)), max(abs(got - r[, 2:1])))))
  }
  for (m in seq_along(q)) {
    got <- as.matrix(q[[m]][match(id, q[[m]]$Label),
                            c("cluster1", "cluster2")])
    # equal to a mode mean, no longer to any single replicate
    expect_lt(near(got, means), 1e-12)
    expect_gt(near(got, reps), 1e-3)
  }
})

test_that("K = 1 after another K is added once [approved 2]", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  sr <- ps_sr(ps_genlight())
  q <- quiet(gl.plot.structure(sr, K = c(2, 1), plot.out = FALSE,
                               verbose = 0))
  expect_equal(unname(sapply(q, function(z) unique(z$K))), c("2", "1"))
  q1 <- quiet(gl.plot.structure(sr, K = 1, plot.out = FALSE, verbose = 0))
  expect_length(q1, 1)
  expect_equal(ncol(q1[[1]]), 5)
  expect_equal(q1[[1]]$cluster1, rep(1, 31))
})

test_that("den = TRUE keeps orig.pop in the returned q-matrices [approved 3]", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- ps_genlight()
  sr <- ps_sr(x)
  q <- quiet(gl.plot.structure(sr, K = 2, den = TRUE, x = x,
                               plot.out = FALSE, verbose = 0))
  q0 <- quiet(gl.plot.structure(sr, K = 2, plot.out = FALSE, verbose = 0))
  expect_equal(as.character(q[[1]]$orig.pop),
               as.character(pop(x))[match(q[[1]]$Label, indNames(x))])
  expect_equal(q, q0)
})

test_that("den = TRUE orders individuals by hclust of dis.mat itself [approved 4]", {
  x <- ps_genlight()
  sr <- ps_sr(x)
  d <- gl.dist.ind(x, method = "Manhattan", plot.display = FALSE,
                   verbose = 0)
  pdf(NULL)
  quiet(gl.plot.structure(sr, K = 2, den = TRUE, x = x, verbose = 0))
  p <- ggplot2::last_plot()
  dev.off()
  dd <- reorder(as.dendrogram(hclust(d)), TRUE, agglo.FUN = mean)
  expected <- labels(d)[order.dendrogram(dd)]
  bars <- p[[1]]$data
  bars <- unique(bars[order(as.integer(as.character(bars$ord))),
                      c("Label", "ord")])
  expect_equal(bars$Label, expected)
})

test_that("den = TRUE needs x or dis.mat; dis.mat alone works [approved 5]", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- ps_genlight()
  sr <- ps_sr(x)
  expect_error(quiet(gl.plot.structure(sr, K = 2, den = TRUE, verbose = 0)),
               "needs the genlight object")
  d <- gl.dist.ind(x, plot.display = FALSE, verbose = 0)
  q <- quiet(gl.plot.structure(sr, K = 2, den = TRUE, dis.mat = d,
                               plot.out = FALSE, verbose = 0))
  expect_length(q, 1)
  d2 <- as.matrix(d)
  dimnames(d2) <- list(paste0("z", 1:31), paste0("z", 1:31))
  expect_error(quiet(gl.plot.structure(sr, K = 2, den = TRUE,
                                       dis.mat = as.dist(d2), verbose = 0)),
               "names in dis.mat do not match")
})

test_that("input errors", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  sr <- ps_sr(ps_genlight())
  expect_error(gl.plot.structure(list(), verbose = 0),
               "not a structure result")
  expect_error(quiet(gl.plot.structure(sr, K = 7, verbose = 0)),
               "No entries for K = 7")
  # [approved 5] only the missing K is named
  expect_error(quiet(gl.plot.structure(sr, K = c(2, 7), verbose = 0)),
               "No entries for K = 7 found")
  # [approved 6] a palette function is accepted
  q <- quiet(gl.plot.structure(sr, K = 2, color_clusters = rainbow,
                               plot.out = FALSE, verbose = 0))
  expect_length(q, 1)
  # [approved 5] argument errors are raised before any work
  expect_error(quiet(gl.plot.structure(sr, K = 3, color_clusters = "red",
                                       verbose = 0)),
               "color_clusters has 1 colours but 3 are needed")
  expect_error(quiet(gl.plot.structure(sr, K = 2, met_clumpp = "foo",
                                       verbose = 0)),
               "met_clumpp must be one of")
})

test_that("verbose = 0 prints nothing; plot.out = FALSE still returns; plot.file saves an RDS", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  sr <- ps_sr(ps_genlight())
  out <- capture.output(q <- quiet(gl.plot.structure(sr, K = 2,
                                                     plot.out = FALSE,
                                                     verbose = 0)))
  expect_length(out, 0)
  expect_length(q, 1)
  d <- withr::local_tempdir()
  quiet(capture.output(gl.plot.structure(sr, K = 2, plot.out = FALSE,
                                         plot.file = "ps", plot.dir = d,
                                         verbose = 0)))
  expect_equal(list.files(d), "ps.RDS")
})

test_that("the plot no longer uses the deprecated aes_() [approved 8]", {
  src <- deparse(body(gl.plot.structure))
  expect_false(any(grepl("aes_(", src, fixed = TRUE)))
})
