# Characterization tests for gl.run.faststructure and gl.plot.faststructure --
# Phase A baseline of the dartR function review (dartR.popgen at a16da26),
# updated in Phase C: assertions tagged [approved n] were flipped with the
# matching approved change; see
# function-review/reports/dartR.popgen/gl.run.faststructure.md (dartR.base).
#
# gl.plot.faststructure is tested on hand-built run objects in the layout
# gl.run.faststructure returns. The gl.run.faststructure tests need the
# fastStructure and PLINK binaries (FASTSTRUCTURE_EXEC, PLINK_DIR) and a
# dartR.base whose gl2plink passes --make-bed (dartR.base dev after 1.2.3).

fs_genlight <- function() {
  x <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1:3],
                   verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  gl.filter.callrate(x, threshold = 0.95, verbose = 0)
}

# q_list[[K]][[replicate]] = data.frame(id, orig.pop, V1..VK)
fs_q <- function(x, k, seed) {
  set.seed(seed)
  m <- matrix(runif(nInd(x) * k), nInd(x), k)
  m <- m / rowSums(m)
  data.frame(id = indNames(x), orig.pop = pop(x), m)
}

fs_sr <- function(x, ks = 2:3, reps = 2) {
  q_list <- lapply(ks, function(k) {
    # identical replicates: one mode per K
    reps_k <- lapply(seq_len(reps), function(r) fs_q(x, k, seed = 100 * k))
    names(reps_k) <- as.character(seq_len(reps))
    reps_k
  })
  names(q_list) <- as.character(ks)
  list(q_list = q_list, plot = NULL)
}

fs_quiet <- function(expr) {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  suppressWarnings(suppressMessages(expr))
}

fs_exec <- function() {
  e <- Sys.getenv("FASTSTRUCTURE_EXEC", unset = "")
  p <- Sys.getenv("PLINK_DIR", unset = "")
  if (!nzchar(e) || !nzchar(p)) skip("FASTSTRUCTURE_EXEC / PLINK_DIR not set")
  if (!any(grepl("make-bed", deparse(dartR.base::gl2plink)))) {
    skip("dartR.base::gl2plink does not pass --make-bed")
  }
  list(exec = path.expand(e), plink = path.expand(p))
}

test_that("plot: one table per K in k.range", {
  x <- fs_genlight()
  q <- fs_quiet(gl.plot.faststructure(fs_sr(x), k.range = 2:3, verbose = 0))
  expect_equal(unname(sapply(q, function(z) unique(z$K))), c("2", "3"))
  expect_named(q[[2]], c("Label", "cluster1", "cluster2", "cluster3", "K",
                         "orig.pop", "ord"))
  # [approved 7] tables come from gl.plot.structure
  expect_s3_class(q[[1]], "data.table")
  expect_equal(nrow(q[[1]]), nInd(x))
})

test_that("plot: k.range = NULL plots every K [approved 7]", {
  q <- fs_quiet(gl.plot.faststructure(fs_sr(fs_genlight()), verbose = 0))
  expect_equal(unname(sapply(q, function(z) unique(z$K))), c("2", "3"))
})

test_that("plot: K = 1 after another K is added once [approved 7]", {
  x <- fs_genlight()
  q <- fs_quiet(gl.plot.faststructure(fs_sr(x, ks = 1:2), k.range = c(2, 1),
                                      verbose = 0))
  expect_equal(unname(sapply(q, function(z) unique(z$K))), c("2", "1"))
})

# K = 2, six replicates in two clearly different modes (3 + 3)
fs_sr_two_modes <- function(x) {
  set.seed(42)
  n <- nInd(x)
  mk <- function(split) {
    base <- cbind(c(rep(1, split), rep(0, n - split)),
                  c(rep(0, split), rep(1, n - split)))
    m <- base + matrix(runif(n * 2, 0, 0.1), n, 2)
    data.frame(id = indNames(x), orig.pop = pop(x), m / rowSums(m))
  }
  reps <- c(lapply(1:3, function(i) mk(11)), lapply(1:3, function(i) mk(21)))
  names(reps) <- as.character(1:6)
  list(q_list = list("2" = reps), plot = NULL)
}

test_that("plot: two modes return the mean of each mode [approved 7]", {
  x <- fs_genlight()
  sr <- fs_sr_two_modes(x)
  set.seed(1)
  q <- fs_quiet(gl.plot.faststructure(sr, k.range = 2, verbose = 0))
  expect_equal(unname(sapply(q, function(z) unique(z$K))), c("2.1", "2.2"))
  reps <- lapply(sr$q_list[["2"]], function(d) as.matrix(d[, 3:4]))
  means <- list(Reduce("+", reps[1:3]) / 3, Reduce("+", reps[4:6]) / 3)
  near <- function(got, L) {
    min(sapply(L, function(r) min(max(abs(got - r)), max(abs(got - r[, 2:1])))))
  }
  for (m in seq_along(q)) {
    got <- as.matrix(q[[m]][match(indNames(x), q[[m]]$Label),
                            c("cluster1", "cluster2")])
    expect_lt(near(got, means), 1e-12)
    expect_gt(near(got, reps), 1e-3)
  }
})

test_that("plot: verbose, plot.out, plot.file, label.size; input check [approved 7]", {
  expect_true(all(c("verbose", "plot.out", "plot.file", "dis.mat") %in%
                    names(formals(gl.plot.faststructure))))
  expect_false(any(grepl("aes_(", deparse(body(gl.plot.faststructure)),
                         fixed = TRUE)))
  x <- fs_genlight()
  out <- capture.output(q <- fs_quiet(gl.plot.faststructure(
    fs_sr(x), k.range = 2, plot.out = FALSE, verbose = 0)))
  expect_length(out, 0)
  d <- withr::local_tempdir()
  fs_quiet(gl.plot.faststructure(fs_sr(x), k.range = 2, plot.out = FALSE,
                                 plot.file = "fs", plot.dir = d,
                                 label.size = 8, verbose = 0))
  expect_equal(list.files(d), "fs.RDS")
  p <- readRDS(file.path(d, "fs.RDS"))
  expect_equal(p$theme$strip.text.x$size, 8)
  expect_error(gl.plot.faststructure(list(), verbose = 0),
               "returned by gl.run.faststructure")
})

test_that("run: k.range = 2:3 returns q_list named by K and a plot", {
  b <- fs_exec()
  x <- fs_genlight()
  d <- withr::local_tempdir()
  r <- fs_quiet(gl.run.faststructure(x, k.range = 2:3, num.k.rep = 2,
                                     exec = b$exec, exec.plink = b$plink,
                                     output = d, verbose = 0,
                                     plot.out = FALSE))
  expect_named(r, c("q_list", "plot"))
  expect_named(r$q_list, c("2", "3"))
  expect_named(r$q_list[["2"]], c("1", "2"))
  expect_named(r$q_list[["3"]][[1]], c("id", "orig.pop", "V1", "V2", "V3"))
  expect_equal(r$q_list[["2"]][[1]]$id, indNames(x))
  # [approved 6] the plot has its data, dartR theme, breaks at the K run
  expect_equal(r$plot$data$K, 2:3)
  expect_equal(r$plot$theme$panel.background, theme_dartR()$panel.background)
  # [approved 2] files go to a new subfolder of output
  expect_length(list.files(d, pattern = "^fastStructure_"), 1)
})

test_that("run: any k.range, K = 1 kept [approved 1]", {
  b <- fs_exec()
  x <- fs_genlight()
  d <- withr::local_tempdir()
  r <- fs_quiet(gl.run.faststructure(x, k.range = c(1, 3), exec = b$exec,
                                     exec.plink = b$plink, output = d,
                                     verbose = 0, plot.out = FALSE))
  expect_named(r$q_list, c("1", "3"))
  expect_equal(unname(vapply(r$q_list, function(z) ncol(z[[1]]) - 2, 0)),
               c(1, 3))
})

test_that("run: repeated calls in one folder never mix [approved 2]", {
  b <- fs_exec()
  x <- fs_genlight()
  d <- withr::local_tempdir()
  fs_quiet(gl.run.faststructure(x, k.range = 2:3, num.k.rep = 2,
                                exec = b$exec, exec.plink = b$plink,
                                output = d, verbose = 0, plot.out = FALSE))
  r <- fs_quiet(gl.run.faststructure(x, k.range = 2, exec = b$exec,
                                     exec.plink = b$plink, output = d,
                                     verbose = 0, plot.out = FALSE))
  expect_named(r$q_list, "2")
  expect_named(r$q_list[["2"]], "1")
  expect_length(list.files(d, pattern = "^fastStructure_"), 2)
})

test_that("run: binaries checked, verbose 0 silent [approved 4, 5]", {
  # on Windows gl.run.faststructure stops before looking for binaries
  skip_on_os("windows")
  x <- fs_genlight()
  d <- withr::local_tempdir()
  expect_error(gl.run.faststructure(x, k.range = 2,
                                    exec = file.path(d, "none"),
                                    exec.plink = d, output = d,
                                    verbose = 0),
               "fastStructure executable was not found")
  b <- fs_exec()
  expect_error(gl.run.faststructure(x, k.range = 2, exec = b$exec,
                                    exec.plink = d, output = d,
                                    verbose = 0),
               "PLINK executable was not found")
  expect_error(gl.run.faststructure(x, k.range = 0, exec = b$exec,
                                    exec.plink = b$plink, output = d,
                                    verbose = 0),
               "positive whole numbers")
  out <- capture.output(r <- suppressWarnings(gl.run.faststructure(
    x, k.range = 2, exec = b$exec, exec.plink = b$plink,
    output = file.path(d, "new_folder"), verbose = 0, plot.out = FALSE)))
  expect_length(out, 0)
  expect_named(r$q_list, "2")
})

test_that("run: a fixed seed gives different, reproducible replicates [approved 3]", {
  b <- fs_exec()
  x <- fs_genlight()
  d <- withr::local_tempdir()
  r <- fs_quiet(gl.run.faststructure(x, k.range = 2, num.k.rep = 2,
                                     seed = 5, exec = b$exec,
                                     exec.plink = b$plink, output = d,
                                     verbose = 0, plot.out = FALSE))
  expect_false(isTRUE(all.equal(r$q_list[["2"]][["1"]],
                                r$q_list[["2"]][["2"]])))
  r2 <- fs_quiet(gl.run.faststructure(x, k.range = 2, num.k.rep = 2,
                                      seed = 5, exec = b$exec,
                                      exec.plink = b$plink, output = d,
                                      verbose = 0, plot.out = FALSE))
  expect_equal(r2$q_list, r$q_list)
})
