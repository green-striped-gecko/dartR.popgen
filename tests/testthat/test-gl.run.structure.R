# Characterization tests for gl.run.structure -- Phase A baseline of the
# dartR function review (dartR.popgen at 31fb40c), updated in Phase C:
# assertions tagged [approved Fn] were flipped with the matching approved
# finding. See function-review/reports/dartR.popgen/gl.run.structure.md
# (dartR.base repo).
#
# STRUCTURE (non-GUI, v2.3.4) is required. The tests run only when the
# environment variable STRUCTURE_EXEC points to the executable, e.g.
#   Sys.setenv(STRUCTURE_EXEC = "~/programs/structure")

structure_exec <- function() {
  p <- Sys.getenv("STRUCTURE_EXEC", unset = "")
  if (!nzchar(p)) skip("STRUCTURE_EXEC not set")
  p <- path.expand(p)
  if (!file.exists(p)) skip(paste("STRUCTURE executable not found at", p))
  p
}

# three populations, polymorphic loci with call rate >= 0.9: 31 individuals,
# 18 loci; small enough for STRUCTURE to run in about a second per K
st_fixture <- function() {
  x <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1:3],
                   verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  gl.filter.callrate(x, threshold = 0.9, verbose = 0)
}

# Run inside a scratch working directory so anything written there is
# detected; STRUCTURE's own output is suppressed by the function below
# verbose 3.
run_in <- function(expr, wd) {
  withr::local_dir(wd)
  suppressWarnings(suppressMessages(expr))
}

test_that("k.range = 1:3 returns a structure.result named by K and replicate, rows in genlight order", {
  exec <- structure_exec()
  wd <- withr::local_tempdir()
  withr::local_options(dartR_wd = NULL)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- st_fixture()
  expect_equal(nInd(x), 31)
  expect_equal(nLoc(x), 18)

  sr <- run_in(gl.run.structure(x, exec = exec, k.range = 1:3,
                                burnin = 100, numreps = 100,
                                randomize = FALSE, seed = 7,
                                plot.out = FALSE, verbose = 0), wd)
  expect_s3_class(sr, "structure.result")
  expect_length(sr, 3)
  # [approved F2] run names no longer carry the timestamp label
  expect_equal(names(sr), c("k1.r1", "k2.r1", "k3.r1"))
  expect_equal(unname(sapply(sr, function(z) z$summary["k"])), c(1, 2, 3))
  expect_equal(unname(sapply(sr, function(z) z$label)), names(sr))
  expect_named(sr[[2]]$summary, c("k", "est.ln.prob", "mean.lnL", "var.lnL"))
  expect_named(sr[[2]]$q.mat, c("id", "pct.miss", "orig.pop", "Group.1", "Group.2"))
  expect_equal(nrow(sr[[2]]$q.mat), 31)
  # [approved F3/F4] rows follow indNames(x), with the full names
  expect_equal(sr[[2]]$q.mat$id, indNames(x))
  expect_equal(as.character(sr[[2]]$q.mat$orig.pop), as.character(pop(x)))
  expect_equal(rowSums(sr[[2]]$q.mat[, 4:5]), rep(1, 31), tolerance = 1e-6)
  expect_null(sr[[2]]$prior.anc)
  expect_null(sr[[2]]$files)
  # nothing in the working directory, no run directory left in tempdir()
  expect_length(list.files(wd, all.files = TRUE, no.. = TRUE), 0)
  expect_length(list.files(tempdir(), pattern = "^structureRun_"), 0)
})

test_that("randomize = FALSE with a seed is reproducible", {
  exec <- structure_exec()
  wd <- withr::local_tempdir()
  withr::local_options(dartR_wd = NULL)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- st_fixture()
  a <- run_in(gl.run.structure(x, exec = exec, k.range = 1:3,
                               burnin = 100, numreps = 100,
                               randomize = FALSE, seed = 7,
                               plot.out = FALSE, verbose = 0), wd)
  b <- run_in(gl.run.structure(x, exec = exec, k.range = 1:3,
                               burnin = 100, numreps = 100,
                               randomize = FALSE, seed = 7,
                               plot.out = FALSE, verbose = 0), wd)
  expect_identical(lapply(a, `[[`, "summary"), lapply(b, `[[`, "summary"))
  expect_identical(lapply(a, `[[`, "q.mat"), lapply(b, `[[`, "q.mat"))
})

test_that("fewer than three K values still returns the runs", {
  exec <- structure_exec()
  wd <- withr::local_tempdir()
  withr::local_options(dartR_wd = NULL)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- st_fixture()
  # [approved F1] the Evanno step used to error here and lose the runs
  out <- capture.output(
    sr1 <- run_in(gl.run.structure(x, exec = exec, k.range = 2, burnin = 100,
                                   numreps = 100, verbose = 1), wd)
  )
  expect_s3_class(sr1, "structure.result")
  expect_length(sr1, 1)
  expect_true(any(grepl("at least three values of K", out)))
  sr2 <- run_in(gl.run.structure(x, exec = exec, k.range = 2:3, burnin = 100,
                                 numreps = 100, plot.out = FALSE, verbose = 0),
                wd)
  expect_equal(names(sr2), c("k2.r1", "k3.r1"))
})

test_that("num.k.rep = 2 with k.range = 1:3 returns six runs and draws the Evanno plot", {
  exec <- structure_exec()
  wd <- withr::local_tempdir()
  withr::local_options(dartR_wd = NULL)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- st_fixture()
  sr <- run_in(gl.run.structure(x, exec = exec, k.range = 1:3, num.k.rep = 2,
                                burnin = 50, numreps = 50, verbose = 0), wd)
  expect_equal(names(sr), c("k1.r1", "k1.r2", "k2.r1", "k2.r2", "k3.r1", "k3.r2"))
})

test_that("Evanno warnings reach the user at verbose >= 1 (follow-up to gl.evanno review)", {
  exec <- structure_exec()
  wd <- withr::local_tempdir()
  withr::local_options(dartR_wd = NULL)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- st_fixture()
  # a fixed seed without randomize gives identical replicates, so sd = 0;
  # the Evanno step runs only when the plot is requested
  run <- function(v) {
    capture.output(sr <- run_in(gl.run.structure(
      x, exec = exec, k.range = 1:3, num.k.rep = 2, burnin = 50,
      numreps = 50, randomize = FALSE, seed = 7, plot.out = TRUE,
      verbose = v), wd))
  }
  expect_true(any(grepl("all replicates agree on LnP\\(K\\) \\(sd = 0\\)",
                        run(1))))
  expect_length(run(0), 0)
})

test_that("delete.files = FALSE keeps the run files under plot.dir, not the working directory", {
  exec <- structure_exec()
  wd <- withr::local_tempdir()
  pd <- withr::local_tempdir()
  withr::local_options(dartR_wd = NULL)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- st_fixture()
  # [approved F2]
  sr <- run_in(gl.run.structure(x, exec = exec, k.range = 1:3,
                                burnin = 100, numreps = 100,
                                plot.out = FALSE, delete.files = FALSE,
                                plot.dir = pd, verbose = 0), wd)
  expect_length(list.files(wd, all.files = TRUE, no.. = TRUE), 0)
  kept <- list.files(pd, pattern = "^structureRun_")
  expect_length(kept, 1)
  expect_named(sr[[1]]$files, c("data", "mainparams", "extraparams", "out"))
  expect_true(all(file.exists(sr[[1]]$files)))
  expect_equal(dirname(sr[[1]]$files[["data"]]), file.path(pd, kept))
  expect_equal(basename(sr[[2]]$files[["out"]]), "k2.r1_out_f")
})

test_that("plot.file follows gl.set.wd() and the Evanno plot is saved", {
  exec <- structure_exec()
  wd <- withr::local_tempdir()
  user_wd <- withr::local_tempdir()
  withr::local_options(dartR_wd = user_wd)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- st_fixture()
  # [approved F6]
  run_in(gl.run.structure(x, exec = exec, k.range = 1:3, burnin = 100,
                          numreps = 100, plot.out = FALSE,
                          plot.file = "evanno_post", verbose = 0), wd)
  expect_true(file.exists(file.path(user_wd, "evanno_post.RDS")))
  pa <- readRDS(file.path(user_wd, "evanno_post.RDS"))
  expect_s3_class(pa, "patchwork")
})

test_that("long individual names and names with spaces are restored in q.mat", {
  exec <- structure_exec()
  wd <- withr::local_tempdir()
  withr::local_options(dartR_wd = NULL)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- st_fixture()
  # [approved F3, F4]
  indNames(x) <- paste0("individual number ", seq_len(nInd(x)))
  sr <- run_in(gl.run.structure(x, exec = exec, k.range = 2, burnin = 100,
                                numreps = 100, plot.out = FALSE, verbose = 0),
               wd)
  expect_equal(sr[[1]]$q.mat$id, indNames(x))
  expect_equal(as.character(sr[[1]]$q.mat$orig.pop), as.character(pop(x)))
})

test_that("an unnamed popflag is matched to individuals in indNames order", {
  exec <- structure_exec()
  wd <- withr::local_tempdir()
  pd <- withr::local_tempdir()
  withr::local_options(dartR_wd = NULL)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- st_fixture()
  pf <- rep(1, nInd(x))
  pf[c(1, 5)] <- 0
  sr <- run_in(gl.run.structure(x, exec = exec, k.range = 3, burnin = 100,
                                numreps = 100, pop.prior = "usepopinfo",
                                popflag = pf, plot.out = FALSE,
                                delete.files = FALSE, plot.dir = pd,
                                verbose = 0), wd)
  # the data file carries id index, population, popflag, then alleles
  rows <- strsplit(readLines(sr[[1]]$files[["data"]])[-1], " +")
  idx <- as.integer(vapply(rows, `[`, "", 1))
  flag <- as.integer(vapply(rows, `[`, "", 3))
  expect_equal(sort(unique(idx)), seq_len(nInd(x)))
  expect_true(all(flag == pf[idx]))
  expect_equal(sort(unique(idx[flag == 0])), c(1L, 5L))
  # results carry the full names; prior ancestry only for flagged individuals
  expect_equal(sr[[1]]$q.mat$id, indNames(x))
  expect_named(sr[[1]]$q.mat, c("id", "pct.miss", "orig.pop", "Group.1", "Group.2", "Group.3"))
  expect_equal(names(sr[[1]]$prior.anc), indNames(x)[pf == 1])
  expect_equal(dim(sr[[1]]$prior.anc[[1]]), c(3L, 3L))
})

test_that("an executable path containing a space works", {
  exec <- structure_exec()
  wd <- withr::local_tempdir()
  withr::local_options(dartR_wd = NULL)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- st_fixture()
  d <- file.path(withr::local_tempdir(), "my tools")
  dir.create(d)
  file.copy(exec, file.path(d, basename(exec)))
  Sys.chmod(file.path(d, basename(exec)), "755")
  # [approved F5]
  sr <- run_in(gl.run.structure(x, exec = file.path(d, basename(exec)),
                                k.range = 2, burnin = 100, numreps = 100,
                                plot.out = FALSE, verbose = 0), wd)
  expect_length(sr, 1)
})

test_that("a STRUCTURE failure is reported with its exit status and leaves nothing behind", {
  wd <- withr::local_tempdir()
  withr::local_options(dartR_wd = NULL)
  x <- st_fixture()
  fake <- file.path(withr::local_tempdir(), "structure")
  writeLines(c("#!/bin/sh", "echo 'fake structure: cannot allocate'", "exit 2"), fake)
  Sys.chmod(fake, "755")
  skip_if(.Platform$OS.type == "windows", "shell script executable")
  expect_error(
    run_in(gl.run.structure(x, exec = fake, k.range = 2, plot.out = FALSE,
                            verbose = 0), wd),
    "exited with status 2"
  )
  expect_length(list.files(wd, all.files = TRUE, no.. = TRUE), 0)
  expect_length(list.files(tempdir(), pattern = "^structureRun_"), 0)
})

test_that("verbose = 0 is silent at the R level", {
  exec <- structure_exec()
  wd <- withr::local_tempdir()
  withr::local_options(dartR_wd = NULL)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  x <- st_fixture()
  # [approved F8] the gl.filter.allna messages used to print here
  out <- capture.output(
    sr <- run_in(gl.run.structure(x, exec = exec, k.range = 2, burnin = 50,
                                  numreps = 50, plot.out = FALSE, verbose = 0),
                 wd)
  )
  expect_length(out, 0)
  expect_length(sr, 1)
})

test_that("input checks stop with a clear message", {
  x <- st_fixture()
  expect_error(
    gl.run.structure(x, exec = "/nonexistent/structure", k.range = 1:3,
                     verbose = 0),
    "Cannot find the STRUCTURE executable"
  )
  exec <- structure_exec()
  expect_error(
    suppressMessages(gl.run.structure(testset.gs, exec = exec, k.range = 1:3,
                                      verbose = 0)),
    "SilicoDArT"
  )
  expect_error(
    gl.run.structure(x, exec = exec, k.range = 0:2, verbose = 0),
    "positive whole numbers"
  )
  expect_error(
    gl.run.structure(x, exec = exec, k.range = 2, popflag = c(1, 0),
                     verbose = 0),
    "one value per individual"
  )
})
