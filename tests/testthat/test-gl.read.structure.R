# Characterization tests for gl.read.structure -- Phase A baseline of the
# dartR function review (dartR.popgen at 063c0c7), updated in Phase C:
# assertions tagged [approved n] were flipped with the matching approved
# change; see
# function-review/reports/dartR.popgen/gl.read.structure.md (dartR.base repo).
#
# Fixtures in fixtures/structure/ are real STRUCTURE 2.3.4 output files
# written by gl.run.structure(delete.files = FALSE) for three populations of
# testset.gl (31 individuals, 18 loci). In these files individuals are
# labelled by their index in indNames(x).
#   plain/      K = 1, 2, 3 (two replicates each), plus a log and a
#               mainparams file from the same folder
#   usepopinfo/ K = 2 and 3 with pop.prior = "usepopinfo"; individuals with
#               an even index had POPFLAG = 1

rs_dir <- function(sub) test_path("fixtures", "structure", sub)

# the genlight the fixtures were run on
rs_genlight <- function() {
  x <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1:3],
                   verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  gl.filter.callrate(x, threshold = 0.9, verbose = 0)
}

rs_read <- function(...) suppressWarnings(gl.read.structure(..., verbose = 0))

# copy the plain output files into a temporary folder under new names
rs_copy <- function(names_out) {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  src <- list.files(rs_dir("plain"), pattern = "_out_f$", full.names = TRUE)
  file.copy(src[seq_along(names_out)], file.path(d, names_out))
  d
}

test_that("plain output files: structure.result in the gl.run.structure layout", {
  r <- rs_read(rs_dir("plain"), pattern = "_out_f$")
  expect_s3_class(r, "structure.result")
  expect_length(r, 6)
  expect_equal(unname(sapply(r, function(z) z$summary[["k"]])),
               c(1, 1, 2, 2, 3, 3))
  expect_named(r[[3]], c("summary", "q.mat", "prior.anc", "files", "label"))
  expect_named(r[[3]]$summary, c("k", "est.ln.prob", "mean.lnL", "var.lnL"))
  expect_named(r[[3]]$q.mat, c("id", "pct.miss", "orig.pop", "Group.1",
                               "Group.2"))
  expect_equal(nrow(r[[3]]$q.mat), 31)
  expect_equal(unname(rowSums(r[[5]]$q.mat[, 4:6])), rep(1, 31),
               tolerance = 1e-8)
  # summary values as printed in the file
  txt <- readLines(file.path(rs_dir("plain"), "k2.r1_out_f"), warn = FALSE)
  lnp <- as.numeric(sub(".*= *", "",
                        grep("Estimated Ln Prob", txt, value = TRUE)))
  expect_equal(r[[3]]$summary[["est.ln.prob"]], lnp)
  # without x, orig.pop is the population number from the file
  expect_type(r[[3]]$q.mat$orig.pop, "double")
  expect_setequal(unique(r[[3]]$q.mat$orig.pop), 1:3)
})

test_that("runs are named k<K>.r<rep>; prefix only when set [approved 4]", {
  r <- rs_read(rs_dir("plain"), pattern = "_out_f$")
  expect_equal(names(r), c("k1.r1", "k1.r2", "k2.r1", "k2.r2",
                           "k3.r1", "k3.r2"))
  expect_equal(r[[3]]$label, "k2.r1")
  rp <- rs_read(rs_dir("plain"), pattern = "_out_f$", prefix = "run")
  expect_equal(names(rp)[1], "run.k1.r1")
})

test_that("replicates are numbered in file-name number order [approved 4]", {
  d <- rs_copy(c("run_K2_rep1_f", "run_K2_rep10_f", "run_K2_rep2_f"))
  # the three copies are K = 1, 1, 2; make them all K = 2 files
  src <- file.path(rs_dir("plain"), "k2.r1_out_f")
  for (f in list.files(d, full.names = TRUE)) file.copy(src, f, overwrite = TRUE)
  r <- rs_read(d)
  files <- basename(vapply(r, `[[`, "", "files"))
  expect_equal(files, c("run_K2_rep1_f", "run_K2_rep2_f", "run_K2_rep10_f"))
  expect_equal(names(r), c("k2.r1", "k2.r2", "k2.r3"))
})

test_that("a gl.run.structure folder is read without pattern [approved 3]", {
  out <- capture.output(r <- suppressWarnings(
    gl.read.structure(rs_dir("plain"), verbose = 2)))
  expect_length(r, 6)
  expect_true(any(grepl("Skipping 2 file\\(s\\).*k2.r1_log", out)))
  d <- withr::local_tempdir()
  file.copy(file.path(rs_dir("plain"), "k2.r1_log"), d)
  expect_error(rs_read(d), "No STRUCTURE output files")
})

test_that("x matches index ids by position and restores names [approved 2]", {
  x <- rs_genlight()
  r <- rs_read(rs_dir("plain"), x = x, pattern = "_out_f$")
  r0 <- rs_read(rs_dir("plain"), pattern = "_out_f$")
  idx <- as.integer(r0[[3]]$q.mat$id)
  expect_equal(r[[3]]$q.mat$id, indNames(x)[idx])
  expect_equal(r[[3]]$q.mat$orig.pop, as.character(pop(x))[idx])
  expect_equal(r[[3]]$q.mat[, 4:5], r0[[3]]$q.mat[, 4:5])
  # the result feeds gl.map.structure with population names
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  q <- suppressWarnings(gl.plot.structure(r, K = 2, plot.out = FALSE,
                                          verbose = 0))
  m <- suppressWarnings(gl.map.structure(q, x, K = 2, plot.out = FALSE,
                                         verbose = 0))
  expect_setequal(names(m$qmats), popNames(x))
})

test_that("ids matching neither way stop with an error [approved 2]", {
  x <- rs_genlight()
  x2 <- gl.keep.ind(x, ind.list = indNames(x)[1:20], verbose = 0)
  expect_error(rs_read(rs_dir("plain"), x = x2, pattern = "_out_f$"),
               "match neither indNames\\(x\\) nor positions")
})

test_that("x matches files labelled with individual names", {
  x <- rs_genlight()
  d <- withr::local_tempdir()
  txt <- readLines(file.path(rs_dir("plain"), "k2.r1_out_f"), warn = FALSE)
  first <- grep("%Miss", txt) + 1
  last <- grep("Estimated Allele", txt) - 2
  for (i in first:last) {
    idx <- as.integer(strsplit(trimws(txt[i]), " +")[[1]][2])
    txt[i] <- sub(paste0("^(\\s*[0-9]+\\s+)", idx, "(\\s)"),
                  paste0("\\1", indNames(x)[idx], "\\2"), txt[i])
  }
  writeLines(txt, file.path(d, "named_K2_f"))
  r <- rs_read(d, x = x)
  q <- r[[1]]$q.mat
  expect_equal(sort(q$id), sort(indNames(x)))
  expect_equal(as.character(q$orig.pop),
               as.character(pop(x))[match(q$id, indNames(x))])
})

test_that("usepopinfo: flagged individuals read from their ancestry blocks [approved 1]", {
  r <- rs_read(rs_dir("usepopinfo"), pattern = "_out_f$")
  q <- r[[1]]$q.mat
  # individual 10 (POPFLAG = 1, population 2); the file line reads
  #   "2  10  (5)  2 :  0.958 | Pop 1: 0.002 0.012 0.028  |"
  # so Group.1 = 0.042 and Group.2 = 0.958
  expect_equal(unlist(q[q$id == "10", c("Group.1", "Group.2")]),
               c(Group.1 = 0.042, Group.2 = 0.958))
  expect_equal(unname(r[[1]]$prior.anc[["10"]]["Pop.1", ]),
               c(0.002, 0.012, 0.028))
  expect_true(all(is.na(r[[1]]$prior.anc[["10"]]["Pop.2", ])))
  expect_equal(unname(rowSums(r[[2]]$q.mat[, 4:6])), rep(1, 31),
               tolerance = 1e-8)
  # one prior.anc matrix per q-matrix line with a '|' in the file
  txt <- readLines(file.path(rs_dir("usepopinfo"), "k2.r1_out_f"),
                   warn = FALSE)
  qtxt <- txt[(grep("%Miss", txt) + 1):(grep("Estimated Allele", txt) - 1)]
  expect_length(r[[1]]$prior.anc, sum(grepl("|", qtxt, fixed = TRUE)))
  # individuals without a prior are read correctly
  expect_equal(unlist(q[q$id == "1", c("Group.1", "Group.2")]),
               c(Group.1 = 0.09, Group.2 = 0.91))
})

test_that("input errors", {
  expect_error(rs_read("no/such/folder"), "doesn't exist")
  expect_error(rs_read(rs_dir("plain"), pattern = "nothing_matches"),
               "No files found")
})

test_that("rename_files refuses to overwrite [approved 5]", {
  d <- withr::local_tempdir()
  file.copy(file.path(rs_dir("plain"), c("k1.r1_out_f", "k2.r1_out_f")), d)
  writeLines("keep me", file.path(d, "k2.r1_out"))
  expect_error(rs_read(d, pattern = "_f$", rename_files = TRUE),
               "No files were renamed")
  expect_true(file.exists(file.path(d, "k1.r1_out_f")))
  expect_equal(readLines(file.path(d, "k2.r1_out")), "keep me")
  unlink(file.path(d, "k2.r1_out"))
  r <- rs_read(d, pattern = "_f$", rename_files = TRUE)
  expect_setequal(list.files(d), c("k1.r1_out", "k2.r1_out"))
})

test_that("verbose 0 is silent; messages single-spaced [approved 6]", {
  out <- capture.output(r <- rs_read(rs_dir("plain"), pattern = "_out_f$"))
  expect_length(out, 0)
  src <- deparse(body(gl.read.structure))
  expect_false(any(grepl("build = ", src, fixed = TRUE)))
  out2 <- capture.output(r <- suppressWarnings(
    gl.read.structure(rs_dir("plain"), pattern = "_out_f$", verbose = 2)))
  expect_true(any(grepl("Processing 6 STRUCTURE output files", out2)))
  expect_error(rs_read(rs_dir("plain"), x = "not a genlight"))
})
