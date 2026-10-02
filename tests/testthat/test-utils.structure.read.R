# utils.structure.read is the single STRUCTURE output reader used by
# gl.read.structure and utils.structure.run (gl.run.structure). Fixtures:
# see test-gl.read.structure.R.

test_that("one reader: gl.read.structure returns what utils.structure.read parses", {
  f <- test_path("fixtures", "structure", "usepopinfo", "k2.r1_out_f")
  a <- utils.structure.read(f)
  b <- suppressWarnings(gl.read.structure(dirname(f), pattern = "k2.r1_out_f",
                                          verbose = 0))[[1]]
  expect_identical(a$summary, b$summary)
  expect_identical(a$q.mat, b$q.mat)
  expect_identical(a$prior.anc, b$prior.anc)
})

test_that("pops labels orig.pop by population number, as gl.run.structure uses it", {
  f <- test_path("fixtures", "structure", "plain", "k2.r1_out_f")
  a <- utils.structure.read(f)
  p <- utils.structure.read(f, pops = c("popA", "popB", "popC"))
  expect_equal(p$q.mat$orig.pop, c("popA", "popB", "popC")[a$q.mat$orig.pop])
  expect_identical(p$q.mat[, -3], a$q.mat[, -3])
})

test_that("utils.structure.run reads its output through utils.structure.read", {
  src <- deparse(body(utils.structure.run))
  expect_true(any(grepl("utils.structure.read(", src, fixed = TRUE)))
  expect_false(any(grepl("structureRead <- function", src, fixed = TRUE)))
})
