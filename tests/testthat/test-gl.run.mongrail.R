# Tests for gl.run.mongrail(). The end-to-end run needs the mongrail2 executable
# and is skipped when it is not available; the validation paths are always tested.

test_that("a missing mongrail.path is a fatal error", {
  old <- getOption("mongrail.path")
  options(mongrail.path = NULL)
  on.exit(options(mongrail.path = old), add = TRUE)
  expect_error(
    gl.run.mongrail(testset.gl, pop.A = popNames(testset.gl)[1],
                    pop.B = popNames(testset.gl)[2], verbose = 0),
    "mongrail.path"
  )
})

test_that("an unknown reference population is a fatal error", {
  expect_error(
    gl.run.mongrail(testset.gl, pop.A = "not_a_pop", pop.B = "also_not",
                    mongrail.path = tempdir(), verbose = 0)
  )
})

test_that("a missing executable is a fatal error", {
  expect_error(
    gl.run.mongrail(testset.gl, pop.A = popNames(testset.gl)[1],
                    pop.B = popNames(testset.gl)[2],
                    mongrail.path = tempdir(), verbose = 0),
    "executable not found"
  )
})

test_that("end-to-end run classifies reference individuals (needs mongrail2)", {
  mp <- getOption("mongrail.path")
  exe <- if (!is.null(mp)) file.path(mp,
    if (Sys.info()[["sysname"]] == "Windows") "mongrail2.exe" else "mongrail2")
  skip_if(is.null(mp) || !file.exists(exe), "mongrail2 executable not available")

  pops <- popNames(testset.gl)[1:2]
  res <- gl.run.mongrail(testset.gl, pop.A = pops[1], pop.B = pops[2],
                         test.pop = pops, mongrail.path = mp, verbose = 0)
  expect_s3_class(res$assignments, "data.frame")
  expect_true(all(c("id", "class", "class.code") %in%
                  colnames(res$assignments)))
  probs <- res$assignments[, grep("^P\\.", colnames(res$assignments))]
  expect_equal(ncol(probs), 6)
  # posteriors sum to ~1 per individual
  expect_true(all(abs(rowSums(probs) - 1) < 0.01))
})
