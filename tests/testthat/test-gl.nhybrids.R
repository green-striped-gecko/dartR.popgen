# Characterization baseline for gl.nhybrids (function-review campaign).
# Snapshots what the function does today; it does not assert correctness.

test_that("input-file-only path (nhyb.directory = NULL) on Unix", {
  skip_on_os("windows")
  gl <- dartR.data::testset.gl
  out <- withr::local_tempdir()
  set.seed(42)
  # Post-fix behaviour (approved change 2, review report F2): with
  # nhyb.directory = NULL the function writes the input file and returns
  # the reduced genlight as documented.
  res <- gl.nhybrids(gl,
    p0 = NULL, p1 = NULL,
    nhyb.directory = NULL,
    outpath = out,
    verbose = 0
  )
  expect_s4_class(res, "genlight")
  expect_equal(nLoc(res), 200)
  expect_equal(nInd(res), nInd(gl))
  # loc.metrics track the sampled loci positionally (approved change 7 / F7)
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
  ids <- sub("\\|.*", "", as.character(res@other$loc.metrics$AlleleID))
  expect_equal(ids, sub("-.*", "", locNames(res)))
  hdr <- readLines(file.path(out, "nhyb.txt"), n = 3)
  expect_match(hdr[1], "^NumIndivs\\s+250")
  expect_match(hdr[2], "^NumLoci\\s+200")
})

test_that("execution path builds aa-PofZ.csv (requires newhybs binary)", {
  bin_dir <- Sys.getenv("NEWHYBS_DIR", "")
  skip_if(bin_dir == "" || !file.exists(file.path(bin_dir, "newhybs")),
          "newhybs binary not available (set NEWHYBS_DIR)")
  gl <- dartR.data::testset.gl
  out <- withr::local_tempdir()
  set.seed(42)
  res <- gl.nhybrids(gl,
    p0 = NULL, p1 = NULL,
    nhyb.directory = bin_dir,
    outpath = out,
    BurnIn = 50, sweeps = 50,
    verbose = 0
  )
  pofz <- read.csv(file.path(out, "aa-PofZ.csv"))
  expect_equal(nrow(pofz), nInd(gl))
  # Post-fix contract (approved change 1, review report F1): id, pop, then
  # the six NewHybrids category probabilities, correctly labelled.
  expect_equal(names(pofz),
               c("id", "pop", "P0", "P1", "F1", "F2", "F1xP0", "F1xP1"))
  probs <- pofz[, c("P0", "P1", "F1", "F2", "F1xP0", "F1xP1")]
  expect_true(all(vapply(probs, is.numeric, logical(1))))
  expect_true(all(abs(rowSums(probs) - 1) < 0.01))
})
