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

# Fake newhybs (Unix shell script) that writes aa-PofZ.txt in the layout of a
# given NewHybrids build. The Unix newhybs binary writes an IndivName column
# ("NoName" when the input file carries no names); NewHybrids_PC_1_1_WOG.exe,
# which the Windows block runs, omits it. Every row carries P0 = 0.9 and
# P1 = 0.1 so column alignment can be asserted.
write_fake_newhybs <- function(dir, indiv_name_column) {
  hdr <- c("50_sweeps", if (indiv_name_column) "IndivName",
           "1.000/0.000/0.000/0.000", "0.000/0.000/0.000/1.000",
           "0.000/0.500/0.500/0.000", "0.250/0.250/0.250/0.250",
           "0.500/0.250/0.250/0.000", "0.000/0.250/0.250/0.500")
  row <- paste(c("$i", if (indiv_name_column) "NoName",
                 "0.90000", "0.10000", "0.00000", "0.00000", "0.00000", "0.00000"),
               collapse = "\t")
  script <- c(
    "#!/bin/sh",
    "while [ $# -gt 0 ]; do",
    "  if [ \"$1\" = \"--data-file\" ]; then DF=\"$2\"; shift; fi",
    "  shift",
    "done",
    "N=$(awk '/^NumIndivs/ {print $2}' \"$DF\")",
    paste0("printf '%s\\n' '", paste(hdr, collapse = "\t"), "' > aa-PofZ.txt"),
    "i=1",
    "while [ \"$i\" -le \"$N\" ]; do",
    paste0("  echo \"", row, "\" >> aa-PofZ.txt"),
    "  i=$((i + 1))",
    "done",
    "for f in aa-LociAndAlleles.txt aa-Pi.hist aa-Theta.hist aa-ThetaAverages.txt aa-EchoedGtypFreqCats.txt EchoedGtypData.txt; do : > \"$f\"; done"
  )
  path <- file.path(dir, "newhybs")
  writeLines(script, path)
  Sys.chmod(path, "755")
  invisible(path)
}

expect_pofz_aligned <- function(indiv_name_column) {
  gl <- dartR.data::testset.gl
  bin_dir <- withr::local_tempdir(.local_envir = parent.frame())
  out <- withr::local_tempdir(.local_envir = parent.frame())
  write_fake_newhybs(bin_dir, indiv_name_column)
  set.seed(42)
  res <- gl.nhybrids(gl,
    p0 = NULL, p1 = NULL,
    nhyb.directory = bin_dir,
    outpath = out,
    BurnIn = 50, sweeps = 50,
    verbose = 0
  )
  pofz <- read.csv(file.path(out, "aa-PofZ.csv"))
  expect_equal(names(pofz),
               c("id", "pop", "P0", "P1", "F1", "F2", "F1xP0", "F1xP1"))
  expect_equal(nrow(pofz), nInd(gl))
  expect_equal(pofz$id, indNames(gl))
  expect_equal(unique(pofz$P0), 0.9)
  expect_equal(unique(pofz$P1), 0.1)
}

test_that("aa-PofZ.csv is aligned when NewHybrids writes an IndivName column (newhybs)", {
  skip_on_os("windows")
  expect_pofz_aligned(indiv_name_column = TRUE)
})

test_that("aa-PofZ.csv is aligned when NewHybrids omits the IndivName column (NewHybrids_PC_1_1_WOG)", {
  skip_on_os("windows")
  expect_pofz_aligned(indiv_name_column = FALSE)
})
