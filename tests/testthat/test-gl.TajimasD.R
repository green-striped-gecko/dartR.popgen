# Characterization tests for gl.TajimasD and utils.get.allele.freq -- Phase A
# baseline of the dartR function review (dartR.popgen at 88b9c22), updated
# in Phase C: assertions tagged [approved n] were changed with the matching
# approved change; see function-review/reports/dartR.popgen/gl.TajimasD.md
# (dartR.base). The neutral simulation no longer uses ms.

tj_genlight <- function() {
  x <- gl.keep.pop(possums.gl, pop.list = popNames(possums.gl)[1:2],
                   verbose = 0)
  gl.filter.callrate(x, threshold = 1, verbose = 0)
}

tj_run <- function(...) {
  out <- NULL
  capture.output(out <- suppressWarnings(gl.TajimasD(..., plot.out = FALSE)))
  out
}

# independent Tajima's D via pegas on pseudo-haplotypes (valid for D, which
# depends only on per-site allele counts)
tj_pegas <- function(x) {
  sapply(levels(pop(x)), function(p) {
    m <- as.matrix(x[pop(x) == p, ])
    m <- m[, colSums(m) > 0 & colSums(m) < 2 * nrow(m), drop = FALSE]
    seqs <- rbind(ifelse(m >= 1, "g", "a"), ifelse(m == 2, "g", "a"))
    unlist(pegas::tajima.test(ape::as.DNAbin(seqs))[c("D", "Pval.beta")])
  })
}

test_that("D matches pegas on complete data [approved 5]", {
  skip_if_not_installed("pegas")
  skip_if_not_installed("ape")
  x <- tj_genlight()
  r <- tj_run(x, verbose = 0)
  ref <- tj_pegas(x)
  expect_equal(r$population, c("A", "B"))
  expect_equal(r$S, c(176, 199))
  expect_equal(r$N, c(30, 30))
  expect_equal(r$D, unname(ref["D", ]), tolerance = 1e-10)
  expect_equal(r$Pval.beta, unname(ref["Pval.beta", ]), tolerance = 1e-4)
})

test_that("verbose = 0 is silent [approved 4]", {
  out <- capture.output(r <- gl.TajimasD(tj_genlight(), plot.out = FALSE,
                                         verbose = 0))
  expect_length(out, 0)
})

test_that("SilicoDArT input is rejected [approved 3]", {
  expect_error(tj_run(testset.gs, verbose = 0), "SilicoDArT")
})

test_that("a genlight without dartR flags works [approved 4]", {
  x <- tj_genlight()
  x@other$loc.metrics.flags <- NULL
  expect_equal(nrow(tj_run(x, verbose = 0)), 2)
})

test_that("missing data: Watterson's term uses per-site sample sizes [approved 5]", {
  x <- gl.keep.pop(platypus.gl, pop.list = popNames(platypus.gl)[1:2],
                   verbose = 0)
  r <- tj_run(x, verbose = 0)
  expect_equal(r$N, c(21.5, 15.766), tolerance = 1e-3)
  expect_equal(r$S, c(487, 452))
  # independent: pi and theta_W from per-locus counts
  m <- as.matrix(x[pop(x) == "SEVERN_BELOW", ])
  n_i <- 2 * colSums(!is.na(m))
  c_i <- colSums(m, na.rm = TRUE)
  ok <- n_i > 1
  p <- c_i[ok] / n_i[ok]
  pi <- sum(n_i[ok] / (n_i[ok] - 1) * 2 * p * (1 - p))
  seg <- ok & c_i > 0 & c_i < n_i
  thetaW <- sum(sapply(n_i[seg], function(n) 1 / sum(1 / (1:(n - 1)))))
  expect_equal(r$pi[2], pi, tolerance = 1e-10)
  expect_equal(sign(r$D[2]), sign(pi - thetaW))
})

test_that("utils.get.allele.freq honours verbose [approved 4]", {
  out <- capture.output(m <- utils.get.allele.freq(tj_genlight(), verbose = 0))
  expect_length(out, 0)
  expect_equal(nrow(m), 2 * 200)
  expect_equal(names(m), c("popn", "locus", "sum", "nobs", "nmissing",
                           "frequency", "n"))
})

test_that("simulated null of unlinked SNPs gives sim_pval [approved 1]", {
  x <- tj_genlight()
  d <- tempfile("simout")
  dir.create(d)
  r <- tj_run(x, rep = 500, seeds = 1, simulation.out = d, verbose = 0)
  expect_equal(r$sim_pval, c(0, 0))
  r2 <- tj_run(x, rep = 500, seeds = 1, verbose = 0)
  expect_identical(r$sim_pval, r2$sim_pval)
  sim <- scan(file.path(d, "Sim_TajimasD_A.txt"), quiet = TRUE)
  expect_length(sim, 500)
  # unlinked null for 60 sequences and 176 sites: sd about 0.2
  expect_true(sd(sim) > 0.15 && sd(sim) < 0.25)
  # the seed does not change the user's random number stream
  set.seed(42)
  a <- runif(1)
  set.seed(42)
  tj_run(x, rep = 10, seeds = 1, verbose = 0)
  expect_equal(runif(1), a)
})

test_that("ms.path without rep, and rep < 1, stop with a message [approved 6]", {
  x <- tj_genlight()
  expect_error(tj_run(x, ms.path = tempdir(), verbose = 0),
               "rep \\(number of simulated replicates\\) is required")
  expect_error(tj_run(x, rep = 0, verbose = 0), "positive number")
})
