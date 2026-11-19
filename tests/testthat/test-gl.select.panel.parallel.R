# Slow file: it starts PSOCK workers twice. Skipped on CRAN to keep the check
# time under 10 minutes; runs with devtools::test() and on GitHub Actions.
skip_on_cran()

# gl.select.panel.parallel keeps the state of each chain on the worker that
# runs it. With fewer cores than chains (n_cores < n_chains) the chains share
# the workers. A worker used to hold one chain only: the last chain
# initialised on it replaced the earlier one, and that chain was then
# advanced once for every chain id sent to the worker. With n_chains = 4 and
# n_cores = 3, chain 4 replaced chain 1, so three chains ran and the result
# listed chain 4 twice.

par_panel <- function(...) {
  out <- NULL
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  utils::capture.output(
    out <- suppressWarnings(suppressMessages(gl.select.panel.parallel(...)))
  )
  out
}

test_that("every chain runs when there are fewer cores than chains", {
  args <- list(x = dartR.data::possums.gl, nl = 10, parameter = "He", n_chains = 3,
               n_iter = 4, burst = 2, plot.out = FALSE, verbose = 0)
  shared <- do.call(par_panel, c(args, n_cores = 2))
  own <- do.call(par_panel, c(args, n_cores = 3))

  chains <- shared@other$sa_all_chains
  expect_length(chains, 3)
  expect_equal(vapply(chains, function(z) z$id, numeric(1)), 1:3)
  # each chain advances n_iter iterations, no more
  expect_equal(vapply(chains, function(z) max(z$history$iteration), numeric(1)),
               rep(4, 3))

  # a chain has its own random stream, so it gives the same result whether
  # it has a worker to itself or shares one
  expect_equal(lapply(chains, function(z) z$history),
               lapply(own@other$sa_all_chains, function(z) z$history))
  expect_equal(adegenet::locNames(shared), adegenet::locNames(own))
})
