# Helpers for the gl.assign.* tests.

# Build a compliant genlight from a 0/1/2/NA matrix (rows = individuals)
make_gl <- function(m, pop, ind = NULL) {
  if (is.null(ind)) ind <- paste0("i", seq_len(nrow(m)))
  gl <- new("genlight", m, ind.names = ind,
            loc.names = paste0("L", seq_len(ncol(m))),
            pop = factor(pop), ploidy = 2L)
  utils::capture.output(gl <- gl.compliance.check(gl, verbose = 0))
  gl
}

# Simulate populations under HWE with population-specific allele
# frequencies scattered around a common frequency
sim_pops <- function(n.pop = 3, n.ind = 20, n.loc = 200, fst = 0.1,
                     seed = 1) {
  set.seed(seed)
  p0 <- stats::runif(n.loc, 0.1, 0.9)
  a <- p0 * (1 - fst) / fst
  b <- (1 - p0) * (1 - fst) / fst
  m <- do.call(rbind, lapply(seq_len(n.pop), function(k) {
    q <- stats::rbeta(n.loc, a, b)
    t(replicate(n.ind, stats::rbinom(n.loc, 2, q)))
  }))
  pop <- rep(LETTERS[seq_len(n.pop)], each = n.ind)
  make_gl(m, pop, ind = paste0(pop, seq_len(nrow(m))))
}

quiet <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressWarnings(expr))
  out
}

kept_pops <- function(gl) sort(setdiff(popNames(gl), "unknown"))
