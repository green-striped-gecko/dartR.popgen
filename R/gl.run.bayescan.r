#' @name gl.run.bayescan
#'
#' @title Identify Fst outlier loci with BayeScan
#'
#' @family selection
#'
#' @description
#' Runs BayeScan (Foll and Gaggiotti 2008) on a SNP genlight object to
#' identify loci under selection, and returns one row of results per locus
#' with its classification as neutral, under diversifying selection or under
#' balancing selection. BayeScan is not distributed with dartR.popgen: the
#' user must download it from \url{https://cmpg.unibe.ch/software/BayeScan/}
#' and give its location in \code{bayescan.path}. Please cite the original
#' publication if you use this function.
#'
#' @details
#' BayeScan splits the Fst of each locus in each population into a locus
#' effect (alpha), shared by all populations, and a population effect
#' (beta), shared by all loci. A locus needs a non-zero alpha to explain its
#' allele frequencies when it is under selection: a positive alpha means
#' more differentiation than expected (diversifying selection), a negative
#' alpha less (balancing selection). BayeScan estimates by reversible-jump
#' MCMC the posterior probability that alpha is non-zero, and from it a
#' q-value, the smallest false discovery rate at which the locus would be
#' called an outlier. Here a locus is classified as under selection when its
#' q-value is at most \code{q.threshold}, and as diversifying or balancing by
#' the sign of alpha.
#'
#' \code{prior.odds} sets how much more likely the neutral model is than the
#' model with selection before the data are seen. The BayeScan default of 10
#' gives many false positives when thousands of loci are tested, or when the
#' populations have a hierarchical structure or a history of range
#' expansion (Lotterhos and Whitlock 2014); larger values, e.g. 100, make
#' the test more conservative.
#'
#' Each population in \code{pop(x)} is a sample in the analysis, so subset
#' \code{x} to the populations to be compared (\code{gl.keep.pop}) before
#' running the function. Loci that are monomorphic in \code{x}, or have no
#' calls, carry no information on differentiation; they are not sent to
#' BayeScan and get NA in the results. A population with no calls at a
#' locus is passed to BayeScan as a sample of zero gene copies.
#'
#' BayeScan has no option to set the random seed, so results vary slightly
#' between runs. Check that the chain has converged before using the
#' results, e.g. with \code{coda::effectiveSize(coda::mcmc(res$trace))} or
#' \code{coda::geweke.diag()}; if it has not, increase \code{burnin} and
#' \code{n.iter}. Runs on tens of thousands of loci take hours.
#'
#' The input file and all BayeScan output files are written to
#' \code{outpath}, with names starting with \code{outfile}: the input file
#' (_input.txt), the results per locus (_fst.txt), the MCMC trace of the
#' log-likelihood and the Fst of each population (.sel), the acceptance
#' rates (_AccRte.txt), a summary of the run settings (_Verif.txt) and the
#' screen output of BayeScan (.log). Files from an earlier run with the same
#' \code{outfile} and \code{outpath} are overwritten.
#'
#' @param x Name of the genlight object containing SNP data, with at least
#'   two populations defined [required].
#' @param bayescan.path Path to the folder holding the BayeScan executable,
#'   or to the executable itself. In a folder, the function looks for the
#'   file names of the BayeScan 2.1 distribution for the operating system
#'   (e.g. BayeScan2.1_macos64bits) and for "bayescan" [required].
#' @param prior.odds Prior odds for the neutral model. Option -pr_odds in
#'   BayeScan [default 10].
#' @param n.iter Number of MCMC iterations kept for the posterior. Option -n
#'   in BayeScan [default 5000].
#' @param thin Thinning interval: number of iterations between two kept
#'   iterations. Option -thin in BayeScan [default 10].
#' @param n.pilot Number of pilot runs used to tune the proposal
#'   distributions. Option -nbp in BayeScan [default 20].
#' @param pilot.length Length of each pilot run. Option -pilot in BayeScan
#'   [default 5000].
#' @param burnin Number of burn-in iterations. Option -burn in BayeScan
#'   [default 50000].
#' @param n.cores Number of threads used by BayeScan. If NULL, BayeScan uses
#'   all available CPUs. Option -threads in BayeScan [default NULL].
#' @param q.threshold Largest q-value (false discovery rate) for a locus to
#'   be classified as under selection [default 0.05].
#' @param outfile Prefix of the names of the input and output files of
#'   BayeScan [default "bayescan"].
#' @param outpath Path where to save the input and output files of BayeScan
#'   [default global working directory or if not specified, tempdir()].
#' @param plot.display Specify if the plot of Fst against log10(q-value) is
#'   to be displayed in the graphics window [default TRUE].
#' @param plot.theme User specified theme [default theme_dartR()].
#' @param plot.colors Vector with three colour names for neutral loci,
#'   loci under diversifying selection and loci under balancing selection
#'   [default c("grey60", "deeppink", "dodgerblue")].
#' @param plot.dir Directory to save the plot RDS files [default as
#'   specified by the global working directory or tempdir()].
#' @param plot.file Filename (minus extension) for the RDS plot file
#'   [Required for plot save].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'   brief progress messages; 3, progress and results summary; 5, full
#'   report [default 2, unless specified using gl.set.verbosity].
#'
#' @return A list with four elements:
#' \itemize{
#' \item results: a data frame with one row per input locus, in input
#'   order: locus (locus name), prob (posterior probability of the model with
#'   selection), log10.PO (log10 of the posterior odds for the model with
#'   selection; BayeScan writes 1000 when prob is 1), q.value, alpha (locus
#'   effect), fst (Fst of the locus averaged over populations) and selection
#'   (neutral, diversifying or balancing). Loci not tested are NA.
#' \item fst.pop: a data frame with the posterior mean and the 2.5\% and
#'   97.5\% quantiles of the Fst of each population.
#' \item trace: a data frame with the MCMC trace of the log-likelihood
#'   (logL) and of the Fst of each population, one row per kept iteration.
#' \item plot: a ggplot of Fst against log10(q-value), with a dashed line at
#'   \code{q.threshold}. q-values below 0.0001 are drawn at 0.0001.
#' }
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#'   \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \dontrun{
#' # BayeScan downloaded from https://cmpg.unibe.ch/software/BayeScan/
#' x <- gl.keep.pop(possums.gl, pop.list = c("A", "B"))
#' res <- gl.run.bayescan(x, bayescan.path = "path/to/BayeScan2.1/binaries",
#'                        prior.odds = 100)
#' # loci under diversifying selection
#' res$results[which(res$results$selection == "diversifying"), ]
#' }
#'
#' @references
#' \itemize{
#' \item Foll, M. and Gaggiotti, O. (2008) A genome-scan method to identify
#'   selected loci appropriate for both dominant and codominant markers: a
#'   Bayesian perspective. Genetics 180: 977-993.
#'   \doi{10.1534/genetics.108.092221}
#' \item Lotterhos, K.E. and Whitlock, M.C. (2014) Evaluation of demographic
#'   history and neutral parameterization on the performance of FST outlier
#'   tests. Molecular Ecology 23: 2178-2192. \doi{10.1111/mec.12725}
#' }
#'
#' @seealso \code{\link{gl.outflank}}, \code{\link[dartR.base]{gl2bayescan}}
#'
#' @importFrom stats quantile setNames
#' @importFrom utils read.table tail
#' @export

gl.run.bayescan <- function(x,
                            bayescan.path,
                            prior.odds = 10,
                            n.iter = 5000,
                            thin = 10,
                            n.pilot = 20,
                            pilot.length = 5000,
                            burnin = 50000,
                            n.cores = NULL,
                            q.threshold = 0.05,
                            outfile = "bayescan",
                            outpath = NULL,
                            plot.display = TRUE,
                            plot.theme = theme_dartR(),
                            plot.colors = c("grey60", "deeppink",
                                            "dodgerblue"),
                            plot.dir = NULL,
                            plot.file = NULL,
                            verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # SET WORKING DIRECTORY
  outpath <- gl.check.wd(outpath, verbose = 0)
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  if (is.null(pop(x)) || any(is.na(pop(x)))) {
    stop(error("Fatal Error: Every individual must be assigned to a",
               "population.\n"))
  }
  pops <- droplevels(pop(x))
  if (nlevels(pops) < 2) {
    stop(error("Fatal Error: BayeScan needs at least two populations;",
               "found", nlevels(pops), "\n"))
  }
  counts <- c(n.iter = n.iter, thin = thin, n.pilot = n.pilot,
              pilot.length = pilot.length, burnin = burnin)
  if (!is.null(n.cores)) counts <- c(counts, n.cores = n.cores)
  bad <- !is.finite(counts) | counts < 0 | counts != round(counts)
  if (any(bad)) {
    stop(error("Fatal Error:", paste(names(counts)[bad], collapse = ", "),
               "must be whole numbers of zero or more\n"))
  }
  if (!is.numeric(prior.odds) || length(prior.odds) != 1 ||
      !(prior.odds > 0)) {
    stop(error("Fatal Error: prior.odds must be a number above 0\n"))
  }
  if (!is.numeric(q.threshold) || length(q.threshold) != 1 ||
      !(q.threshold > 0 && q.threshold < 1)) {
    stop(error("Fatal Error: q.threshold must be a number between 0 and 1\n"))
  }
  if (length(plot.colors) != 3) {
    stop(error("Fatal Error: plot.colors must have three colours (neutral,",
               "diversifying, balancing)\n"))
  }

  # Locate the executable: a folder is searched for the names used in the
  # BayeScan 2.1 distribution for this operating system, then "bayescan"
  os <- tolower(Sys.info()[["sysname"]])
  if (dir.exists(bayescan.path)) {
    candidates <- switch(os,
      windows = c("BayeScan2.1_win32bits_cmd_line.exe", "bayescan.exe"),
      darwin = c("BayeScan2.1_macos64bits", "bayescan"),
      c("BayeScan2.1_linux64bits", "bayescan_linux64", "bayescan")
    )
    exec <- file.path(bayescan.path, candidates)
    exec <- exec[file.exists(exec) & !dir.exists(exec)]
    if (length(exec) == 0) {
      stop(error(
        "Fatal Error: Cannot find a BayeScan executable (",
        paste(candidates, collapse = ", "),
        ") in the folder given by bayescan.path:", bayescan.path,
        "\n  Download BayeScan from https://cmpg.unibe.ch/software/BayeScan/",
        "and set bayescan.path to the folder holding the executable, or to",
        "the executable itself.\n"
      ))
    }
    exec <- exec[1]
  } else if (file.exists(bayescan.path)) {
    exec <- bayescan.path
  } else {
    stop(error(
      "Fatal Error: bayescan.path does not exist:", bayescan.path,
      "\n  Download BayeScan from https://cmpg.unibe.ch/software/BayeScan/",
      "and set bayescan.path to the folder holding the executable, or to",
      "the executable itself.\n"
    ))
  }
  exec <- normalizePath(exec)
  if (os != "windows" && file.access(exec, mode = 1) != 0) {
    stop(error("Fatal Error: BayeScan executable", exec, "is not executable.",
               "Make it executable with: chmod +x", exec, "\n"))
  }

  # DO THE JOB
  # Monomorphic loci and loci without calls are not sent to BayeScan
  alf <- gl.alf(x, verbose = 0)$alf2
  tested <- which(!is.na(alf) & alf > 0 & alf < 1)
  if (length(tested) == 0) {
    stop(error("Fatal Error: No polymorphic loci to test\n"))
  }
  if (length(tested) < nLoc(x) && verbose >= 2) {
    cat(warn("  Warning:", nLoc(x) - length(tested), "loci are monomorphic",
             "or have no calls; they are not tested and get NA\n"))
  }
  xt <- x[, tested]

  # Population x locus combinations without calls are written by gl2bayescan
  # as samples of zero gene copies
  zero.cells <- sum(vapply(levels(pops), function(p) {
    xp <- xt[pops == p, ]
    sum(glNA(xp, alleleAsUnit = FALSE) == nInd(xp))
  }, numeric(1)))
  if (zero.cells > 0 && verbose >= 1) {
    cat(warn("  Warning:", zero.cells, "population x locus combinations",
             "have no calls and are passed to BayeScan as samples of zero",
             "gene copies. Consider filtering on call rate first.\n"))
  }

  infile <- paste0(outfile, "_input.txt")
  gl2bayescan(xt, outfile = infile, outpath = outpath, verbose = 0)
  # gl2bayescan numbers the populations in their order of appearance in
  # pop(x); BayeScan reports the Fst of each population in that order
  pop.names <- as.character(unique(pop(xt)))

  # Output files of an earlier run would be read if this run fails
  out.files <- file.path(outpath, paste0(outfile, c("_fst.txt", ".sel",
                                                    "_AccRte.txt",
                                                    "_Verif.txt")))
  unlink(out.files)
  logfile <- file.path(outpath, paste0(outfile, ".log"))
  args <- c(shQuote(file.path(outpath, infile)),
            "-od", shQuote(outpath),
            "-o", shQuote(outfile),
            "-n", n.iter,
            "-thin", thin,
            "-nbp", n.pilot,
            "-pilot", pilot.length,
            "-burn", burnin,
            "-pr_odds", prior.odds,
            if (!is.null(n.cores)) c("-threads", n.cores))

  if (verbose >= 2) {
    cat(report("  Running BayeScan on", length(tested), "loci and",
               length(pop.names), "populations. Progress is written to",
               logfile, "\n"))
  }
  timing <- proc.time()
  status <- system2(exec, args, stdout = logfile, stderr = logfile)
  if (status != 0 || !file.exists(out.files[1])) {
    msg <- if (file.exists(logfile)) {
      tail(readLines(logfile, warn = FALSE), 5)
    } else {
      character(0)
    }
    stop(error("Fatal Error: BayeScan failed (exit status", status, "):",
               paste(msg, collapse = " "), "\n  See", logfile, "\n"))
  }
  if (verbose >= 3) {
    cat(report("  BayeScan finished in",
               round((proc.time() - timing)[["elapsed"]] / 60, 1),
               "minutes\n"))
  }

  # Results per locus. The header of _fst.txt has one name fewer than the
  # rows, which start with the index of the locus in the input file
  bs <- read.table(out.files[1], skip = 1,
                   col.names = c("index", "prob", "log10.PO", "q.value",
                                 "alpha", "fst"))
  if (nrow(bs) != length(tested)) {
    stop(error("Fatal Error: BayeScan returned", nrow(bs), "loci but",
               length(tested), "were sent. See", logfile, "\n"))
  }
  bs <- bs[order(bs$index), ]
  res <- data.frame(locus = locNames(x), prob = NA_real_,
                    log10.PO = NA_real_, q.value = NA_real_,
                    alpha = NA_real_, fst = NA_real_)
  res[tested, c("prob", "log10.PO", "q.value", "alpha", "fst")] <-
    bs[, c("prob", "log10.PO", "q.value", "alpha", "fst")]
  sel <- ifelse(res$q.value <= q.threshold & res$alpha > 0, "diversifying",
                ifelse(res$q.value <= q.threshold & res$alpha < 0,
                       "balancing", "neutral"))
  res$selection <- factor(sel, levels = c("neutral", "diversifying",
                                          "balancing"))

  # MCMC trace; the first column of each row is the iteration, which
  # read.table takes as row names because the header is one name short
  trace <- read.table(out.files[2], header = TRUE)
  rownames(trace) <- NULL
  colnames(trace) <- c("logL", pop.names)
  fst.trace <- trace[, pop.names, drop = FALSE]
  fst.pop <- data.frame(
    population = pop.names,
    fst.mean = colMeans(fst.trace),
    fst.lower = apply(fst.trace, 2, quantile, probs = 0.025),
    fst.upper = apply(fst.trace, 2, quantile, probs = 0.975),
    row.names = NULL
  )

  if (verbose >= 3) {
    n.sel <- table(res$selection)
    cat(report("  Loci tested:", length(tested), "of", nLoc(x), "\n"))
    cat(report("  Loci with q-value <=", q.threshold, ": diversifying",
               n.sel[["diversifying"]], "; balancing",
               n.sel[["balancing"]], "\n"))
    cat(report("  Posterior mean Fst per population:\n"))
    print(fst.pop, row.names = FALSE)
  }

  # PLOT
  log10.q <- fst <- selection <- NULL
  pdat <- res[tested, ]
  # BayeScan writes a q-value of 0 when every kept iteration includes
  # selection; it is drawn at 0.0001 as in the plot script of BayeScan
  pdat$log10.q <- log10(pmax(pdat$q.value, 1e-4))
  p1 <- ggplot(pdat, aes(x = log10.q, y = fst, colour = selection)) +
    geom_point(alpha = 0.7) +
    geom_vline(xintercept = log10(q.threshold), linetype = "dashed") +
    scale_x_reverse() +
    scale_colour_manual(values = setNames(plot.colors,
                                          levels(res$selection)),
                        drop = FALSE) +
    labs(x = "log10(q-value)", y = "Fst", colour = "Selection") +
    plot.theme

  # PRINTING OUTPUTS
  if (plot.display) {
    print(p1)
  }
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p1, dir = plot.dir, file = plot.file,
                           verbose = verbose)
  }
  if (verbose >= 2) {
    cat(report("  BayeScan files written to", outpath, "\n"))
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN
  return(list(results = res, fst.pop = fst.pop, trace = trace, plot = p1))
}
