#' @name gl.assign.mahalanobis
#' @title Calculate probabilities of assignment of an individual of unknown
#' provenance to populations based on Mahalanobis distance
#' @family assignment
#'
#' @description
#' Assigns an individual of unknown provenance to one or more target
#' populations based on the unknown individual's distance from each
#' population centroid in a PCA, measured as a squared Mahalanobis distance
#' and tested against the spread of the population.
#'
#' @details
#' The following process is followed:
#' \enumerate{
#' \item Loci scored as missing in the unknown are removed, and remaining
#' missing values are imputed (method "neighbour" of \code{gl.impute()}).
#' \item A PCA is undertaken on the populations to yield a series of
#' orthogonal axes.
#' \item A workable subset of dimensions is chosen: dim.limit, or the number
#' of leading eigenvalues larger than expected under the broken-stick model,
#' whichever is the smaller (at least one).
#' \item The squared Mahalanobis distance (D2) of the unknown from the
#' centroid of each population is calculated in those dimensions, using the
#' population's own mean and covariance matrix.
#' \item D2 is tested with Hotelling's T-squared test for a new observation:
#' with n individuals in the population and p dimensions,
#' F = n (n - p) D2 / ((n + 1) p (n - 1)) follows an F distribution with p
#' and n - p degrees of freedom. This test allows for the mean and covariance
#' being estimated from a small sample; the chi-square approximation does not
#' and rejects true source populations too often when n is small.
#' }
#'
#' There are three considerations to assignment. First, consider only those
#' populations for which the unknown has no excess of private alleles. This
#' can be evaluated with \code{gl.assign.pa()}. Next, consider the PCA plot
#' for the remaining populations and the position of the unknown in relation
#' to their confidence ellipses, as produced by \code{gl.assign.pca()}. The
#' third step (delivered by this function) considers the position of the
#' unknown in relation to the confidence envelope in all selected dimensions
#' of the ordination. The larger the probability, the greater the confidence
#' in the assignment. If the unknown is an outlier with respect to a
#' population, say at less than 0.001 probability, that population can be
#' eliminated from further consideration.
#'
#' Warning: gl.assign.mahalanobis() treats each selected dimension equally,
#' without regard to the percentage variation explained. If the unknown is an
#' outlier in a lower dimension with an explanatory variance of, say, 0.1\%,
#' the putative population will be eliminated. This is why the function only
#' uses dimensions retained by the broken-stick criterion.
#'
#' Populations with no more individuals than selected dimensions, or with a
#' singular covariance matrix (for example when imputation makes individuals
#' identical), cannot be tested. They are reported and retained.
#'
#' Each of these approaches provides evidence, none are 100\% definitive.
#' They need to be interpreted cautiously.
#'
#' @param x Name of the input genlight object [required].
#' @param unknown Identity label of the focal individual whose provenance is
#' unknown [required].
#' @param nmin Minimum sample size for a target population to be included in
#' the analysis [default 10].
#' @param dim.limit Maximum number of dimensions to consider for the
#' confidence envelope [default nPop(x) - 1].
#' @param plevel Probability level below which a population is eliminated as
#' a putative source [default 0.001].
#' @param n.best If given a value, dictates the best n = n.best populations to
#' retain for consideration (or more if there are ties). If not specified,
#' the populations for which the probability is >= plevel are retained
#' [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @return A genlight object containing the unknown individual (assigned to
#' population "unknown") and the putative source populations, with their
#' original (not imputed) genotypes. Loci scored as missing in the unknown and
#' monomorphic loci are removed.
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \donttest{
#' # Focal individual from the Macleay River (EmmacMaclGeor)
#' test <- gl.assign.mahalanobis(testset.gl, unknown = "UC_00146", verbose = 3)
#' }
#'
#' @seealso \code{\link{gl.assign.pa}}, \code{\link{gl.assign.on.genotype}},
#' \code{\link{gl.assign.pca}}
#'
#' @importFrom stats cov mahalanobis pf
#' @export

gl.assign.mahalanobis <- function(x,
                                  unknown,
                                  nmin = 10,
                                  dim.limit = NULL,
                                  plevel = 0.001,
                                  n.best = NULL,
                                  verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  if (any(is.na(indNames(x))) || any(duplicated(indNames(x)))) {
    stop(error(
      "Fatal Error: Individual names must be unique and non-missing.\n"
    ))
  }
  if (length(unknown) != 1) {
    stop(error("Fatal Error: unknown must name a single individual.\n"))
  }
  if (!(unknown %in% indNames(x))) {
    stop(error(
      "Fatal Error: Unknown must be listed among the individuals in the",
      "genlight object!\n"
    ))
  }
  if (is.null(pop(x)) || any(is.na(pop(x)))) {
    stop(error(
      "Fatal Error: Population assignments are missing or incomplete.\n"
    ))
  }
  if (nLoc(x) < nPop(x)) {
    stop(error("Fatal Error: Number of loci less than number of populations\n"))
  }
  if (plevel > 1 || plevel < 0) {
    if (verbose >= 1) {
      cat(warn(
        "  Warning: Value of plevel must be between 0 and 1, set to 0.001\n"
      ))
    }
    plevel <- 0.001
  }
  if (nmin <= 0) {
    if (verbose >= 1) {
      cat(warn(
        "  Warning: the minimum size of the target population must be",
        "greater than zero, set to 10\n"
      ))
    }
    nmin <- 10
  }
  if (!is.null(n.best) && n.best < 1) {
    if (verbose >= 1) {
      cat(warn(
        "  Warning: the n.best parameter for retention of best match",
        "populations must be a positive integer, set to NULL\n"
      ))
    }
    n.best <- NULL
  }
  if (is.null(dim.limit)) {
    dim.limit <- nPop(x) - 1
  }

  hard.min <- 10
  if (nmin < hard.min && verbose >= 1) {
    cat(warn(
      "  Warning: The specified minimum sample size is less than", hard.min,
      "individuals. Population covariance matrices are poorly estimated in",
      "small samples.\n"
    ))
  }

  # DO THE JOB

  # Remove all known populations with fewer than nmin individuals
  pop.sizes <- table(pop(x))
  pop.keep <- names(pop.sizes)[pop.sizes >= nmin]
  # Input chained from another gl.assign.* function holds the unknown as
  # population "unknown"; it is not a discarded reference population
  pop.toss <- setdiff(names(pop.sizes)[pop.sizes < nmin], "unknown")
  if (verbose >= 2 && length(pop.toss) > 0) {
    cat(report(
      "  Discarding", length(pop.toss), "populations with sample size <",
      nmin, "\n"
    ))
    if (verbose >= 3) {
      cat("   ", paste(pop.toss, collapse = ", "), "\n")
    }
  }

  # Assign the unknown individual to its own population and keep it with the
  # populations of adequate size
  tmp <- as.character(pop(x))
  tmp[indNames(x) == unknown] <- "unknown"
  pop(x) <- as.factor(tmp)
  pop.keep <- setdiff(pop.keep, "unknown")
  pop.keep <- pop.keep[pop.keep %in% levels(pop(x))]
  if (length(pop.keep) == 0) {
    stop(error(
      "Fatal Error: All target populations excluded based on minimum",
      "sample size.\n"
    ))
  }
  x <- gl.keep.pop(x, pop.list = c(pop.keep, "unknown"), verbose = 0)

  # Remove loci that are missing in the unknown
  na.loc <- which(is.na(as.matrix(x[indNames(x) == unknown, ])[1, ]))
  if (length(na.loc) > 0) {
    x <- gl.drop.loc(x, loc.list = locNames(x)[na.loc], verbose = 0)
  }

  # Keep the observed genotypes for the returned object
  x.obs <- x

  # Mahalanobis distances need a dense matrix
  if (verbose >= 3) {
    cat(report("  Rendering the data matrix dense by imputation\n"))
  }
  x <- gl.impute(x, method = "neighbour", verbose = 0)

  if (verbose >= 3) {
    cat(report("  Undertaking a PCA\n"))
  }
  pcoa <- gl.pcoa(x, nfactors = dim.limit, verbose = 0)

  if (any(!is.finite(pcoa$scores))) {
    stop(error("Fatal Error: PCA ordination returned non-finite values.\n"))
  }
  if (any(pcoa$eig <= 0)) {
    stop(error(
      "Fatal Error: PCA ordination returned non-positive eigenvalues, which",
      "would give negative Mahalanobis distances. Run again with a smaller",
      "dim.limit.\n"
    ))
  }

  # Broken stick: number of leading eigenvalues larger than their expected
  # value under the broken-stick model
  eig <- pcoa$eig
  n.eig <- length(eig)
  broken.stick <- sum(eig) * rev(cumsum(1 / n.eig:1)) / n.eig
  above <- eig > broken.stick
  first.est <- if (above[1]) {
    if (all(above)) n.eig else which(!above)[1] - 1
  } else {
    0
  }
  if (first.est == 0 && verbose >= 2) {
    cat(warn(
      "  Warning: no eigenvalue exceeds the broken-stick expectation; using",
      "one dimension\n"
    ))
  }
  dim <- max(1, min(first.est, dim.limit, ncol(pcoa$scores)))
  if (verbose >= 2) {
    cat(report(
      "  Number of leading dimensions with substantial eigenvalues",
      "(broken-stick criterion):", first.est, ". Hardwired limit", dim.limit,
      "\n"
    ))
    cat(report("    Dimension of confidence envelope set at", dim, "\n"))
  }

  scores <- pcoa$scores[, seq_len(dim), drop = FALSE]
  pops <- as.character(pop(x))
  unknown.score <- scores[pops == "unknown", , drop = FALSE]

  df <- data.frame(
    pop = pop.keep,
    MahalD = NA_real_,
    pval = NA_real_,
    assign = NA_character_,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(pop.keep)) {
    A <- scores[pops == pop.keep[i], , drop = FALSE]
    n <- nrow(A)
    # Hotelling's test needs n > dim and an invertible covariance matrix
    if (n <= dim) {
      next
    }
    covariance <- stats::cov(A)
    if (any(!is.finite(covariance)) || rcond(covariance) < 1e-10) {
      next
    }
    D2 <- stats::mahalanobis(unknown.score, colMeans(A), covariance)
    f.stat <- n * (n - dim) * D2 / ((n + 1) * dim * (n - 1))
    df$MahalD[i] <- D2
    df$pval[i] <- stats::pf(f.stat, dim, n - dim, lower.tail = FALSE)
    df$assign[i] <- if (df$pval[i] >= plevel) "yes" else "no"
  }

  untested <- df$pop[is.na(df$pval)]
  if (length(untested) > 0 && verbose >= 1) {
    cat(warn(
      "  Warning: populations not tested (no more individuals than",
      "dimensions, or a singular covariance matrix), retained:",
      paste(untested, collapse = ", "), "\n"
    ))
  }

  # Order the dataframe in descending order on pval; untested last
  df <- df[order(df$pval, decreasing = TRUE), ]
  row.names(df) <- NULL

  if (verbose >= 3) {
    cat("  Assignment of unknown individual:", unknown, "\n")
    cat("  Alpha level of significance:", plevel, "\n")
    print(df)
    best <- df$pop[df$assign %in% "yes"][1]
    if (!is.na(best)) {
      cat(report(
        "  Best assignment is the population with the largest probability",
        "of assignment, in this case", best, "\n"
      ))
    }
  }

  # Retain the n.best populations (with ties), or the populations not
  # eliminated
  if (!is.null(n.best)) {
    k <- min(n.best, sum(!is.na(df$pval)))
    if (k > 0) {
      pop.keep <- df$pop[!is.na(df$pval) & df$pval >= df$pval[k]]
    } else {
      pop.keep <- character(0)
    }
  } else {
    pop.keep <- df$pop[!(df$assign %in% "no")]
  }

  gl.out <- gl.keep.pop(x.obs, pop.list = c(pop.keep, "unknown"), verbose = 0)
  gl.out <- gl.filter.monomorphs(gl.out, verbose = 0)

  if (nInd(gl.out) == 1 && verbose >= 1) {
    cat(warn(
      "  Warning: Final genlight object contains only the unknown",
      "individual, no populations assigned.\n"
    ))
  }
  if (verbose >= 2) {
    cat(report(
      "  Returning a genlight object with the putative source populations",
      "and the unknown\n"
    ))
  }

  # ADD TO HISTORY
  nh <- length(gl.out@other$history)
  gl.out@other$history[[nh + 1]] <- match.call()

  # FLAG SCRIPT END
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }

  return(gl.out)
}
