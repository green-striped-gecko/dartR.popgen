#' @name gl.assign.pa
#' @title Use private alleles to identify populations as possible source
#' populations for an individual of unknown provenance
#' @family assignment
#'
#' @description
#' Identifies as putative source populations those for which the individual
#' of unknown provenance has no more private alleles than expected. The
#' putative source populations are retained and returned in a genlight object.
#'
#' @details
#' A private allele of the unknown with respect to a population is an allele
#' (reference or alternate) carried by the unknown at a locus where no member
#' of that population carries it. Loci with no genotype in the population, and
#' loci scored as missing in the unknown, are not counted.
#'
#' The expectation for each population is built from its own members: each
#' member is compared with the remaining members of its population, and the
#' counts of private alleles are log10(count + 1) transformed to approximate
#' a normal distribution with a mean and standard deviation. The count for the
#' unknown is compared with this expectation, and populations for which the
#' unknown has significantly more private alleles (one-tailed, p < alpha) are
#' unlikely to be its source. If every member of a population has the same
#' count (standard deviation zero or undefined), the population is retained
#' when the unknown's count does not exceed that common value.
#'
#' An excess of private alleles is evidence that the unknown does not belong
#' to a population provided the sample size is adequate, say >= 10.
#'
#' A population named "unknowns" in the input is taken to hold other
#' individuals of unknown provenance and is excluded from the reference
#' populations.
#'
#' WARNING: If a putative population is not in Hardy-Weinberg equilibrium, as
#' might occur if it includes F1 hybrids and backcrosses, then the standard
#' deviation for the expectation will be inflated. This inflation may result
#' in false identification of the population as a putative source for the
#' focal unknown individual. For this reason, you may wish to remove
#' populations that contain individuals likely to be subject to contemporary
#' hybridization or admixture.
#'
#' @param x Name of the input genlight object [required].
#' @param unknown SpecimenID label (indName) of the focal individual whose
#' provenance is unknown [required].
#' @param nmin Minimum sample size for a target population to be included in
#' the analysis [default 10].
#' @param n.best If given a value, dictates the best n = n.best populations to
#' retain for consideration (or more if there are ties) based on private
#' alleles. If not specified, then the putative source populations for which
#' the count is within expectation (p >= alpha) are retained [default NULL].
#' @param alpha The critical value used to select populations for which the
#' unknown individual has a count of private alleles within expectation
#' [default 0.01].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @return A genlight object containing the focal individual (assigned to
#' population "unknown") and the putative source populations, with loci
#' scored as missing in the unknown and monomorphic loci removed. If no
#' population qualifies, the genlight object contains only the unknown and a
#' warning is given.
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # Focal individual from the Macleay River (EmmacMaclGeor)
#' test <- gl.assign.pa(testset.gl, unknown = "UC_00146", nmin = 10,
#'                      verbose = 3)
#'
#' @seealso \code{\link{gl.assign.on.genotype}},
#' \code{\link{gl.assign.pca}}, \code{\link{gl.assign.mahalanobis}}
#'
#' @importFrom stats pnorm sd
#' @export

gl.assign.pa <- function(x,
                         unknown,
                         nmin = 10,
                         n.best = NULL,
                         alpha = 0.01,
                         verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  if (any(duplicated(indNames(x)))) {
    stop(error("Fatal Error: Duplicate individual names in genlight object.\n"))
  }
  if (any(is.na(indNames(x)))) {
    stop(error("Fatal Error: NA found in individual names.\n"))
  }
  if (length(unknown) != 1) {
    stop(error("Fatal Error: unknown must name a single individual.\n"))
  }
  if (!(unknown %in% indNames(x))) {
    stop(error(
      "Fatal Error: nominated focal individual (of unknown provenance) is",
      "not present in the dataset!\n"
    ))
  }
  if (is.null(pop(x))) {
    stop(error("Fatal Error: Population assignments (pop(x)) are NULL.\n"))
  }
  if (any(is.na(pop(x)))) {
    stop(error("Fatal Error: NA values found in population assignments.\n"))
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

  hard.min <- 10
  if (nmin < hard.min && verbose >= 1) {
    cat(warn(
      "  Warning: The specified minimum sample size is less than", hard.min,
      "individuals. Alleles present in the unknown risk being missed when",
      "sampling populations of fewer than", hard.min, "individuals.\n"
    ))
  }

  # DO THE JOB

  # Separate the unknown individual from x
  unknown.ind <- gl.keep.ind(x, ind.list = unknown, verbose = 0)
  pop(unknown.ind) <- rep("unknown", nInd(unknown.ind))
  knowns <- gl.drop.ind(x, ind.list = unknown, verbose = 0)
  # A population called "unknowns" holds other unknowns, not a reference
  if (any(popNames(knowns) == "unknowns")) {
    knowns <- gl.drop.pop(knowns, pop.list = "unknowns", verbose = 0)
  }

  # Remove all known populations with fewer than nmin individuals
  pop.sizes <- table(pop(knowns))
  pop.keep <- names(pop.sizes)[pop.sizes >= nmin]
  pop.toss <- names(pop.sizes)[pop.sizes < nmin]
  if (length(pop.keep) == 0) {
    stop(error(
      "Fatal Error: All target populations excluded based on minimum",
      "sample size.\n"
    ))
  }
  if (verbose >= 2 && length(pop.toss) > 0) {
    cat(report(
      "  Discarding", length(pop.toss), "populations with sample size <",
      nmin, "\n"
    ))
    if (verbose >= 3) {
      cat("   ", paste(pop.toss, collapse = ", "), "\n")
    }
  }
  knowns <- gl.keep.pop(knowns, pop.list = pop.keep, verbose = 0)

  # Remove loci scored as NA for the unknown. Index by position, because
  # data.frame() would rewrite hyphens in locus names.
  na.loc <- which(is.na(as.matrix(unknown.ind)[1, ]))
  if (length(na.loc) > 0) {
    knowns <- gl.drop.loc(knowns, loc.list = locNames(knowns)[na.loc],
                          verbose = 0)
    unknown.ind <- gl.drop.loc(unknown.ind,
                               loc.list = locNames(unknown.ind)[na.loc],
                               verbose = 0)
  }
  focal <- as.matrix(unknown.ind)[1, ]

  # One genotype matrix per population
  pop.list <- seppop(knowns)
  matrix.list <- lapply(pop.list, as.matrix)
  names(matrix.list) <- names(pop.list)

  # Private alleles of a focal genotype g with respect to a group summarised
  # by, per locus, the number of members carrying the alternate allele
  # (n.alt), carrying the reference allele (n.ref) and scored (n.scored).
  # The alternate allele is private when g >= 1 and n.alt == 0; the reference
  # allele is private when g <= 1 and n.ref == 0. Loci not scored in the
  # group or in the focal individual are not counted.
  count.pa <- function(g, n.alt, n.ref, n.scored) {
    ok <- !is.na(g) & n.scored > 0
    sum(ok & ((g >= 1 & n.alt == 0) | (g <= 1 & n.ref == 0)))
  }

  result <- data.frame(
    population = names(matrix.list),
    count = NA_integer_,
    z = NA_real_,
    p = NA_real_,
    flag = "",
    stringsAsFactors = FALSE
  )

  for (i in seq_along(matrix.list)) {
    pop.mat <- matrix.list[[i]]
    scored <- !is.na(pop.mat)
    has.alt <- scored & pop.mat >= 1
    has.ref <- scored & pop.mat <= 1
    n.alt <- colSums(has.alt)
    n.ref <- colSums(has.ref)
    n.scored <- colSums(scored)

    # Each member against the other members (leave one out)
    n.pa <- vapply(seq_len(nrow(pop.mat)), function(j) {
      count.pa(pop.mat[j, ],
               n.alt - has.alt[j, ],
               n.ref - has.ref[j, ],
               n.scored - scored[j, ])
    }, numeric(1))
    log.n.pa <- log10(n.pa + 1)
    mu <- mean(log.n.pa)
    sigma <- if (length(log.n.pa) > 1) sd(log.n.pa) else NA_real_

    # The unknown against all members
    result$count[i] <- count.pa(focal, n.alt, n.ref, n.scored)
    count.log <- log10(result$count[i] + 1)

    if (is.na(sigma) || sigma == 0) {
      # All members have the same count; the z-score is undefined
      result$flag[i] <- if (count.log <= mu) "yes" else "no"
    } else {
      result$z[i] <- (count.log - mu) / sigma
      result$p[i] <- round(1 - pnorm(result$z[i]), 6)
      result$flag[i] <- if (result$p[i] >= alpha) "yes" else "no"
    }
  }

  result <- result[order(result$p, decreasing = TRUE), ]
  row.names(result) <- NULL
  names(result) <- c("pop", "count", "Z-score", "p-value", "assign")

  if (verbose >= 3) {
    print(result)
  }

  # Retain the n.best populations (with ties), or the populations whose
  # count is within expectation
  if (!is.null(n.best)) {
    k <- min(n.best, nrow(result))
    p.k <- result[["p-value"]][k]
    if (is.na(p.k)) {
      pop.keep <- result$pop[seq_len(k)]
    } else {
      pop.keep <- result$pop[!is.na(result[["p-value"]]) &
                               result[["p-value"]] >= p.k]
    }
  } else {
    pop.keep <- result$pop[result$assign == "yes"]
  }

  if (length(pop.keep) > 0) {
    gl.out <- gl.keep.pop(knowns, pop.list = pop.keep, verbose = 0)
    gl.out <- gl.join(gl.out, unknown.ind, method = "end2end", verbose = 0)
  } else {
    gl.out <- unknown.ind
  }
  gl.out <- gl.filter.monomorphs(gl.out, verbose = 0)

  if (nInd(gl.out) == 1 && verbose >= 1) {
    cat(warn(
      "  Warning: Final genlight object contains only the unknown",
      "individual, no populations assigned.\n"
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
