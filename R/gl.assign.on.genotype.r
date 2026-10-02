#' @name gl.assign.on.genotype
#' @title Use genotype to identify populations as possible source populations
#' for an individual of unknown provenance
#' @family assignment
#'
#' @description
#' Identifies populations from which the individual of unknown provenance
#' could have been drawn, given its genotype and the allele frequencies in
#' the putative source populations. The putative source populations that
#' survive are retained and returned in a genlight object.
#'
#' @details
#' For each putative source population, the function computes the
#' log-likelihood of the unknown's multilocus genotype under Hardy-Weinberg
#' equilibrium, treating loci as independent. Allele frequencies are
#' estimated with the Rannala & Mountain (1997) prior, (x + 1/2) / (n + 1),
#' where x is the count of the alternate allele and n the number of gene
#' copies scored in the population. An allele absent from a population's
#' sample therefore has a small non-zero frequency that depends on the sample
#' size.
#'
#' All populations are scored on the same loci: loci scored as missing in the
#' unknown, or with no genotype in one or more of the retained populations,
#' are excluded. Otherwise a population with more missing loci would sum
#' fewer log-likelihood terms and appear more likely.
#'
#' The log-likelihoods are converted to AIC values (-2 log-likelihood) and AIC
#' weights. Because no parameters are fitted, the AIC weights equal the
#' likelihoods normalised across populations, that is, the posterior
#' probability of origin under equal prior probabilities. Populations with an
#' AIC weight below aic.threshold are unlikely to be the source. The weights
#' assume independent loci; with many linked loci they concentrate on the
#' single best population more than the data justify.
#'
#' A suitable estimate of the allele frequencies requires an adequate sample
#' size, say >= 10.
#'
#' A population named "unknowns" in the input is taken to hold other
#' individuals of unknown provenance and is excluded from the reference
#' populations.
#'
#' WARNING: If a putative population is not in Hardy-Weinberg equilibrium, as
#' might occur if it includes F1 hybrids and backcrosses, the genotype
#' likelihoods are misspecified. For this reason, you may wish to remove
#' populations that contain individuals likely to be subject to contemporary
#' hybridization or admixture.
#'
#' @param x Name of the input genlight object [required].
#' @param unknown SpecimenID label (indName) of the focal individual whose
#' provenance is unknown [required].
#' @param nmin Minimum sample size for a target population to be included in
#' the analysis [default 10].
#' @param n.best If given a value, dictates the best n = n.best populations to
#' retain for consideration (or more if there are ties) based on AIC weight.
#' If not specified, then the putative source populations with
#' AIC.wt >= aic.threshold are retained [default NULL].
#' @param aic.threshold The critical value used to select populations for
#' which there is some support as a putative source based on AIC weights
#' [default 0.05].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @return A genlight object containing the focal individual (assigned to
#' population "unknown") and the putative source populations based on AIC
#' weights, with monomorphic loci removed. If no population qualifies, the
#' genlight object contains only the unknown and a warning is given.
#'
#' @references
#' Rannala, B., & Mountain, J. L. (1997). Detecting immigration by using
#' multilocus genotypes. Proceedings of the National Academy of Sciences,
#' 94(17), 9197-9201.
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # Focal individual from the Macleay River (EmmacMaclGeor)
#' test <- gl.assign.on.genotype(testset.gl, unknown = "UC_00146",
#'                               nmin = 10, verbose = 3)
#'
#' @seealso \code{\link{gl.assign.pa}}, \code{\link{gl.assign.pca}},
#' \code{\link{gl.assign.mahalanobis}}
#'
#' @export

gl.assign.on.genotype <- function(x,
                                  unknown,
                                  nmin = 10,
                                  n.best = NULL,
                                  aic.threshold = 0.05,
                                  verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

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
  if (aic.threshold < 0 || aic.threshold > 1) {
    if (verbose >= 1) {
      cat(warn(
        "  Warning: the aic.threshold must be between 0 and 1, set to",
        "default 0.05\n"
      ))
    }
    aic.threshold <- 0.05
  }

  hard.min <- 10
  if (nmin < hard.min && verbose >= 1) {
    cat(warn(
      "  Warning: The specified minimum sample size is less than", hard.min,
      "individuals. Allele frequencies in the putative source populations",
      "may be poorly estimated.\n"
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

  # One genotype matrix per population
  pop.list <- seppop(knowns)
  matrix.list <- lapply(pop.list, as.matrix)
  names(matrix.list) <- names(pop.list)
  focal <- as.matrix(unknown.ind)[1, ]

  # Score every population on the same loci: scored in the unknown and in
  # every retained population
  n.scored <- sapply(matrix.list, function(m) colSums(!is.na(m)))
  n.scored <- matrix(n.scored, ncol = length(matrix.list))
  use <- !is.na(focal) & apply(n.scored > 0, 1, all)
  if (sum(use) == 0) {
    stop(error(
      "Fatal Error: no locus is scored in the unknown and in every retained",
      "population.\n"
    ))
  }
  if (verbose >= 2) {
    cat(report(
      "  Scoring", sum(use), "of", length(use), "loci (scored in the unknown",
      "and in every retained population)\n"
    ))
  }
  g <- focal[use]

  # Log-likelihood of the focal genotype under HWE, Rannala & Mountain (1997)
  # allele frequencies
  log.lik <- function(pop.mat) {
    pop.mat <- pop.mat[, use, drop = FALSE]
    n.copies <- 2 * colSums(!is.na(pop.mat))
    n.alt <- colSums(pop.mat, na.rm = TRUE)
    q <- (n.alt + 0.5) / (n.copies + 1)
    prob <- ifelse(g == 0, (1 - q)^2, ifelse(g == 1, 2 * q * (1 - q), q^2))
    sum(log(prob))
  }

  result <- data.frame(
    population = names(matrix.list),
    logL = vapply(matrix.list, log.lik, numeric(1)),
    stringsAsFactors = FALSE
  )
  result$aic <- -2 * result$logL
  result$delta.aic <- result$aic - min(result$aic)
  result$aic.wt <- exp(-0.5 * result$delta.aic) /
    sum(exp(-0.5 * result$delta.aic))
  result$flag <- ifelse(result$aic.wt >= aic.threshold, "yes", "no")

  result <- result[order(result$aic.wt, decreasing = TRUE), ]
  row.names(result) <- NULL
  names(result) <- c("population", "Log Likelihood", "AIC", "dAIC", "AIC.wt",
                     "assign")

  if (verbose >= 3) {
    print(result)
    cat(report(
      "\n  Best prospect for source population is", result$population[1],
      "\n\n"
    ))
  }

  # Retain the n.best populations (with ties), or the populations with
  # AIC.wt >= aic.threshold
  if (!is.null(n.best)) {
    k <- min(n.best, nrow(result))
    pop.keep <- result$population[result$AIC.wt >= result$AIC.wt[k]]
  } else {
    pop.keep <- result$population[result$assign == "yes"]
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
