#' @name gl.assign.pca
#' @title Eliminate from consideration putative source populations for a
#' specified individual of unknown provenance using PCA
#' @family assignment
#'
#' @description
#' Eliminates from consideration putative source populations for a specified
#' individual of unknown provenance based on its position relative to the
#' confidence ellipse of each putative source population in the top two
#' dimensions of a PCA. Populations for which the unknown lies outside the
#' specified confidence ellipse are eliminated.
#'
#' @details
#' There are four approaches to population assignment of which this is one.
#' \enumerate{
#' \item Eliminate those populations from consideration where the genotype of
#' the unknown is not consistent with their allelic profiles. This can be
#' evaluated with \code{gl.assign.on.genotype()}.
#' \item Eliminate those populations for which the unknown has substantial
#' private alleles. Substantial numbers of private alleles are an indication
#' that the unknown does not belong to a target population (provided that the
#' sample size is adequate, say >= 10). This can be evaluated with
#' \code{gl.assign.pa()}.
#' \item Consider the assignment probabilities using
#' \code{gl.assign.mahalanobis()}. This approach calculates the squared
#' Generalised Distance (Mahalanobis distance) of the unknown from the
#' centroid of each remaining putative source population and tests it
#' against the population's spread. This index takes into account the
#' position of the unknown in relation to the confidence envelope in all
#' selected dimensions of the ordination.
#' \item Consider the PCA plot for populations and the position of the unknown
#' in relation to confidence ellipses as produced by \code{gl.assign.pca()}.
#' }
#' Each of these approaches is useful for decisions on which populations to
#' eliminate as putative sources. They provide evidence for a decision on
#' population assignment, but none are 100\% definitive, and they need to be
#' interpreted cautiously.
#'
#' This function considers only the top two dimensions of the ordination. An
#' unknown outside the ellipse in two dimensions has a two-dimensional squared
#' Mahalanobis distance above the chi-square (2 df) critical value. Adding
#' dimensions never reduces that distance, but the critical value also rises
#' with the number of dimensions, so an unknown outside the ellipse in two
#' dimensions may still lie inside the envelope in more dimensions.
#' Conversely, an unknown inside the ellipse in two dimensions may lie outside
#' the envelope in deeper dimensions. Use a stringent plevel (e.g. the
#' default of 0.001).
#'
#' Loci scored as missing in the unknown are removed, and the remaining
#' missing values are imputed (method "neighbour" of \code{gl.impute()})
#' before the PCA. Populations whose covariance matrix in the two dimensions
#' is singular, or that have fewer than three individuals, cannot be tested;
#' they are reported and retained.
#'
#' @param x Name of the input genlight object [required].
#' @param unknown Identity label of the focal individual whose provenance is
#' unknown [required].
#' @param nmin Minimum sample size for a target population to be included in
#' the analysis [default 10].
#' @param plevel Alpha level for the bounding ellipses in the PCA plot
#' [default 0.001].
#' @param plot.out If TRUE, plot the 2D PCA showing the position of the
#' unknown [default TRUE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @return A genlight object containing the unknown individual (assigned to
#' population "unknown") and the populations that are not eliminated, with
#' their original (not imputed) genotypes. Loci scored as missing in the
#' unknown are removed.
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("SIBER", quietly = TRUE)) {
#' # Focal individual from the Macleay River (EmmacMaclGeor)
#' test <- gl.assign.pca(testset.gl, unknown = "UC_00146", verbose = 3)
#' }
#' }
#'
#' @seealso \code{\link{gl.assign.pa}}, \code{\link{gl.assign.on.genotype}},
#' \code{\link{gl.assign.mahalanobis}}
#'
#' @importFrom stats cov
#' @export

gl.assign.pca <- function(x,
                          unknown,
                          nmin = 10,
                          plevel = 0.001,
                          plot.out = TRUE,
                          verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  if (verbose == 0) {
    plot.out <- FALSE
  }

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  pkg <- "SIBER"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop(error(
      "Package", pkg, "needed for this function to work. Please install it.\n"
    ))
  }

  if (any(duplicated(indNames(x)))) {
    stop(error("Fatal Error: Duplicate individual names in genlight object.\n"))
  }
  if (any(is.na(indNames(x)))) {
    stop(error("Fatal Error: NA values found in individual names.\n"))
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
  if (is.null(pop(x))) {
    stop(error("Fatal Error: Population assignments are NULL.\n"))
  }
  if (any(is.na(pop(x)))) {
    stop(error("Fatal Error: NA values found in population assignments.\n"))
  }
  if (nPop(x) < 2) {
    stop(error(
      "Fatal Error: Only one population, including the unknown, no putative",
      "source\n"
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
  if (plevel > 0.001 && verbose >= 1) {
    cat(warn(
      "  Warning: Value of plevel greater than 0.001 may result in",
      "unacceptable false elimination of putative source populations\n"
    ))
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

  hard.min <- 10
  if (nmin < hard.min && verbose >= 1) {
    cat(warn(
      "  Warning: The specified minimum sample size is less than", hard.min,
      "individuals. Confidence ellipses of small populations are poorly",
      "estimated.\n"
    ))
  }

  # DO THE JOB

  plevel <- 1 - plevel

  # Assign the unknown individual to its own population
  tmp <- as.character(pop(x))
  tmp[indNames(x) == unknown] <- "unknown"
  pop(x) <- as.factor(tmp)

  # Remove all known populations with fewer than nmin individuals (the
  # unknown is not counted)
  pop.sizes <- table(pop(x))
  pop.keep <- setdiff(names(pop.sizes)[pop.sizes >= nmin], "unknown")
  pop.toss <- setdiff(names(pop.sizes)[pop.sizes < nmin], "unknown")
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
  x <- gl.keep.pop(x, pop.list = c(pop.keep, "unknown"), verbose = 0)

  # Remove loci scored as NA for the unknown. Index by position, because
  # data.frame() would rewrite hyphens in locus names.
  na.loc <- which(is.na(as.matrix(x[indNames(x) == unknown, ])[1, ]))
  if (length(na.loc) > 0) {
    x <- gl.drop.loc(x, loc.list = locNames(x)[na.loc], verbose = 0)
  }

  # Keep the observed genotypes for the returned object
  x.obs <- x

  # Impute remaining missing values for the PCA
  x <- gl.impute(x, method = "neighbour", verbose = 0)

  if (verbose >= 2) {
    cat(report(
      "  Calculating a PCA to represent the unknown in the context of",
      "putative sources\n"
    ))
  }
  pcoa <- gl.pcoa(x, nfactors = 2, verbose = 0)
  if (any(!is.finite(pcoa$scores))) {
    stop(error(
      "Fatal Error: PCA ordination failed. Check for missing data or",
      "insufficient variation.\n"
    ))
  }

  if (plot.out) {
    # The single-individual "unknown" population cannot have an ellipse
    suppressWarnings(suppressMessages(gl.pcoa.plot(
      pcoa, x, ellipse = TRUE, plevel = plevel, verbose = 0
    )))
  }

  # Determine whether the unknown lies within the confidence ellipse of each
  # population
  scores <- pcoa$scores[, 1:2, drop = FALSE]
  pops <- as.character(pop(x))
  unknown.score <- scores[pops == "unknown", , drop = FALSE]
  result <- data.frame(pop = pop.keep, hit = NA, stringsAsFactors = FALSE)

  for (i in seq_along(pop.keep)) {
    A <- scores[pops == pop.keep[i], , drop = FALSE]
    if (nrow(A) < 3) {
      next
    }
    mu <- colMeans(A)
    sigma <- stats::cov(A)
    # A singular covariance has no ellipse; the population cannot be tested
    if (any(!is.finite(sigma)) || rcond(sigma) < 1e-10) {
      next
    }
    testset <- rbind(unknown.score, A)
    transform <- SIBER::pointsToEllipsoid(testset, sigma, mu)
    inside.or.out <- SIBER::ellipseInOut(transform, p = plevel)
    result$hit[i] <- inside.or.out[1]
  }

  untested <- result$pop[is.na(result$hit)]
  eliminated <- result$pop[result$hit %in% FALSE]
  retained <- result$pop[!(result$hit %in% FALSE)]

  if (length(untested) > 0 && verbose >= 1) {
    cat(warn(
      "  Warning: populations not tested (fewer than three individuals or a",
      "singular covariance matrix), retained:",
      paste(untested, collapse = ", "), "\n"
    ))
  }
  if (verbose >= 2) {
    cat(report(
      "  Eliminating populations for which the unknown is outside their",
      "confidence envelope\n"
    ))
  }
  if (verbose >= 3) {
    if (any(result$hit %in% TRUE)) {
      cat("  Putative source populations:",
          paste(result$pop[result$hit %in% TRUE], collapse = ", "), "\n")
    } else {
      cat("  No putative source populations identified\n")
    }
    if (length(eliminated) == 0) {
      cat("  No populations eliminated from consideration\n")
    } else {
      cat("  Populations eliminated from consideration:",
          paste(eliminated, collapse = ", "), "\n")
    }
  }

  # Return the observed genotypes of the retained populations and the unknown
  x2 <- gl.keep.pop(x.obs, pop.list = c(retained, "unknown"), verbose = 0)

  if (verbose >= 2) {
    cat(report(
      "  Returning a genlight object with remaining putative source",
      "populations plus the unknown\n"
    ))
  }
  if (nPop(x2) <= 1 && verbose >= 1) {
    cat(warn("  Warning: No putative source populations remain after",
             "filtering.\n"))
  }

  # ADD TO HISTORY
  nh <- length(x2@other$history)
  x2@other$history[[nh + 1]] <- match.call()

  # FLAG SCRIPT END
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }

  return(x2)
}
