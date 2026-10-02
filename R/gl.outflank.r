#' @name gl.outflank
#'
#' @title Identify Fst outlier loci with the OutFLANK method of Whitlock and
#' Lotterhos (2015)
#'
#' @family selection
#'
#' @description
#' Runs OutFLANK on a SNP genlight (or genind) object: a single Fst outlier
#' scan across all populations that flags loci whose Fst is higher (or
#' lower) than expected under a neutral chi-square distribution inferred from
#' the trimmed centre of the Fst distribution.
#'
#' @details
#' The code of OutFLANK (Whitlock and Lotterhos) is bundled in dartR.popgen,
#' so the OutFLANK package itself is not needed; the Bioconductor package
#' qvalue is. Results are identical to running OutFLANK::MakeDiploidFSTMat()
#' and OutFLANK::OutFLANK() on the same genotypes.
#'
#' Each SNP enters the analysis once, as the number of copies of the
#' alternate allele (0, 1, 2) per individual. Every input locus is kept in
#' the output, in input order. Loci that are monomorphic or missing in every
#' individual are not tested and get NA. Loci with expected heterozygosity
#' below \code{Hmin} are not tested either, but OutFLANK reports them with
#' OutlierFlag FALSE (GoodH = "lowH" in the results), so their index is
#' TRUE.
#'
#' Negative Fst estimates cannot be tested against the chi-square
#' distribution; they are flagged as outliers only when the lowest positive
#' Fst values are outliers. Low-Fst outliers are not reliable (Whitlock and
#' Lotterhos 2015).
#'
#' @param gi A genlight object with SNP data, or a genind object, with at
#'   least two populations defined [required].
#' @param plot Logical: plot the histogram of Fst (uncorrected for sample
#'   size) with the inferred neutral distribution [default TRUE].
#' @param LeftTrimFraction The proportion of loci that are trimmed from the
#'   lower end of the range of Fst before the likelihood function is applied
#'   [default 0.05].
#' @param RightTrimFraction The proportion of loci that are trimmed from the
#'   upper end of the range of Fst before the likelihood function is applied
#'   [default 0.05].
#' @param Hmin The minimum expected heterozygosity required for a locus to be
#'   tested [default 0.1].
#' @param qthreshold The false discovery rate threshold used to flag outliers
#'   from the q-values [default 0.05].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'   brief progress messages; 3, progress and results summary; 5, full report
#'   [default 2, unless specified using gl.set.verbosity].
#'
#' @return A list with two elements:
#' \itemize{
#' \item index: a logical vector with one value per input locus, in input
#'   order: TRUE for loci that are NOT outliers (including loci below
#'   \code{Hmin}), FALSE for outliers, NA for monomorphic loci and loci
#'   without calls. Outliers are selected with \code{which(!res$index)}.
#' \item outflank: the OutFLANK output, a list with FSTbar (mean Fst of
#'   non-outlier loci), FSTNoCorrbar (the same without sample size
#'   correction), dfInferred (inferred degrees of freedom of the neutral
#'   chi-square distribution), numberLowFstOutliers, numberHighFstOutliers,
#'   and results, a data frame with one row per input locus (LocusName, He,
#'   FST, T1, T2, FSTNoCorr, T1NoCorr, T2NoCorr, meanAlleleFreq (frequency of
#'   the reference allele), indexOrder, GoodH, qvalues, pvalues,
#'   pvaluesRightTail, OutlierFlag).
#' }
#'
#' @author Author(s): Bernd Gruber; OutFLANK code by Whitlock and Lotterhos.
#'   Custodian: Bernd Gruber -- Post to
#'   \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("qvalue", quietly = TRUE)) {
#'   res <- gl.outflank(bandicoot.gl, plot = FALSE)
#'   # names of the outlier loci
#'   locNames(bandicoot.gl)[which(!res$index)]
#' }
#' }
#'
#' @references
#' Whitlock, M.C. and Lotterhos K.J. (2015) Reliable detection of loci
#' responsible for local adaptation: inference of a neutral model through
#' trimming the distribution of Fst. The American Naturalist 186: 24 - 36.
#'
#' Github repository: Whitlock & Lotterhos:
#'  \url{https://github.com/whitlock/OutFLANK}
#'
#' @seealso \code{\link{utils.outflank}}, \code{\link{utils.outflank.plotter}},
#'  \code{\link{utils.outflank.MakeDiploidFSTMat}}
#'
#' @importFrom stats optim pgamma quantile
#' @export

gl.outflank <- function(gi,
                        plot = TRUE,
                        LeftTrimFraction = 0.05,
                        RightTrimFraction = 0.05,
                        Hmin = 0.1,
                        qthreshold = 0.05,
                        verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK IF PACKAGES ARE INSTALLED
  pkg <- "qvalue"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop(error(
      "Package ", pkg, " needed for this function to work. Please install it",
      " from Bioconductor (BiocManager::install('qvalue')).\n"
    ))
  }

  # CHECK DATATYPE
  if (is(gi, "genind")) {
    gi <- gi2gl(gi, verbose = 0)
  }
  datatype <- utils.check.datatype(gi, accept = "SNP", verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  if (is.null(pop(gi)) || any(is.na(pop(gi)))) {
    stop(error("Every individual must be assigned to a population.\n"))
  }
  pops <- droplevels(pop(gi))
  if (nlevels(pops) < 2) {
    stop(error(paste0("OutFLANK needs at least two populations; found ",
                      nlevels(pops), ".\n")))
  }

  # DO THE JOB
  # OutFLANK expects one column per SNP with the number of copies of the
  # focal allele (0, 1, 2) and 9 for missing calls.
  snpmat <- as.matrix(gi)
  snpmat[is.na(snpmat)] <- 9
  loc_names <- locNames(gi)

  # Loci missing in every individual cannot be computed; keep them as NA rows
  # so that the output has one row per input locus, in input order.
  has_data <- colSums(snpmat != 9) > 0
  if (verbose >= 2) {
    cat(report("  Calculating Fst for", sum(has_data), "loci\n"))
    if (any(!has_data)) {
      cat(warn("  ", sum(!has_data), "loci have no calls and are not tested\n"))
    }
  }
  mdfm_ok <- utils.outflank.MakeDiploidFSTMat(
    SNPmat = snpmat[, has_data, drop = FALSE],
    locusNames = list(loc_names[has_data]),
    popNames = list(as.character(pops)),
    verbose = verbose
  )
  mdfm <- data.frame(LocusName = loc_names,
                     matrix(NA_real_, nrow = length(loc_names),
                            ncol = ncol(mdfm_ok) - 1,
                            dimnames = list(NULL, names(mdfm_ok)[-1])))
  mdfm[has_data, -1] <- mdfm_ok[, -1]

  # run outflank
  if (verbose >= 2) cat(report("  Fitting the neutral Fst distribution\n"))
  outf <- utils.outflank(
    FstDataFrame = mdfm,
    LeftTrimFraction = LeftTrimFraction,
    RightTrimFraction = RightTrimFraction,
    Hmin = Hmin,
    NumberOfSamples = nlevels(pops),
    qthreshold = qthreshold,
    verbose = verbose
  )

  index.outflank <- !(outf$results$OutlierFlag)

  if (verbose >= 3) {
    n_tested <- sum(!is.na(outf$results$FSTNoCorr) &
                      outf$results$GoodH == "goodH", na.rm = TRUE)
    cat(report("  Loci tested:", n_tested, "of", length(index.outflank),
               paste0("(He >= ", Hmin, ")\n")))
    cat(report("  Outliers: high Fst", outf$numberHighFstOutliers,
               "; low Fst", outf$numberLowFstOutliers, "\n"))
    cat(report("  Inferred df:", signif(outf$dfInferred, 4),
               "; mean Fst of non-outliers:", signif(outf$FSTbar, 4), "\n"))
  }

  if (plot) {
    utils.outflank.plotter(outf, Hmin = Hmin)
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(list(index = index.outflank, outflank = outf))
}
