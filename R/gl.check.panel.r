#' @name gl.check.panel
#'
#' @title Check how well a SNP panel reproduces a population parameter
#'
#' @family panel selection
#'
#' @description
#' This function checks how well a panel of loci (for example from
#' \code{\link{gl.select.panel}}) reproduces a parameter of conservation
#' concern (Fst, He, Ho, Fis, allelic richness or Ne) estimated from the full
#' data, by computing the parameter per population (per population pair for
#' Fst) in both data sets and plotting panel against full data with a linear
#' regression.
#'
#' @param x A genlight object with the SNP panel [required].
#' @param xorig A genlight object with the full SNP data; it must hold the
#'   same individuals (by name) and populations as \code{x} [required].
#' @param parameter The parameter to check: "Fst" (pairwise, from
#'   gl.fst.pop), "He", "Ho", "Fis" (from gl.report.heterozygosity), "Na"
#'   (mean corrected allelic richness, from gl.report.allelerich; "Nall" is
#'   accepted as an alias) or "Ne" (from gl.LDNe, lowest allele frequency
#'   0.05) [default "Fst"].
#' @param neest.path Path to the folder with the NeEstimator executable (see
#'   gl.LDNe); needed only for "Ne" [default NULL].
#' @param plot.out Logical. If TRUE, the plot is printed [default TRUE].
#' @param plot.file Name for the RDS binary file to save the plot (base name
#'   only, exclude extension). If NULL, the plot is not saved [default NULL].
#' @param plot.dir Directory in which to save the plot [default = tempdir(),
#'   unless specified using gl.set.wd].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'   brief progress messages; 3, progress and results summary; 5, full report
#'   [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' Each point in the plot is one population (one population pair for Fst);
#' the line is a least-squares fit with its equation and R2. Populations for
#' which the parameter cannot be estimated (NA, or infinite Ne) are dropped.
#'
#' @return A data frame with two columns, the parameter from the full data
#'   and from the panel, one row per population (per population pair for
#'   Fst).
#'
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#'   \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # Select 50 loci randomly
#' selected <- gl.select.panel(possums.gl, method = "random", nl = 50)
#' gl.check.panel(selected, possums.gl, parameter = "Fst")
#'
#' @export
#' @importFrom ggpmisc stat_poly_eq use_label

gl.check.panel <- function(x, xorig, parameter="Fst", neest.path = NULL,
                           plot.out = TRUE, plot.file = NULL, plot.dir = NULL,
                           verbose = NULL) {

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = 0)
  datatype <- utils.check.datatype(xorig, accept = "SNP", verbose = 0)

  # FUNCTION SPECIFIC ERROR CHECKING
  parameters <- c("Fst", "He", "Ho", "Fis", "Na", "Ne")
  if (!is.character(parameter) || length(parameter) != 1) {
    stop(error("'parameter' must be one of:",
               paste(parameters, collapse = ", "), "\n"))
  }
  if (parameter == "Nall") parameter <- "Na"
  if (!parameter %in% parameters) {
    stop(error("'parameter' must be one of:",
               paste(parameters, collapse = ", "), "\n"))
  }

  # The panel and the full data must hold the same individuals; align the
  # panel to xorig by individual name.
  if (!setequal(indNames(x), indNames(xorig)) ||
      anyDuplicated(indNames(x)) > 0) {
    stop(error("'x' and 'xorig' must hold the same individuals (by name).",
               "Please ensure that the panel is a subset of loci of the",
               "original data.\n"))
  }
  x <- x[match(indNames(xorig), indNames(x)), ]
  if (!identical(as.character(pop(x)), as.character(pop(xorig)))) {
    stop(error("Individuals are assigned to different populations in 'x'",
               "and 'xorig'.\n"))
  }

  # DO THE JOB
  if (verbose >= 2) {
    cat(report("  Computing", parameter, "for the full data and the panel\n"))
  }

#check fsts

if (parameter == "Fst") {


fst_x <- as.numeric(gl.fst.pop(xorig, verbose = 0, nboots=1))

fst_test <- as.numeric(gl.fst.pop(x, verbose = 0, nboots = 1))

fsts <- data.frame(fst_orig=fst_x, fst_panel=fst_test)
fsts <- fsts[complete.cases(fsts),]

res <- fsts
}
#check expected heterozygosity

if (parameter=="He") {
het_x <- gl.report.heterozygosity(xorig, verbose=0)
het_test <- gl.report.heterozygosity(x, verbose=0)
hets <- data.frame(het_orig=het_x$He, het_panel=het_test$He)
res <- hets
}

#check number of alleles
if (parameter=="Na") {
nall_x <- gl.report.allelerich(xorig, verbose=0,plot.display = F)
nall_test <-  gl.report.allelerich(x, verbose=0,plot.display = F)
nalls <- data.frame(nall_orig=nall_x$`Allelic Richness per population`$`mean_corrected_richness`, nall_panel= nall_test$`Allelic Richness per population`$`mean_corrected_richness`)

res <- nalls
}


#check FIS
if (parameter=="Fis") {


fis_x <- gl.report.heterozygosity(xorig, verbose=0)
fis_test <- gl.report.heterozygosity(x, verbose=0)

fiss <- data.frame(fis_orig=fis_x$FIS, fis_panel= fis_test$FIS)

res <- fiss
}


#check observed heterozygosity
if (parameter=="Ho") {
hetso_x <- gl.report.heterozygosity(xorig, verbose=0)
hetso_test <- gl.report.heterozygosity(x, verbose=0)

hetso <- data.frame(hetso_orig=hetso_x$Ho, hetso_panel =hetso_test$Ho)

res <- hetso
}

#check Ne
if (parameter=="Ne") {
ne_x <- gl.LDNe(xorig, neest.path = neest.path, critical = 0.05, verbose=0,mating = "random")

ne_test <- gl.LDNe(x, neest.path = neest.path, critical = c(0.05), verbose=0,mating = "random", singleton.rm = F)

# Ne estimate at the lowest allele frequency (first column), by row label
ne_at <- function(d) {
  d[trimws(d$Statistic) == "Estimated Ne^", "Frequency 1"]
}
nes_x <- suppressWarnings(as.numeric(unlist(lapply(ne_x, ne_at))))

nes_test <- suppressWarnings(as.numeric(unlist(lapply(ne_test, ne_at))))

nes<- data.frame(nes_orig=nes_x,  nes_panel = nes_test)

nes[sapply(nes, function(x) x=="Inf")] <- NA  # replace Inf with NA]

res <- nes
}

  # PLOT: panel against full data, with the regression line
  xcol <- names(res)[1]
  ycol <- names(res)[2]
  gg <- ggplot(res, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x) +
    stat_poly_eq(use_label(c("eq", "R2")), formula = y ~ x)

  if (verbose >= 3) {
    ok <- stats::complete.cases(res)
    if (sum(ok) >= 3) {
      cat(report("  Correlation of panel with full data (", sum(ok),
                 "points):", signif(stats::cor(res[ok, 1], res[ok, 2]), 3),
                 "\n"))
    }
  }

  if (plot.out) {
    suppressWarnings(print(gg))
  }

  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(gg,
      dir = plot.dir,
      file = plot.file,
      verbose = verbose
    )
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(res)
}
