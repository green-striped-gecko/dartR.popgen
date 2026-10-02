#' @name gl.evanno
#' @title Creates an Evanno plot from a STRUCTURE run object
#' @family population structure
#' @description
#' This function takes a structure run object (output from
#' \code{\link{gl.run.structure}}), computes the Evanno et al. (2005)
#' statistics used to choose the number of clusters K, and plots them.
#' @param sr Structure run object from \code{\link{gl.run.structure}}
#' [required].
#' @param plot.out TRUE: the plots are shown. FALSE: the plots are returned as
#' ggplot objects but not shown [default TRUE].
#' @param plot.theme Theme for the plots. See Details for options
#' [default theme_dartR()].
#' @param plot.dir Directory in which to save files [default = tempdir(),
#'  unless specified using gl.set.wd].
#' @param plot.file Name for the RDS binary file to save the combined plot
#' (base name only, exclude extension) [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  brief progress messages; 3, progress and results summary; 5, full report
#'  [default 2, unless specified using gl.set.verbosity].
#' @details
#' For each K the function takes the mean and standard deviation of the
#' estimated Ln probability of the data, LnP(K), over the replicate runs,
#' and computes LnP'(K) = LnP(K) - LnP(K - 1),
#' |LnP''(K)| = |LnP'(K + 1) - LnP'(K)| and
#' delta K = |LnP''(K)| / sd(LnP(K)). The differences are taken on the mean
#' LnP(K) of each K, as in the strataG package (Archer et al. 2016) from
#' which the code was adapted. The K with the largest delta K is the usual
#' choice.
#'
#' The method needs at least three consecutive values of K and at least two
#' replicates per K. Delta K is NA for the smallest and largest K, for K
#' values whose neighbours K - 1 or K + 1 were not run, for K with a single
#' replicate, and for K whose replicates all report the same LnP(K)
#' (sd = 0); a warning names these K at verbose >= 1.
#'
#' Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#' @return A list with two elements: df, a data frame with one row per K and
#' columns k, reps (number of replicates), mean.ln.k, sd.ln.k, ln.pk
#' (LnP'(K)), ln.ppk (|LnP''(K)|) and delta.k; and plots, a list of ggplot
#' objects mean.ln.k, ln.pk, ln.ppk, delta.k (only when delta K could be
#' computed for at least one K) and combined (all panels in one figure).
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # examples need structure to be installed on the system (see above)
#' \dontrun{
#'  bc <- bandicoot.gl[,1:100]
#'  sr <- gl.run.structure(bc, k.range = 2:5, num.k.rep = 3, exec = './structure.exe')
#'  ev <- gl.evanno(sr)
#'  ev$df
#'  qmat <- gl.plot.structure(sr, K=3)
#'  head(qmat)
#'  gl.map.structure(qmat, bc, K=3, scalex=1, scaley=0.5)
#' }
#' @import patchwork
#' @export
#' @seealso \code{\link{gl.run.structure}}, \code{\link{gl.plot.structure}}
#' @references
#' \itemize{
#' \item Pritchard, J.K., Stephens, M., Donnelly, P. (2000) Inference of
#' population structure using multilocus genotype data. Genetics 155, 945-959.
#' \item Archer, F. I., Adams, P. E. and Schneiders, B. B. (2016) strataG: An R
#' package for manipulating, summarizing and analysing population genetic data.
#' Mol Ecol Resour. doi:10.1111/1755-0998.12559
#' \item Evanno, G., Regnaut, S., and J. Goudet. 2005. Detecting the number of
#'  clusters of individuals using the software STRUCTURE: a simulation study.
#'   Molecular Ecology 14:2611-2620.
#' }

gl.evanno <- function(sr,
                      plot.out = TRUE,
                      plot.theme = theme_dartR(),
                      plot.dir = NULL,
                      plot.file = NULL,
                      verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  if (!is(sr, "structure.result")) {
    stop(error(
      "sr is not a structure.result object returned by gl.run.structure.\n"
    ))
  }

  n_k <- length(unique(sapply(sr, function(r) r$summary[["k"]])))
  if (n_k < 3) {
    stop(error(
      "The Evanno method needs at least three values of K;", n_k,
      "found.\n"
    ))
  }

  # DO THE JOB

  evno <- utils.structure.evanno(sr,
                                 plot = FALSE,
                                 plot.theme = plot.theme,
                                 verbose = verbose)

  if (verbose >= 3 && any(!is.na(evno$df$delta.k))) {
    best <- evno$df$k[which.max(evno$df$delta.k)]
    cat(report("  Largest delta K at K =", best, "\n"))
  }

  # PRINTING OUTPUTS
  if (plot.out) {
    print(evno$plots$combined)
  }

  # Optionally save the plot ---------------------
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(evno$plots$combined,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(evno)
}
