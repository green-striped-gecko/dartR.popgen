#' @name gl.ld.distance
#' @title Plots linkage disequilibrium against distance by population
#' @description
#' The function bins SNP pairs by the distance between them, in base pairs
#' pooled over all chromosomes, and plots the mean pairwise LD of each bin
#' for each population, with a red line at R.squared = 0.2, the threshold
#' commonly used to imply that two loci are unlinked (Delourme et al., 2013;
#' Li et al., 2014).
#' @details
#' Bins are (1, 1 + ld.resolution], (1 + ld.resolution, 1 + 2 *
#' ld.resolution], ..., with the first bin closed on the left and the last
#' bin ending at the largest distance in \code{ld.report}. Each point is the
#' mean LD of the pairs in a bin and is plotted at the upper edge of the bin.
#'
#' The red threshold line applies to R.squared, the default LD measure of
#' \code{\link[dartR.base]{gl.report.ld.map}}. If the report was computed
#' with another measure (e.g. D.prime), the line has no meaning.
#' @param ld.report Output from function \code{\link[dartR.base]{gl.report.ld.map}}
#' [required].
#' @param ld.resolution Width of the distance bins in number of base pairs
#' [default 100000].
#' @param pop.colors A palette function (e.g. \code{rainbow}) or a vector
#'  with at least as many colors as there are populations in the dataset
#' [default NULL].
#' @param plot.title Title of the plot [default " "].
#' @param plot.theme User specified theme [default NULL].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot.dir Directory in which to save files [default = tempdir(),
#'   unless specified using gl.set.wd].
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @references
#' \itemize{
#' \item Delourme, R., Falentin, C., Fomeju, B. F., Boillot, M., Lassalle, G.,
#' Andre, I., . . . Marty, A. (2013). High-density SNP-based genetic map
#' development and linkage disequilibrium assessment in Brassica napusL. BMC
#' genomics, 14(1), 120.
#' \item Li, X., Han, Y., Wei, Y., Acharya, A., Farmer, A. D., Ho, J., . . .
#' Brummer, E. C. (2014). Development of an alfalfa SNP array and its use to
#' evaluate patterns of population structure and linkage disequilibrium. PLoS
#' One, 9(1), e84329.
#'  }
#' @return A data.table with one row per population and bin: pop, distance
#'   (upper edge of the bin in base pairs), ld.stat (mean LD of the pairs in
#'   the bin; NA for an empty bin) and n.pairs (number of pairs with an LD
#'   value in the bin).
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # gl.report.ld.map (dartR.base) needs snpStats and fields
#' if (requireNamespace("snpStats", quietly = TRUE) &&
#'     requireNamespace("fields", quietly = TRUE)) {
#'   require("dartR.data")
#'   x <- platypus.gl
#'   x <- gl.filter.callrate(x, threshold = 1)
#'   x <- gl.filter.monomorphs(x)
#'   x$position <- x$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
#'   x$chromosome <- as.factor(x$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1)
#'   ld_res <- gl.report.ld.map(x, ld.max.pairwise = 10000000)
#'   ld_res_2 <- gl.ld.distance(ld_res, ld.resolution = 1000000)
#' }
#' @family ld functions
#' @export

gl.ld.distance <- function(ld.report,
                           ld.resolution = 100000,
                           pop.colors = NULL,
                           plot.title = " ",
                           plot.theme = NULL,
                           plot.out = TRUE,
                           plot.file = NULL,
                           plot.dir = NULL,
                           verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  if (!is.data.frame(ld.report) ||
      !all(c("pop", "distance", "ld.stat") %in% names(ld.report))) {
    stop(error("'ld.report' must be the output of gl.report.ld.map, with",
               "columns pop, distance and ld.stat.\n"))
  }
  if (nrow(ld.report) == 0) {
    stop(error("'ld.report' has no SNP pairs.\n"))
  }
  if (!is.numeric(ld.resolution) || length(ld.resolution) != 1 ||
      is.na(ld.resolution) || ld.resolution <= 0) {
    stop(error("'ld.resolution' must be a single positive number of base",
               "pairs.\n"))
  }

  # DO THE JOB
  ld_max_pairwise <- max(ld.report$distance, na.rm = TRUE)
  # unique() removes the repeated last break when the sequence already ends
  # on the largest distance
  break_bins <- unique(c(seq(1, ld_max_pairwise, ld.resolution),
                         ld_max_pairwise))
  if (length(break_bins) < 2) {
    break_bins <- c(0, ld_max_pairwise)
  }

  # Mean LD and number of pairs per bin: bins (a, b], first bin closed
  bin_pop <- function(x) {
    bins <- cut(x$distance, breaks = break_bins, include.lowest = TRUE)
    ok <- !is.na(x$ld.stat)
    n <- as.vector(table(bins[ok]))
    m <- as.vector(tapply(x$ld.stat[ok], bins[ok], mean))
    data.table::data.table(distance = break_bins[-1], ld.stat = m, n.pairs = n)
  }

  split_df <- split(ld.report, f = ld.report$pop, drop = TRUE)
  bins_ld <- data.table::rbindlist(lapply(split_df, bin_pop), idcol = "pop")
  bins_ld$pop <- as.factor(bins_ld$pop)

  # pairwise LD by population
  npops <- nlevels(bins_ld$pop)

  if (is.null(plot.theme)) {
    plot.theme <- theme_dartR()
  }

  if (is.null(pop.colors)) {
    pop.colors <-
      gl.select.colors(
        library = "gr.palette", palette = "Alphabet",
        ncolors = npops, verbose = 0
      )
  } else if (is.function(pop.colors)) {
    pop.colors <- pop.colors(npops)
  }
  if (length(pop.colors) < npops) {
    stop(error(paste0("'pop.colors' has ", length(pop.colors),
                      " colors but there are ", npops, " populations.\n")))
  }

  distance <- ld.stat <- NULL
  p3 <-
    ggplot(bins_ld, aes(x = distance, y = ld.stat, colour = pop)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(
      aes(yintercept = 0.2, linetype = "R.squared = 0.2 (unlinked threshold)"),
      colour = "red",
      linewidth = 1
    ) +
    scale_linetype_manual(name = NULL, values = "solid") +
    xlab("Base pairs") +
    ylab("Linkage disequilibrium") +
    labs(color = "") +
    scale_color_manual(values = pop.colors) +
    plot.theme +
    theme(legend.position = "bottom") +
    ggtitle(plot.title)


  # PRINTING OUTPUTS
  if (plot.out) {
    print(p3)
  }
  if (verbose >= 3) {
    print(bins_ld, row.names = FALSE)
  }

  # Optionally save the plot ---------------------

  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p3,
      dir = plot.dir,
      file = plot.file,
      verbose = verbose
    )
  }


  # FLAG SCRIPT END

  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN

  return(invisible(bins_ld))
}
