#' @name gl.map.snmf
#' @title Maps a snmf plot using a genlight object
#' @family population structure
#' @description
#' This function takes the Q matrix of a snmf run (as returned by
#' \code{\link{gl.plot.snmf}}) and maps it using the population centres from
#' the genlight object that was used to run the snmf analysis via
#' \code{\link{gl.run.snmf}}. It plots the typical snmf bar plots on a spatial
#' map, providing a barplot for each subpopulation. Therefore it requires
#' coordinates from a genlight object. This kind of plot should support the
#' interpretation of the spatial structure of a population, but in principle
#' is not different from \code{\link{gl.plot.snmf}}.
#' @param x Name of the genlight object containing the coordinates in the
#'  \code{\@other$latlon} slot (columns named lon and lat) to calculate the
#'  population centres [required].
#' @param qmat Q-matrix returned by \code{\link{gl.plot.snmf}} (or an element
#'  of the matrix list returned by \code{\link{gl.run.snmf}}) [required].
#' @param color_clusters A colour palette function (for example
#'  \code{rainbow}), which is called with the number of clusters, or a vector
#'  with at least as many colours as clusters [default NULL, which uses
#'  gl.select.colors()].
#' @param provider Provider passed to leaflet. Check \link[leaflet]{providers}
#' for a list of possible backgrounds [default "Esri.NatGeoWorldMap"].
#' @param scalex Scaling factor to determine the size of the bars in x direction
#' [default 1].
#' @param scaley Scaling factor to determine the size of the bars in y direction
#'  [default 1].
#' @param movepops A data frame with two columns, lon and lat, that moves the
#' centre of the barplots manually in case they overlap. Rows are matched to
#' populations by row name when the row names are population names;
#' otherwise there must be one row per population of x, in the order of
#' levels(pop(x)) (see example) [default NULL].
#' @param pop.labels Switch for population labels below the barplots
#' [default TRUE].
#' @param pop.labels.cex Size of population labels [default 12].
#' @param plot.out Specify if the map is to be printed [default TRUE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  brief progress messages; 3, progress and results summary; 5, full report
#'  [default 2, unless specified using gl.set.verbosity].
#' @return A list with the map (element map) and the Q matrix split into
#'  tables per population, in the order the bars are drawn (element Q_name,
#'  with the population in column Pop_name). This can be used to create your
#'  own map.
#' @author Author(s): Ching Ching Lau. Custodian: Ching Ching Lau -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @details
#' Creates a mapped version of snmf plots with
#' \code{\link{gl.map.structure}}, which places each population's bars at its
#' own centre. Individuals are matched to x by name (column Label); only the
#' populations present in qmat are mapped. For possible background maps
#' check as specified via the provider:
#' \url{http://leaflet-extras.github.io/leaflet-providers/preview/index.html}.
#' You may need to adjust scalex and scaley values [default 1], as the size
#' depends on the scale of the map and the position of the populations.
#' @examples
#' \dontrun{
#' m <- gl.run.snmf(x = bandicoot.gl, minK = 1, maxK = 5, rep = 10)
#' Q <- gl.plot.snmf(snmf.result = m, plot.K = 3, ind.name = TRUE)
#' gl.map.snmf(bandicoot.gl, qmat = Q)
#' # move population 4 (out of 5) 0.5 degrees to the right and population 1
#' # 0.3 degrees to the south of the map.
#' mp <- data.frame(lon = c(0, 0, 0, 0.5, 0), lat = c(-0.3, 0, 0, 0, 0))
#' gl.map.snmf(bandicoot.gl, qmat = Q, movepops = mp)
#' }
#' @export
#' @seealso \code{\link{gl.run.snmf}}, \code{\link{gl.plot.snmf}},
#' \code{\link{gl.map.structure}}
#' @references
#' \itemize{
#' \item Frichot E, Mathieu F, Trouillon T, Bouchard G, Francois O. (2014). Fast and Efficient 
#' Estimation of Individual Ancestry Coefficients. Genetics, 194(4): 973--983.
#' 
#' }

gl.map.snmf <- function(x,
                        qmat,
                        color_clusters = NULL,
                        provider = "Esri.NatGeoWorldMap",
                        scalex = 1,
                        scaley = 1,
                        movepops = NULL,
                        pop.labels = TRUE,
                        pop.labels.cex = 12,
                        plot.out = TRUE,
                        verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  if (!is.data.frame(qmat) || !"Label" %in% names(qmat) ||
      !any(grepl("^Pop_[0-9]+$", names(qmat)))) {
    stop(error(
      "qmat must be the Q matrix returned by gl.plot.snmf (columns Pop_1,",
      "..., Label).\n"
    ))
  }

  # DO THE JOB

  q_cols <- grep("^Pop_[0-9]+$", names(qmat), value = TRUE)
  labels <- as.character(qmat$Label)
  pops <- as.character(pop(x))[match(labels, indNames(x))]

  if (all(is.na(pops))) {
    stop(error(
      "None of the individuals in qmat (column Label) is in x.\n"
    ))
  }
  if (anyNA(pops) && verbose >= 1) {
    cat(warn(
      "  Warning:", sum(is.na(pops)), "individual(s) in qmat are not in x",
      "and are not mapped.\n"
    ))
  }

  if (verbose >= 2) {
    extra <- setdiff(levels(pop(x)), pops)
    if (length(extra) > 0) {
      cat(report(
        "  Populations in x that are not in qmat are not mapped:",
        paste(extra, collapse = ", "), "\n"
      ))
    }
  }

  # the layout gl.map.structure takes (as returned by gl.plot.structure)
  q_struct <- data.frame(Label = labels, stringsAsFactors = FALSE)
  for (j in seq_along(q_cols)) {
    q_struct[[paste0("cluster", j)]] <- qmat[[q_cols[j]]]
  }
  q_struct$K <- as.character(length(q_cols))
  q_struct$orig.pop <- pops
  q_struct$ord <- seq_len(nrow(q_struct))

  res <- gl.map.structure(list("1" = q_struct),
                          x = x,
                          K = length(q_cols),
                          provider = provider,
                          scalex = scalex,
                          scaley = scaley,
                          movepops = movepops,
                          pop.labels = pop.labels,
                          pop.labels.cex = pop.labels.cex,
                          plot.colors = color_clusters,
                          plot.out = plot.out,
                          verbose = 0)

  # per-population tables of the original rows, in the order of the bars
  qmat$Pop_name <- pops
  Q_name <- lapply(res$qmats, function(tab) {
    qmat[match(tab$Label, labels), , drop = FALSE]
  })

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(list(Q_name = Q_name, map = res$map))
}
