#' @name gl.map.popcluster
#' @title Maps a PopCluster plot using a genlight object
#' @family population structure
#' @description
#' This function takes the Q matrix of a PopCluster run (as returned by
#' \code{\link{gl.plot.popcluster}}) and maps it using the population centres
#' from the genlight object that was used to run PopCluster via
#' \code{\link{gl.run.popcluster}}. It plots the typical structure bar plots
#' on a spatial map, providing a barplot for each subpopulation. Therefore it
#' requires coordinates from a genlight object. This kind of plot should
#' support the interpretation of the spatial structure of a population, but
#' in principle is not different from \code{\link{gl.plot.popcluster}}.
#' @param x Name of the genlight object containing the coordinates in the
#'  \code{\@other$latlon} slot (columns named lon and lat) to calculate the
#'  population centres [required].
#' @param qmat Q-matrix returned by \code{\link{gl.plot.popcluster}}
#'  [required].
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
#'  with the population in column Pop_name).
#' @author Author(s): Ching Ching Lau. Custodian: Ching Ching Lau -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @details
#' The map is drawn by \code{\link{gl.map.snmf}} (and so by
#' \code{\link{gl.map.structure}}), which places each population's bars at
#' its own centre; see there for details. For possible background maps
#' check as specified via the provider:
#' \url{http://leaflet-extras.github.io/leaflet-providers/preview/index.html}.
#' You may need to adjust scalex and scaley values [default 1], as the size
#' depends on the scale of the map and the position of the populations.
#' @examples
#' \dontrun{
#' m <- gl.run.popcluster(x = bandicoot.gl,
#'   popcluster.path = "/User/PopCluster/Bin/", minK = 1, maxK = 3, rep = 2)
#' Q <- gl.plot.popcluster(pop_cluster_result = m, plot.K = 3, ind_name = TRUE)
#' gl.map.popcluster(x = bandicoot.gl, qmat = Q)
#' # move population 4 (out of 5) 0.5 degrees to the right and population 1
#' # 0.3 degrees to the south of the map.
#' mp <- data.frame(lon = c(0, 0, 0, 0.5, 0), lat = c(-0.3, 0, 0, 0, 0))
#' gl.map.popcluster(bandicoot.gl, qmat = Q, movepops = mp)
#' }
#' @export
#' @importFrom dplyr left_join
#' @seealso \code{\link{gl.run.popcluster}}, \code{\link{gl.plot.popcluster}}
#' @references
#' \itemize{
#' \item Wang, J. (2022). Fast and accurate population admixture inference
#' from genotype data from a few microsatellites to millions of SNPs.
#' Heredity, 129(2), 79-92.
#' }

gl.map.popcluster <- function(x,
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

  # DO THE JOB

  # the PopCluster Q matrix has the Pop_1..K and Label columns gl.map.snmf
  # reads
  res <- gl.map.snmf(x,
                     qmat = qmat,
                     color_clusters = color_clusters,
                     provider = provider,
                     scalex = scalex,
                     scaley = scaley,
                     movepops = movepops,
                     pop.labels = pop.labels,
                     pop.labels.cex = pop.labels.cex,
                     plot.out = plot.out,
                     verbose = verbose)

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(res)
}
