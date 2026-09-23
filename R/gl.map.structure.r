#' @name gl.map.structure
#' @title Maps a STRUCTURE plot using a genlight object
#' @family population structure
#' @description
#' This function takes the output of \code{\link{gl.plot.structure}} (the q
#' matrix) and maps the q-matrix using the population centres from the
#' genlight object that was used to run the structure analysis via
#' \code{\link{gl.run.structure}}. It plots the typical structure bar plots on
#' a spatial map, providing a barplot for each subpopulation. Therefore it
#' requires coordinates from a genlight object. This kind of plot should
#' support the interpretation of the spatial structure of a population, but
#' in principle is not different from \code{\link{gl.plot.structure}}.
#'
#' @param qmat List of q-matrices returned by \code{\link{gl.plot.structure}}
#'  [required].
#' @param x Name of the genlight object containing the coordinates in the
#'  \code{\@other$latlon} slot (columns named lon and lat) to calculate the
#'  population centres [required].
#' @param K The q-matrix to be mapped: either the number of clusters (e.g. 3),
#'  which maps the first q-matrix with that many clusters, or a panel label
#'  from the K column of the q-matrices (e.g. "2.2" to map the second mode of
#'  K = 2) [required].
#' @param provider Provider passed to leaflet. Check \link[leaflet]{providers}
#' for a list of possible backgrounds [default "Esri.NatGeoWorldMap"].
#' @param scalex Scaling factor to determine the size of the bars in x direction
#' [default 1].
#' @param scaley Scaling factor to determine the size of the bars in y direction
#'  [default 1].
#' @param movepops A data frame with two columns, lon and lat, that moves the
#' centre of the barplots manually in case they overlap (often when
#' populations are horizontally close to each other). Each row gives the
#' longitude and latitude units to move one population. Rows are matched to
#' populations by row name when the row names are population names;
#' otherwise there must be one row per population of x, in the order of
#' levels(pop(x)). Columns not named lon and lat are read as lon, lat (see
#' example) [default NULL].
#' @param pop.labels Switch for population labels below the barplots
#' [default TRUE].
#' @param pop.labels.cex Size of population labels [default 12].
#' @param plot.colors A colour palette function (for example \code{rainbow}),
#'  which is called with the number of clusters, or a vector with at least as
#'  many colours as clusters [default NULL, which uses gl.select.colors(), the
#'  default of gl.plot.structure].
#' @param plot.out Specify if the map is to be printed [default TRUE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  brief progress messages; 3, progress and results summary; 5, full report
#'  [default 2, unless specified using gl.set.verbosity].
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @details
#' Creates a mapped version of structure plots. For possible background maps
#' check as specified via the provider:
#' \url{http://leaflet-extras.github.io/leaflet-providers/preview/index.html}.
#' You may need to adjust scalex and scaley values [default 1], as the size
#' depends on the scale of the map and the position of the populations.
#'
#' Only the populations present in qmat are mapped; other populations of x
#' are ignored. Each population in qmat needs coordinates in x.
#' @return A list with the map (element map) and the q-matrix split into
#'  sorted tables per population (element qmats). This can be used to create
#'  your own map.
#' @examples
#' # examples need structure to be installed on the system (see above)
#' \dontrun{
#' bc <- bandicoot.gl[,1:100]
#' sr <- gl.run.structure(bc, k.range = 2:5, num.k.rep = 3, exec = './structure.exe')
#' ev <- gl.evanno(sr)
#' ev
#' qmat <- gl.plot.structure(sr, K = 2:4)
#' head(qmat)
#' gl.map.structure(qmat, bc, K = 3)
#' gl.map.structure(qmat, bc, K = 4)
#' # move population 4 (out of 5) 0.5 degrees to the right and population 1
#' # 0.3 degrees to the south of the map.
#' mp <- data.frame(lon = c(0, 0, 0, 0.5, 0), lat = c(-0.3, 0, 0, 0, 0))
#' gl.map.structure(qmat, bc, K = 4, movepops = mp)
#' }
#' @export
#' @seealso \code{\link{gl.run.structure}},  \code{clumpp},
#' \code{\link{gl.plot.structure}}
#' @references
#' \itemize{
#' \item Pritchard, J.K., Stephens, M., Donnelly, P. (2000) Inference of
#' population structure using multilocus genotype data. Genetics 155, 945-959.
#' \item Archer, F. I., Adams, P. E. and Schneiders, B. B. (2016) strataG: An R
#'  package for manipulating, summarizing and analysing population genetic data.
#'   Mol Ecol Resour. doi:10.1111/1755-0998.12559
#' \item Evanno, G., Regnaut, S., and J. Goudet. 2005. Detecting the number of
#' clusters of individuals using the software STRUCTURE: a simulation study.
#' Molecular Ecology 14:2611-2620.
#' \item Mattias Jakobsson and Noah A. Rosenberg. 2007. CLUMPP: a cluster
#' matching and permutation program for dealing with label switching and
#' multimodality in analysis of population structure. Bioinformatics
#' 23(14):1801-1806. Available at
#' \href{http://web.stanford.edu/group/rosenberglab/clumppDownload.html}{clumpp}
#' }

gl.map.structure <- function(qmat,
                             x,
                             K,
                             provider = "Esri.NatGeoWorldMap",
                             scalex = 1,
                             scaley = 1,
                             movepops = NULL,
                             pop.labels = TRUE,
                             pop.labels.cex = 12,
                             plot.colors = NULL,
                             plot.out = TRUE,
                             verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  if (!(requireNamespace("leaflet", quietly = TRUE))) {
    stop(error(
      "Package leaflet needed for this function to work. Please install it.\n"
    ))
  }

  if (length(K) != 1) {
    stop(error("K must be a single value, e.g. 3 or \"2.2\".\n"))
  }

  is_qmat <- function(y) {
    is.data.frame(y) && all(c("K", "orig.pop") %in% names(y)) &&
      any(grepl("^cluster", names(y)))
  }
  if (!is.list(qmat) || length(qmat) == 0 ||
      !all(vapply(qmat, is_qmat, logical(1)))) {
    stop(error(
      "qmat must be the list of q-matrices returned by gl.plot.structure.\n"
    ))
  }

  latlon <- x@other$latlon
  if (is.null(latlon) || !all(c("lon", "lat") %in% colnames(latlon))) {
    stop(error(
      "x@other$latlon must hold the coordinates in columns named lon and",
      "lat.\n"
    ))
  }

  # choose the q-matrix: by panel label when K is character, otherwise by
  # its number of clusters
  k_label <- vapply(qmat, function(y) {
    as.character(unique(y$K))[1]
  }, character(1))
  n_clusters <- vapply(qmat, function(y) {
    sum(grepl("^cluster", names(y)))
  }, integer(1))
  if (is.character(K)) {
    eq.k <- k_label == K
  } else {
    eq.k <- n_clusters == K
  }

  if (sum(eq.k) == 0) {
    stop(error("No entries for K =", K, "found in 'qmat'.\n"))
  }

  if (sum(eq.k) > 1 && verbose >= 2) {
    cat(warn(
      "  Warning: K =", K, "matches", sum(eq.k), "modes (",
      paste(k_label[eq.k], collapse = ", "), "); mapping",
      k_label[eq.k][1], ". Set K to one of these labels to map another",
      "mode.\n"
    ))
  }

  qmat <- as.data.frame(qmat[eq.k][[1]])
  qmat$orig.pop <- factor(qmat$orig.pop)
  pops <- levels(qmat$orig.pop)

  nk <- sum(grepl("^cluster", colnames(qmat)))
  # a palette function is called with the number of clusters
  if (is.function(plot.colors)) {
    plot.colors <- plot.colors(nk)
  }
  if (is.null(plot.colors)) {
    plot.colors <- gl.select.colors(ncolors = nk, verbose = 0)
  }
  if (length(plot.colors) < nk) {
    stop(error(
      "plot.colors has", length(plot.colors), "colours but", nk,
      "are needed (one per cluster).\n"
    ))
  }

  # DO THE JOB

  # population centres, one row per level of pop(x)
  centers <- cbind(
    lon = tapply(latlon[, "lon"], pop(x), mean, na.rm = TRUE),
    lat = tapply(latlon[, "lat"], pop(x), mean, na.rm = TRUE)
  )

  if (!is.null(movepops)) {
    movepops <- as.data.frame(movepops)
    if (ncol(movepops) != 2) {
      stop(error("movepops needs two columns, lon and lat.\n"))
    }
    if (!all(c("lon", "lat") %in% names(movepops))) {
      names(movepops) <- c("lon", "lat")
    }
    # rows matched by population name when named, otherwise by position
    named_rows <- all(rownames(movepops) %in% rownames(centers))
    if (named_rows) {
      idx <- match(rownames(movepops), rownames(centers))
    } else {
      if (nrow(movepops) != nrow(centers)) {
        stop(error(
          "movepops needs one row per population of x (", nrow(centers),
          "), in the order of levels(pop(x)), or row names that are",
          "population names.\n"
        ))
      }
      idx <- seq_len(nrow(centers))
    }
    centers[idx, "lon"] <- centers[idx, "lon"] + movepops$lon
    centers[idx, "lat"] <- centers[idx, "lat"] + movepops$lat
  }

  # only the populations in qmat are mapped, in the order of the bars
  has_coords <- rownames(centers)[stats::complete.cases(centers)]
  missing_pops <- setdiff(pops, has_coords)
  if (length(missing_pops) > 0) {
    stop(error(
      "No coordinates in x for these populations of qmat:",
      paste(missing_pops, collapse = ", "), "\n"
    ))
  }
  extra_pops <- setdiff(rownames(centers), pops)
  if (length(extra_pops) > 0 && verbose >= 2) {
    cat(report(
      "  Populations in x that are not in qmat are not mapped:",
      paste(extra_pops, collapse = ", "), "\n"
    ))
  }
  centers <- centers[match(pops, rownames(centers)), , drop = FALSE]

  cx <- centers[, "lon"]
  cy <- centers[, "lat"]
  lon_range <- abs(diff(range(centers[, "lon"])))
  # all centres at one longitude: bars 0.01 degrees wide (times scalex)
  if (lon_range == 0) {
    lon_range <- 1
  }
  sx <- lon_range / (100) * scalex
  sy <- 20 * sx * scaley
  #
  npops <- length(pops)
  ff <- qmat[, grepl("^cluster", colnames(qmat)), drop = FALSE]
  ll <- data.frame(cbind(as.numeric(qmat$orig.pop), ff))
  zz <- do.call(order, unname(as.list(ll)))
  bb <- qmat[zz, ]
  bb$orig.pop <- factor(bb$orig.pop)

  ff <- bb[, grepl("^cluster", colnames(bb)), drop = FALSE]

  out <- list()
  m1 <- leaflet::leaflet() %>%
    leaflet::addProviderTiles(provider = provider)
  for (p in 1:npops) {
    qmi <- ff[bb$orig.pop == levels(bb$orig.pop)[p], , drop = FALSE]
    out[[p]] <- bb[bb$orig.pop == levels(bb$orig.pop)[p], ]
    names(out)[p] <- levels(bb$orig.pop)[p]
    qmi1 <- cbind(rep(0, nrow(qmi)), qmi)
    for (xx in 1:nrow(qmi1)) {
      qmi1[xx, ] <- cumsum(as.numeric(qmi1[xx, ]))
    }


    for (ii in 1:nrow(qmi1)) {
      for (i in 1:(ncol(qmi1) - 1)) {
        oo <- (ii - nrow(qmi) / 2) * sx

        m1 <- m1 %>%
          leaflet::addRectangles(
            cx[p] + oo,
            cy[p] + qmi1[ii, i] * sy,
            cx[p] + oo + sx,
            cy[p] + qmi1[ii, i + 1] * sy,
            opacity = 0,
            color = plot.colors[i],
            fillOpacity = 0.8
          )
      }
    }
  }
  if (pop.labels) {
    m1 <- m1 %>%
      leaflet::addLabelOnlyMarkers(
        lng = centers[, "lon"],
        lat = centers[, "lat"] - sy * 0.1,
        label = rownames(centers),
        labelOptions = leaflet::labelOptions(
          noHide = T,
          direction = "center",
          textOnly = T,
          textsize = paste0(pop.labels.cex, "px")
        )
      )
  }

  if (plot.out) {
    print(m1)
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(list(qmats = out, map = m1))
}
