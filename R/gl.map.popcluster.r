#' @name gl.map.popcluster
#' @title Maps a PopCluster plot using a genlight object
#' @description
#' This function takes the output of gl.plot.popcluster (the Q matrix) and maps the
#' Q-matrix across using the population centers from the genlight object that
#' was used to run the PopCluster analysis via \code{\link{gl.run.popcluster}})
#' and plots the typical PopCluster bar plots on a spatial map, providing a
#' barplot for each subpopulation. Therefore it requires coordinates from a
#'  genlight object. This kind of plots should support the interpretation of the
#'   spatial PopCluster of a population, but in principle is not different from
#'   \code{\link{gl.plot.popcluster}}
#' @param gl genlight object for PopCluster [required].
#' @param qmat Q-matrix from a gl.plot.popcluster
#' [from \code{\link{gl.run.popcluster}} and \code{\link{gl.plot.popcluster}}]
#'  [required].
#' @param x Name of the genlight object containing the coordinates in the
#'  \code{\@other$latlon} slot to calculate the population centers [required].
#' @param plot.K The number for K to be plotted [required].
#' @param provider Provider	passed to leaflet. Check \link[leaflet]{providers}
#' for a list of possible backgrounds [default "Esri.NatGeoWorldMap"].
#' @param scalex Scaling factor to determine the size of the bars in x direction
#' [default 1].
#' @param scaley Scaling factor to determine the size of the bars in y direction
#'  [default 1].
#' @param movepops A two-dimensional data frame that allows to move the center of
#' the barplots manually in case they overlap. Often if populations are
#' horizontally close to each other. This needs to be a data.frame of the
#' dimensions [rows=number of populations, columns = 2 (lon/lat)]. For each
#' population you have to specify the x and y (lon and lat) units you want to
#' move the center of the plot, (see example for details) [default NULL].
#' @param pop.labels Switch for population labels below the parplots
#' [default TRUE].
#' @param pop.labels.cex Size of population labels [default 12].
#' @return An interactive map that shows the PopCluster plots broken down by
#' population.
#' @author Ching Ching Lau (Post to \url{https://groups.google.com/d/forum/dartr})
#' @details
#' Creates a mapped version of PopCluster plots. For possible background maps
#' check as specified via the provider:
#' \url{http://leaflet-extras.github.io/leaflet-providers/preview/index.html}.
#' You may need to adjust scalex and scaley values [default 1], as the size
#' depends on the scale of the map and the position of the populations.
#' @return returns the map and a list of the qmat split into sorted matrices per
#'  population. This can be used to create your own map.
#' @examples
#' # examples need popcluster to be installed on the system
#' \dontrun{
#' m <- gl.run.popcluster(testset.gl, popcluster.path,output.path,
#'   filename, minK, maxK, 
#'   rep, search_relate, allele_freq,PopData,PopFlag,
#'   model, location, loc_admixture, relatedness,
#'   kinship, pr_allele_freq, cleanup=TRUE, verbose=NULL
#' )
#' qmat <- gl.plot.popcluster(gl=testset.gl,
#'                            filename=NULL,
#'                            ind_name=F,
#'                            bestK_file=NULL,
#'                            border_ind=0.25,
#'                            plot.K = NULL,
#'                            plot_theme=NULL,
#'                            k_name=NULL,
#'                            plot.out=TRUE,
#'                            plot.file=NULL,
#'                            plot.dir=NULL,
#'                              cleanup=TRUE)
#'                               
#' gl.map.popcluster(testset.gl, qmat,K=4)
#' # move population 4 (out of 5) 0.5 degrees to the right and populations 1
#' # 0.3 degree to the north of the map.
#' mp <- data.frame(lon=c(0,0,0,0.5,0), lat=c(-0.3,0,0,0,0))
#' gl.map.popcluster(testset.gl,qmat,K=4, movepops=mp)
#' }
#' @export
#' @seealso \code{\link{gl.run.popcluster}},  \code{clumpp},
#' \code{\link{gl.plot.popcluster}}
#' @references
#' \itemize{
#' \item Wang, J. (2022). Fast and accurate population admixture inference 
#' from genotype data from a few microsatellites to millions of SNPs. Heredity, 129(2), 79-92.
#' 
#' }

gl.map.popcluster <- function(gl,
                             qmat,
                             x,
                             K,
                             provider = "Esri.NatGeoWorldMap",
                             scalex = 1,
                             scaley = 1,
                             movepops = NULL,
                             pop.labels = TRUE,
                             pop.labels.cex = 12) {

  ff <- Q[, which(grepl("Pop_", colnames(Q)))]

  df <- gl@other$latlon
  centers <-
    apply(df, 2, function(xx) {
      tapply(xx, pop(gl), mean, na.rm = TRUE)
    })

  if (!is.null(movepops)) {
    if (nrow(movepops) != nrow(centers)) {
      stop(
        error(
          "The provided movepops data.frame has not the corret number of rows, please check. It needs to have the same numbers of rows as the number populations in your genlight object."
        )
      )
    }
    centers[, 1] <- centers[, 1] + movepops[, 1]
    centers[, 2] <- centers[, 2] + movepops[, 2]
  }
  
  Q_name <- left_join(Q, data.frame(Label=gl$ind.names, Pop_name=gl$pop))
  
  sc <-
    match(rownames(centers), levels(factor(Q_name$Pop_name)))
  if (any(is.na(sc))) {
    message(
      error(
        "Population names (coordinates) in the genlight object do not match population in your q-matrix. Please check both."
      )
    )
  }
  
  centers <- centers[sc, ]
  cx <- centers[, "lon"]
  cy <- centers[, "lat"]
  sx <- abs(diff(range(centers[, "lon"]))) / (100) * scalex
  sy <- 20 * sx * scaley
  #
  Q_name$Pop_name <- factor(Q_name$Pop_name)
  npops <- length(levels(Q_name$Pop_name))
  ll <- data.frame(cbind(as.numeric(Q_name$Pop_name), ff))
  zz <- do.call(order, unname(as.list(ll)))
  bb <- Q_name[zz, ]
  bb$Pop_name <- factor(bb$Pop_name)
  # ff <- bb[, 4:(ncol(bb))]

  ff <- bb[, which(grepl("Pop_", colnames(bb)))]

  out <- list()
  m1 <- leaflet::leaflet() %>%
    leaflet::addProviderTiles(provider = provider)
  for (p in 1:npops) {
    qmi <- ff[bb$Pop_name == levels(bb$Pop_name)[p], ]
    out[[p]] <- bb[bb$Pop_name == levels(bb$Pop_name)[p], ]
    names(out)[p] <- levels(bb$Pop_name)[p]
    qmi1 <- cbind(rep(0, nrow(qmi)), qmi)
    for (xx in 1:nrow(qmi1)) {
      qmi1[xx, -c(ncol(qmi1))] <- cumsum(as.numeric(qmi1[xx, -c(ncol(qmi1))]))
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
            color = rainbow(ncol(ff))[i],
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


  print(m1)
  
  # mapshot(m1, file='./Rplot.png', remove_controls = TRUE)
  return(list(Q_name=out,map=m1))
  # %>% addLegend(labels=paste('Group',1:ncol(ff)), colors=rainbow(ncol(ff)),position ='topright' )

  # if (save) mapshot(m1, file='./Rplot.png', remove_controls = TRUE)
}
