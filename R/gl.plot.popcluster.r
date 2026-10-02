#' @name gl.plot.popcluster
#'
#' @title Plots PopCluster analysis results (Admixture Model)
#'
#' @family population structure
#'
#' @description
#' This function takes a PopCluster run object (output from
#'  \code{\link{gl.run.popcluster}}) and plots the typical structure bar
#'   plot that visualises the Q matrix of the best run for one K.
#'   
#' @param pop_cluster_result run object from \code{\link{gl.run.popcluster}} [required].
#' @param border_ind The width of the border line between individuals
#' [default 0.25].
#' @param plot.K The K of the Q matrix to be plotted: a single value among
#'  the K values in the PopCluster run object [required].
#' @param plot_theme Theme for the plot. See Details for options
#' [default NULL].
#' @param color_clusters A colour palette function (for example
#'  \code{rainbow}), which is called with plot.K, or a vector with at least
#'  plot.K colours [default NULL, which uses gl.select.colors()].
#' @param ind_name Whether to plot individual names [default TRUE].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot.file Name for the RDS binary file to save (base name only, exclude
#' extension) [default NULL]
#' @param plot.dir Directory in which to save files [default = tempdir(),
#'  unless specified using gl.set.wd].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  brief progress messages; 3, progress and results summary; 5, full report
#'  [default 2, unless specified using gl.set.verbosity].
#'
#' @details The function outputs a barplot and Q matrix which is the typical
#'  output of PopCluster, for the best run of plot.K chosen by PopCluster
#'  (replicates are not averaged). The plot is saved as an RDS file in
#'  plot.dir if plot.file is set.
#'
#' Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'The Q matrices can be input to other R packages for plotting ancestry proportion, e.g. FSTruct
#'\url{https://github.com/MaikeMorrison/FSTruct}
#' @return The Q matrix of plot.K (invisible), as in the matrix element of the
#' gl.run.popcluster result: Index, Order, Label, PercentMiss, Cluster,
#' Pop_1 ... Pop_K, Pop.
#'
#' @author Author(s): Ching Ching Lau. Custodian: Ching Ching Lau -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # examples need popcluster to be installed on the system (see above)
#' \dontrun{
#' m <- gl.run.popcluster(x = bandicoot.gl,
#'   popcluster.path = "/User/PopCluster/Bin/", minK = 1, maxK = 3, rep = 10)
#' Q <- gl.plot.popcluster(pop_cluster_result = m, plot.K = 3, ind_name = TRUE)
#' gl.map.popcluster(x = bandicoot.gl, qmat = Q)
#' }
#' @export
#' @seealso \code{\link{gl.run.popcluster}}, \code{\link{gl.map.popcluster}}
#' @references
#' \itemize{
#' \item Wang, J. (2022). Fast and accurate population admixture inference 
#' from genotype data from a few microsatellites to millions of SNPs. Heredity, 129(2), 79-92.
#' 
#' }

gl.plot.popcluster <- function(pop_cluster_result,
                              border_ind=0.25,
                              plot.K,
                              plot_theme=NULL,
                              color_clusters=NULL,
                              ind_name=T,
                              plot.out=TRUE,
                              plot.file=NULL,
                              plot.dir=NULL,
                              verbose=NULL) {
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(
    func = funname,
    verbose = verbose
  )

  # FUNCTION SPECIFIC ERROR CHECKING

  run_K <- as.numeric(pop_cluster_result$best_run$K)
  if (length(plot.K) != 1 || !plot.K %in% run_K) {
    stop(error(
      "plot.K must be one of the K values in pop_cluster_result:",
      paste(run_K, collapse = ", "), "\n"
    ))
  }

  if (is.function(color_clusters)) {
    color_clusters <- color_clusters(plot.K)
  }
  if (is.null(color_clusters)) {
    color_clusters <- gl.select.colors(ncolors = plot.K, verbose = 0)
  }
  if (length(color_clusters) < plot.K) {
    stop(error(
      "color_clusters has", length(color_clusters), "colours but", plot.K,
      "are needed.\n"
    ))
  }
  
  if (is.null(plot_theme)) {
    plot_theme <- theme_dartR()
  }
  
  # extract admixture analysis from best run
  best_run <- pop_cluster_result$best_run[which(run_K == plot.K), "BestRun"]
  Q <- pop_cluster_result$matrix[best_run][[1]]
  Q_long <- tidyr::pivot_longer(Q, cols = starts_with("Pop_"), names_to = "K", values_to = "values")
  
  p3 <- ggplot(Q_long, aes(x = factor(.data$Order), y = .data$values,
                           fill = .data$K)) +
    geom_col(color = "black", linewidth = border_ind, width = 1) +
    facet_grid( ~ factor(.data$Pop, levels = unique(Q_long$Pop)), 
                scales = "free", 
                space = "free") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(
      breaks = unique(Q_long$Order),
      labels = unique(Q_long$Label),
      expand = c(0, 0)
    ) +
    scale_fill_manual(values = color_clusters) +
    plot_theme +
    theme(
      panel.spacing = unit(0, "lines"),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      ),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 12, angle = 90),
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        size = 8,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    )
  
  if (ind_name == FALSE) {
    p3 <- p3 + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  
  if (plot.out) {
    print(p3)
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
  
  return(invisible(Q))
}
