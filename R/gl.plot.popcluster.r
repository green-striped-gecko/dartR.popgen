#' @name gl.plot.popcluster
#'
#' @title Plots PopCluster analysis results (Admixture Model)
#'
#' @description
#' This function takes a Q matrix (output from
#'  \code{\link{gl.run.popcluster}}) and plots the typical structure bar
#'   plot that visualize the Q matrix of a structure run.
#'   
#' @param pop_cluster_result run object from \code{\link{gl.run.popcluster}} [required].
#' @param border_ind The width of the border line between individuals
#' [default 0.25].
#' @param plot.K The number for K of the q matrix that should be plotted. Needs to
#'  be within you simulated range of K's in your PopCluster run object. [required]
#' @param plot_theme Theme for the plot. See Details for options
#' [default NULL].
#' @param color_clusters A color palette for clusters (K) or a list with
#' as many colors as there are clusters (K) [default NULL].
#' @param ind_name Whether to plot individual names [default TRUE].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot.file Name for the RDS binary file to save (base name only, exclude
#' extension) [default NULL]
#' @param plot.dir Directory in which to save files [default = working directory]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log ; 3, progress and results summary; 5, full report [default
#'   2, unless specified using gl.set.verbosity]
#'
#' @details The function outputs a barplot and Q matrix which is the typical output of
#'  PopCluster.
#'  Plots and table are saved to the working directory specified in plot.dir
#'  if plot.file is set.
#'
#' Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'The Q matrices can be input to other R packages for plotting ancestry proportion, e.g. FSTruct
#'\url{https://github.com/MaikeMorrison/FSTruct}
#' @return Q-matrix and structure plot
#'
#' @author Ching Ching Lau (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' # examples need popcluster to be installed on the system (see above)
#' \dontrun{
#' m <- gl.run.popcluster(x=bandicoot.gl, popcluster.path="/User/PopCluster/Bin/",
#' output.path="/User/Documents/Output/",
#' minK=1, maxK=3, 
#' rep=10, PopData=1, location=1)
#' Q <- gl.plot.popcluster(pop_cluster_result=m, plot.K = 3, ind_name=T)
#' gl.map.popcluster(x = bandicoot.gl, qmat = Q)
#' # move population 4 (out of 5) 0.5 degrees to the right and populations 1
#' # 0.3 degree to the north of the map.
#' mp <- data.frame(lon=c(0,0,0,0.5,0), lat=c(-0.3,0,0,0,0))
#' gl.map.popcluster(bandicoot.gl, qmat=Q, movepops=mp)}
#' @export
#' @seealso \code{gl.run.popcluster}, \code{gl.plot.popcluster}
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
                              verbose=2) {
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(
    func = funname,
    build = "Jody",
    verbose = verbose
  )
  
  if (is.null(plot_theme)) {
    plot_theme <- theme_dartR()
  }
  
  # extract admixture analysis from best run
  best_run <- pop_cluster_result$best_run[which(pop_cluster_result$best_run$K == plot.K),'BestRun']
  Q <- pop_cluster_result$matrix[best_run][[1]]
  Q_long <- tidyr::pivot_longer(Q, cols = starts_with("Pop_"), names_to = "K", values_to = "values")
  
  if (is.null(color_clusters)) {
    color_clusters <- gl.select.colors(ncolors = max(plot.K), verbose = 0)
  }
  
  p3 <- ggplot(Q_long, aes_(x = ~ factor(Order), y = ~values, fill = ~K)) +
    geom_col(color = "black", size = border_ind, width = 1) +
    facet_grid( ~ factor(Pop, levels=unique(Q_long$Pop)), 
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
        size = 1
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
