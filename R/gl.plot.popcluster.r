#' @name gl.plot.popcluster
#'
#' @title Plots PopCluster analysis results (Admixture Model)
#'
#' @description
#' This function takes a structure run object (output from
#'  \code{\link{gl.plot.popcluster}}) and plots the typical structure bar
#'   plot that visualize the q matrix of a structure run.
#'   
#' @param gl genlight object for PopCluster [required].
#' @param filename prefixe of run object from \code{\link{gl.run.popcluster}} [required].
#' @param plot.K The number for K of the q matrix that should be plotted. Needs to
#'  be within you simulated range of K's in your sr structure run object. If
#'  NULL, all the K's are plotted [default NULL].
#' @param plot_theme Theme for the plot. See Details for options
#' [default NULL].
#' @param color_clusters A color palette for clusters (K) or a list with
#' as many colors as there are clusters (K) [default NULL].
#' @param ind_name Whether to plot individual names [default TRUE].
#' @param input.dir Directory with PopCluster output and best run summary
#' @param k_name Name of the structure plot to plot. It should be character
#'  [default NULL].
#' @param border_ind The width of the border line between individuals
#' [default 0.25].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot.dir Directory in which to save files [default = working directory]
#' @param plot.file Name for the RDS binary file to save (base name only, exclude
#' extension) [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log ; 3, progress and results summary; 5, full report [default
#'   NULL, unless specified using gl.set.verbosity]
#'
#' @details The function outputs a barplot which is the typical output of
#'  PopCluster For a Evanno plot use gl.evanno.
#'  Plots and table are saved to the working directory specified in plot.dir (tempdir)
#'  if plot.file is set.
#'
#' Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return List of Q-matrices
#'
#' @author Ching Ching Lau (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' # examples need structure to be installed on the system (see above)
#' \dontrun{
#' m <- gl.run.popcluster(testset.gl, popcluster.path,output.path,
#'   filename, minK, maxK, 
#'   rep, search_relate, allele_freq,PopData,PopFlag,
#'   model, location, loc_admixture, relatedness,
#'   kinship, pr_allele_freq, cleanup=TRUE, verbose=NULL
#' )
#'gl.plot.popcluster <- function(testset.gl, 
#'                               filename=NULL,
#'                               ind_name=F,
#'                               input.dir="/Users/Documents/PopCluster/Results"
#'                               border_ind=0.25,
#'                               plot.K = NULL,
#'                               plot_theme=NULL,
#'                               k_name=NULL,
#'                               plot.out=TRUE,
#'                               plot.file=NULL,
#'                               plot.dir=NULL,
#'                               cleanup=TRUE)}
#' @export
#' @seealso \code{gl.run.popcluster}, \code{gl.plot.popcluster}
#' @references
#' \itemize{
#' \item Wang, J. (2022). Fast and accurate population admixture inference 
#' from genotype data from a few microsatellites to millions of SNPs. Heredity, 129(2), 79-92.
#' 
#' }

gl.plot.popcluster <- function(gl,
                              filename=NULL,
                              ind_name=F,
                              input.dir=NULL,
                              border_ind=0.25,
                              plot.K = NULL,
                              plot_theme=NULL,
                              k_name=NULL,
                              plot.out=TRUE,
                              plot.file=NULL,
                              plot.dir=NULL,
                              cleanup=TRUE) {
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
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
  summary <- read.table(paste0(input.dir, paste0(filename, ".popcluster.best_run_summary")), header = T)
  best_run_file <- summary[which(summary$K == plot.K),'BestRun']
  best <- readLines(paste0(input.dir, best_run_file))
  
  tempd <-  tempfile(pattern = "dir")
  dir.create(tempd, showWarnings = FALSE)
  
  write.table(best[(which(startsWith(best, "Inferred ancestry of individuals"))+2):
                     (which(startsWith(best, "Inferred ancestry of individuals"))+1+adegenet::nInd(gl))], 
              paste0(tempd, "/", filename,".popcluster.Qmatrix"), 
              quote = F, row.names = F, col.names = F)
  
  Q <- read.table(paste0(tempd, "/", filename,".popcluster.Qmatrix"))[,-6]
  colnames(Q) <- c("Index", "Order", "Label", "PercentMiss", "Pop", paste0("Pop_", seq(1, plot.K, by=1)))
  #Q$Pop <- "Species"
  Q_long <- tidyr::pivot_longer(Q, cols = starts_with("Pop_"), names_to = "K", values_to = "values")
  
  
  if (is.null(color_clusters)) {
    color_clusters <- gl.select.colors(ncolors = max(plot.K), verbose = 0)
  }
  
  p3 <- ggplot(Q_long, aes_(x = ~ factor(Order), y = ~values, fill = ~K)) +
    geom_col(color = "black", size = 0.15, width = 1, position = "fill") +
    facet_grid( ~ Pop, scales = "free", space = "free") +
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
  
  if (cleanup) unlink(tempd, recursive = T)
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
