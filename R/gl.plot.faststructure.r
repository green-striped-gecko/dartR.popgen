#' @name gl.plot.faststructure
#'
#' @title Plots fastStructure analysis results (Q-matrix)
#'
#' @family population structure
#'
#' @description
#' This function takes a fastStructure run object (output from
#'  \code{\link{gl.run.faststructure}}) and plots the typical structure bar
#'   plot that visualises the q matrix of a fastStructure run. The replicate
#'   runs are aligned, grouped into modes and averaged exactly as in
#'   \code{\link{gl.plot.structure}}, which draws the plot.
#'
#' @param sr fastStructure run object from \code{\link{gl.run.faststructure}}
#'  [required].
#' @param k.range The values of K to be plotted. Need to be among the K values
#'  in sr. If NULL, all the K's are plotted [default NULL].
#' @param met_clumpp The algorithm to use to infer the correct permutations.
#' One of 'greedy' or 'greedyLargeK' or 'stephens' [default "greedyLargeK"].
#' @param iter_clumpp The number of iterations to use if running either
#'  'greedy' or 'greedyLargeK' [default 100].
#' @param clumpak Whether use the Clumpak method (see details) [default TRUE].
#' @param plot_theme Theme for the plot. See Details for options
#' [default NULL, which uses theme_dartR()].
#' @param colors_clusters A colour palette function (for example
#'  \code{rainbow}), which is called with the largest K, or a vector with at
#'  least as many colours as clusters in the largest K [default NULL, which
#'  uses gl.select.colors()].
#' @param ind_name Whether to plot individual names [default TRUE].
#' @param k_name Label of the K panel to plot, as shown in the K column of the
#'  returned tables: "3" for K = 3, or "2.1", "2.2", ... when a K has more
#'  than one mode. It should be character [default NULL, all panels].
#' @param label.size Size of the population labels [default 12].
#' @param border_ind The width of the border line between individuals
#' [default 0.15].
#' @param den Whether to include a dendrogram. It needs either the genlight
#'  object used in gl.run.faststructure (parameter x) or a distance matrix
#'  (parameter dis.mat) [default FALSE].
#' @param x The genlight object used in gl.run.faststructure; needed for the
#'  dendrogram when dis.mat is not supplied [default NULL].
#' @param dis.mat A dist object (distance matrix) among the individuals, used
#'  to build the dendrogram; its labels must be the individual names
#'  [default NULL].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot.dir Directory in which to save files [default = tempdir(),
#'  unless specified using gl.set.wd].
#' @param plot.file Name for the RDS binary file to save (base name only,
#'  exclude extension) [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  brief progress messages; 3, progress and results summary; 5, full report
#'  [default 2, unless specified using gl.set.verbosity].
#'
#' @details The function outputs a barplot which is the typical output of
#'  fastStructure. The fastStructure run object is converted to the layout of
#'  a structure run object and passed to \code{\link{gl.plot.structure}}; see
#'  its help for the CLUMPP and Clumpak methods, the averaging of the
#'  replicates within each mode, and the dendrogram.
#'
#' Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return A list (invisible) with one data.table per plotted K panel, as
#' returned by \code{\link{gl.plot.structure}}: one row per individual, with
#' columns Label, cluster1 ... clusterK, K, orig.pop and ord.
#'
#' @author Author(s): Bernd Gruber & Luis Mijangos. Custodian: Luis Mijangos
#'  -- Post to \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \dontrun{
#' t1 <- gl.filter.callrate(platypus.gl, threshold = 1)
#' res <- gl.run.faststructure(t1,
#'   exec = "./fastStructure", k.range = 2:3,
#'   num.k.rep = 2, exec.plink = getwd()
#' )
#' qmat <- gl.plot.faststructure(res, k.range = 2:3)
#' gl.map.structure(qmat, K = 2, t1, scalex = 1, scaley = 0.5)
#' }
#' @export
#' @importFrom stats as.dendrogram dist hclust order.dendrogram reorder runif
#' @import ggdendro
#' @seealso \code{\link{gl.run.faststructure}}, \code{\link{gl.plot.structure}}
#' @references
#' \itemize{
#' \item Raj, A., Stephens, M., & Pritchard, J. K. (2014). fastSTRUCTURE:
#' variational inference of population structure in large SNP data sets.
#' Genetics, 197(2), 573-589.
#' \item Kopelman, Naama M., et al. "Clumpak: a program for identifying
#' clustering modes and packaging population structure inferences across K."
#' Molecular ecology resources 15.5 (2015): 1179-1191.
#' \item Mattias Jakobsson and Noah A. Rosenberg. 2007. CLUMPP: a cluster
#' matching and permutation program for dealing with label switching and
#' multimodality in analysis of population structure. Bioinformatics
#' 23(14):1801-1806. Available at
#' \href{http://web.stanford.edu/group/rosenberglab/clumppDownload.html}{clumpp}
#' }

gl.plot.faststructure <- function(sr,
                                  k.range = NULL,
                                  met_clumpp = "greedyLargeK",
                                  iter_clumpp = 100,
                                  clumpak = TRUE,
                                  plot_theme = NULL,
                                  colors_clusters = NULL,
                                  ind_name = TRUE,
                                  k_name = NULL,
                                  label.size = 12,
                                  border_ind = 0.15,
                                  den = FALSE,
                                  x = NULL,
                                  dis.mat = NULL,
                                  plot.out = TRUE,
                                  plot.dir = NULL,
                                  plot.file = NULL,
                                  verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  if (!is.list(sr) || is.null(sr$q_list) || length(sr$q_list) == 0) {
    stop(error(
      "sr must be the object returned by gl.run.faststructure.\n"
    ))
  }

  # DO THE JOB

  # one structure-like run per fastStructure replicate: q.mat with id,
  # pct.miss, orig.pop and Group.1..K, as gl.plot.structure expects
  runs <- list()
  for (k in names(sr$q_list)) {
    for (r in names(sr$q_list[[k]])) {
      q <- sr$q_list[[k]][[r]]
      if (!is.data.frame(q)) next
      grp <- as.matrix(q[, -(1:2), drop = FALSE])
      colnames(grp) <- paste0("Group.", seq_len(ncol(grp)))
      q.mat <- data.frame(
        id = as.character(q$id),
        pct.miss = 0,
        orig.pop = as.character(q$orig.pop),
        grp,
        stringsAsFactors = FALSE
      )
      lab <- paste0("k", k, ".r", r)
      runs[[lab]] <- list(
        summary = c(k = ncol(grp), est.ln.prob = NA, mean.lnL = NA,
                    var.lnL = NA),
        q.mat = q.mat,
        prior.anc = NULL,
        label = lab
      )
    }
  }
  if (length(runs) == 0) {
    stop(error("sr contains no q-matrices.\n"))
  }
  class(runs) <- c("structure.result", "list")

  Q_list <- gl.plot.structure(runs,
                              K = k.range,
                              met_clumpp = met_clumpp,
                              iter_clumpp = iter_clumpp,
                              clumpak = clumpak,
                              plot_theme = plot_theme,
                              color_clusters = colors_clusters,
                              ind_name = ind_name,
                              k_name = k_name,
                              border_ind = border_ind,
                              den = den,
                              dis.mat = dis.mat,
                              x = x,
                              plot.out = plot.out,
                              plot.file = plot.file,
                              plot.dir = plot.dir,
                              verbose = 0,
                              label.size = label.size)

  if (!is.null(plot.file) && verbose >= 2) {
    cat(report(
      "  Plot saved as", file.path(gl.check.wd(plot.dir, verbose = 0),
                                   paste0(plot.file, ".RDS")), "\n"
    ))
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(invisible(Q_list))
}
