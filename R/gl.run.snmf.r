#' @name gl.run.snmf
#' @title Runs a snmf analysis (package LEA) using a genlight object
#' @family population structure
#' @description
#' Creates an input file for the function \link[LEA]{snmf} (package LEA)
#' and runs it. Refer to the documentation of function \link[LEA]{snmf} for
#' further information on the method and its parameters.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param filename File name of output data [default "output"].
#' @param minK Minimum K [default 1].
#' @param maxK Maximum K [default 2].
#' @param rep Number of replicates runs per K [default 1].
#' @param regularization Alpha value for regularization when analyzing small
#' dataset [default 10].
#' @param ploidy_lv Ploidy level of dataset [default 2].
#' @param ncores How many cores should be used [default 1].
#' @param cleanup Remove the LEA run files from the temporary folder when
#' finished [default TRUE].
#' @param plot.out Specify if the cross-entropy plot is to be shown; it is
#' returned either way [default TRUE].
#' @param plot.dir Directory in which to save files [default = tempdir(),
#'  unless specified using gl.set.wd].
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  brief progress messages; 3, progress and results summary, including the
#'  output of LEA; 5, full report
#'  [default 2, unless specified using gl.set.verbosity].
#' @param ... Parameters passed to function \link[LEA]{snmf} (package LEA),
#' e.g. seed.
#' @details
#' For each K, only the replicate run with the lowest cross-entropy (the
#' best run) is kept; replicates are not averaged.
#' @return A list with three elements: best_run, the best run for each K
#' (the path of its folder when cleanup = FALSE, otherwise its name, e.g.
#' "K2/run1"); cross_entropy, a list with the cross-entropy plot (a
#' recordedplot); and matrix, a list named K1, K2, ... of data frames with one
#' row per individual: columns Pop_1 ... Pop_K (ancestry coefficients of the
#' best run), Cluster (the cluster with the largest coefficient), Pop
#' (population), Label (individual name) and Order (position in the bar
#' plot, sorted by population and cluster).
#' @export
#' @importFrom LEA snmf
#' @importFrom LEA cross.entropy
#' @importFrom utils capture.output
#' @importFrom grDevices recordPlot
#' @references
#' \itemize{
#' \item Frichot E, Mathieu F, Trouillon T, Bouchard G, Francois O. (2014).
#' Fast and Efficient Estimation of Individual Ancestry Coefficients. Genetics,
#' 194(4): 973--983.
#' }
#' @author Author(s): Ching Ching Lau. Custodian: Ching Ching Lau -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
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

gl.run.snmf <- function(x,
                        filename = "output",
                        minK = 1,
                        maxK = 2,
                        rep = 1,
                        regularization = 10,
                        ploidy_lv = 2,
                        ncores = 1,
                        cleanup = TRUE,
                        plot.out = TRUE,
                        plot.dir = NULL,
                        plot.file = NULL,
                        verbose = NULL,
                        ...) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- dartR.base::utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  whole <- function(v) {
    is.numeric(v) && length(v) == 1 && !is.na(v) && v >= 1 && v == round(v)
  }
  if (!whole(minK) || !whole(maxK) || minK > maxK) {
    stop(error(
      "minK and maxK must be whole numbers of at least 1, with minK <= maxK.\n"
    ))
  }
  if (!whole(rep)) {
    stop(error("rep must be a whole number of at least 1.\n"))
  }

  # DO THE JOB
  
  #create tempdir
  tempd <-  tempfile(pattern = "dir")
  dir.create(tempd, showWarnings = FALSE)

  # LEA prints its own log; show it only at verbose >= 3
  run_quiet <- function(expr) {
    if (verbose >= 3) {
      expr
    } else {
      utils::capture.output(res <- expr)
      res
    }
  }
  
  #convert genlight object to geno file
  run_quiet(gl2geno(x, outpath = tempd, outfile = filename, verbose = 0))
  
  ss <- run_quiet(LEA::snmf(
    input.file = file.path(tempd, paste0(filename, ".geno")),
    K = minK:maxK,
    entropy = TRUE,
    repetitions = rep,
    CPU = ncores,
    project = "new",
    alpha = regularization,
    ploidy = ploidy_lv,
    ...
  ))
  
  #Choose best run
  best_run <- NULL
  best_run_path <- NULL
  K_range <- minK:maxK
  for (i in 1:length(K_range)) {
    ce <- run_quiet(LEA::cross.entropy(ss, K = K_range[i]))
    best_run <- c(best_run, paste0("run", which.min(ce)))
    best_run_path <- c(best_run_path, (file.path(
      tempd,
      paste0(filename, ".snmf"),
      paste0("K", K_range[i]),
      best_run[i]
    )))
  }
  
  #extract Q matrices from best run
  Q_matrices <- NULL
  for (i in 1:length(K_range)) {
    Q <- read.table(list.files(
      best_run_path[i],
      pattern = ".Q",
      full.names = T
    ))

    colnames(Q) <- paste0("Pop_", seq(1, K_range[i]))
    Q$Cluster <- apply(Q, 1, which.max)
    Q$Pop <- as.character(x$pop)
    Q$Label <- as.character(x$ind.names)
    Q <- Q[with(Q, order(Q$Pop, as.numeric(Q$Cluster))), ]
    Q$Order <- 1:nrow(Q)
    Q_matrices[[i]] <- Q
  }
  names(Q_matrices) <- paste0("K", K_range)
  
  # plot cross-entropy; drawn off-screen when plot.out = FALSE so that the
  # recorded plot is still returned
  plot.list = list()
  if (!plot.out) {
    grDevices::pdf(NULL)
  }
  plot(ss,
       cex = 1.2,
       col = "lightblue",
       pch = 19)
  plot.list[[1]] = recordPlot()
  if (!plot.out) {
    grDevices::dev.off()
  }
  names(plot.list) <- "cross-entropy"
  
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(plot.list,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }

  if (cleanup) {
    unlink(tempd, recursive = TRUE)
    best_run_out <- paste0("K", K_range, "/", best_run)
  } else {
    best_run_out <- best_run_path
  }
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # return all Q matrices and best run
  return(list(
    best_run = best_run_out,
    cross_entropy = plot.list,
    matrix = Q_matrices
  ))
}
