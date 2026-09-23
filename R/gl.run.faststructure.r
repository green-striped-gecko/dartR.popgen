#' @name gl.run.faststructure
#'
#' @title Runs a faststructure analysis using a genlight object
#'
#' @family population structure
#'
#' @description
#' This function takes a genlight object and runs a faststructure analysis.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param k.range Values of K (number of populations) to run, e.g. 2:5 or
#' c(1, 3, 5) [required].
#' @param num.k.rep Number of replicates for each K [default 1].
#' @param exec Full path and name of the fastStructure executable
#' [default "./fastStructure"].
#' @param exec.plink Path of the folder that holds the PLINK executable,
#' named plink [default getwd()].
#' @param output Folder in which a new subfolder is created for the files of
#' this run (PLINK input, fastStructure output) [default tempdir()].
#' @param tol Convergence criterion [default 10e-6].
#' @param prior Choice of prior: simple or logistic [default "simple"].
#' @param cv Number of test sets for cross-validation, 0 implies no CV step
#'  [default 0].
#' @param seed Seed for the random number generator; replicate r uses
#' seed + r - 1, so replicates differ and can be repeated [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  one line per run; 3, progress and results summary, including the output
#'  of PLINK and fastStructure; 5, full report
#'  [default 2, unless specified using gl.set.verbosity].
#' @param plot.out Specify if the marginal likelihood plot is to be produced
#' [default TRUE].
#' @param plot.theme Theme for the plot [default theme_dartR()].
#' @param plot.dir Directory in which to save the plot [default = tempdir(),
#'  unless specified using gl.set.wd].
#' @param plot.file Name for the RDS binary file to save the plot (base name
#' only, exclude extension) [default NULL].
#'
#' @details
#' fastStructure does not run on Windows. Download the faststructure binary
#' for Mac or Linux from here:
#'
#' https://github.com/StuntsPT/Structure_threader/tree/master/structure_threader/bins
#'
#' and make it executable, e.g.
#' \code{system(paste0("chmod u+x ", getwd(), "/fastStructure"))}.
#'
#' Download the PLINK 1.9 binary for your system from here:
#'
#' https://www.cog-genomics.org/plink/
#'
#' and pass the folder that holds it as exec.plink. The genlight object is
#' converted to PLINK binary files with gl2plink.
#'
#' To install fastStructure dependencies follow these directions:
#' https://github.com/rajanil/fastStructure
#'
#' fastStructure performs inference for the simplest, independent-loci,
#' admixture model, with two choices of priors that can be specified using
#' the --prior parameter. Thus, unlike Structure, fastStructure does not require
#' the mainparams and extraparam files. The inference algorithm used by
#'  fastStructure is fundamentally different from that of Structure and
#'  requires the setting of far fewer options.
#'
#'  To identify the number of populations that best approximates the marginal
#'  likelihood of the data, the marginal likelihood is extracted from each run
#'  of K, averaged across replicates and plotted.
#'
#'  Each call writes its files to a new subfolder of output (named
#'  fastStructure_<date>_<time>), so files from earlier runs are never read.
#'
#' @return A list with two elements: q_list, a list named by K, each element
#' a list named by replicate ("1", "2", ...) of data frames with one row per
#' individual (columns id, orig.pop, V1 ... VK: the fastStructure
#' q-matrix); and plot, the ggplot of the mean marginal likelihood per K.
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \dontrun{
#' # Please note: faststructure needs to be installed
#' # Please note: faststructure is not available for windows
#' t1 <- gl.filter.callrate(platypus.gl, threshold = 1)
#' res <- gl.run.faststructure(t1,
#'   exec = "./fastStructure", k.range = 2:3,
#'   num.k.rep = 2, exec.plink = getwd()
#' )
#' qmat <- gl.plot.faststructure(res, k.range = 2:3)
#' gl.map.structure(qmat, K = 2, t1, scalex = 1, scaley = 0.5)
#' }
#' @export
#' @seealso \code{\link{gl.plot.faststructure}}
#' @references
#' \itemize{
#' \item Raj, A., Stephens, M., & Pritchard, J. K. (2014). fastSTRUCTURE:
#' variational inference of population structure in large SNP data sets.
#' Genetics, 197(2), 573-589.
#' }

gl.run.faststructure <- function(x,
                                 k.range,
                                 num.k.rep = 1,
                                 exec = "./fastStructure",
                                 exec.plink = getwd(),
                                 output = tempdir(),
                                 tol = 10e-6,
                                 prior = "simple",
                                 cv = 0,
                                 seed = NULL,
                                 verbose = NULL,
                                 plot.out = TRUE,
                                 plot.theme = theme_dartR(),
                                 plot.dir = NULL,
                                 plot.file = NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  if (.Platform$OS.type == "windows") {
    stop(error(
      "fastStructure is not available for Windows; run it on Mac or Linux.\n"
    ))
  }

  if (!file.exists(exec)) {
    stop(error(
      "The fastStructure executable was not found at", exec, ". Download it",
      "from https://github.com/StuntsPT/Structure_threader/tree/master/",
      "structure_threader/bins and set exec to its path.\n"
    ))
  }

  if (!file.exists(file.path(exec.plink, "plink"))) {
    stop(error(
      "The PLINK executable was not found in", exec.plink, ". Download PLINK",
      "1.9 from https://www.cog-genomics.org/plink/ and set exec.plink to",
      "the folder that holds it.\n"
    ))
  }

  if (!is.numeric(k.range) || any(k.range < 1) ||
      any(k.range != round(k.range))) {
    stop(error("k.range must hold positive whole numbers.\n"))
  }
  k.range <- unique(k.range)

  if (!is.numeric(num.k.rep) || length(num.k.rep) != 1 || num.k.rep < 1) {
    stop(error("num.k.rep must be a whole number of at least 1.\n"))
  }

  # DO THE JOB

  # a new folder for this call, so files from earlier runs are never read
  if (!dir.exists(output)) {
    dir.create(output, recursive = TRUE)
  }
  run_dir <- file.path(
    output, paste0("fastStructure_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  i <- 1
  base_dir <- run_dir
  while (dir.exists(run_dir)) {
    run_dir <- paste0(base_dir, "_", i)
    i <- i + 1
  }
  dir.create(run_dir, recursive = TRUE)
  if (verbose >= 2) {
    cat(report("  Writing the run files to", run_dir, "\n"))
  }

  plink_run <- function() {
    dartR.base::gl2plink(x,
      bed.files = TRUE,
      outpath = run_dir,
      verbose = 0,
      plink.bin.path = exec.plink
    )
  }
  if (verbose >= 3) {
    plink_run()
  } else {
    # gl2plink reports PLINK's output with message() at every verbosity
    utils::capture.output(suppressMessages(plink_run()))
  }
  plink_in <- file.path(run_dir, "gl_plink")
  if (!file.exists(paste0(plink_in, ".bed"))) {
    stop(error(
      "PLINK did not write", paste0(plink_in, ".bed"), ". Rerun with",
      "verbose = 3 to see its output; gl2plink needs a dartR.base version",
      "that runs PLINK with --make-bed.\n"
    ))
  }

  show <- if (verbose >= 3) "" else FALSE
  q_list <- list()
  lik <- data.frame(K = numeric(0), rep = numeric(0), ml = numeric(0))

  for (k_n in k.range) {
    for (rep_n in seq_len(num.k.rep)) {
      if (verbose >= 2) {
        cat(report("  Running K =", k_n, "replicate", rep_n, "\n"))
      }
      out_prefix <- file.path(run_dir, paste0("k", k_n, ".r", rep_n))
      args <- c(
        "-K", k_n,
        paste0("--input=", shQuote(plink_in)),
        paste0("--output=", shQuote(out_prefix)),
        paste0("--tol=", tol),
        paste0("--prior=", prior),
        paste0("--cv=", cv)
      )
      if (!is.null(seed)) {
        args <- c(args, paste0("--seed=", seed + rep_n - 1))
      }
      system2(exec, args, stdout = show, stderr = show)

      # fastStructure appends .<K>.meanQ and .<K>.log to the output name
      q_file <- paste0(out_prefix, ".", k_n, ".meanQ")
      log_file <- paste0(out_prefix, ".", k_n, ".log")
      if (!file.exists(q_file) || !file.exists(log_file)) {
        stop(error(
          "fastStructure failed for K =", k_n, "replicate", rep_n,
          "; rerun with verbose = 3 to see its output.\n"
        ))
      }

      log_lines <- readLines(log_file, warn = FALSE)
      ml_line <- grep("^Marginal Likelihood = ", log_lines, value = TRUE)
      ml <- as.numeric(sub("^Marginal Likelihood = ", "",
                           ml_line[length(ml_line)]))
      lik <- rbind(lik, data.frame(K = k_n, rep = rep_n, ml = ml))

      q_df <- utils::read.table(q_file)
      q_df <- cbind(id = indNames(x), orig.pop = pop(x), q_df)
      q_list[[as.character(k_n)]][[as.character(rep_n)]] <- q_df
    }
  }

  lik_k <- data.frame(
    K = k.range,
    ml = vapply(k.range, function(k) mean(lik$ml[lik$K == k]), numeric(1))
  )

  p3 <- ggplot(lik_k, aes(x = .data$K, y = .data$ml)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2, color = "blue") +
    plot.theme +
    xlab("K") +
    ylab("Marginal Likelihood") +
    scale_x_continuous(breaks = sort(k.range))

  if (verbose >= 3) {
    cat(report("  Mean marginal likelihood per K:\n"))
    print(lik_k, row.names = FALSE)
  }

  # PRINTING OUTPUTS
  if (plot.out) {
    print(p3)
  }

  # Optionally save the plot ---------------------
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p3,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(list(q_list = q_list, plot = p3))
}
