#' @name gl.run.nQuack
#' 
#' @title Run nQuack ploidy inference pipeline on BAM files
#'
#' @description
#' Wrapper function to execute a complete \pkg{nQuack} workflow for
#' ploidy inference from BAM files, including data preparation,
#' filtering, mixture model fitting, ploidy interpretation, and
#' accuracy summarisation.
#'
#' @param outputPath Character. Path to the directory where all output
#'   folders will be created [required].
#' @param bamInput Character. Path to a directory containing BAM files
#'   to be analysed [required].
#' @param parallel Logical. Run computationally intensive steps in
#'   parallel using \pkg{foreach} [default TRUE].
#' @param numCores Integer. Number of CPU cores to use. If set to 0,
#'   all available cores minus one are used [default 0].
#' @param min.depth Integer. Minimum total read depth threshold applied
#'   during locus filtering [default 2].
#' @param max.depth.quantile.prob Numeric. Upper quantile probability
#'   used to filter excessively high depth loci [default 0.9].
#' @param error Numeric. Expected sequencing error rate supplied to
#'   \pkg{nQuack} filtering routines [default 0.01].
#' @param trunc Numeric vector of length two. Truncation limits applied
#'   to allele frequency filtering [default c(0,0)].
#' @param verbose Verbosity: 0 (silent), 1 (start/end), 2 (progress),
#'   3 (summary), 5 (full report) [default as set by gl.set.verbosity].
#'
#' @details
#' The function creates a structured output directory
#' \code{runnQuack_output} inside \code{outputPath}. Intermediate and
#' final outputs are written to disk and are not returned as R objects.
#'
#' This function is intended for batch processing and HPC environments.
#'
#' @return Invisibly returns \code{NULL}. All outputs are written to disk.
#'
#' @author
#' Ethan Halford  
#' (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' \dontrun{
#' gl.runnQuack(
#'   outputPath = "results",
#'   bamInput = "bam_files",
#'   parallel = TRUE,
#'   numCores = 8
#' )
#' }
#'
#' @seealso
#' \code{\link[nQuack]{prepare_data}},
#' \code{\link[nQuack]{process_data}},
#' \code{\link[nQuack]{quackit}}
#'
#' @family pipelines
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom dplyr left_join mutate group_by summarise
#' @importFrom kableExtra kbl
#'
#' @export
gl.run.nQuack <- function(outputPath,
                         bamInput,
                         parallel = TRUE,
                         numCores = 0,
                         min.depth = 2,
                         max.depth.quantile.prob = 0.9,
                         error = 0.01,
                         trunc = c(0, 0),
                         verbose = NULL) {
  
  #---------------------------------------------------------------
  # SET VERBOSITY
  #---------------------------------------------------------------
  verbose <- gl.check.verbosity(verbose)
  
  #---------------------------------------------------------------
  # FLAG SCRIPT START
  #---------------------------------------------------------------
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, v = verbose)
  
  #---------------------------------------------------------------
  # FUNCTION-SPECIFIC ERROR CHECKING
  #---------------------------------------------------------------
  if (!dir.exists(bamInput)) {
    stop("Input BAM directory does not exist:", bamInput, "\n")
  }
  
  pkgs <- c("nQuack", "foreach", "doParallel", "kableExtra", "dplyr")
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package", pkg, "is required but not installed.\n")
    }
  }
  
  #---------------------------------------------------------------
  # CREATE OUTPUT DIRECTORIES
  #---------------------------------------------------------------
  quackOutput <- paste0(outputPath, "/runnQuack_output")
  
  dir.create(quackOutput, recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(quackOutput, "/02_prepared/temp"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(quackOutput, "/03_processed"), showWarnings = FALSE)
  dir.create(paste0(quackOutput, "/04_output"), showWarnings = FALSE)
  dir.create(paste0(quackOutput, "/05_interpret"), showWarnings = FALSE)
  dir.create(paste0(quackOutput, "/06_quackitTable"), showWarnings = FALSE)
  
  if (verbose >= 2) {
    cat(report("Output directories created in", quackOutput, "\n"))
  }
  
  #---------------------------------------------------------------
  # PARALLEL BACKEND
  #---------------------------------------------------------------
  if (parallel) {
    cores <- if (numCores == 0) {
      parallel::detectCores() - 1
    } else {
      numCores
    }
    
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    
    if (verbose >= 2) {
      cat(report("Parallel backend registered using", cores, "cores\n"))
    }
  }
  
  #---------------------------------------------------------------
  # 1. PREPARE DATA
  #---------------------------------------------------------------
  bamList <- list.files(bamInput, pattern = "\\.bam$")
  bamList <- gsub("\\.bam$", "", bamList)
  
  foreach::foreach(i = seq_along(bamList), .packages = "nQuack") %dopar% {
    
    nQuack::prepare_data(
      name = bamList[i],
      inpath = paste0(bamInput, "/"),
      outpath = paste0(quackOutput, "/02_prepared/"),
      tempfolder = paste0(quackOutput, "/02_prepared/temp")
    )
  }
  
  #---------------------------------------------------------------
  # 2. PROCESS PREPARED FILES
  #---------------------------------------------------------------
  textFiles <- list.files(paste0(quackOutput, "/02_prepared"), pattern = "\\.txt$")
  
  foreach::foreach(i = seq_along(textFiles), .packages = "nQuack") %dopar% {
    samp <- textFiles[i]
    df <- nQuack::process_data(
      paste0(quackOutput, "/02_prepared/", samp),
      min.depth = min.depth,
      max.depth.quantile.prob = max.depth.quantile.prob,
      error = error,
      trunc = trunc
    )
    
    write.csv(
      df,
      file = paste0(quackOutput, "/03_processed/",
                    gsub("\\.txt$", "", samp), ".csv"),
      row.names = FALSE
    )
  }
  
  #---------------------------------------------------------------
  # 3. FIT MIXTURE MODELS
  #---------------------------------------------------------------
  samples <- list.files(paste0(quackOutput, "/03_processed"),
                        pattern = "\\.csv$")
  
  foreach::foreach(i = seq_along(samples), .packages = "nQuack") %dopar% {
    
    xm <- as.matrix(
      read.csv(paste0(quackOutput, "/03_processed/", samples[i]))
    )
    
    out <- rbind(
      nQuack::quackNormal(xm, samples[i], cores = 10, parallel = FALSE),
      nQuack::quackBeta(xm, samples[i], cores = 10, parallel = FALSE),
      nQuack::quackBetaBinom(xm, samples[i], cores = 10, parallel = FALSE)
    )
    
    write.csv(out,
              file = paste0(quackOutput, "/04_output/", samples[i]),
              row.names = FALSE)
  }
  
  #---------------------------------------------------------------
  # 4. MODEL INTERPRETATION
  #---------------------------------------------------------------
  modelexplore <- list.files(paste0(quackOutput, "/04_output"),
                             pattern = "\\.csv$")
  
  for (f in modelexplore) {
    
    df <- read.csv(paste0(quackOutput, "/04_output/", f))
    
    summary <- nQuack::quackit(
      model_out = df,
      summary_statistic = "BIC",
      mixtures = c("diploid", "triploid", "tetraploid")
    )
    
    write.csv(summary,
              file = paste0(quackOutput, "/05_interpret/", f),
              row.names = FALSE)
  }
  
  #---------------------------------------------------------------
  # 5. ACCURACY SUMMARY
  #---------------------------------------------------------------
  quackitFiles <- list.files(paste0(quackOutput, "/05_interpret"),
                             pattern = "\\.csv$")
  
  ploidy <- sapply(quackitFiles, function(f) {
    read.csv(paste0(quackOutput, "/05_interpret/", f))[1, 2]
  })
  
  key <- data.frame(sample = quackitFiles,
                    ploidal.level = ploidy)
  
  dfs <- lapply(paste0(quackOutput, "/05_interpret/", quackitFiles),
                read.csv)
  
  alloutput <- dplyr::left_join(
    do.call(rbind, dfs),
    key
  )
  
  sumcheck <- alloutput %>%
    dplyr::mutate(accuracy = ifelse(winnerBIC == ploidal.level, 1, 0)) %>%
    dplyr::group_by(Distribution, Type) %>%
    dplyr::summarise(
      total = n(),
      correct = sum(accuracy),
      accuracy_rate = correct / total,
      .groups = "drop"
    )
  
  write.csv(sumcheck,
            file = paste0(quackOutput,
                          "/06_quackitTable/quackit_accuracy_summary.csv"),
            row.names = FALSE)
  
  #---------------------------------------------------------------
  # CLEAN UP
  #---------------------------------------------------------------
  if (parallel) {
    parallel::stopCluster(cl)
  }
  
  #---------------------------------------------------------------
  # FLAG SCRIPT END
  #---------------------------------------------------------------
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  invisible(NULL)
}


























