#' @name gl.run.popcluster
#' 
#' @title Runs a PopCluster analysis using a genlight object
#'
#' @family population structure
#' 
#' @description
#' Creates an input file for the program PopCluster and runs it if
#' PopCluster is installed (can be installed at:
#'  https://www.zsl.org/about-zsl/resources/software/popcluster)
#'  
#' If you specify a directory for the PopCluster executable file, then the
#' script will create the input file (DataForm=0) from the SNP data then run 
#' PopCluster.
#'
#' PopCluster infers population admixture by coupling a clustering stage with
#' a subsequent admixture-analysis stage. First, it uses simulated annealing to
#' assign individuals to clusters under a mixture model, thus identifying
#' discrete populations and estimating allele frequencies without prematurely
#' converging to local optima. In the second step, these results provide
#' starting points for an expectation–maximization (EM) algorithm under an
#' admixture model, where each individual’s genetic contributions from multiple
#' populations are refined.
#' 
#' Refer to the PopCluster manual for further information on the parameters to
#' set. 
#' 
#' @param x Name of the genlight object containing the SNP data [required].
#' @param popcluster.path Path to the directory that contains the PopCluster
#' program (PopClusterMac, PopClusterLnx or PopClusterWin.exe)
#' [default getwd()].
#' @param output.path Folder in which the PopCluster parameter file and input
#' data are written; created if missing [default tempdir()].
#' @param filename Prefix of all the files that will be produced
#'  [default “output”].
#' @param minK Minimum K [default 1].
#' @param maxK Maximum K [default 2].
#' @param rep Number of replicates runs per K [default 1].
#' @param Scaling Scaling to be applied in the clustering analysis: none (0), 
#' weak (1), medium (2), strong (3) and very strong (4), see details 
#' section [default 0].
#' @param search_relate Method for proposing a configuration in clustering 
#' analysis. 0 for the assignment probability method and 1 for relatedness 
#' method. [default 0].
#' @param allele_freq Output allele frequency: 0=N, 1=Y [default 1].
#' @param ISeed  Seed for random number generator [default 333].
#' @param PopFlag Whether to use population information stored in the genlight 
#' object in the slot "pop" in structure analysis. 0=No and 1=Yes [default 0].
#' @param model 1=Clustering, 2=Admixture, 3=Hybridyzation, 4=Migration 
#' model [default 2].
#' @param loc_admixture Whether to estimate and output the admixture 
#' proportions for each individual at each locus (=1) or not (=0) [default 0].
#' @param relatedness Compute relatedness = 0=No, 1=Wang, 2=LynchRitland
#'  [default 0].
#' @param kinship Estimate kinship: 0=N, 1=Y [default 0].
#' @param pr_allele_freq Whether allele frequency prior should be determined 
#' by the program (0), the Equal Frequency prior (1) or Unequal Frequency 
#' prior (2) [default 2].
#' @param parallel Use parallelisation (implemented only in LINUX for the
#'  moment) [default FALSE].
#' @param ncores How many cores should be used [default 1].
#' @param cleanup Remove the temporary folder in which PopCluster runs (a
#' copy of the program and all its output files) when finished
#' [default TRUE].
#' @param plot.dir Directory in which to save files [default = tempdir(),
#'  unless specified using gl.set.wd].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot.file Name for the RDS binary file to save (base name only, 
#' exclude extension) [default NULL].
#' @param plot_theme Theme of the plots [default theme_dartR()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  brief progress messages; 3, progress and results summary, including the
#'  output of PopCluster; 5, full report
#'  [default 2, unless specified using gl.set.verbosity].
#' 
#' @details
#'
#' For best results, run multiple replicates with different starting seeds to
#' verify convergence and consistency. 
#' 
#' Use scaling when your sampling is highly unbalanced (e.g., one population 
#' with few individuals vs. another with many). Applying an appropriate 
#' scaling level (1, 2, 3, or 4) can substantially improve structure 
#' inference in these cases.
#' 
#' If your sample has many closely related individuals, using the Equal 
#' Frequency Prior (pr_allele_freq = 1)  gives better admixture results. If 
#' your sample doesn't include many relatives, the Unequal Frequency Prior 
#' (pr_allele_freq = 2) is more accurate. If you're unsure about how related 
#' the individuals in your sample are, set pr_allele_freq = 0. This 
#' will let the program check for relatedness and automatically choose the 
#' best prior (Equal or Unequal) based on the results.
#'
#' For each K, the Q matrix of the best run chosen by PopCluster is returned;
#' replicates are not averaged.
#'
#' PopCluster builds carry an expiry date, printed when the program starts
#' (seen with verbose = 3); download a current build when it has passed.
#'
#' @return A list with: output_path, the folder with the parameter and input
#' files; best_run, a data frame with one row per K and the numeric columns
#' K, LogL_Mean, LogL_Min, LogL_Max, DLK1, DLK2, FST.FIS (NA where PopCluster
#' reports "-") plus BestRun (the name of the best run); plots, a list of four
#' ggplots (LogL_Mean, DLK1, DLK2, FST.FIS against K); and matrix, a list
#' named by best run of data frames with one row per individual: Index,
#' Order (position in the bar plot), Label, PercentMiss, Cluster,
#' Pop_1 ... Pop_K (ancestry proportions) and Pop (population).
#'
#' @importFrom pillar align
#' @importFrom stringr str_split
#' @importFrom stats na.omit
#' @importFrom utils capture.output
#' @references
#' \itemize{
#' \item Wang, J. (2022). Fast and accurate population admixture inference
#' from genotype data from a few microsatellites to millions of SNPs. 
#' Heredity, 129(2), 79-92.
#' }
#' @author Author(s): Ching Ching Lau. Custodian: Ching Ching Lau -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
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
#'
#' @export

gl.run.popcluster <- function(x,
                              popcluster.path = getwd(),
                              output.path = tempdir(),
                              filename = "output",
                              minK = 1,
                              maxK = 2,
                              rep = 1,
                              Scaling = 0,
                              search_relate = 0,
                              allele_freq = 1,
                              ISeed = 333,
                              PopFlag = 0,
                              model = 2,
                              loc_admixture = 0,
                              relatedness = 0,
                              kinship = 0,
                              pr_allele_freq = 2,
                              parallel = FALSE,
                              ncores = 1,
                              cleanup = TRUE,
                              plot.dir = NULL,
                              plot.out = TRUE,
                              plot.file = NULL,
                              plot_theme = theme_dartR(),
                              verbose = NULL) {

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

  if (model == 4 &&
      (minK != maxK || maxK != nPop(x) || PopFlag == 0)) {
    stop(error(
      "For the migration model (model = 4), minK and maxK must be equal to",
      "each other and to the number of populations in x, and PopFlag must",
      "be 1.\n"
    ))
  }

  if (PopFlag == 1 && minK < nPop(x)) {
    stop(error(
      "If population information is used (PopFlag = 1), minK must be at",
      "least the number of populations in x.\n"
    ))
  }
  
  # check OS
  os <- Sys.info()['sysname']
  if (os == "Windows") {
    popcluster_version <- c(paste0("PopCluster", "Win.exe"),
                            "impi.dll",
                            "libiomp5md.dll")
  } else if (os == "Darwin") {
    popcluster_version <- paste0("PopCluster", "Mac")
  } else if (os == "Linux") {
    if (parallel) {
      popcluster_version <- "PopClusterLnx_impi"
    } else {
      popcluster_version <- paste0("PopCluster", "Lnx")
    }
  } else {
    stop(error(
      "PopCluster runs on Windows, macOS and Linux; this system is", os, "\n"
    ))
  }

  fex <- file.exists(file.path(popcluster.path, popcluster_version))
  if (!all(fex)) {
    stop(error(
      "Cannot find", paste(popcluster_version[!fex], collapse = ", "),
      "in popcluster.path:", popcluster.path, ". Download PopCluster from",
      "https://www.zsl.org/about-zsl/resources/software/popcluster\n"
    ))
  }
  
  # DO THE JOB

  if (!dir.exists(output.path)) {
    dir.create(output.path, recursive = TRUE)
  }

  #create tempdir
  tempd <- tempfile(pattern = "dir")
  dir.create(tempd, showWarnings = FALSE)
  
  # create INPUT FILE
  if(PopFlag == 0){
    PopData <- 0
  }else if(PopFlag == 1){
    PopData <- 1
  }
  
  genotype <- as.matrix.genlight(x)
  genotype[is.na(genotype)] <- 3
  sample_name <- x@ind.names
  ind_numbers <- seq_len(nInd(x))
  family <- x@pop
  rownames(genotype) <- NULL
  # IndivLoc is not used in structure inference. It is used solely for
  # visualizing population structuring in relation to individual geographic
  # locations in PopCluster’s GUI.
  names <- data.frame(id = paste0(ind_numbers,
                                  " ", family,
                                  " ", PopFlag))
  
  names2 <- apply(names, 1, paste0, collapse = " ")
  genotype2 <- apply(genotype, 1, paste0, collapse = "")
  
  writeLines(capture.output(for (i in 1:nInd(x)) {
    cat(names2[i], genotype2[i], sep = "\n")
  }), con = file.path(output.path, paste0(filename, ".popcluster.dat")))
  
  # PARAMETER from user input
  parameter <- c(
    nInd(x),
    nLoc(x),
    1,
    0,
    ISeed,
    paste0(filename, ".popcluster.dat") ,
    paste0(filename, ".popcluster"),
    Scaling,
    minK,
    maxK,
    rep,
    search_relate,
    allele_freq,
    PopData,
    PopFlag,
    model,
    0,
    1,
    0,
    1,
    loc_admixture,
    relatedness,
    kinship,
    pr_allele_freq
  )
  
  ## default parameter name
  parameter_name <- c(
    "Integer, #Individuals",
    "Integer, #Loci",
    "Boolean, All loci SNP (1/0=Y/N)",
    "String, Missing allele",
    "Integer, Random number seed",
    "String, Genotypefilename",
    "String, Outputfilename",
    "integer, 3/2/1/0 = strong/medium/weak/no scaling",
    "Integer, Minimum K",
    "Integer, Maximum K",
    "Integer, Num replicate runs per K",
    "Integer, 0/1=Search using assignment_prob/relatedness",
    "Boolean, 1/0=Output allele frequency:YES/NO",
    "Boolean, 1/0=PopData available:YES/NO",
    "Boolean, 1/0=PopFlag available:YES/NO",
    "Integer, 1/2/3/4=Mixture/Admixture/Hybridyzation/Migration model",
    "Boolean, 1/0=Estimate locus-specific F-Statistics=Y/N",
    "Boolean, 1/0=Use K-Means clustering method=Y/N",
    "Boolean, 1/0=Individual location data available=Y/N",
    "Integer, 0/1/2=Individual data in 1-row/2-rows/1-column",
    "Boolean, 1/0=Infer locus admixture=Y/N",
    "Boolean, 0/1/2: Compute relatedness = No/Wang/LynchRitland",
    "Boolean, 0/1 estimate kinship = No/Yes",
    "Integer, 0/1/2=Undefined/equal/unequal prior allele freq"
  )
  
  #create PARAMETER FILE
  write.table(
    cbind(
      pillar::align(parameter, align = "left"),
      paste0("!", parameter_name)
    ),
    file.path(output.path, paste0(filename, ".popcluster", ".PcPjt")),
    sep = " ",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
  )
  
  input_file <- c(paste0(filename, ".popcluster.PcPjt"),
                  paste0(filename, ".popcluster.dat"))
  file.copy(
    file.path(popcluster.path, popcluster_version),
    to = tempd,
    overwrite = TRUE,
    recursive = TRUE
  )
  file.copy(
    file.path(output.path, input_file),
    to = tempd,
    overwrite = TRUE,
    recursive = TRUE
  )
  
  old.path <- getwd()
  setwd(tempd)
  on.exit(setwd(old.path), add = TRUE)

  # PopCluster's own output only at verbose >= 3
  quiet <- verbose < 3
  if (os == "Linux" | os == "Darwin") {
    for (f in c(popcluster_version, input_file)) {
      system(paste0("chmod 777", " ", f),
             ignore.stdout = quiet, ignore.stderr = quiet)
    }
  }
  
  # RUN POPCLUSTER
  if (verbose >= 2) {
    cat(report("  Running PopCluster for K =", minK, "to", maxK, "\n"))
  }
  if (os == "Linux" & parallel) {
    system(paste0(
      "mpirun -np ", ncores, " --use-hwthread-cpus ",
      file.path(tempd, popcluster_version),
      " INP:",
      paste0(filename, ".popcluster.PcPjt MPI:1")
    ), ignore.stdout = quiet, ignore.stderr = quiet)
  } else {
    system(paste0(
      file.path(tempd, popcluster_version[1]),
      " INP:",
      paste0(filename, ".popcluster.PcPjt")
    ), ignore.stdout = quiet, ignore.stderr = quiet)
  }

  k_file <- file.path(tempd, paste0(filename, ".popcluster.K"))
  if (!file.exists(k_file)) {
    stop(error(
      "PopCluster did not produce its results; rerun with verbose = 3 to see",
      "its output.\n"
    ))
  }

  # Summarise best run and likelihood
  res <- readLines(con <- file(k_file), n = ((maxK - minK) + 2))[-1]
  close(con)
  res2 <- stringr::str_split(gsub('\"', "", res), " ")
  for (i in 1:length(res2)) {
    res2[[i]][which(res2[[i]] == "")] <- NA
    res2[[i]] <- na.omit(res2[[i]])
  }
  best_run_file <- NULL
  for (j in 1:length(res2)) {
    best_run_file <- data.frame(rbind(best_run_file, res2[[j]]))
  }
  colnames(best_run_file) <- c("K",
                               "BestRun",
                               "LogL_Mean",
                               "LogL_Min",
                               "LogL_Max",
                               "DLK1",
                               "DLK2",
                               "FST.FIS")
  # numbers are read as text; PopCluster writes "-" where a value is undefined
  num_cols <- setdiff(colnames(best_run_file), "BestRun")
  best_run_file[num_cols] <- lapply(best_run_file[num_cols], function(v) {
    v[v == "-"] <- NA
    as.numeric(v)
  })
  
  # plot likelihood and related statistics against K
  plot_stat <- function(col) {
    d <- best_run_file[!is.na(best_run_file[[col]]), , drop = FALSE]
    ggplot2::ggplot(d, aes(x = .data$K, y = .data[[col]])) +
      geom_line() +
      geom_point(fill = "white", shape = 21, size = 3) +
      scale_x_continuous(breaks = best_run_file$K) +
      ylab(col) +
      plot_theme +
      theme(axis.title.x = element_blank())
  }
  plot.list <- lapply(c("LogL_Mean", "DLK1", "DLK2", "FST.FIS"), plot_stat)
  names(plot.list) <- c("LogL_Mean", "DLK1", "DLK2", "FST.FIS")

  if (plot.out) {
    print(patchwork::wrap_plots(plot.list, ncol = 2) +
            patchwork::plot_annotation(
              caption = "K",
              theme = theme(plot.caption = element_text(hjust = 0.5))
            ))
  }
  
  #extract admixture analysis from best run
  Q_matrices <- NULL
  Q <- NULL
  
  for (i in best_run_file$BestRun) {
     
    if(abs(minK - maxK)== 0){
      i <- best_run_file$BestRun[1]
    }
    
    best <- readLines(con <- file(file.path(tempd, i)))
    close(con)
    hdr <- which(startsWith(best, "Inferred ancestry of individuals"))
    anc_lines <- best[(hdr + 2):(hdr + 1 + nInd(x))]

    # Each ancestry line has the layout
    #   Index Order <PopName...> %Miss Cluster : prop_1 prop_2 ... prop_K
    # The <PopName> column is the genlight population label and may contain
    # spaces (e.g. "Upper Murray") or wide multi-byte characters (e.g. the
    # en-dash in "Macquarie-Castlereagh"), so the number of whitespace tokens
    # varies between individuals. Anchor on the ":" and read Index/Order from
    # the front and %Miss/Cluster from the back, so parsing never depends on
    # the (variable) token count of the middle. The previous approach split on
    # spaces and used fixed positional columns, which misaligned these rows and
    # recycled the Index into the last cluster, producing ancestry values > 1.
    Q_rows <- lapply(anc_lines, function(line) {
      sp    <- strsplit(line, ":", fixed = TRUE)[[1]]
      left  <- strsplit(trimws(sp[1]), "\\s+")[[1]]
      props <- as.numeric(strsplit(trimws(sp[2]), "\\s+")[[1]])
      n <- length(left)
      data.frame(
        Index       = as.integer(left[1]),
        Order       = as.integer(left[2]),
        PercentMiss = left[n - 1],
        Cluster     = left[n],
        t(props),
        stringsAsFactors = FALSE
      )
    })
    Q <- do.call(rbind, Q_rows)
    colnames(Q) <- c("Index",
                     "Order",
                     "PercentMiss",
                     "Cluster",
                     paste0("Pop_", seq_len(ncol(Q) - 4)))
    Q$Label <- sample_name[Q$Index]
    Q$Cluster <- as.character(Q$Cluster)
    # population by the same individual index as the label
    Q$Pop <- as.character(x$pop)[Q$Index]
    # restore the documented column order
    Q <- Q[, c("Index", "Order", "Label", "PercentMiss", "Cluster",
               paste0("Pop_", seq_len(sum(startsWith(names(Q), "Pop_")))),
               "Pop")]
    # change the Order
    Q <- Q[with(Q, order(Q$Pop, as.numeric(Q$Cluster))), ]
    Q$Order <- 1:nrow(Q)
    Q_matrices[[i]] <- Q
    Q <- NULL
  }
  
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(plot.list,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }

  setwd(old.path)
  if (cleanup) {
    unlink(tempd, recursive = TRUE)
  }
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # return all Q matrices and best run summary
  return(
    list(
      output_path = output.path,
      best_run = best_run_file,
      plots = plot.list,
      matrix = Q_matrices
    )
  )
}
