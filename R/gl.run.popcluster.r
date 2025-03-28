#' @name gl.run.popcluster
#' 
#' @title Runs a PopCluster analysis using a genlight object
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
#' @param popcluster.path Path to the directory that contain the PopCluster
#' program [default getwd()].
#' @param output.path Path to store the parameter file and input
#' data [default getwd()].
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
#' @param cleanup clean data in tmp [default  TRUE].
#' @param plot.dir Directory in which to save files [default getwd()].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot.file Name for the RDS binary file to save (base name only, 
#' exclude extension) [default NULL].
#' @param plot_theme Theme of the plot [default theme_dartR()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
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
#' @return The plot of likelihood, DLK1, DLK2, FST.FIS, best run, Q-matrices 
#' of PopCluster.
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
#' @author Custodian: Ching Ching Lau -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' m <- gl.run.popcluster(x=bandicoot.gl, 
#' popcluster.path="/User/PopCluster/Bin/",
#' output.path="/User/Documents/Output/",
#' minK=1, maxK=3,
#' rep=10, PopData=1, location=1)
#' Q <- gl.plot.popcluster(pop_cluster_result=m, plot.K = 3, ind_name=T)
#' gl.map.popcluster(x = bandicoot.gl, qmat = Q)
#' # move population 4 (out of 5) 0.5 degrees to the right and populations 1
#' # 0.3 degree to the north of the map.
#' mp <- data.frame(lon=c(0,0,0,0.5,0), lat=c(-0.3,0,0,0,0))
#' gl.map.popcluster(bandicoot.gl, qmat=Q, movepops=mp)
#' }
#'
#' @export

gl.run.popcluster <- function(x,
                              popcluster.path = getwd(),
                              output.path = getwd(),
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
                              cleanup = TRUE,
                              plot.dir = NULL,
                              plot.out = TRUE,
                              plot.file = NULL,
                              plot_theme = theme_dartR(),
                              verbose = NULL) {

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(
    func = funname,
    build = "Jody",
    verbose = verbose
  )
  
  # CHECK DATATYPE
  datatype <- dartR.base::utils.check.datatype(x, verbose = verbose)
  
  pkg <- "stringr"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  pkg <- "pillar"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  #create tempdir
  tempd <- tempfile(pattern = "dir")
  dir.create(tempd, showWarnings = FALSE)
  
  if(model == 4 &&
     nPop(x) != maxK && 
     maxK != minK ){
   cat(error(
     "For migration model, K must be fixed (i.e. maxK = minK) and equal to sampling locations or known populations\n"
   )) 
    stop()
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
    popcluster_version <- paste0("PopCluster", "Lnx")
  }
  
  # create INPUT FILE
  genotype <- as.matrix.genlight(x)
  genotype[is.na(genotype)] <- 3
  sample_name <- x@ind.names
  family <- x@pop
  rownames(genotype) <- NULL
  # IndivLoc is not used in structure inference. It is used solely for 
  # visualizing population structuring in relation to individual geographic 
  # locations in PopCluster’s GUI.
  # if (location == 1) {
  #   lat <- x@other$latlon$lat
  #   lon <- x@other$latlon$lon
  #   names <- data.frame(id = paste0(sample_name, 
  #                                   " ", family, 
  #                                   " ", PopFlag,
  #                                   " ", lat, 
  #                                   " ", lon))
  # } else if (location == 0) {
  #   lat <- NULL
  #   lon <- NULL
    names <- data.frame(id = paste0(sample_name,
                                    " ", family,
                                    " ", PopFlag))
  # }
  
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
    1,
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
  create_parameter <- file(file.path(output.path, 
                                     paste0(filename, ".popcluster", ".PcPjt")))
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
  close(create_parameter)
  
  # check input files existence
  input_file <- c(paste0(filename, ".popcluster.PcPjt"),
                  paste0(filename, ".popcluster.dat"))
  fex <- file.exists(file.path(popcluster.path, popcluster_version))
  fex2 <- file.exists(file.path(output.path, input_file))
  
  if (all(fex)) {
    file.copy(
      file.path(popcluster.path, popcluster_version),
      to = tempd,
      overwrite = TRUE,
      recursive = TRUE
    )
  } else {
    cat(error(
      "  Cannot find",
      popcluster_version[!fex],
      "in the specified folder given by popcluster.path:",
      popcluster.path,
      "\n"
    ))
    stop()
  }
  
  if (all(fex2)) {
    file.copy(
      file.path(output.path, input_file),
      to = tempd,
      overwrite = TRUE,
      recursive = TRUE
    )
  } else {
    cat(error(
      "  Cannot find",
      input_file[!fex2],
      "in the specified folder given by output.path:",
      output.path,
      "\n"
    ))
    stop()
  }
  
  old.path <- getwd()
  setwd(tempd)
  on.exit(setwd(old.path))
  if (os == "Linux" | os == "Darwin") {
    system(paste0("chmod 777", " ", popcluster_version))
    system(paste0("chmod 777", " ", paste0(filename, ".popcluster.PcPjt")))
    system(paste0("chmod 777", " ", paste0(filename, ".popcluster.dat")))
  }
  
  # RUN POPCLUSTER
  system(paste0(
    file.path(tempd, popcluster_version[1]),
    " INP:",
    paste0(filename, ".popcluster.PcPjt")
  ))
  
  # Summarise best run and likelihood
  res <- readLines(con <- file(file.path(
    tempd, paste0(filename, ".popcluster.K")
  )), n = maxK + 1)[-1]
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
  K <- LogL_Mean <- LogL_Min <- LogL_Max <- DLK1 <- DLK2 <- FST.FIS <- NA
  colnames(best_run_file) <- c("K",
                               "BestRun",
                               "LogL_Mean",
                               "LogL_Min",
                               "LogL_Max",
                               "DLK1",
                               "DLK2",
                               "FST.FIS")
  
  # plot likelihood
  plot.list <- list()
  plot.list[[1]] <- ggplot2::ggplot(best_run_file, 
                                    aes(K, LogL_Mean, group = 1)) +
    geom_line() + 
    geom_point(fill = "white",shape = 21, size = 3) + 
    theme(axis.title.x = element_blank()) +
    theme_dartR()
  
  plot.list[[2]] <- ggplot2::ggplot(best_run_file, 
                                    aes(K, DLK1, group = 1)) + 
    geom_line() +
    geom_point(fill = "white",shape = 21, size = 3) + 
    theme(axis.title.x = element_blank()) +
    theme_dartR()
  
  plot.list[[3]] <- ggplot2::ggplot(best_run_file, 
                                    aes(K, DLK2, group =1)) + 
    geom_line() + geom_point(fill = "white",shape = 21, size = 3) + 
    theme(axis.title.x = element_blank()) +
    theme_dartR()
  
  plot.list[[4]] <- ggplot2::ggplot(best_run_file,
                                    aes(K, FST.FIS, group = 1)) + 
    geom_line() + 
    geom_point(fill = "white", shape = 21, size = 3) + 
    theme(axis.title.x = element_blank()) +
    theme_dartR()
  
  names(plot.list) <- c("LogL_Mean", "DLK1", "DLK2", "FST.FIS")
  
  p <- plot.list %>% purrr::map(function(x) {
    ggplot2::ggplot_gtable(ggplot2::ggplot_build(x))
  })
  maxWidth <- do.call(grid::unit.pmax, purrr::map(p, function(x)
    x$widths[2:3]))
  for (i in 1:length(p))
    p[[i]]$widths[2:3] <- maxWidth
  p$bottom <- "K"
  p$ncol <- 2
  if (plot.out) {
    do.call(gridExtra::grid.arrange, p)
  }
  names(plot.list) <- c("LogL_Mean", "DLK1", "DLK2", "FST.FIS")
  
  #extract admixture analysis from best run
  Q_matrices <- NULL
  Q <- NULL
  
  for (i in best_run_file$BestRun) {
    
    if(abs(minK - maxK)== 0){
      i <- best_run_file$BestRun[1]
    }
    best <- readLines(con <- file(file.path(tempd, i)))
    close(con)
    Q_raw <- stringr::str_split(best[(which(startsWith(
      best, "Inferred ancestry of individuals"
    )) + 2):(which(startsWith(
      best, "Inferred ancestry of individuals"
    )) + 1 + nInd(x))], " ")
    for (j in 1:length(Q_raw)) {
      Q_raw[[j]][which(Q_raw[[j]] == "")] <- NA
      Q_raw[[j]][which(Q_raw[[j]] == ":")] <- NA
      Q_raw[[j]] <- na.omit(Q_raw[[j]])
      Q <- data.frame(rbind(Q, Q_raw[[j]]))
    }
    colnames(Q) <- c("Index",
                     "Order",
                     "Label",
                     "PercentMiss",
                     "Cluster",
                     paste0("Pop_", seq(1, (ncol(
                       Q
                     ) - 5), by = 1)))
    Q$Label <- as.character(Q$Label)
    Q$Cluster <- as.character(Q$Cluster)
    Q$Pop <- as.character(x$pop)
    # change the Order
    Q <- Q[with(Q, order(Q$Pop, as.numeric(Q$Cluster))), ]
    Q$Order <- 1:nrow(Q)
    Q <- Q %>%
      mutate_at(paste0("Pop_", seq(1, (ncol(
        Q
      ) - 6), by = 1)), as.numeric)
    Q_matrices[[i]] <- Q
    Q <- NULL
  }
  
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(plot.list,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
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
  setwd(old.path)
  if (cleanup)
    unlink(tempd, recursive = T)
}
