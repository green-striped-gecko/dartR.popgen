#' @name gl.run.snmf
#' @title Creates an input file for the R package LEA and runs it if
#' LEA is installed
#' @description
#' Refer to the LEA manual for further information on the parameters to
#' set
#' -- ##########################################
#'
#'
#'@param x Name of the genlight object containing the SNP data [required].
#' @param filename file name of output data  [default = "output"]
#' @param minK Minimum K
#' @param maxK Maximum K 
#' @param rep Number of replicates runs per K
#' @param regularization alpha value for regularization when analyzing small dataset [default=10]
#' @param ploidy_lv ploidy level of dataset [default=2]
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot.dir Directory in which to save files [default = working directory]
#' @param plot.file Name for the RDS binary file to save (base name only, exclude
#' extension) [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return The plot of likelihood, DLK1, DLK2, FST.FIS, best run, Q-matrices of PopCluster, 
#' @export
#' @importFrom LEA snmf
#' @importFrom LEA cross.entropy
#' @references
#' \itemize{
#' \item Frichot E, Mathieu F, Trouillon T, Bouchard G, Francois O. (2014). Fast and Efficient 
#' Estimation of Individual Ancestry Coefficients. Genetics, 194(4): 973--983.
#' }
#' @author Custodian: Ching Ching Lau -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' m <- gl.run.snmf(x=testset.gl, minK=1, 
#' maxK=5, rep=10)
#' }
#' Wrapper function to run snmf from LEA
#' 
#' @export 


gl.run.snmf <- function(x=NULL, filename="output", minK=NULL, maxK=NULL, rep=NULL, regularization=10,
                        ploidy_lv=2, plot.out=TRUE, plot.dir=NULL,
                        plot.file=NULL, verbose =2) 
{
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  #utils.flag.start(func = funname,
  #                 build = "Jody",
  #                 verbosity = verbose)
  
  # CHECK DATATYPE
  datatype <- dartR.base::utils.check.datatype(x, verbose = verbose)
  
  
  #create tempdir
  tempd <-  tempfile(pattern = "dir")
  dir.create(tempd, showWarnings = FALSE)
  
  #TODO: convert genlight object to geno file
  gl2geno(x, outpath=tempd, outfile = filename)
  
  old.path <- getwd()
  
  ss <- LEA::snmf(input.file = file.path(tempd, paste0(filename, ".geno")),
             K = minK:maxK,
             entropy = TRUE,
             repetitions = rep,
             project = "new",
             alpha = regularization, ploidy = ploidy_lv)
  
  #Choose best run
  best_run <- NULL
  best_run_path <- NULL
  K_range <- minK:maxK
  for (i in 1:length(K_range)) {
  ce <-  LEA::cross.entropy(ss, K = K_range[i])
  best_run <- c(best_run, paste0("run", which.min(ce)))
  best_run_path <- c(best_run_path, (file.path(tempd, paste0(filename,".snmf"), 
                     paste0("K", K_range[i]), best_run[i])))
   }
  
  #extract Q matrices from best run
  Q_matrices <- NULL
  for (i in 1:length(K_range)){
  Q <- read.table(list.files(best_run_path[i], pattern = ".Q", full.names = T))
  colnames(Q) <- paste0("Pop_",seq(1, K_range[i]))
  Q$Label <- as.character(x$ind.names)
  Q$Pop <- as.character(x$pop)
  for (j in 1:nrow(Q)){
    Q$Cluster[j] <- sub("Pop_", "", names(which.max(Q[j,paste0("Pop_",seq(1, K_range[i]))])))
  }
  Q <- Q[with(Q, order(Q$Pop,as.numeric(Q$Cluster))),]
  Q$Order <- 1:nrow(Q)
  Q_matrices[[i]] <- Q
  }
  names(Q_matrices) <- paste0("K", K_range)
  
  # plot cross-entropy
  plot.list=list()
  plot(ss, cex = 1.2, col = "lightblue", pch = 19)
  plot.list[[1]]=recordPlot()
  names(plot.list) <- "cross-entropy"
    
    if (!is.null(plot.file)) {
      tmp <- utils.plot.save(plot.list,
                             dir = plot.dir,
                             file = plot.file,
                             verbose = verbose
      )
    }
  # reture all Q matrices and best run  
  return(list(best_run = best_run_path, cross_entropy = plot.list, matrix=Q_matrices))
  setwd(old.path)
  if (cleanup) unlink(tempd, recursive = T)
}
