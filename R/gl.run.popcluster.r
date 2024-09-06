#' @name gl.run.popcluster
#' @title Creates an input file for the program PopCluster and runs it if
#' PopCluster is installed (can be installed at: https://www.zsl.org/about-zsl/resources/software/popcluster)
#' @description
#' This function infers the population admixture
#'
#' If you specify a directory for the PopCluster executable file, then the
#' script will create the input file (DataForm=0) from the SNP data then run PopCluster.
#'
#' Refer to the PopCluster manual for further information on the parameters to
#' set
#' -- ##########################################
#'
#'
#'@param x Name of the genlight object containing the SNP data [required].
# @param outfile Name of the file that will be the input file for PopCluster
# [default popcluster.txt].
#' @param popcluster.path absolute path to the directory that contain the PopCluster program 
#' (eg: c:/User/bin/PopCluster/)
#' @param filename file name of output data 
#' @param minK Minimum K
#' @param maxK Maximum K 
#' @param rep Number of replicates runs per K
#' @param search_relate Search using assignment_prob/relatedness: 0=N, 1=Y [default 0]
#' @param allele_freq Output allele frequency: 0=N, 1=Y [default 1]
#' @param PopData PopData available: 0=N, 1=Y
#' @param PopFlag PopFlag available: 0=N, 1=Y. [default 0]
#' @param model 1=Mixture, 2=Admixture, 3=Hybridyzation, 4=Migration model [default 2]
#' @param location Individual location data available: 0=N, 1=Y 
#' @param loc_admixture Infer locus admixture: 0=N, 1=Y [default 1]
#' @param relatedness Compute relatedness = 0=No, 1=Wang, 2=LynchRitland [default 0]
#' @param kinship Estimate kinship: 0=N, 1=Y [default 0]
#' @param pr_allele_freq 0=Undefined, 1=equal, 2=unequal prior allele freq [default 2]
#' @return The results of PopCluster
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @export
#' @importFrom pillar align
#' @references
#' \itemize{
#' \item Wang, J. (2022). Fast and accurate population admixture inference 
#' from genotype data from a few microsatellites to millions of SNPs. Heredity, 129(2), 79-92.
#' }
#' @author Custodian: Ching Ching Lau -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' m <- gl.run.popcluster(x=testset.gl, popcluster.path="/User/PopCluster/Bin/",
#' output.path=/User/Documents/Output/,
#'   filename="prefix", minK=1, maxK=3, 
#'   rep=10, PopData=1, location=1
#' )
#' }
#' Wrapper function to run PopCluster
#' 
#' @export 


gl.run.popcluster <- function(x, popcluster.path, output.path, filename, minK, maxK, 
                              rep, search_relate=0, allele_freq=1,PopData, PopFlag=0,
                              model=2, location, loc_admixture=1, relatedness=0, 
                              kinship=0, pr_allele_freq=2, cleanup=TRUE, 
                             # plot.display=TRUE,
                             # plot.out = TRUE,
                             # plot_theme = theme_dartR(),
                             # plot.file = NULL,
                              verbose=NULL) 
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
  
  
  # check OS
  os <- Sys.info()['sysname'] 
  if (os == "Windows") {
    
      popcluster_version <- paste0("PopCluster", "Win.exe")} else if (os == "Darwin"){
    
          popcluster_version <- paste0("PopCluster", "Mac")} else if (os == "Linux") {
              popcluster_version <- paste0("PopCluster", "Lnx") }
  
  
  #TODO: create INPUT FILE
  genotype <- as.matrix.genlight(x)
  genotype[is.na(genotype)] <- 3
  sample_name <- x@ind.names
  family <- x@pop
  rownames(genotype) <- NULL
  if (location==1){
    lat <- x@other$latlon$lat
    lon <- x@other$latlon$lon 
    names <- data.frame(id = paste0(sample_name, " ",family, " ",PopFlag, " ", lat, " ", lon))} else if (location==0) {
      lat <- NULL
      lon <- NULL
      names <- data.frame(id = paste0(sample_name, " ",family, " ",PopFlag))
    }
  names2 <- apply(names, 1, paste0, collapse = " ")
  genotype2 <- apply(genotype, 1, paste0, collapse = "")
 
   writeLines(capture.output(
    for (i in 1:nInd(x)){
    cat(names2[i],genotype2[i], sep = "\n")}),
    con=file.path(output.path, paste0(filename,".popcluster.dat")))
  
  # parameter from user input
  parameter <- c(nInd(x), nLoc(x),1, 0, 333,paste0(filename,".popcluster.dat") , 
                 paste0(filename,".popcluster"), 
                 0, minK, maxK, 
                 rep, search_relate, allele_freq, PopData, PopFlag, 
                 model, 1, 1, location, 1, loc_admixture, relatedness, kinship, pr_allele_freq)
  
  ## default parameter name
  parameter_name <- c("Integer, #Individuals", 
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
                      "Integer, 0/1/2=Undefined/equal/unequal prior allele freq")
  
  
  #create parameter file
  create_parameter<-file(file.path(output.path, paste0(filename,".popcluster",".PcPjt")))
  write.table(cbind(pillar::align(parameter,align="left"), paste0("!",parameter_name)), 
              file.path(output.path, paste0(filename,".popcluster",".PcPjt")),sep=" ",
              quote=F, col.names = F, row.names = F)
  close(create_parameter)
  
  # check file existence
  input_file <- c(paste0(filename,".popcluster.PcPjt"), paste0(filename,".popcluster.dat"))
  fex <- file.exists(file.path(popcluster.path, popcluster_version))
  fex2 <- file.exists(file.path(output.path, input_file))
  
  if (all(fex) & all(fex2)) {
    file.copy(file.path(popcluster.path, popcluster_version),
              to = tempd,
              overwrite = TRUE, recursive = TRUE)
    file.copy(file.path(output.path, input_file),
             to = tempd,
          overwrite = TRUE, recursive = TRUE)
  } else{
    cat("  Cannot find",
        progs[!fex],
        "in the specified folder given by popcluster.path:",
        popcluster.path,
        "\n")
    stop()
  }
  
  old.path <- getwd()
  setwd(tempd)
  on.exit(setwd(old.path))
  if (os=="Linux"|os=="Darwin") {
    system(paste0("chmod 777", " ", popcluster_version))
    system(paste0("chmod 777", " ", paste0(filename,".popcluster.Pcpjt")))
           system(paste0("chmod 777", " ", paste0(filename,".popcluster.dat")))
  }
  
  #SET PATH TO RUN POPCLUSTER
  #system("export DYLD_LIBRARY_PATH=/usr/local/opt/gcc/lib/gcc/11:/usr/local/homebrew/lib/gcc/14")
  system(paste0(file.path(popcluster.path,popcluster_version), " INP:", paste0(filename,".popcluster.PcPjt")))
  
  # SET WORKING DIRECTORY
  # Select file to save and plot later
  list_of_files <- list.files(tempd, filename) 
  file.copy(file.path(tempd, list_of_files),
            to = output.path,
            overwrite = F, recursive = TRUE)
  setwd(old.path)
  res <- readLines(file.path(output.path, paste0(filename,".popcluster.K")))[0:maxK+1]
  write.table(res, file.path(output.path, paste0(filename,".popcluster.best_run_summary")), quote = F, row.names = F, col.names = F)
 
  best_run <- read.table(file.path(output.path,paste0(filename,".popcluster.best_run_summary")), header = T)
  plot.list <- list()
  plot.list[[1]] <- ggplot2::ggplot(best_run, aes(K, LogL_Mean, group=1)) + geom_line() + geom_point(fill = "white", shape = 21,
                    size = 3) + theme(axis.title.x = element_blank())
  plot.list[[2]] <- ggplot2::ggplot(best_run, aes(K, DLK1, group=1)) + geom_line() + geom_point(fill = "white", shape = 21,
                    size = 3) + theme(axis.title.x = element_blank())
  plot.list[[3]] <- ggplot2::ggplot(best_run, aes(K, DLK2, group=1)) + geom_line() + geom_point(fill = "white", shape = 21,
                    size = 3) + theme(axis.title.x = element_blank())
  plot.list[[4]] <- ggplot2::ggplot(best_run, aes(K, FST.FIS, group=1)) + geom_line() + geom_point(fill = "white", shape = 21,
                    size = 3) + theme(axis.title.x = element_blank())
  names(plot.list) <- c("LogL_Mean", "DLK1", "DLK2", "FST.FIS")
  
  p <- plot.list %>% purrr::map(function(x) {
    ggplot2::ggplot_gtable(ggplot2::ggplot_build(x))
  })
  maxWidth <- do.call(grid::unit.pmax, purrr::map(p, function(x) x$widths[2:3]))
  for (i in 1:length(p)) p[[i]]$widths[2:3] <- maxWidth
  p$bottom <- "K"
  p$ncol <- 2
  do.call(gridExtra::grid.arrange, p)
  return(list(df = best_run, plots = plot.list))

  if (cleanup) unlink(tempd, recursive = T)
}
