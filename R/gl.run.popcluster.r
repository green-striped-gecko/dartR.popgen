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
#'@param gl Name of the genlight object containing the SNP data [required].
# @param outfile Name of the file that will be the input file for PopCluster
# [default popcluster.txt].
#' @param popcluster.path absolute path to the directory that contain the PopCluster program 
#' (eg: c:/User/bin/PopCluster/)
#' @param filename file name of output data 
#' @param minK Minimum K
#' @param maxK Maximum K 
#' @param rep Number of replicates runs per K
#' @param search_relate Search using assignment_prob/relatedness: 0=N, 1=Y 
#' @param allele_freq Output allele frequency: 0=N, 1=Y
#' @param PopData PopData available: 0=N, 1=Y
#' @param PopFlag PopFlag available: 0=N, 1=Y
#' @param model 1=Mixture, 2=Admixture, 3=Hybridyzation, 4=Migration model
#' @param location Individual location data available: 0=N, 1=Y 
#' @param loc_admixture Infer locus admixture: 0=N, 1=Y 
#' @param relatedness Compute relatedness = 0=No, 1=Wang, 2=LynchRitland
#' @param kinship Estimate kinship: 0=N, 1=Y 
#' @param pr_allele_freq 0=Undefined, 1=equal, 2=unequal prior allele freq
#' @param K number of K to be saved and plotted
#' @return The reduced genlight object, if parentals are provided; output of
#'  PopCluster is saved to the working directory.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @export
#' @importFrom pillar align
#' @references Wang, J. (2022). Fast and accurate population admixture 
#' inference from genotype data from a few microsatellites to millions of SNPs. Heredity, 129(2), 79-92.
#' @author Custodian: Ching Ching Lau -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' m <- gl.run.popcluster(testset.gl, popcluster.path,
#'   filename, minK, maxK, 
#'   rep, search_relate, allele_freq,PopData,PopFlag,
#'   model, location, loc_admixture, relatedness,
#'   kinship, pr_allele_freq, K, cleanup=TRUE, verbose=NULL
#' )
#' }
#' Wrapper function to run PopCluster
#' 
#' @export 


gl.run.popcluster <- function(gl, popcluster.path, filename, minK, maxK, 
                              rep, search_relate, allele_freq,PopData,PopFlag,
                              model, location, loc_admixture, relatedness, 
                              kinship, pr_allele_freq, K, cleanup=TRUE, 
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
  datatype <- utils.check.datatype(gl, verbose = verbose)
  
  
  #create tempdir
  tempd <-  tempfile(pattern = "dir")
  dir.create(tempd, showWarnings = FALSE)
  
  
  # check OS
  os <- Sys.info()['sysname'] 
  if (os == "Windows") {
    if (verbose == 2) {
      popcluster_version <- paste0("PopCluster", "Win")}} else if (os == "Darwin"){
        if (verbose == 2) {
          popcluster_version <- paste0("PopCluster", "Mac")}} else if (os == "Linux") {
            if (verbose == 2) {
              popcluster_version <- paste0("PopCluster", "Lnx") }} 
  
  
  #TODO: create INPUT FILE
  genotype <- as.matrix(gl)
  genotype[is.na(genotype)] <- 3
  sample_name <- gl@ind.names
  family <- gl@pop
  rownames(genotype) <- NULL
  if (location==1){
    lat <- gl@other$latlon$lat
    lon <- gl@other$latlon$lon } else if (location==0) {
      lat <- NULL
      lon <- NULL
    }
  names <- data.frame(id = paste0(sample_name, " ",family, " ", lat, " ", lon))
  names2 <- apply(names, 1, paste0, collapse = " ")
  genotype2 <- apply(genotype, 1, paste0, collapse = "")
 
   writeLines(capture.output(
    for (i in 1:nInd(gl)){
    cat(names2[i],genotype2[i], sep = "\n")}),
    con=paste0(popcluster.path, filename,".popcluster.dat"))
  
  # parameter from user input
  parameter <- c(nInd(platypus.gl), nLoc(platypus.gl),1, 3, 333,paste0(filename,".popcluster.dat") , 
                 paste0(filename,".popcluster"), 
                 0, minK, maxK, 
                 rep, search_relate, allele_freq, PopData, PopFlag, 
                 model, 1, 1, location, 1, loc_admixture, relatedness, kinship, pr_allele_freq)
  
  ## default parameter name
  parameter_name <- c("Integer, #Individuals", "Integer, #Loci", "Boolean, All loci SNP (1/0=Y/N)", 
                      "String, Missing allele", "Integer, Random number seed", "String, Genotypefilename", 
                      "String, Outputfilename", "integer, 3/2/1/0 = strong/medium/weak/no scaling",                
                      "Integer, Minimum K", "Integer, Maximum K", "Integer, Num replicate runs per K", 
                      "Integer, 0/1=Search using assignment_prob/relatedness", "Boolean, 1/0=Output allele frequency:YES/NO",                     
                      "Boolean, 1/0=PopData available:YES/NO","Boolean, 1/0=PopFlag available:YES/NO",                           
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
  create_parameter<-file(paste0(popcluster.path, "/", filename,".popcluster",".PcPjt"))
  write.table(cbind(pillar::align(parameter,align="left"), paste0("!",parameter_name)), 
              paste0(popcluster.path, "/", filename,".popcluster",".PcPjt"),sep=" ",
              quote=F, col.names = F, row.names = F)
  close(create_parameter)
  
  # check file existence
  progs <- c(popcluster_version, paste0(filename,".popcluster.PcPjt"), paste0(filename,".popcluster.dat"))
  fex <- file.exists(file.path(popcluster.path, progs))
  if (all(fex)) {
    file.copy(file.path(popcluster.path, progs),
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
#  popcluster.path <- tempd
  
  old.path <- getwd()
  setwd(tempd)
  on.exit(setwd(old.path))
  if (os=="Linux"|os=="Darwin") {
    system(paste0("chmod 777", " ", popcluster_version))
    system(paste0("chmod 777", " ", paste0(filename,".popcluster.Pcpjt")))
           system(paste0("chmod 777", " ", paste0(filename,".popcluster.dat")))
  }
  
  #SET PATH TO RUN POPCLUSTER
  #Sys.setenv(DYLD_LIBRARY_PATH="/usr/local/opt/gcc/lib/gcc/11:/usr/local/homebrew/lib/gcc/14")
  #system("export DYLD_LIBRARY_PATH=/usr/local/opt/gcc/lib/gcc/11:/usr/local/homebrew/lib/gcc/14")
  system(paste0("/Users/chingchinglau/Documents/dartR/dartR_popgen_cc/PopCluster/Bin/",popcluster_version, " INP:", paste0(filename,".popcluster.PcPjt")))
  
  # SET WORKING DIRECTORY
  #plot.dir <- gl.check.wd(tempd,verbose=0)
  #res <- read.delim(file.path(tempd,paste0(filename,".popcluster_K_", maxK, "_R_", rep)))
  
  # Select file to save and plot later
  file.copy(file.path(tempd, paste0(filename,".popcluster_K_", maxK, "_R_", rep)),
            to = popcluster.path,
            overwrite = F, recursive = TRUE)
  
  setwd(old.path)
  if (cleanup) unlink(tempd, recursive = T)
  #RETURN
  #return(res)
}
