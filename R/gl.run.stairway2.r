#' Run Stairway Plot 2 for Demographic History Inference
#'
#' This function runs Stairway Plot 2 to infer demographic history using folded SNP frequency spectra.
#' Stairway Plot 2 is a method for inferring demographic history using folded SNP frequency spectra. The key features and methodology of Stairway Plot 2 include:
#' \itemize{
#'   \item \strong{Folded SNP Frequency Spectra}: The method uses folded SNP frequency spectra, which are less sensitive to errors in ancestral state inference compared to unfolded spectra.
#'   \item \strong{Demographic Inference}: By analyzing the SNP frequency spectra, Stairway Plot 2 can infer changes in population size over time, providing insights into historical demographic events.
#'   \item \strong{Bootstrap Replicates}: The method employs bootstrap replicates to estimate confidence intervals for the inferred demographic history, ensuring robust and reliable results.
#'   \item \strong{Flexible Modeling}: Stairway Plot 2 allows for flexible modeling of demographic history without assuming a specific parametric form for population size changes.
#' }
#' To be able to run Stairway Plot 2, the binaries need to be provided in a single folder and can be downloaded via the \link[dartRverse]{gl.download.binary} function. In this case your system need to have Java installed as well. for more details on the method and how to install on your system refer to the githubh repository: \url{https://github.com/xiaoming-liu/stairway-plot-v2}. Please also refer to the original publication for more details on the method: \url{https://doi.org/10.1186/s13059-020-02196-9}. **Also if you use this method, make sure you cite the original publication in your work.**
#' This function implements the theoretical and computational procedures described by Liu and Fu (2020), making it suitable for a wide range of population-genomic datasets to uncover historical demographic patterns.
#' Please note: There is currently not really a good way to estimate L, the length 
#' of all sequences. Often users of dart data use the number of loci multiplied 
#' by 69, but this is definitely an underestimate as monomorphic loci need to be 
#' included (also the length of the restriction site should be added for each loci).
#' For mutation rate u, the default value is set to 5e-9, but should be adapted 
#' to the species of interest. The good news is, that settings of L and mu affects 
#' only the axis of the inferred history, but not the shape of the history. 
#' So users can infer the shape, but need to be careful with a temporal interpretation 
#' as both x and y axis are affected by the mutation rate and L.
#' @param x A genlight/dartR object containing SNP data.
#' @param L the length of the sequence in base pairs. (see notes below)
#' @param mu the mutation rate per base pair per generation. (see notes below)
#' @param stairway2.path the path to the Stairway Plot 2 executable. (check the example)
#' @param minbinsize the minumum bin size for the SFS that should be used. (default=1)
#' @param maxbinsize the maximum bin size for the SFS that should be used. (default=NULL, 
#' so the maximum bin size is set to the number of samples in the dataset)
#' @param gentime the generation time in years. (default=1)
#' @param sfs the folded site frequency spectrum (SFS) to be used for the analysis. 
#' If not provided the SFS is created from the genlight/dartR object (default=NULL)
#' @param parallel the number of parallel processes to use for the analysis. (default=1)
#' @param run logical. If TRUE, the analysis is run immediately. Otherwise only the
#' blueprint files are created [might be useful to run on a cluster]. (default=FALSE)
#' @param blueprint the name of the blueprint file. (default="blueprint")
#' @param filename the name of the filename. Also used for the plot. (default="sample")
#' @param pct_training the percentage of the data to use for training. (default=0.67)
#' @param nrand the number of breakpoint to use for the analysis. (default=NULL)
#' @param stairway_plot_dir the name of the directory where the stairway plot is saved. 
#' (default="stairway_plot_es")
#' @param nreps the number of bootstrap replicates to use for the analysis. (default=200)
#' @param seed the random seed to use for the analysis. (default=NULL)
#' @param plot_title the title of the plot. (default="Ne"+filename)
#' @param xmin minimum x value for the plot. (default=0)
#' @param xmax maximum x value for the plot. (default=0)
#' @param ymin minimum y value for the plot. (default=0)
#' @param ymax maximum y value for the plot. (default=0)
#' @param xspacing spacing between x values for the plot. (default=2)
#' @param yspacing spacing between y values for the plot. (default=2)
#' @param fontsize the font size for the plot. (default=12)
#' @param cleanup logical. If TRUE, the stairway 2 plot output files are removed. 
#' (default=TRUE)
#' @param plot.display Specify if plot is to be produced [default TRUE].
#' @param plot.theme User specified theme [default theme_dartR()].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]
#' @param plot.file Filename (minus extension) for the RDS plot file [Required for plot save]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return returns a list with two components: 
#' \itemize{
#' \item{history: Ne estimates of over generations (generation, median, low and high)} 
#' \item{plot: a ggplot of history }
#' }
#'
#' @references Liu, X., & Fu, Y. X. (2020). Stairway Plot 2: demographic history inference with folded SNP frequency spectra. Genome Biology, 21(1), 280.
#' @importFrom parallel detectCores
#' @importFrom future plan
#' @importFrom furrr future_map
#' @references Liu, X., Fu, YX. Stairway Plot 2: demographic history inference with folded SNP frequency spectra. Genome Biol 21, 280 (2020). \doi{10.1186/s13059-020-02196-9}

#' @export
#' @examples
#' \dontrun{
#' #download binary, if not already installed, to tempdir()
#' gl.download.binary(software="stairway2",os="windows")
#' require(dartR.data)
#' sw<- gl.run.stairway2(possums.gl[1:50,1:100], L=1e5, mu = 1e-9, 
#'            stairway2.path = file.path(tempdir(),"stairway2"), 
#'            parallel=5, nreps = 10)
#' head(sw$history)
#' }

gl.run.stairway2 <- 
  function(x, 
           L = NULL,
           mu=NULL,
           stairway2.path,
           minbinsize=1,
           maxbinsize=NULL,
           gentime=1,
           sfs=NULL, 
           parallel=1,
           run=TRUE,
           blueprint="blueprint", 
           filename="sample",
           pct_training=0.67,
           nrand=NULL,
           stairway_plot_dir="stairway_plot_es",
           nreps=200,
           seed=NULL,
           plot_title="Ne",
           xmin=0, xmax=0, ymin=0, ymax=0,
           xspacing=2,
           yspacing=2,
           fontsize=12,
           cleanup=TRUE,
           plot.display=TRUE,
           plot.theme = theme_dartR(),
           plot.dir=NULL,
           plot.file=NULL,
           verbose=NULL) {
    
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # SET WORKING DIRECTORY
    plot.dir <- gl.check.wd(plot.dir,verbose=0)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    #utils.flag.start(func = funname,
    #                 build = "Jody",
    #                 verbosity = verbose)
    
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    
    # RUN STAIRWAY PLOT 2
    tempd <-  tempdir()
    dir.create(tempd, showWarnings = FALSE)
    #blueprint <- paste0("blueprint")
    outfilespec <- file.path(tempd, blueprint)
    
    # check OS
    os <- tolower(Sys.info()['sysname']) 
    #create individual tempdir
    
    progs <- c("stairway_plot_es")
    fex <- file.exists(file.path(stairway2.path, progs))
    if (all(fex)) {
      file.copy(file.path(stairway2.path, progs),
                to = tempd,
                overwrite = TRUE, recursive = TRUE)
      stairway2.path <- tempd
    } else{
      cat("  Cannot find",
          progs[!fex],
          "in the specified folder given by stairway2.path:",
          stairway2.path,
          "\n"      )
      stop()
    }
    
    
    
    
    # SET VERBOSITY
    
    if (is.null(verbose)){ 
      if(!is.null(x@other$verbose)){ 
        verbose <- x@other$verbose
      } else { 
        verbose <- 2
      }
    } 
    
    if (verbose < 0 | verbose > 5){
      cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
      verbose <- 2
    }
    
    # FLAG SCRIPT START
    
    if (verbose >= 1){
      if(verbose==5){
        cat("Starting",funname,"[ Build =",build,"]\n")
      } else {
        cat("Starting",funname,"\n")
      }
    }
    
    # STANDARD ERROR CHECKING
    
    if(!is(x,"genlight")) {
      stop("  Fatal Error: genlight object required!\n")
    }
    
    if (verbose >= 2){
      if (all(x@ploidy == 1)){
        stop("Fatal Error: Detected Presence/Absence (SilicoDArT) data. Please provide a SNP dataset\n")
      } else if (all(x@ploidy == 2)){
        cat("  Processing a SNP dataset\n")
      } else {
        stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
      }
    }
    
    # SCRIPT SPECIFIC ERROR CHECKING
    
    if(is.null(maxbinsize)){
      maxbinsize <- nInd(x)
      if(verbose >= 3){cat("  Max Bin Size not specified, set to",nInd(x),"\n")}
    }
    if(is.null(nrand))  nrand <- c(round((nInd(x)-1)/2), round(nInd(x)-1), round((nInd(x)-1)*3/2), round(2*(nInd(x)-1))) else nrand <- round(seq(0.5,2,len=nrand)*(nInd(x)-1))
    
    if(verbose >= 3){cat("  No. of break points, set to:",paste0(nrand, collapse=" "),"\n")}
    
    if(is.null(stairway_plot_dir)){
      stop("Fatal Error: Directory path for the Stairway Plot 2 executables not specified\n")
    }
    if(is.null(mu)){
      stop("Fatal Error: Mutation rate per site per generation not specified\n")
    }
    if(is.null(gentime)){
      stop("Fatal Error: Generation time (years) not specified\n")
    }
    whether_folded <- "true"
    nseq <- 2*nInd(x)
    
    
    
    #usfs <- colSums(as.matrix(x), na.rm=T)
    
    #fold usfs (in case of uneven individuals, drop the 0s and minbinsize)
    #usfs <- ifelse(usfs>nInd(x), nInd(x)*2-usfs, usfs)
    #sfs <- table(usfs)[-(1:minbinsize)]
    
    #fill in empty bins in sfs with zeros
    #if (length(sfs) != (nInd(x)-(minbinsize-1) )) {
    
    #sfs2 <- rep(0,(nInd(x)-(minbinsize-1) ))
    #sfs2[as.numeric(names(sfs))]  <- sfs
    #sfs <- sfs2
    
    #}
    if (is.null(sfs)) sfs <- gl.sfs(x, minbinsize=1, plot.out=FALSE, singlepop = TRUE)
    
    #cut at maxbinsize
    #sfs <- sfs[1:maxbinsize] (not necessary)
    
    sfs <- paste(sfs,collapse = " ")
    
    #if total length of sequence is not specificed simply assume nLoc*69 (standard from dart)
    if (is.null(L)) L = nLoc(x)*69 
    
    if (is.null(seed)) {set.seed= as.numeric(Sys.time()) ;seed = round(runif(1)*1e6)}
    
    # DO THE JOB
    
    if (verbose >= 2) {cat(paste("  Extacting SNP data and creating records for each individual\n"))}
    
    # Output the results
    
    write.table(paste("#Ne analysis from SNPs for ",as.character(substitute(x))),
                file=outfilespec,row.names=FALSE,col.names=FALSE,
                quote=FALSE)
    write.table(paste("popid:",filename,"# id of the population (no white space)"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("nseq:",nseq,"# number of haploid sequences = 2n"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("L:",L,"# total number of nucleic sites, including polymorphic and monomorphic"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("whether_folded:", whether_folded, "# whethr the SFS is folded (true or false)"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("SFS:", sfs, "# snp frequency spectrum: number of singleton, number of doubleton, etc. (separated by white space)"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("smallest_size_of_SFS_bin_used_for_estimation:", minbinsize, "# default is 1; to ignore singletons, change this number to 2"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("largest_size_of_SFS_bin_used_for_estimation:", maxbinsize, "# default is nseq/2 for folded SFS"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("pct_training:", pct_training, "# proportion of sites for training"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("nrand:", paste0(nrand, collapse=" "), "# number of random break points for each try (separated by white space)"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("project_dir:", "files", "# project directory"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("stairway_plot_dir:", stairway_plot_dir, "# directory to the stairway plot files"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("ninput:", nreps, "# number of input files to be created for each estimation"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("random_seed:", seed),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("mu:", mu, "# assumed mutation rate per site per generation"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("year_per_generation:", gentime,"# assumed generation time (in years)"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    plottitle <- plot_title
    write.table(paste("plot_title:",plottitle , "# title of the plot"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("xrange:", paste0(xmin,',',xmax), "# Time (1k year) range; format: xmin,xmax; 0,0 for default"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("yrange:", paste0(ymin,",",ymax), "# Ne (1k individual) range; format: ymin,ymax; 0,0 for default"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("xspacing:", xspacing, "# X axis spacing"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("yspacing:",yspacing, "# Y axis spacing"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    write.table(paste("fontsize:", fontsize, "# Font size"),
                file=outfilespec,
                row.names=FALSE,col.names=FALSE,
                quote=FALSE, sep=" ",append=TRUE)
    
    if (verbose > 2) {cat(paste("    Stairway Plot 2 blueprint written to",outfilespec,"\n"))}
    
    #set blueprint.sh to executable under linux
    if (os=="linux" | os=="darwin") system(paste0("chmod 777 ",outfilespec)) 
    # FLAG SCRIPT END
    
    if (verbose > 0) {
      cat("Completed:",funname,"\n")
    }
    
    if (run==TRUE)
    {
      oldpath <- getwd()
      setwd(stairway2.path)
      on.exit(setwd(oldpath))
      
      system(paste0("java -cp stairway_plot_es Stairbuilder ",blueprint ))
      #make .sh file executable...
      if (os=="linux" | os=="darwin") system(paste0("chmod 777 ",blueprint, ".sh"))
      
      
      #run on multiple cores
      if (parallel>1)
      {
        if (os=="windows") ff <- readLines(paste0(blueprint, ".bat"))
        if (os=="linux" | os=="darwin") ff <- readLines(paste0(blueprint, ".sh"))
        no_cores <- min(parallel, detectCores()-1)
        future::plan(future::multisession, workers = no_cores)
        
        index = grep("Stairway_fold_training_testing7", ff)
        
        runs <- ff[index]
        
        runstair <- function(i) {
          system(runs[i])
          return(i)
        }
        
        temp <- furrr::future_map(1:length(runs), function(x) runstair(x))
        
        
        if (os=="windows") index = grep("MOVE", ff)
        if (os=="linux" | os=="darwin") index = grep("mv -f",ff)
        runs <- ff[index]
        if (os=="windows") runs <- gsub("MOVE /y", "cp " ,runs)
        
        for (i in 1:length(runs)) system(runs[i])
        index = grep("Stairpainter", ff)
        system( ff[index])
        er <- {
          if (os=="linux" | os=="darwin") system(paste0("bash ",blueprint, ".plot.sh"))
          if (os=="windows") system(paste0(blueprint, ".plot.bat"))
        }
        
      } else er <- {
        if (os=="linux" | os=="darwin") {
          system(paste0("./",blueprint, ".sh"))
        }
        if (os=="windows") system(paste0(blueprint, ".bat"))
      }
      
      
      
      if (length(er)>0)
      {
        cat("Attempt to rerun last step with different settings (lower memory allocation fo the Java Virtual Machine")
        ff <- readLines(paste0(blueprint, ".plot.bat"))
        ff <-  gsub("-Xmx4g", "-Xmx1g", ff)
        con <- file(paste0(blueprint, ".plot.bat"),"w")
        writeLines(ff, con)
        close(con)
        
        if (os=="linux" | os=="darwin") system(paste0("chmod 777 ",blueprint, ".plot.sh"))
        if (os=="linux" | os=="darwin") system(paste0("./",blueprint, ".plot.sh"))
        if (os=="windows") system(paste0(blueprint, ".plot.bat"))
      }
      cat("Check plots (pdf and png files) in folder:", stairway2.path,".\n")
      
     res <- read.csv(file.path(tempd,"files",paste0(plottitle,".final.summary")),sep="\t") 
      setwd(oldpath)
      
      # PLOT
      if (!is.null(res)) 
      {
        
      mutation_per_site <- n_estimation <- theta_per_site_median <- theta_per_site_2.5 <- theta_per_site_97.5 <- year <- Ne_median <- low95 <- high95 <- low75 <- high75 <- NULL
  
      colnames(res) <- c("mutation_per_site" ,"n_estimation", "theta_per_site_median", "theta_per_site_2.5","theta_per_site_97.5", "year" , "Ne_median" ,"low95" ,"high95"," low75" ,"high75")
      
      p1 <- ggplot(res, aes(x=year, y=Ne_median))+geom_point()+geom_line()+geom_ribbon(aes(ymin=low95, ymax=high95), alpha=0.2)+ylab("Effective population size")+xlab("year")#+plot.theme
      
      
      if (cleanup) unlink(file.path(tempd), recursive = TRUE)    
      # PRINTING OUTPUTS
      if (plot.display) {print(p1)}
      #print(epp)
      
      if(!is.null(plot.file)){
        tmp <- utils.plot.save(p1,
                               dir=plot.dir,
                               file=plot.file,
                               verbose=verbose)
      }
      } else {res<-NULL;p1=NULL}  #if only blueprints are created, return NULL
      
      
    }
   
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    return(list(history=res, plot=p1))
    
  }
