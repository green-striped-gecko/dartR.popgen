#' @name gl.run.stairway2
#' @title Run Stairway Plot 2 for demographic history inference
#' @family demographic history
#' @description
#' This function runs Stairway Plot 2 to infer demographic history using folded SNP frequency spectra.
#' Stairway Plot 2 is a method for inferring demographic history using folded SNP frequency spectra. The key features and methodology of Stairway Plot 2 include:
#' \itemize{
#'   \item \strong{Folded SNP Frequency Spectra}: The method uses folded SNP frequency spectra, which are less sensitive to errors in ancestral state inference compared to unfolded spectra.
#'   \item \strong{Demographic Inference}: By analyzing the SNP frequency spectra, Stairway Plot 2 can infer changes in population size over time, providing insights into historical demographic events.
#'   \item \strong{Bootstrap Replicates}: The method employs bootstrap replicates to estimate confidence intervals for the inferred demographic history.
#'   \item \strong{Flexible Modeling}: Stairway Plot 2 allows for flexible modeling of demographic history without assuming a specific parametric form for population size changes.
#' }
#' To be able to run Stairway Plot 2, the binaries need to be provided in a single folder and can be downloaded via the \code{gl.download.binary} function. In this case your system needs to have Java installed as well. For more details on the method and how to install on your system refer to the github repository: \url{https://github.com/xiaoming-liu/stairway-plot-v2}. Please also refer to the original publication for more details on the method: \doi{10.1186/s13059-020-02196-9}. **Also if you use this method, make sure you cite the original publication in your work.**
#' This function implements the theoretical and computational procedures described by Liu and Fu (2020), making it suitable for a wide range of population-genomic datasets to uncover historical demographic patterns.
#' @details
#' Please note: There is currently not really a good way to estimate L, the length
#' of all sequences. Often users of dart data use the number of loci multiplied
#' by 69, but this is definitely an underestimate as monomorphic loci need to be
#' included (also the length of the restriction site should be added for each loci).
#' For mutation rate u, the default value is set to 5e-9, but should be adapted
#' to the species of interest. The good news is, that settings of L and mu affects
#' only the axis of the inferred history, but not the shape of the history.
#' So users can infer the shape, but need to be careful with a temporal interpretation
#' as both x and y axis are affected by the mutation rate and L.
#'
#' The SFS assumes that every locus is scored in all 2n sequences. A locus with
#' missing calls is counted in a lower frequency class than its true one, so
#' filter loci to full call rate (\code{gl.filter.callrate(x, threshold = 1)})
#' or impute them (\code{gl.impute}) before running the analysis.
#'
#' Stairway Plot 2 runs in a new subfolder of \code{tempdir()}. With
#' \code{cleanup = TRUE} only that subfolder is deleted at the end; with
#' \code{cleanup = FALSE} or \code{run = FALSE} it is kept and its path is
#' returned as \code{run.dir}. With \code{run = FALSE} the subfolder holds the
#' blueprint file and the Stairway Plot 2 folder, ready to be copied to another
#' machine (e.g. a cluster) and run there with
#' \code{java -cp stairway_plot_es Stairbuilder blueprint}.
#' @param x A genlight/dartR object containing SNP data [required].
#' @param L the length of the sequence in base pairs (see details)
#' [default NULL, number of loci x 69].
#' @param mu the mutation rate per base pair per generation (see details) [required].
#' @param stairway2.path the path to the folder that contains the Stairway
#' Plot 2 folder \code{stairway_plot_es} (check the example) [required].
#' @param minbinsize the minimum bin size for the SFS that should be used [default 1].
#' @param maxbinsize the maximum bin size for the SFS that should be used
#' [default NULL, the number of individuals in the dataset].
#' @param gentime the generation time in years [default 1].
#' @param sfs the folded site frequency spectrum (SFS) to be used for the analysis.
#' If not provided the SFS is created from the genlight/dartR object [default NULL].
#' @param parallel the number of parallel processes to use for the analysis [default 1].
#' @param run logical. If TRUE, the analysis is run immediately. Otherwise only the
#' blueprint file is created [might be useful to run on a cluster] [default TRUE].
#' @param blueprint the name of the blueprint file [default "blueprint"].
#' @param filename the name of the population, written as \code{popid} in the
#' blueprint and used for the Stairway Plot 2 input file names (no white space)
#' [default "sample"].
#' @param pct_training the proportion of sites to use for training [default 0.67].
#' @param nrand the number of break point settings to try, spread between 0.5
#' and 2 times (number of individuals - 1). If NULL, the four values
#' recommended by Stairway Plot 2, (nseq-2)/4, (nseq-2)/2, (nseq-2)*3/4 and
#' nseq-2, are used [default NULL].
#' @param stairway_plot_dir the name of the Stairway Plot 2 folder, as written
#' in the blueprint. It must stay "stairway_plot_es", the name of the folder
#' that is copied from \code{stairway2.path} [default "stairway_plot_es"].
#' @param nreps the number of bootstrap replicates to use for the analysis [default 200].
#' @param seed the random seed to use for the analysis [default NULL].
#' @param plot_title the title of the plot; also the name of the Stairway Plot 2
#' output files [default "Ne"].
#' @param xmin minimum x value for the plot [default 0].
#' @param xmax maximum x value for the plot [default 0].
#' @param ymin minimum y value for the plot [default 0].
#' @param ymax maximum y value for the plot [default 0].
#' @param xspacing spacing between x values for the plot [default 2].
#' @param yspacing spacing between y values for the plot [default 2].
#' @param fontsize the font size for the plot [default 12].
#' @param cleanup logical. If TRUE, the folder with the Stairway Plot 2 run
#' files is removed [default TRUE].
#' @param plot.display Specify if plot is to be produced [default TRUE].
#' @param plot.theme User specified theme [default theme_dartR()].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()]
#' @param plot.file Filename (minus extension) for the RDS plot file [Required for plot save]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary, including the
#' output of Stairway Plot 2; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return returns a list with three components:
#' \itemize{
#' \item{history: the Stairway Plot 2 summary table: time in mutations per
#' site and in years (\code{year}), the median effective population size
#' (\code{Ne_median}) and its 95\% (\code{low95}, \code{high95}) and 75\%
#' (\code{low75}, \code{high75}) limits; NULL when \code{run = FALSE}}
#' \item{plot: a ggplot of history; NULL when \code{run = FALSE}}
#' \item{run.dir: the folder with the Stairway Plot 2 run files; NULL when
#' it was removed (\code{cleanup = TRUE} and \code{run = TRUE})}
#' }
#'
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @references Liu, X., Fu, YX. Stairway Plot 2: demographic history inference with folded SNP frequency spectra. Genome Biol 21, 280 (2020). \doi{10.1186/s13059-020-02196-9}
#' @importFrom parallel detectCores
#' @importFrom future plan
#' @importFrom furrr future_map
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
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbose = verbose)

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)

    # FUNCTION SPECIFIC ERROR CHECKING

    if (datatype == "SilicoDArT") {
      stop(error("Fatal Error: Detected Presence/Absence (SilicoDArT) data. Please provide a SNP dataset\n"))
    }

    if(is.null(stairway_plot_dir)){
      stop(error("Fatal Error: Directory path for the Stairway Plot 2 executables not specified\n"))
    }
    if(is.null(mu)){
      stop(error("Fatal Error: Mutation rate per site per generation not specified\n"))
    }
    if(is.null(gentime)){
      stop(error("Fatal Error: Generation time (years) not specified\n"))
    }

    progs <- c("stairway_plot_es")
    if (!all(dir.exists(file.path(stairway2.path, progs)))) {
      stop(error(
        "Fatal Error: Cannot find the folder", progs,
        "in the folder given by stairway2.path:", stairway2.path,
        "\n  Download it with gl.download.binary(software = \"stairway2\")",
        "and set stairway2.path to the folder it was saved in.\n"
      ))
    }
    if (run && !nzchar(Sys.which("java"))) {
      stop(error(
        "Fatal Error: Java was not found. Stairway Plot 2 needs Java;",
        "install it and make sure the 'java' command is on the PATH.\n"
      ))
    }

    # check OS
    os <- tolower(Sys.info()['sysname'])

    # RUN STAIRWAY PLOT 2
    # a folder of its own, so that cleanup never touches other files in
    # the session tempdir
    tempd <- tempfile("stairway2_")
    dir.create(tempd, showWarnings = FALSE)
    outfilespec <- file.path(tempd, blueprint)

    file.copy(file.path(stairway2.path, progs),
              to = tempd,
              overwrite = TRUE, recursive = TRUE)

    if(is.null(maxbinsize)){
      maxbinsize <- nInd(x)
      if(verbose >= 3){cat(report("  Max Bin Size not specified, set to",nInd(x),"\n"))}
    }
    if(is.null(nrand))  nrand <- c(round((nInd(x)-1)/2), round(nInd(x)-1), round((nInd(x)-1)*3/2), round(2*(nInd(x)-1))) else nrand <- round(seq(0.5,2,len=nrand)*(nInd(x)-1))

    if(verbose >= 3){cat(report("  No. of break points, set to:",paste0(nrand, collapse=" "),"\n"))}

    whether_folded <- "true"
    nseq <- 2*nInd(x)

    if (is.null(sfs)) sfs <- gl.sfs(x, minbinsize=1, plot.out=FALSE, singlepop = TRUE, verbose = verbose)

    sfs <- paste(sfs,collapse = " ")

    #if total length of sequence is not specificed simply assume nLoc*69 (standard from dart)
    if (is.null(L)) L = nLoc(x)*69

    if (is.null(seed)) seed = round(runif(1)*1e6)

    # DO THE JOB

    if (verbose >= 2) {cat(report("  Writing the Stairway Plot 2 blueprint\n"))}

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
    

    if (verbose >= 3) {cat(report("  Stairway Plot 2 blueprint written to",outfilespec,"\n"))}

    res <- NULL
    p1 <- NULL

    if (run==TRUE)
    {
      oldpath <- getwd()
      setwd(tempd)
      on.exit(setwd(oldpath), add = TRUE)

      # Stairway Plot 2 output goes to the console only from verbose 3
      quiet <- verbose < 3
      run_cmd <- function(cmd) {
        system(cmd, ignore.stdout = quiet, ignore.stderr = quiet)
      }
      script <- paste0(blueprint, if (os == "windows") ".bat" else ".sh")
      plot_script <- paste0(blueprint, if (os == "windows") ".plot.bat" else ".plot.sh")
      run_script <- function(s) {
        if (os == "windows") run_cmd(s) else run_cmd(paste("bash", s))
      }

      if (verbose >= 2) {cat(report("  Running Stairway Plot 2\n"))}

      status <- run_cmd(paste0("java -cp stairway_plot_es Stairbuilder ",blueprint ))
      if (status != 0 || !file.exists(script)) {
        stop(error("Fatal Error: Stairway Plot 2 could not create the run script from the blueprint (Stairbuilder). Run files kept in",
                   tempd, "; rerun with verbose = 3 to see the output of Stairway Plot 2.\n"))
      }
      if (os != "windows") Sys.chmod(script, mode = "0755")

      #run on multiple cores
      if (parallel>1)
      {
        ff <- readLines(script)
        no_cores <- max(1, min(parallel, parallel::detectCores() - 1))
        oplan <- future::plan(future::multisession, workers = no_cores)
        on.exit(future::plan(oplan), add = TRUE)

        runs <- ff[grep("Stairway_fold_training_testing7", ff)]

        runstair <- function(i, wd) {
          owd <- setwd(wd)
          on.exit(setwd(owd))
          system(runs[i], ignore.stdout = quiet, ignore.stderr = quiet)
        }

        status <- unlist(furrr::future_map(seq_along(runs), runstair, wd = tempd))
        # restore the user's plan now, so the workers shut down before
        # their working folder is removed by cleanup
        future::plan(oplan)
        if (any(status != 0) && verbose >= 1) {
          cat(warn("  Warning:", sum(status != 0), "of", length(status),
                   "Stairway Plot 2 estimation runs failed\n"))
        }

        # move the training/testing files into the folders of each number of
        # break points: last two words of each "mv -f"/"MOVE /y" line
        moves <- ff[grep("^(mv -f|MOVE /y) ", ff)]
        for (mv in moves) {
          w <- strsplit(trimws(mv), "[[:space:]]+")[[1]]
          from <- w[length(w) - 1]
          to <- file.path(w[length(w)], basename(from))
          file.rename(from, to)
        }
        run_cmd(ff[grep("Stairpainter", ff)])
        run_script(plot_script)

      } else {
        run_script(script)
      }

      summary_file <- file.path(tempd,"files",paste0(plot_title,".final.summary"))

      # the summary step asks Java for 4 GB; retry once with 1 GB
      if (!file.exists(summary_file) && file.exists(plot_script)) {
        if (verbose >= 2) {
          cat(warn("  Attempt to rerun last step with different settings (lower memory allocation for the Java Virtual Machine)\n"))
        }
        ff <- readLines(plot_script)
        ff <-  gsub("-Xmx4g", "-Xmx1g", ff)
        writeLines(ff, plot_script)
        run_script(plot_script)
      }

      if (!file.exists(summary_file)) {
        stop(error("Fatal Error: Stairway Plot 2 did not produce the summary file", summary_file,
                   "\n  Run files kept in", tempd, "; rerun with verbose = 3 to see the output of Stairway Plot 2.\n"))
      }

      res <- read.csv(summary_file,sep="\t")
      setwd(oldpath)

      # PLOT
      mutation_per_site <- n_estimation <- theta_per_site_median <- theta_per_site_2.5 <- theta_per_site_97.5 <- year <- Ne_median <- low95 <- high95 <- low75 <- high75 <- NULL

      colnames(res) <- c("mutation_per_site" ,"n_estimation", "theta_per_site_median", "theta_per_site_2.5","theta_per_site_97.5", "year" , "Ne_median" ,"low95" ,"high95","low75" ,"high75")

      p1 <- ggplot(res, aes(x=year, y=Ne_median))+geom_point()+geom_line()+geom_ribbon(aes(ymin=low95, ymax=high95), alpha=0.2)+ylab("Effective population size")+xlab("Years ago")+plot.theme

      if (cleanup) {
        unlink(tempd, recursive = TRUE)
        tempd <- NULL
      } else if (verbose >= 2) {
        cat(report("  Check plots (pdf and png files) in folder:", tempd, "\n"))
      }

      # PRINTING OUTPUTS
      if (plot.display) {print(p1)}

      if(!is.null(plot.file)){
        tmp <- utils.plot.save(p1,
                               dir=plot.dir,
                               file=plot.file,
                               verbose=verbose)
      }
    } else if (verbose >= 2) {
      cat(report("  Blueprint and Stairway Plot 2 files written to:", tempd, "\n"))
    }

    # FLAG SCRIPT END

    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }

    # RETURN
    return(list(history=res, plot=p1, run.dir=tempd))

  }
