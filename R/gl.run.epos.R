#' @name gl.run.epos
#' @title Run EPOS for inference of historical population-size changes
#' @family demographic history
#' @description
#' This function runs EPOS (based on Lynch et al. 2019) to estimate historical population-size
#' \url{https://github.com/EvolBioInf/epos}. It relies on a compiled version of the software
#'  epos, epos2plot and if a bootstrap output is required bootSfs. For more information on the
#'  approach check the publication (Lynch at al. 2019), the github repository
#'  \url{https://github.com/EvolBioInf/epos} and look out for the manual epos.pdf
#'  (\url{https://github.com/EvolBioInf/epos/blob/master/doc/epos.pdf}).
#' The binaries need to be provided in a single folder and can be downloaded via the
#' \code{gl.download.binary} function (including the necessary dlls for windows; under Linux gsl, blas need to be installed on your system). Please note: if you use this method, make sure you cite the original publication in your work.
#' @details
#' EPOS (Estimation of Population Size changes) is a software tool that infers historical
#' changes in population size using allele-frequency data from population-genomic surveys.
#' The method relies on the site-frequency spectrum (SFS) of nearly neutral polymorphisms.
#' The underlying theory uses coalescence models, which describe how gene sequences have
#' originated from a common ancestor. By analyzing the probability distributions of the
#' starting and ending points of branch segments over all possible coalescence trees,
#' EPOS can estimate historic population sizes.\cr
#' The function uses a model-flexible approach, meaning it estimates historic population
#' sizes, without the necessity to provide a candidate scenario. An efficient statistical
#' procedure is employed, to estimate historic effective population sizes.\cr
#'  For all the possible settings, please refer to the manual of EPOS. \cr
#'  The main parameters that are necessary to run the function are a genlight/dartR object,
#' L (length of sequences), u (mutation rate), and the path to the epos binaries.
#' For details check the example below.\cr
#' Please note: There is currently not really a good way to estimate L, the length
#' of all sequences. Often users of dart data use the number of loci multiplied
#' by 69, but this is definitely an underestimate as monomorphic loci need to be
#' included (also the length of the restriction site should be added for each loci).
#' For mutation rate u, the default value is set to 5e-9, but should be adapted
#' to the species of interest. The good news is, that settings of L and mu affects
#' only the axis of the inferred history, but not the shape of the history.
#' So users can infer the shape, but need to be careful with a temporal interpretation
#' as both x and y axis are affected by the mutation rate and L.\cr
#' EPOS needs the number of monomorphic sites, either as the sequence length L
#' (option -l) or as the zero class of the SFS (\code{minbinsize = 0}). The zero
#' class counts only the monomorphic loci kept in the genlight object, which for
#' DArT data is usually far below the true number of monomorphic sites, so L is
#' usually the better choice.\cr
#' The SFS is sent to EPOS with its true class numbers: classes below
#' \code{minbinsize} are excluded with the EPOS option -x, and an unfolded SFS
#' (\code{folded = FALSE}) is sent as classes 1 to 2n-1 with option -U.\cr
#' The SFS assumes that every locus is scored in all 2n sequences. A locus with
#' missing calls is counted in a lower frequency class than its true one, so
#' filter loci to full call rate (\code{gl.filter.callrate(x, threshold = 1)})
#' or impute them (\code{gl.impute}) before running the analysis.
#' @param x dartR/genlight object with SNP data [required].
#' @param epos.path path to the folder with epos and epos2plot (always required)
#' and bootSfs (required if boot > 0) [required].
#' @param sfs the SFS to use instead of calculating it from x. Either named as
#' returned by \code{gl.sfs} ("d0", "d1", ...; the names give the classes), or
#' an unnamed vector of counts for the classes minbinsize, minbinsize + 1, ...
#' [default NULL, calculated with gl.sfs(x, singlepop = TRUE)].
#' @param minbinsize the smallest SFS class used. 0 keeps the zero class
#' (monomorphic sites), which EPOS then uses instead of L; 1 (default) drops it;
#' 2 also excludes singletons, and so on (option -x in epos). If your
#' genlight object has more than one population the sfs is calculated with
#' singlepop set to TRUE (one sfs for all individuals), as epos does not work
#' with a multidimensional sfs [default 1].
#' @param folded if set to TRUE (default) a folded sfs (minor allele frequency sfs) is used. If set to FALSE an unfolded (derived allele frequency sfs) is used. It is assumed that 0 is homozygote for the reference and 2 is homozygote for the derived allele. So you need to make sure your coding is correct. Option -U in epos [default TRUE].
#' @param L length of sequences (including monomorphic and polymorphic sites).
#' Required unless minbinsize = 0. Option -l in epos [default NULL].
#' @param u mutation rate. If not provided the default value of epos is used (5e-9).
#' Option -u in epos [default NULL].
#' @param boot if set to a value >0 the program bootSfs is used to provide multiple
#' bootstrapped sfs, which allows to calculate confidence intervals of the historic Ne
#' sizes. Be aware the runtime can be extended. 0: no bootstrapped simulations are
#' run, otherwise boot number of bootstraps are run (option -i in bootSfs) [default 0].
#' @param upper upper quantile of the bootstrap (only used if boot>0).
#' Option -u in epos2plot [default 0.975].
#' @param lower lower quantile of the bootstrap (only used if boot>0).
#' Option -l in epos2plot [default 0.025].
#' @param method either "exhaustive" or "greedy". Check the epos manual for details. If method="exhaustive" then the parameter depth is used [default "greedy"].
#' @param depth if method="exhaustive" then this parameter is used to set the search depth. If method is set to greedy this setting is ignored [default 2].
#' @param other.options additional options for epos, as typed on the command
#' line (e.g "-m 10") [default ""].
#' @param seed seed for the random number generators of epos and bootSfs
#' (option -s). If NULL, both programs use a seed from the system clock and
#' results differ between runs [default NULL].
#' @param outfile File name of the output file (the raw epos output) [default 'epos.out'].
#' @param outpath Path where to save the output file [default global working
#' directory or if not specified, tempdir()].
#' @param cleanup if set to true intermediate tempfiles are deleted after the run [default TRUE].
#' @param plot.display Specify if plot is to be produced [default TRUE].
#' @param plot.theme User specified theme [default theme_dartR()].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()]
#' @param plot.file Filename (minus extension) for the RDS plot file [Required for plot save]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary, including the
#' messages of epos; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return returns a list with four components:
#' \itemize{
#' \item{history: Ne estimates over generations (generation, low, median and high)}
#' \item{plot: a ggplot of history }
#' \item{sfs: the sfs used for the analysis (class r and count fr)}
#' \item{diagnostics: a list with the several diagnostics and a plot of observed and expected sfs}
#' }
#' @export
#' @examples
#' \dontrun{
#' #gl.download.binary("epos",os="windows")
#' require(dartR.data)
#' epos <- gl.run.epos(possums.gl, epos.path = file.path(tempdir(),"epos"), L=1e5, u = 1e-8)
#' epos$history
#' }
#'
#'
#' @references Lynch, Michael, Bernhard Haubold, Peter Pfaffelhuber, and Takahiro Maruki. 2019. Inference of Historical Population-Size Changes with Allele-Frequency Data. G3: Genes|Genomes|Genetics 10, no. 1: 211-23. \doi{10.1534/g3.119.400854}.
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}


gl.run.epos <- function(x,
                        epos.path,
                        sfs=NULL,
                        minbinsize=1,
                        folded=TRUE,
                        L=NULL,
                        u=NULL,
                        boot=0,
                        upper=0.975,
                        lower=0.025,
                        method="greedy",
                        depth=2,
                        other.options="",
                        seed=NULL,
                        outfile="epos.out",
                        outpath=NULL,
                        cleanup=TRUE,
                        plot.display=TRUE,
                        plot.theme = theme_dartR(),
                        plot.dir=NULL,
                        plot.file=NULL,
                        verbose=NULL)
{

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir,verbose=0)
  # SET WORKING DIRECTORY for file
  outpath <- gl.check.wd(outpath,verbose=0)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  if (datatype == "SilicoDArT") {
    stop(error("Fatal Error: Detected Presence/Absence (SilicoDArT) data. Please provide a SNP dataset\n"))
  }
  methods <- c("exhaustive","greedy")
  if (is.na(pmatch(method, table = methods))) stop(error(paste0("method must be one of ",paste(methods,collapse=", "))))
  if (minbinsize > 0 && is.null(L)) {
    stop(error("Fatal Error: L (sequence length) is required unless minbinsize = 0, in which case the zero class of the sfs is used\n"))
  }

  # check OS
  os <- tolower(Sys.info()['sysname'] )

  # check if epos epos2plot and [bootSfs are there]
  progs <- c("epos", "epos2plot")
  if (boot>0) progs <- c(progs,"bootSfs")
  if (os=="windows") progs <- paste0(progs,".exe")
  if (os=="windows") progs <- c(progs, "libblas.dll","libgsl.dll","libgslcblas.dll")
  fex <- file.exists(file.path(epos.path, progs))
  if (!all(fex)) {
    stop(error(
      "Fatal Error: Cannot find", paste(progs[!fex], collapse = ", "),
      "in the folder given by epos.path:", epos.path,
      "\n  Download them with gl.download.binary(software = \"epos\")",
      "and set epos.path to the folder they were saved in.\n"
    ))
  }

  tempd <-  tempfile(pattern = "dir")
  dir.create(tempd, showWarnings = FALSE)
  file.copy(file.path(epos.path, progs),
            to = tempd,
            overwrite = TRUE)
  if (os != "windows") Sys.chmod(file.path(tempd, progs), mode = "0755")

  #initialise variables in data.frames
  e <- o <- r <- sde <- fr <- meane <- NULL

  # SFS with its true class numbers
  if (is.null(sfs)) {
    sfs <- gl.sfs(x, minbinsize = 0, folded = folded, singlepop = TRUE, plot.out = FALSE, verbose = 0)
    cls <- as.numeric(substr(names(sfs), 2, 100))
    # unfolded: the last class holds loci fixed for the derived allele
    if (!folded) {
      keep <- cls < max(cls)
      sfs <- sfs[keep]
      cls <- cls[keep]
    }
  } else if (!is.null(names(sfs)) && all(grepl("^d[0-9]+$", names(sfs)))) {
    cls <- as.numeric(substr(names(sfs), 2, 100))
  } else {
    cls <- minbinsize - 1 + seq_along(sfs)
  }
  sfs <- as.numeric(sfs)
  # classes 1..max are all written (missing ones as 0) so that epos reads the
  # right sample size; classes below minbinsize are excluded with -x
  maxcls <- max(cls)
  send <- data.frame(r = 1:maxcls, fr = 0)
  send$fr[match(cls[cls >= 1], send$r)] <- sfs[cls >= 1]
  if (minbinsize == 0) {
    send <- rbind(data.frame(r = 0, fr = sum(sfs[cls == 0])), send)
  }
  xcls <- seq_len(max(minbinsize - 1, 0))
  write.table(send, file = file.path(tempd, "dummy.sfs"), row.names = F, sep = "\t", col.names = TRUE, quote = FALSE)

  # run a program from the run folder; stop with its message if it fails
  run_prog <- function(prog, args) {
    if (os == "windows") prog <- paste0(prog, ".exe")
    errf <- tempfile(tmpdir = tempd)
    out <- suppressWarnings(system2(file.path(tempd, prog), args, stdout = TRUE, stderr = errf))
    msg <- if (file.exists(errf)) readLines(errf, warn = FALSE) else character(0)
    status <- attr(out, "status")
    if (!is.null(status) && status != 0) {
      stop(error("Fatal Error:", prog, "failed (exit status", paste0(status, "):"),
                 paste(msg, collapse = " "), "\n"))
    }
    if (verbose >= 3 && length(msg) > 0) cat(report(paste0("  ", msg, "\n")))
    out
  }
  seedarg <- if (is.null(seed)) character(0) else c("-s", seed)

  eposargs <- c(if (minbinsize > 0) c("-l", L),
                if (!is.null(u)) c("-u", u),
                if (length(xcls) > 0) c("-x", paste(xcls, collapse = ",")),
                if (!folded) "-U",
                if (pmatch(method, table = methods) == 1) c("-E", depth),
                seedarg,
                other.options,
                "-o",
                if (boot > 0) "bs.sfs" else "dummy.sfs")

  # DO THE JOB
  old.path <- getwd()
  setwd(tempd)
  on.exit(setwd(old.path))
  if (verbose >= 2) cat(report("  Running epos\n"))
  if (boot>0) {
    bsdummy <- run_prog("bootSfs", c("-i", boot, seedarg, "dummy.sfs"))
    writeLines(bsdummy,file.path(tempd,"bs.sfs"))
  }
  epdummy <- run_prog("epos", eposargs)
  writeLines(epdummy,file.path(tempd,"ep.dat"))
  eposout <- run_prog("epos2plot", c("-l", lower, "-u", upper, "ep.dat"))
  setwd(old.path)
  ep2 <- (do.call(rbind,(strsplit(eposout,split = "\t"))))
  epp <- data.frame(ep2[-1,])
  colnames(epp)<- ep2[1,]
  epp <- data.frame(apply(epp,2, as.numeric))

  # PLOT
  generation <-low <- high <- median <- NULL
  colnames(epp) <- c("generation", "low", "median", "high")

  p1 <- ggplot(epp, aes(x=generation, y=median))+geom_point()+geom_line()+geom_ribbon(aes(ymin=low, ymax=high), alpha=0.2)+ylab("Effective population size")+xlab("Generation")+plot.theme


  #parse ep.dat
  con <- file(file.path(tempd,"ep.dat"), "r")
  ep <- readLines(con)
  close(con)

  #poly
  ll <- which(substr(ep,1,12)=="#Polymorphic")
  polymorphic_sites <- as.numeric(gsub(".*#Polymorphic sites surveyed:\\s*([0-9]+).*", "\\1",ep[ll]))
  #mono
  ll <- which(substr(ep,1,12)=="#Monomorphic")
 monomorphic_sites <- as.numeric(gsub(".*#Monomorphic sites surveyed:\\s*([0-9]+).*", "\\1",ep[ll]))


 #likelihood
 ll <- which(substr(ep,1,10)=="#Final Log")
 fll<- as.numeric(gsub(".*Log\\(Likelihood\\):\\s*([-0-9.]+).*", "\\1", ep[ll]))

 #d2
 ll <- which(substr(ep,1,5)=="#d^2:")
 d2<- as.numeric(gsub(".*#d\\^2:\\s*([-0-9.]+).*", "\\1", ep[ll]))

 #find sfs(s)

 ll <- which(substr(ep,1,4)=="#sfs")
 if (length(ll)>0) {
 ll <- c(min(ll)-1,ll) #find the header
 sfss <- read.csv(text=ep[ll], header =T, sep = "\t")
 if (boot==0) boot=1
 sfsl <-list()
 ff <- rep(1:boot, each=(length(ll)-1)/boot)
 for (i in 1:boot) sfsl[[i]] <- sfss[ff==i,]

 } else sfsl <- NULL

  # OUTPUT
  #if outpath not null copy to outpath
  if (!is.null(outfile)) {
    file.copy(file.path(tempd,"ep.dat"), file.path(outpath,outfile), overwrite = TRUE)
    if (verbose >= 2) cat(report(paste("  Output written to", file.path(outpath,outfile), "\n")))

  }

  #sfs format: the classes used
  used <- cls >= minbinsize
  dfsfs <- data.frame(r=cls[used], fr=sfs[used])

  #PLOT EXPECTED VS OBSERVED SFS
  if (!is.null(sfsl)) {
  xx <- do.call(rbind,sfsl)

  sfsm <- plyr::ddply(xx,.variables =  "r", plyr::summarise, meane=mean(e), sde=sd(e), meano=mean(o), sdo=sd(o))

  p2 <- ggplot(dfsfs, aes(x=r, y=fr) )+ geom_bar(stat="identity",color="darkgrey", fill="darkgrey")  + geom_errorbar(data=sfsm, aes(x=r-0.1, y=meane, ymin = meane-1.96*sde, ymax = meane+1.96*sde), color="purple") +  labs(x="bin", y="frequency")+plot.theme

  if (boot==1)
  p2 <- p2+geom_point(data=data.frame(sfsl[[1]]),aes( x=r, y=e), color="purple", size=1)
 }
  else p2 <- NULL

  out <- list(
    polymorphic_sites = polymorphic_sites,
    monomorphic_sites = monomorphic_sites,
    likelihood = fll,
    d2 = d2,
    sfs = sfsl,
    sfs_plot = p2
  )


  if (cleanup) unlink(tempd, recursive = T)
  # PRINTING OUTPUTS
  if (plot.display) {print(p1)}

  if(!is.null(plot.file)){
    tmp <- utils.plot.save(p1,
                           dir=plot.dir,
                           file=plot.file,
                           verbose=verbose)
  }

  # FLAG SCRIPT END

  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN


  return(list(history=epp, plot=p1, sfs=dfsfs, diagnostics=out))
}
