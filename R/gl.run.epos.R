#' Run EPOS for Inference of Historical Population-Size Changes
#'
#' This function runs EPOS (based on Lynch et al. 2019) to estimate historical population-size
#' \url{https://github.com/EvolBioInf/epos}. It relies on a compiled version of the software
#'  epos, epos2plot and if a bootstrap output is required bootSfs. For more information on the 
#'  approach check the publication (Lynch at al. 2019), the github repository 
#'  \url{https://github.com/EvolBioInf/epos} and look out for the manual epos.pdf 
#'  (\url{https://github.com/EvolBioInf/epos/blob/master/doc/epos.pdf}. 
#' The binaries need to be provided in a single folder and can be downloaded via the 
#' \code{gl.download.binary} function (including the necessary dlls for windows; under Linux gls, blas need to be installed on your system). Please note: if you use this method, make sure you cite the original publication in your work.
#' EPOS (Estimation of Population Size changes) is a software tool developed based on the theoretical framework outlined by Lynch et al. (2019). It is designed to infer historical changes in population size using allele-frequency data obtained from population-genomic surveys. Below is a brief summary of the main concepts of EPOS:\cr\cr
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
#' as both x and y axis are affected by the mutation rate and L.
#' @param x dartR/genlight object
#' @param epos.path path to epos and other required programs (epos, epos2plot are always required and bootSfs in case a bootstrap and confidence estimate is required )
#' @param sfs if no sfs is provided function gl.sfs(x, minbinsize=1, singlepop=TRUE) is used to calculate the sfs that is provided to epos 
#' @param minbinsize remove bins from the left of the sfs. if you run epos from a genlight object the sfs is calculated by the function (using gl.sfs) and as default minbinsize is set to 1 (the monomorphic loci of the sfs are removed). This parameter is ignored if sfs is provide via the sfs parameter (see below). Be aware even if you genlight object has more than one population the sfs is calculated with singlepop set to true (one sfs for all individuals) as epos does not work with multidimensional sfs)
#' @param folded if set to TRUE (default) a folded sfs (minor allele frequency sfs) is returned. If set to FALSE then an unfolded (derived allele frequency sfs) is returned. It is assumed that 0 is homozygote for the reference and 2 is homozygote for the derived allele. So you need to make sure your coding is correct. option -U in epos.
#' @param L length of sequences (including monomorphic and polymorphic sites). If the sfs is provided with minbinsize=1 (default) then L needs to be specified. 
#' option -l in epos
#' @param u mutation rate. If not provided the default value of epos is used (5e-9).
#' option -u in epos
#' @param boot if set to a value >0 the programm bootSfs is used to provide multiple
#' bootstrapped sfs, which allows to calculate confidence intervals of the historic Ne
#' sizes. Be aware the runtime can be extended. default:0 no bootstrapped simulations are
#' run, otherwise boot number of bootstraps are run (option -i in bootSfs)
#' @param upper upper quantile of the bootstrap (only used if boot>0). default 0.975.
#' (option -u in epos2plot)
#' @param lower lower quantile of the bootstrap (only used if boot>0). default 0.025.
#' (option -l in epos2plot)
#' @param method either "exhaustive" or "greedy". check the epos manual for details. If method="exhaustive" then the paramter depth is used. default: "greedy".
#' @param depth if method="exhaustive" then this parameter is used to set the search depth, default is 2. If method is set to greedy this is setting is ignored.
#' @param other.options additional options for epos (e.g -m, -x etc.)
#' @param outfile File name of the output file [default 'genepop.gen'].
#' @param outpath Path where to save the output file [default global working 
#' directory or if not specified, tempdir()].
#' @param cleanup if set to true intermediate tempfiles are deleted after the run
#' @param plot.display Specify if plot is to be produced [default TRUE].
#' @param plot.theme User specified theme [default theme_dartR()].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]
#' @param plot.file Filename (minus extension) for the RDS plot file [Required for plot save]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return returns a list with four components: 
#' \itemize{
#' \item{history: Ne estimates of over generations (generation, median, low and high)} 
#' \item{plot: a ggplot of history }
#' \item{sfs: the sfs used for the analysis}
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
#'@references Lynch, Michael, Bernhard Haubold, Peter Pfaffelhuber, and Takahiro Maruki. 2019. Inference of Historical Population-Size Changes with Allele-Frequency Data. G3: Genes|Genomes|Genetics 10, no. 1: 211–23. \doi{10.1534/g3.119.400854}.
#'@author Custodian: Bernd Gruber -- Post to 
#'\url{https://groups.google.com/d/forum/dartr}


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
  #utils.flag.start(func = funname,
  #                 build = "Jody",
  #                 verbosity = verbose)
  

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  # all sanity checks (e.g. if folded length(sfs)=nInd etc.)
  methods <- c("exhaustive","greedy")
  if (is.na(pmatch(method, table = methods))) stop(paste0("method must be one of ",paste(methods,collapse=", ")))
  
  # check OS
  os <- tolower(Sys.info()['sysname'] )
  
  tempd <-  tempfile(pattern = "dir")
  dir.create(tempd, showWarnings = FALSE)
  #if (os=="Linux" | os=="Darwin" )  os=="Windows"
  # check if epos epos2plot and [bootSfs are there]
  progs <- c("epos", "epos2plot")
  if (boot>0) progs <- c(progs,"bootSfs")
  if (os=="windows") progs <- paste0(progs,".exe") 
  if (os=="windows") progs <- c(progs, "libblas.dll","libgsl.dll","libgslcblas.dll")
  fex <- file.exists(file.path(epos.path, progs))
  if (all(fex)) {
    file.copy(file.path(epos.path, progs),
              to = tempd,
              overwrite = TRUE)
  } else{
    cat(
      "  Cannot find",
      progs[!fex],
      "in the specified folder given by epos.path:",
      epos.path,
      "\n"
      
    )
    stop()
  }
  
  if (is.null(sfs)) sfs <- gl.sfs(x, minbinsize=minbinsize,folded=folded,singlepop = TRUE, plot.out = FALSE, verbose = 0 )
  df <- data.frame(r=1:length(sfs), fr=sfs)
  write.table(df,file=file.path(tempd,"dummy.sfs"), row.names = F, sep="\t", col.names = TRUE, quote = FALSE)
  bootcmd<-""
  
  if (boot>0) bootcmd <- paste("bootSfs -i",boot, "dummy.sfs ")
  if (os!="windows") bootcmd <- paste0("./", bootcmd)
  if (minbinsize==0) lcmd="" else lcmd=paste0(" -l ", L)
  ucmd <- paste0(" -u ",u)
  if (boot>0) sfsfile <-"bs.sfs" else sfsfile <- " dummy.sfs"
  if (pmatch(method, table=methods)==1) depthcmd <- paste0(" -E ",depth) else depthcmd <- ""
  eposcmd <- paste0("epos ",lcmd, ucmd,depthcmd, other.options, " -o ", sfsfile)
  if (os!="windows") eposcmd <- paste0("./", eposcmd)
  # DO THE JOB
  old.path <- getwd()
  setwd(tempd)
  on.exit(setwd(old.path))
  if (boot>0) {
    if (os=="linux") system("chmod 777 bootSfs")
    bsdummy <- system(bootcmd, intern = TRUE)
    writeLines(bsdummy,file.path(tempd,"bs.sfs"))
  }
  if (os=="linux") system("chmod 777 epos")
  epdummy <- system(eposcmd, intern=TRUE)
  writeLines(epdummy,file.path(tempd,"ep.dat"))
  eposplotcmd <-"epos2plot ep.dat"
  if (os=="linux") system("chmod 777 epos2plot")
  if (os!="windows") eposplotcmd <- paste0("./", eposplotcmd)
  eposout <- system(eposplotcmd, intern = TRUE)
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
    cat(report(paste("  Output written to", file.path(outpath,outfile), "\n")))
    
  }
  
  #sfs format
  bins <- as.numeric(substr(names(sfs),2,100))
  dfsfs <- data.frame(r=bins, fr=sfs)

  #PLOT EXPECTED VS OBSERVED SFS
  if (!is.null(sfsl)) {
  xx <- do.call(rbind,sfsl)
  
  sfsm <- plyr::ddply(xx,.variables =  "r", plyr::summarise, meane=mean(e), sde=sd(e), meano=mean(o), sdo=sd(o))
  
  
  p2 <- ggplot(dfsfs, aes(x=r, y=fr) )+ geom_bar(stat="identity",color="darkgrey", fill="darkgrey")  + geom_errorbar(data=sfsm, aes(x=r-0.1, y=meane, ymin = meane-1.96*sde, ymax = meane+1.96*sde), color="purple") +  labs(x="bin", y="frequency")+theme_bw() #+ geom_errorbar(data=sfsm, aes(x=r+0.1, y=meano, ymin = meano-1.96*sdo, ymax = meano+1.96*sdo), color="orange") + geom_point(data=sfsm, aes(x=r-0.1, y=meane), color="blue", size=0.8)
  
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
  #print(epp)
  
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

