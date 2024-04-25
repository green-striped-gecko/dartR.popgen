#' @name gl.run.structure
#'
#' @title Runs a STRUCTURE analysis using a genlight object
#'
#' @description
#' This function takes a genlight object and runs a STRUCTURE analysis based on
#' functions from \code{strataG}
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param exec Full path and name+extension where the structure executable is
#' located. E.g. \code{'c:/structure/structure.exe'} under Windows. For Mac and
#' Linux it might be something like \code{'./structure/structure'} if the
#' executable is in a subfolder 'structure' in your home directory
#' [default working directory "."].
#' @param k.range Range of the number of populations [required].
#' @param num.k.rep Number of replicates [default 1].
#' @param burnin Number of iterations for MCMC burnin [default 1000].
#' @param numreps Number of MCMC replicates [default 1000].
#' @param noadmix Logical. No admixture? [default TRUE].
#' @param freqscorr Logical. Correlated frequencies? [default FALSE].
#' @param randomize Randomize [default TRUE].
#' @param seed Set random seed [default 0].
#' @param pop.prior A character specifying which population prior model to use:
#'  "locprior" or "usepopinfo" [default NULL].
#' @param locpriorinit Parameterizes locprior parameter r - how informative the
#'  populations are. Only used when pop.prior = "locprior" [default 1].
#' @param maxlocprior Specifies range of locprior parameter r. Only used when 
#' pop.prior = "locprior" [default 20].
#' @param gensback Integer defining the number of generations back to test for 
#' immigrant ancestry. Only used when pop.prior = "usepopinfo" [default 2].
#' @param migrprior Numeric between 0 and 1 listing migration prior. Only used 
#' when pop.prior = "usepopinfo" [default 0.05].
#' @param pfrompopflagonly Logical. update allele frequencies from individuals 
#' specified by popflag. Only used when pop.prior = "usepopinfo" [default TRUE].
#' @param popflag A vector of integers (0, 1) or logicals identifiying whether 
#' or not to use strata information. Only used when pop.prior = "usepopinfo"
#'  [default NULL].
#' @param inferalpha Logical. Infer the value of the model parameter # from the 
#' data; otherwise is fixed at the value alpha which is chosen by the user. 
#' This option is ignored under the NOADMIX model. Small alpha implies that 
#' most individuals are essentially from one population or another, while 
#' alpha > 1 implies that most individuals are admixed [default FALSE].
#' @param alpha Dirichlet parameter for degree of admixture. This is the 
#' initial value if inferalpha = TRUE [default 1].
#' @param unifprioralpha Logical. Assume a uniform prior for alpha which runs 
#' between 0 and alphamax. This model seems to work fine; the alternative model
#'  (when unfprioralpha = 0) is to take alpha as having a Gamma prior, with 
#'  mean alphapriora × alphapriorb, and variance alphapriora × alphapriorb^2 
#'  [default TRUE].
#' @param alphamax Maximum for uniform prior on alpha when 
#' unifprioralpha = TRUE [default 20].
#' @param alphapriora Parameters of Gamma prior on alpha when 
#' unifprioralpha = FALSE [default 0.05].
#' @param alphapriorb Parameters of Gamma prior on alpha when 
#' unifprioralpha = FALSE [default 0.001].
#' @param plot.out Create an Evanno plot once finished. Be aware k.range needs
#' to be at least three different k steps [default TRUE].
#' @param plot_theme Theme for the plot. See details for options
#' [default theme_dartR()].
#' @param plot.dir Directory to save the plot RDS files [default as specified 
#' by the global working directory or tempdir()]
#' @param plot.file Name for the RDS binary file to save (base name only, 
#' exclude extension) [default NULL]
#' @param verbose Set verbosity for this function (though structure output
#' cannot be switched off currently) [default NULL].
#' @param ... Parameters to specify the STRUCTURE run (check 
#' \code{structureRun}
#'  within strataG.
#' @details The function is basically a convenient wrapper around the beautiful
#' strataG function \code{structureRun} (Archer et al. 2016). For a detailed
#' description please refer to this package (see references below).
#' 
#' Before running STRUCTURE, we suggest reading its manual (see link below) and 
#' the literature in mentioned in the references section. 
#'  
#' https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/structure_doc.pdf
#' To make use of this function you need to download STRUCTURE for you system
#' (\bold{non GUI version}) from here
#' \href{https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html}{STRUCTURE}.
#'
#' \bold{Format note}
#'
#' For this function to work, make sure that individual and population names
#' have no spaces. To substitute spaces by underscores you could use the R
#' function \code{gsub} as below.
#'
#' \code{
#' popNames(gl) <- gsub(" ","_",popNames(gl));
#' }
#' \code{
#' indNames(gl) <- gsub(" ","_",indNames(gl))
#' }
#'
#' It's also worth noting that Structure truncates individual names at 11
#' characters. The function will fail if the names of individuals are not unique
#'  after truncation. To avoid this possible problem, a number sequence, as
#'  shown in the code below, might be used instead of individual names.
#' \code{
#' indNames(gl) <- as.character(1:length(indNames(gl)))
#' }
#' @return An sr object (structure.result list output). Each list entry is a
#' single structurerun output (there are k.range * num.k.rep number of runs).
#' For example the summary output of the first run can be accessed via
#' \code{sr[[1]]$summary} or the q-matrix of the third run via
#' \code{sr[[3]]$q.mat}. To conveniently summarise the outputs across runs
#' (clumpp) you need to run gl.plot.structure on the returned sr object. For
#' Evanno plots run gl.evanno on your sr object.
#'
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' # examples need structure to be installed on the system (see above)
#' \dontrun{
#' bc <- bandicoot.gl[,1:100]
#' sr <- gl.run.structure(bc, k.range = 2:5, num.k.rep = 3,
#' exec = './structure.exe')
#' ev <- gl.evanno(sr)
#' ev
#' qmat <- gl.plot.structure(sr, K=3)
#' head(qmat)
#' gl.map.structure(qmat, bc, scalex=1, scaley=0.5)
#' }
#' @import patchwork
#' @importFrom dplyr bind_rows mutate_at vars starts_with mutate group_by
#' ungroup arrange n rename select everything n_distinct bind_rows starts_with
#' @export
#' @references
#' \itemize{
#' \item Pritchard, J.K., Stephens, M., Donnelly, P. (2000) Inference of
#' population structure using multilocus genotype data. Genetics 155, 945-959.
#' \item Archer, F. I., Adams, P. E. and Schneiders, B. B. (2016) strataG: An R
#' package for manipulating, summarizing and analysing population genetic data.
#'  Mol Ecol Resour. doi:10.1111/1755-0998.12559
#' \item Wang, Jinliang. "The computer program structure for assigning 
#' individuals to populations: easy to use but easier to misuse." Molecular
#'  ecology resources 17.5 (2017): 981-990.
#' \item Lawson, Daniel J., Lucy Van Dorp, and Daniel Falush. "A tutorial on 
#' how not to over-interpret STRUCTURE and ADMIXTURE bar plots." Nature
#'  communications 9.1 (2018): 3258.
#' \item Porras-Hurtado, Liliana, et al. "An overview of STRUCTURE: 
#' applications, parameter settings, and supporting software." Frontiers in
#'  genetics 4 (2013): 98.
#' }

gl.run.structure <- function(x,
                             exec = "./structure",
                             k.range = NULL,
                             num.k.rep = 1,
                             burnin = 1000,
                             numreps = 1000,
                             noadmix = TRUE,
                             freqscorr = FALSE,
                             randomize = TRUE,
                             seed = 0,
                             pop.prior = NULL,
                             locpriorinit = 1,
                             maxlocprior = 20,
                             gensback = 2,
                             migrprior = 0.05,
                             pfrompopflagonly = TRUE,
                             popflag = NULL,
                             inferalpha = FALSE,
                             alpha = 1,
                             unifprioralpha = TRUE,
                             alphamax = 20,
                             alphapriora = 0.05,
                             alphapriorb = 0.001,
                             plot.out = TRUE,
                             plot_theme = theme_dartR(),
                             plot.dir = tempdir(),
                             plot.file = NULL,
                             verbose = NULL
                             # ,
                             # ...
                             ) {
  pkg <- "purrr"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  # check that Structure is installed
  structure <- file.exists(exec)

  if (!structure) {
    stop(error(
      paste(
        "Cannot find Structure executable in the exex path provided:\n",
        exec,
        "\nCheck the help page of ?gl.run.structure on how to download and the exec
parameter to locate it."
      )
    ))
  }
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir,verbose=0)
  

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(
    func = funname,
    build = "Jody",
    verbose = verbose
  )

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  if (datatype != "SNP") {
    stop(error(
      "You need to provide a SNP genlight object (ploidy=2)!"
    ))
  }

  # DO THE JOB
  gg <- utils.structure.genind2gtypes(gl2gi(x, verbose = 0))

  sr <- utils.structure.run(g = gg, 
                            exec = exec,
                            k.range = k.range,
                            num.k.rep = num.k.rep,
                            burnin = burnin,
                            numreps = numreps,
                            noadmix = noadmix,
                            freqscorr = freqscorr,
                            randomize = randomize,
                            seed = seed,
                            pop.prior = pop.prior,
                            locpriorinit = locpriorinit,
                            maxlocprior = maxlocprior,
                            gensback = gensback,
                            migrprior = migrprior,
                            pfrompopflagonly = pfrompopflagonly,
                            popflag = popflag,
                            inferalpha = inferalpha,
                            alpha = alpha,
                            unifprioralpha = unifprioralpha,
                            alphamax = alphamax,
                            alphapriora = alphapriora,
                            alphapriorb = alphapriorb
                            # ...
                            )

  ev <- utils.structure.evanno(sr)

  if(num.k.rep==1){
    pa <- ((ev$plots$mean.ln.k + ev$plots$mean.ln.k) /
             (ev$plots$ln.ppk)) + plot_theme
  }else{
    pa <- ((ev$plots$mean.ln.k + ev$plots$mean.ln.k) /
             (ev$plots$ln.ppk + ev$plots$delta.k)) + plot_theme
  }

  # PRINTING OUTPUTS
  if (plot.out) {
    suppressMessages(print(pa))
  }

  if(!is.null(plot.file)){
    tmp <- utils.plot.save(pa,
                           dir=plot.dir,
                           file=plot.file,
                           verbose=verbose)
  }
  
  
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n\n"))
  }

  # RETURN
  return(sr)
}
