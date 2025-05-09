#' @name gl.TajimasD
#' @title Calculation of Tajima's D
#' simulation results can only be output if ms and sample_stats are compiled (can be installed at: https://uchicago.app.box.com/s/l3e5uf13tikfjm7e1il1eujitlsjdx13)
#' @description
#' This function calculate Tajima's D, with p-values from beta distribution, standard normal distribution and simulation from ms (Hudson, 2002)
#'
#' Refer to the ms manual for further information on the parameters to set
#' ##########################################
#'
#'
#'@param x Name of the genlight object containing the SNP data [required]
#' @param ms.path absolute path that stores the ms program 
#' (eg: /User/msdir/) [default  NULL]
#' @param simulation.out Directory in which to save simulated summary statistics from MS, given ms.path is provided  [default  NULL]
#' @param rep Number of replicates in ms [default  NULL, required for simulation]
#' @param seeds Seeds for the random number generator in ms [default  NULL,
#'  seeds are randomly generated]
#' @param cleanup clean data in tmp [default  TRUE]
#' @param plot.dir Directory in which to save files [default = working directory]
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot.file Name for the RDS binary file to save (base name only, exclude
#' extension) [default NULL]
#' @param plot_theme Theme of the plot [default theme_dartR()]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return A plot and table of Tajima's D (results from MS and plot of simulated Tajima's D can be returned if ms.path is provided) 
#' @export
#' @importFrom stringr str_split_fixed
#' @importFrom terra split
#' @references
#' \itemize{
#' \item Tajima, F. (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. Genetics, 123(3), 585-595.
#' \item Hudson, R. R. (2002). Generating samples under a Wright–Fisher neutral model of genetic variation. Bioinformatics, 18(2), 337-338.
#' \item Paradis, E. (2010). pegas: an R package for population genetics with an integrated–modular approach. Bioinformatics, 26(3), 419-420.
#' }
#' @author Renee Catullo (Custodian: Ching Ching Lau) -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # To run without ms simulation
#' Tajima <- gl.TajimasD(x=bandicoot.gl)
#' #' # To run with ms simulation
#' \dontrun{
#' Tajima <- gl.TajimasD(x=bandicoot.gl, rep=10, ms.path="/User/msdir/")
#' }
#' @export 
#' @importFrom stats D pbeta pnorm sd


gl.TajimasD <- function(x, ms.path=NULL,
                        simulation.out=NULL,
                        rep=NULL, 
                        seeds=NULL,
                        cleanup=TRUE, 
                        plot.dir=NULL,
                        plot.out = TRUE,
                        plot.file = NULL,
                        plot_theme = NULL,
                        verbose=2) 
{
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  
  #create tempdir
  tempd <-  tempfile(pattern = "dir")
  dir.create(tempd, showWarnings = FALSE)
  
  if (is.null(plot_theme)) {
    plot_theme <- theme_dartR()
  }
  
  population <- D  <- meane <- NULL
  
  
  # get Tajima's D (code adopted from Paradis, E. (2010) -- pegas package -- written by Renee Catullo)
  get_tajima_D <- function(x){
    # Find allele frequencies (p1 and p2) for every locus in every population
    allele_freqs <- utils.get.allele.freq(x)
    names(allele_freqs)[names(allele_freqs) == "frequency"] <- "p1"
    allele_freqs$p1 <- allele_freqs$p1 / 100
    allele_freqs$p2 <- 1 - allele_freqs$p1
    
    # Get the names of all the populations
    pops <- unique(allele_freqs$popn)
    
    #split each population
    allele_freqs_by_pop <- split(allele_freqs, allele_freqs$popn)
    
    # Internal function to calculate pi
    calc_pi <- function(allele_freqs) {
      n = allele_freqs$nobs * 2  # vector of n values
      pi_sqr <- allele_freqs$p1 ^ 2 + allele_freqs$p2 ^ 2
      h = (n / (n - 1)) * (1 - pi_sqr) # vector of values of h
      sum(h,na.rm = T) # return pi, which is the sum of h across loci
    }
    
    get_tajima_D_for_one_pop <- function(allele_freqs_by_pop) {
      pi <- calc_pi(allele_freqs_by_pop)
      
      #Calculate number of segregating sites, ignoring missing data (missing data will not appear in the allele freq calculations)
      #S <- sum(!(allele_freqs_by_pop$p1 == 0 | allele_freqs_by_pop$p1 == 1))
      S <- sum(allele_freqs_by_pop$p1 >0 & allele_freqs_by_pop$p1 <1,na.rm = T)
      if(S == 0) {
        warning("No segregating sites")
        data.frame(pi = NaN, 
                   S = NaN, 
                   D = NaN, 
                   Pval.normal = NaN, 
                   Pval.beta = NaN)
      }
      
      n <- mean(allele_freqs_by_pop$nobs * 2 )
      
      tmp <- 1:(n - 1)
      a1 <- sum(1/tmp)
      a2 <- sum(1/tmp^2)
      b1 <- (n + 1)/(3 * (n - 1))
      b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
      c1 <- b1 - 1/a1
      c2 <- b2 - (n + 2)/(a1 * n) + a2/a1^2
      e1 <- c1/a1
      e2 <- c2/(a1^2 + a2)
      
      # calculate D and do beta testing
      D <- (pi - S/a1) / sqrt(e1 * S + e2 * S * (S - 1))
      Dmin <- (2/n - 1/a1)/sqrt(e2)
      Dmax <- ((n/(2*(n - 1))) - 1/a1)/sqrt(e2)
      tmp1 <- 1 + Dmin * Dmax
      tmp2 <- Dmax - Dmin
      a <- -tmp1 * Dmax/tmp2
      b <- tmp1 * Dmin/tmp2
      p <- pbeta((D - Dmin)/tmp2, b, a)
      p <- ifelse(p < 0.5, 2 * p, 2 * (1 - p))
      
      data.frame(pi = pi, 
                 S = S, 
                 D = D, 
                 Pval.normal = 2 * pnorm(-abs(D)), 
                 Pval.beta = p,
                 N=n/2)
    }
    
    output <- do.call("rbind", lapply(allele_freqs_by_pop, 
                                      get_tajima_D_for_one_pop))
    data.frame(population = rownames(output), output, row.names = NULL)
  }
  
  tmp_tajD <- get_tajima_D(x)
  tmp_tajD$theta_per_site <- tmp_tajD$pi/x@n.loc
  
  
  
  
  
  
  old.path <- getwd()
  setwd(tempd)
  on.exit(setwd(old.path))
  
  # RUN MS in tempd when ms.path is provided
  
  if (!is.null(ms.path)) {
    tmp_tajD$sim_pval <- NA
    sim_sum <- NULL
    # check OS
    os <- tolower(Sys.info()['sysname'] )
    
    progs<- c("ms", "sample_stats")
    if (os=="windows") progs <- paste0(progs,".exe") 
    fex <- file.exists(file.path(ms.path, progs))
    
    if (all(fex)) {
      file.copy(file.path(ms.path,progs),
                to = tempd,
                overwrite = TRUE, recursive = TRUE)
    } else {
      cat("  Cannot find ",progs[!fex],
          "in the specified folder given by ms.path:",
          ms.path,
          "\n, please download it e.g using gl.download.binary('ms')")
      stop()
    }
    
    
    #simulation summary and p-value
    for (p in unique(x$pop)){
      
      
      if (is.null(seeds)) {
        seeds <- sample(1:10000, 3)
      }
      # Run the command depending on the operating system
      if (.Platform$OS.type == "windows") {
        # On Windows, use shell() which handles piping through the Windows command interpreter
        cmd <- paste0( "ms ",
                       tmp_tajD[which(tmp_tajD$population==p), 'N'], 
                       " ", rep, " ", "-t ", 
                       tmp_tajD[which(tmp_tajD$population==p), 'theta_per_site'], 
                       " -seed ", paste(seeds, collapse = " "), " -s ", 
                       tmp_tajD[which(tmp_tajD$population==p), 'S'], " | ", 
                       "sample_stats")
        assign(paste0("sim_", p),shell(cmd, intern = TRUE))
      } else {
        cmd <- paste0( "./ms ",
                       tmp_tajD[which(tmp_tajD$population==p), 'N'], 
                       " ", rep, " ", "-t ", 
                       tmp_tajD[which(tmp_tajD$population==p), 'theta_per_site'],
                       " -seed ", paste(seeds, collapse = " "), " -s ", 
                       tmp_tajD[which(tmp_tajD$population==p), 'S'], " | ", 
                       "./sample_stats")
        
        # On Linux, use system() which passes the command to the shell
        assign(paste0("sim_", p),system(cmd, intern = TRUE))
      }
      
      write.table(get(paste0("sim_", p)), file.path(tempd, paste0("MS_sim_", p, ".txt")), row.names = F, col.names = F, quote = F)
      assign(paste0("sim_taj_", p), as.numeric(unlist(str_split_fixed(get(paste0("sim_", p)), "\t",12))[,c(6)]))
      tmp_tajD[which(tmp_tajD$population==p), 'sim_pval'] <- 
        mean(abs(get(paste0("sim_taj_", p))) >= abs(tmp_tajD[which(tmp_tajD$population==p), 'D']))
      sim_sum <- cbind(sim_sum, as.numeric(unlist(str_split_fixed(get(paste0("sim_", p)), "\t",12))[,6]))
    }
    
    colnames(sim_sum) <- unique(x$pop)
    
    
    # plot Tajima's D distribution
    plot.list <- list()
    for (p in unique(x$pop)){
      plot.list[[p]] <- ggplot(sim_sum, aes(x=get(p))) + geom_histogram(bins = 20) + 
        geom_vline(xintercept = tmp_tajD[which(tmp_tajD$population==p), 'D'], colour="red") +
        plot_theme +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
        labs(subtitle = p)
    }
    
    sim_plot <- plot.list %>% purrr::map(function(x) {
      ggplot2::ggplot_gtable(ggplot2::ggplot_build(x))
    })
    maxWidth <- do.call(grid::unit.pmax, purrr::map(sim_plot, function(x) x$widths[2:3]))
    for (i in 1:length(sim_plot)) sim_plot[[i]]$widths[2:3] <- maxWidth
    sim_plot$bottom <- "simulated Tajima's D"
    sim_plot$left <- "count"
    
    
    # if choose to output the simulation from MS
    if (!is.null(simulation.out)) {
      for (p in unique(x@pop)) {
        file.copy(file.path(tempd, paste0("MS_sim_", p, ".txt")),
                  to = simulation.out,
                  overwrite = FALSE, recursive = TRUE)
      }
    }
  }
  
  p1 <-
    ggplot(tmp_tajD, aes(
      x = population,
      y = D,
      fill = population
    )) + geom_bar(position = "dodge",
                  stat = "identity",
                  color = "black") + plot_theme + theme(
                    axis.ticks.y = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position = "none"
                  ) +
    labs(fill = "Population", x="Populations") +
    ggtitle("Tajima's D by Population")
  
  if (plot.out & exists("sim_plot")) {
    print(p1)
    do.call(gridExtra::grid.arrange, sim_plot)
  } else if (plot.out) {
    print(p1)
  }
  
  if (!is.null(plot.file) & exists("plot.list")) {
    tmp <- utils.plot.save(plot.list,
                           dir = plot.dir,
                           file = paste0(plot.file, "_distribution"),
                           verbose = verbose
    )
    tmp2 <- utils.plot.save(p1,
                            dir = plot.dir,
                            file = paste0(plot.file, "_TajimasD"),
                            verbose = verbose
    )
  } else if (!is.null(plot.file) & exists("plot.list")==F) {
    tmp <- utils.plot.save(p1,
                           dir = plot.dir,
                           file = paste0(plot.file, "_TajimasD"),
                           verbose = verbose
    )
  }
  
  return(invisible(tmp_tajD))
  setwd(old.path)
  if (cleanup) unlink(tempd, recursive = T)
}
