#' @name gl.LDNe
#' @title Estimates effective population size using the Linkage Disequilibrium
#' method based on NeEstimator (V2)
#' @description
#' This function is basically a convenience function that runs the LD Ne
#'  estimator using Neestimator2
#'  (\url{http://www.molecularfisherieslaboratory.com.au/neestimator-software/})
#'  within R using the provided genlight object. To be able to do so, the
#'  software has to be downloaded from their website and the appropriate
#'  executable Ne2-1 has to be copied into the path as specified in the function
#'  (see example below).
#'  
#' @references \itemize{
#'  \item Waples, R. S. (2006). "A bias correction for estimates of effective 
#'  population size based on linkage disequilibrium at unlinked gene loci*." 
#'  Conservation Genetics 7(2): 167-184.
#'
#'  \item Waples, R. K., et al. (2016). "Estimating contemporary effective 
#'  population size in non-model species using linkage disequilibrium across 
#'  thousands of loci." Heredity 117(4): 233-240.
#'  }
#'  
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file  with
#' all results from Neestimator 2 [default 'genepopLD.txt'].
#' @param outpath Path where to save the output file. Use outpath=getwd() or
#' outpath='.' when calling this function to direct output files to your working
#'  directory [default tempdir(), mandated by CRAN].
#' @param neest.path Path to the folder of the NE2-1 file.
#'  Please note there are 3 different executables depending on your OS:
#'  Ne2-1.exe (=Windows), Ne2-1M (=Mac), Ne2-1L (=Linux). You only need to point
#'  to the folder (the function will recognise which OS you are running)
#'  [default getwd()].
#' @param critical (vector of) Critical values that are used to remove alleles
#' based on their minor allele frequency. This can be done before using the
#' gl.filter.maf function, therefore the default is set to 0 (no loci are
#' removed). To run for MAF 0 and MAF 0.05 at the same time specify: critical =
#' c(0,0.05) [default 0].
#' @param singleton.rm Whether to remove singleton alleles [default TRUE].
#' @param mating Formula for Random mating='random' or monogamy= 'monogamy'
#' [default 'random'].
#' @param pairing 'all' [default] if all possible loci should be paired, or 'separate'
#'    if only loci on different chromosomes should be used.
#' @param Waples.correction The type of Waples et al 2016 correction to apply. 
#'    This is ignored if \code{pairing} is set to 'separate'.
#'    Options are 'nChromosomes', for eq 1a, or 'genomeLength' for eq 1b. 
#'    NULL if none should be applied [default NULL]. 
#' @param Waples.correction.value The number of chromosomes or the genome length 
#'    in cM. See Waples et al 2016 for details.
#' @param naive Whether the naive (uncorrected for samples size - see 
#'    eq 7 and eq 8 in Waples 2006) should also be reported. This is mostly 
#'    to diagnose the source of Inf estimate.
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot_theme User specified theme [default theme_dartR()].
#' @param plot_colors_pop  population colors with as many colors as there are populations in the dataset
#' [default discrete_palette].
#' @param plot.dir Directory in which to save files [default = working directory]
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return Dataframe with the results as table
#' @author Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @importFrom stats weighted.mean
#' @examples
#' \dontrun{
#' # SNP data (use two populations and only the first 100 SNPs)
#' pops <- possums.gl[1:60, 1:100]
#' nes <- gl.LDNe(pops,
#'   outfile = "popsLD.txt", outpath = tempdir(),
#'   neest.path = "./path_to Ne-21",
#'   critical = c(0, 0.05), singleton.rm = TRUE, mating = "random"
#' )
#' nes
#'
#' # Using only pairs of loci on different chromosomes

#' made up some chromosome location
#' pops@chromosome <- as.factor(sample(1:10, size = nLoc(pops), replace = TRUE))
#' nessep <- gl.LDNe(pops,
#'               outfile = "popsLD.txt", outpath = "./TestNe", pairing="separate",
#'               neest.path = "./path_to Ne-21",
#'               critical = c(0, 0.05), singleton.rm = TRUE, mating = "random"
#' nessep
#' }
#' @export

gl.LDNe <- function(x,
                    outfile = "genepopLD.txt",
                    outpath = tempdir(),
                    neest.path = getwd(),
                    critical = 0,
                    singleton.rm = TRUE,
                    mating = "random",
                    pairing = "all",
                    Waples.correction=NULL,
                    Waples.correction.value=NULL, 
                    naive=FALSE,
                    plot.out = TRUE,
                    plot_theme = theme_dartR(),
                    plot_colors_pop = gl.select.colors(x, verbose = 0),
                    plot.file = NULL,
                    plot.dir = NULL,
                    verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(
    func = funname,
    build = "Jody",
    verbose = verbose
  )

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  # works only with SNP data
  if (datatype != "SNP") {
    message(error(
      "  Only SNPs (diploid data can be transformed into genepop format!\n"
    ))
  }
  
  # Correct arg options?
  if(!pairing %in% c("all", "separate")) {
    message(error(
      "  'pairing' can only be either 'all' or 'separate'!\n"
    ))
  }
  
  if(pairing == "separate") Waples.correction <- NULL
  
  if(!is.null(Waples.correction)) {
    if(!Waples.correction %in% c('nChromosomes', 'genomeLength'))
    message(error(
      "  'Waples.correction' can only be either 'nChromosomes' or 'genomeLength'!\n"
    ))
    if(!(is.numeric(Waples.correction.value) & length(Waples.correction.value == 1)))
      message(error(
        "  'Waples.correction.value' should be a numeric vector of length == 1!\n"
      ))
  }

  # DO THE JOB
  
  # Helper FUN to obtain naive Ne estimate
  naiveNe <- function(tmp, bur="dummyBur.txt", matingsys) {
    # Figure out where the data are
    fnBur <- file.path(tmp, bur)
    rl <- readLines(fnBur)
    headings <- grep("Loc1   Loc2   LowP1   LowP2  Samp.Size", rl)
    starts <- headings + 2
    ends <- grep("Total locus pairs investigated", rl) - 2
    headr <-  c("Loc1", "Loc2", "LowP1", "LowP2", "Samp.Size", "Mean_rsq", "rsq_drift")
    
    # do the actual calculations. This return a vector with estimates for each 
    # frequencies threshold, for each pop, in this order
    Ne <- function(i, fn, starts=starts, ends=ends, headr=headr, ms=matingsys) {
      d<-data.table::fread(file = fn, skip = starts[i], 
               nrows = ends[i] - starts[i], col.names = headr)
      
      d[, rsq_sample := 1/Samp.Size]
      d[, pc.rsq_drift := Mean_rsq - rsq_sample]
      
      wmean.rsq_drift <- d[, weighted.mean(x = pc.rsq_drift, w = Samp.Size)]
      num <- ifelse(ms == "random", 1, 2)
      Ne <- num/(3*wmean.rsq_drift)
      return(Ne)
    }
    
    Ne <- sapply(seq_along(starts), Ne, fn=fnBur,
                 starts=starts, ends=ends, headr=headr)
    return(Ne)
    
  }
  
  #-------End helper FUN----------------#
  # Set NULL to variables to pass CRAN checks
  "Lowest Allele Frequency Used" <- "CI high Parametric" <- "CI low Parametric" <- "Estimated Ne^" <- NULL
  rsq_sample <- Samp.Size <- pc.rsq_drift <- Mean_rsq <- naiveNe <- Mean_rsq <- Samp.Size <- pc.rsq_drift <- rsq_sample <- NULL

  xx <- gl2genepop(x, outfile = "dummy.gen", outpath = tempdir())

  if (singleton.rm == TRUE) {
    critical[length(critical) + 1] <- 1
  }

  # copy info file to tempdir
  info <- NA
  info[1] <- "1"
  info[2] <- "./" # path of input file
  info[3] <- "dummy.gen" # input file
  info[4] <- 2 # Genepop format
  info[5] <- "./" # path of output file
  info[6] <- outfile # output file
  info[7] <- length(critical)
  info[8] <- paste(critical, collapse = " ")

  mm <- pmatch(mating, c("random", "mono")) - 1
  if (mm == 0 | mm == 1) {
    info[9] <- mm
  } else {
    cat(error("  Mating is not either 'random' or 'monogamy'. Please check\n"))
    stop()
  }

  con <- file(file.path(tempdir(), "infodummy"), "w")
  writeLines(info, con)
  close(con)
  
  # set the pairing option and generate the map file if necessary 
  if(pairing == "all") {
    setPairs <- 0
  } else {
    if(pairing == "separate") {
      setPairs <- "2 ChrMap"
      write.table(
      data.frame(x@chromosome, locNames(x)), 
      file = file.path(tempdir(), "ChrMap"), 
      row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
  # copy option file to tempdir
  option <- NA
  option[1] <- paste(c(1,0, length(critical), 0), collapse = " ")
  option[2] <- 0 # Maximum individuals/pop. If 0: no limit
  option[3] <- -1 # -1: Freq. output up to population 50
  option[4] <- ifelse(naive, -1, 0) # Burrow output 0: No output; -1: first 50 pop
  option[5] <- 1 # Parameter CI: 1 for Yes
  option[6] <- 1 # Jackknife CI: 1 for Yes
  option[7] <- 0 # Up to population. 0: no restriction
  option[8] <- 0 # All loci accepted
  option[9] <- 0 # No file with missing data summary
  option[10] <- setPairs  # 0: no pairing restriction; 1: loci within same chrs; 2: loci on separate Chrs
  
  con <- file(file.path(tempdir(), "option"), "w")
  writeLines(option, con)
  close(con)

  if (Sys.info()["sysname"] == "Windows") {
    prog <- "Ne2-1.exe"
    cmd <- "Ne2-1.exe i:infodummy o:option"
  }

  if (Sys.info()["sysname"] == "Linux") {
    prog <- "Ne2-1L"
    cmd <- "./Ne2-1L i:infodummy o:option"
  }

  if (Sys.info()["sysname"] == "Darwin") {
    prog <- "Ne2-1M"
    cmd <- "./Ne2-1M i:infodummy o:option"
  }

  # check if file program can be found
  if (file.exists(file.path(neest.path, prog))) {
    file.copy(file.path(neest.path, prog),
      to = tempdir(),
      overwrite = TRUE
    )
  } else {
    cat(
      error(
        "  Cannot find",
        prog,
        "in the specified folder given by neest.path:",
        neest.path,
        "\n"
      )
    )
    stop()
  }

  # change into tempdir (run it there)
  old.path <- getwd()
  setwd(tempdir())
  on.exit(setwd(old.path))
  system(cmd)
  res <- read.delim(outfile)
  res <-
    unlist(lapply(res[, 1], function(x) {
      x <- gsub(pattern = "Infinite", replacement = "Inf", x)
    }))

  pops <-
    sapply(res[res %like% "Population"], function(x) {
      str_extract(x, "(?<=\\[).*(?=\\])")
    }, USE.NAMES = F)
  pops <- sub("_[^_]+$", "", pops)
  freq <- str_split(res[res %like% "Lowest"], "\\s{3,}")[[1]][-1]
  Estimated_Ne <-
    lapply(res[res %like% "Estimated"], function(x) {
      strsplit(x, "\\s{3,}")[[1]][-1]
    })
  CI_low_Parametric <-
    lapply(res[res %like% "* Parametric"], function(x) {
      strsplit(x, "\\s{3,}")[[1]][-1]
    })
  CI_high_Parametric <-
    lapply(res[grep("^\\* Parametric", res) + 1], function(x) {
      strsplit(x, "\\s{3,}")[[1]][-1]
    })
  CI_low_JackKnife <-
    lapply(res[res %like% "* JackKnife"], function(x) {
      strsplit(x, "\\s{3,}")[[1]][-1]
    })
  CI_high_JackKnife <-
    lapply(res[grep("^\\* JackKnife", res) + 1], function(x) {
      strsplit(x, "\\s{3,}")[[1]][-1]
    })

  # JackKnife works only with more than 2 individuals

  ind_threshold <- which(table(pop(x)) < 3)

  if (length(ind_threshold) > 0) {
    for (r in unname(ind_threshold)) {
      CI_low_JackKnife <- c(CI_low_JackKnife[1:(r - 1)], NA, CI_low_JackKnife[r:length(CI_low_JackKnife)])
      CI_high_JackKnife <- c(CI_high_JackKnife[1:(r - 1)], NA, CI_high_JackKnife[r:length(CI_high_JackKnife)])
    }
  }

  harmonic_mean <-
    lapply(res[res %like% "Harmonic"], function(x) {
      strsplit(x, "\\s{3,}")[[1]][-1]
    })
  comparisons <-
    lapply(res[res %like% "Independent"], function(x) {
      strsplit(x, "\\s{3,}")[[1]][-1]
    })
  overall_r2 <-
    lapply(res[res %like% "OverAll"], function(x) {
      strsplit(x, "\\s{3,}")[[1]][-1]
    })
  expected_r2 <-
    lapply(res[res %like% "Expected"], function(x) {
      strsplit(x, "\\s{3,}")[[1]][-1]
    })

  pop_list <- lapply(1:length(pops), function(i) {
    df_temp <- as.data.frame(cbind(
      c(
        "Lowest Allele Frequency Used",
        "Harmonic Mean Sample Size",
        "Independent Comparisons",
        "OverAll r^2",
        "Expected r^2 Sample",
        "Estimated Ne^",
        "CI low Parametric",
        "CI high Parametric",
        "CI low JackKnife",
        "CI high JackKnife"
      ),
      rbind(
        freq,
        as.numeric(harmonic_mean[[i]]),
        as.numeric(comparisons[[i]]),
        as.numeric(overall_r2[[i]]),
        as.numeric(expected_r2[[i]]),
        as.numeric(Estimated_Ne[[i]]),
        as.numeric(CI_low_Parametric[[i]]),
        as.numeric(CI_high_Parametric[[i]]),
        as.numeric(CI_low_JackKnife[[i]]),
        as.numeric(CI_high_JackKnife[[i]])
      )
    ))

    df_temp <- df_temp[!duplicated(as.list(df_temp))]

    colnames(df_temp) <-
      c("Statistic", paste("Frequency", 1:sum(!duplicated(freq))))
    rownames(df_temp) <- 1:nrow(df_temp)
    return(df_temp)
  })

  names(pop_list) <- pops

  file.copy(outfile, file.path(outpath, outfile))
  setwd(old.path)
  
  # Apply correction if relevant
  cr.est <- function(pop, crtn, fr) {
    # Pull out the numerical values and apply correction
    m <- round(matrix(as.numeric(as.matrix(pop[6:nrow(pop), -1])) / crtn, 
             nrow = nrow(pop) - 5), digits = 1)
    df <- data.frame(m)
    names(df) <- paste("Frequency", 1:sum(!duplicated(fr)))
    # labels
    lbs <- c("Waples' corrected Ne",
             "Waples' corrected CI low Parametric",
             "Waples' corrected CI high Parametric",
             "Waples' corrected CI low JackKnife",
             "Waples' corrected CI high JackKnife")
    # append to existing df
    res <- rbind(pop,
      cbind(Statistic=lbs, df))
    return(res)
  } 
  
  if(!is.null(Waples.correction)) {
    if(Waples.correction == "nChromosomes") {
      crc <- 0.098 + 0.219 * log(Waples.correction.value)
    } else {
      crc <- -0.910 + 0.219 * log(Waples.correction.value)
    }
    pop_list <- lapply(pop_list, cr.est, crtn=crc, fr=freq)
  }
  # --- Apply correction if relevant END --- #
  #------------------------------------------#
  
  # Naive Ne estimates
  if(naive) {
    nNe <- naiveNe(tmp = tempdir(), matingsys = mating)
    pop_list <- lapply(seq_along(pop_list), function(i, fr=freq, nPops=length(pop_list)){
    nValuesPop <- length(nNe) / nPops  
    v <- nNe[((i - 1)*nValuesPop + 1):(i*nValuesPop)]
    v <- round(v[!duplicated(fr)], 1)
    m <- matrix(v, ncol=sum(!duplicated(fr)), byrow = TRUE)
    tmpdf <- data.frame(m)
    names(tmpdf) <- paste("Frequency", 1:sum(!duplicated(fr)))
    tmpdf2 <- data.frame(Statistic="Naive Estimated Ne^")
    updated <- rbind(pop_list[[i]], cbind(tmpdf2, tmpdf))
      return(updated)
    })
  }
  
  # PLOTS
  if (plot.out) {
    # printing plots and reports assigning colors to populations

    pop_list_plot <- lapply(pop_list, function(x) {
      stats::setNames(data.frame(t(x[, -1])), x[, 1])
    })

    pop_list_plot <- lapply(1:length(pops), function(i) {
      pop_temp <- pop_list_plot[[i]]
      pop_temp$pop <- pops[i]
      return(pop_temp)
    })

    pop_list_plot <- as.data.frame(rbindlist(pop_list_plot))
    pop_list_plot$pop <- as.factor(pop_list_plot$pop)
    pop_list_plot[pop_list_plot == Inf] <- NA
    pop_list_plot$color <- rep(plot_colors_pop, each = sum(!duplicated(freq)))
    pop_list_plot$`CI low Parametric` <-
      as.numeric(pop_list_plot$`CI low Parametric`)
    pop_list_plot$`CI high Parametric` <-
      as.numeric(pop_list_plot$`CI high Parametric`)
    pop_list_plot$`Estimated Ne^` <-
      as.numeric(pop_list_plot$`Estimated Ne^`)

    pop_size <- unlist(lapply(harmonic_mean, function(x) {
      mean(as.numeric(x), na.rm = T)
    }))

    p3 <-
      ggplot(data = pop_list_plot, aes(
        x = pop,
        y = `Estimated Ne^`
      )) +
      geom_bar(
        position = "dodge2",
        stat = "identity",
        color = "black",
        fill = pop_list_plot$color
      ) +
      scale_x_discrete(labels = paste(unique(pop_list_plot$pop),
        round(
          pop_size,
          0
        ),
        sep = " | "
      )) +
      geom_errorbar(
        aes(
          ymin = `CI low Parametric`,
          ymax = `CI high Parametric`
        ),
        position = position_dodge2(padding = 0.5)
      ) +
      geom_text(
        aes(label = `Lowest Allele Frequency Used`),
        position = position_dodge2(width = 0.9),
        stat = "identity",
        vjust = -0.50,
        size = 4,
        fontface = "bold"
      ) +
      plot_theme +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          face = "bold",
          size = 12
        ),
        axis.title.x = element_blank(),
        legend.position = "none"
      ) +
      ylab("Estimated Ne") +
      ggtitle("Effective population size (Ne) by Population")
  }

  # PRINTING OUTPUTS
  if (plot.out) {
    print(p3)
  }

  print(pop_list, row.names = FALSE)

  # Optionally save the plot ---------------------

  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p3,
      dir = plot.dir,
      file = plot.file,
      verbose = verbose
    )
  }
  # save also table (automatically if plot is not null)
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(pop_list,
      dir = plot.dir,
      file = paste0(plot.file, "_tab"),
      verbose = verbose
    )
  }



  if (verbose >= 1) {
    cat(report(
      "  The results are saved in:",
      file.path(outpath, outfile),
      "\n"
    ))
  }

  # FLAG SCRIPT END

  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN
  return(invisible(pop_list))
}
