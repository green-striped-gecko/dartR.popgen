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
#' @param mating Mating system assumed by the estimator: 'random' or
#' 'monogamy'. The abbreviation 'mono' is also accepted
#' [default 'random'].
#' @param pairing 'all' [default] if all possible loci should be paired, or 'separate'
#'    if only loci on different chromosomes should be used.
#' @param Waples.correction The type of Waples et al 2016 correction to apply. 
#'    This is ignored if \code{pairing} is set to 'separate'.
#'    Options are 'nChromosomes', for eq 1a, or 'genomeLength' for eq 1b. 
#'    NULL if none should be applied [default NULL]. 
#' @param Waples.correction.value A single positive number: the number of
#'    chromosomes for 'nChromosomes', or the genome length in cM for
#'    'genomeLength'. Required when \code{Waples.correction} is set. See
#'    Waples et al 2016 for details [default NULL].
#' @param naive Whether the naive (uncorrected for samples size - see 
#'    eq 7 and eq 8 in Waples 2006) should also be reported. This is mostly 
#'    to diagnose the source of Inf estimate.
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot_theme User specified theme [default theme_dartR()].
#' @param plot_colors_pop Either a vector with at least one colour per
#' population, or a palette function taking the number of populations
#' [default gl.select.colors(x)].
#' @param plot.dir Directory in which to save the RDS files
#' [default tempdir(), mandated by CRAN; see gl.check.wd].
#' @param plot.file Name for the RDS binary files to save (base name only,
#' exclude extension). The plot is saved under this name and the table under
#' the same name with the suffix '_tab'. Nothing is saved when NULL
#' [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return Invisibly, a named list with one data frame per population. Each
#' data frame holds one row per statistic (lowest allele frequency used,
#' harmonic mean sample size, independent comparisons, overall r^2, expected
#' r^2, estimated Ne and its parametric and jackknife confidence limits) and
#' one column per allele frequency threshold. Five rows are added when
#' \code{Waples.correction} is set and one when \code{naive} is TRUE. The
#' full NeEstimator output is written to \code{file.path(outpath, outfile)}.
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @importFrom stats weighted.mean
#' @examples
#' \dontrun{
#' # SNP data (use two populations and only the first 100 SNPs)
#' pops <- possums.gl[1:60, 1:100]
#' nes <- gl.LDNe(pops,
#'   outfile = "popsLD.txt", outpath = tempdir(),
#'   neest.path = "./path_to_Ne2-1",
#'   critical = c(0, 0.05), singleton.rm = TRUE, mating = "random"
#' )
#' nes
#'
#' # Using only pairs of loci on different chromosomes
#' # make up some chromosome locations
#' pops@chromosome <- as.factor(sample(1:10, size = nLoc(pops), replace = TRUE))
#' nessep <- gl.LDNe(pops,
#'   outfile = "popsLD.txt", outpath = tempdir(), pairing = "separate",
#'   neest.path = "./path_to_Ne2-1",
#'   critical = c(0, 0.05), singleton.rm = TRUE, mating = "random"
#' )
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
    verbose = verbose
  )

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  # works only with SNP data
  if (datatype != "SNP") {
    stop(error(
      "  Only SNP (diploid) data can be transformed into genepop format!\n"
    ))
  }

  # Correct arg options? An invalid value used to be reported with a message
  # and the run continued, failing later or returning a silently wrong table.
  if (!is.character(pairing) || length(pairing) != 1) {
    stop(error("  'pairing' can only be either 'all' or 'separate'!\n"))
  }
  pairing <- match.arg(pairing, c("all", "separate"))

  # 'mono' is accepted as an abbreviation of the documented 'monogamy'
  if (identical(mating, "mono")) {
    mating <- "monogamy"
  }
  if (!is.character(mating) || length(mating) != 1) {
    stop(error("  'mating' can only be either 'random' or 'monogamy'!\n"))
  }
  mating <- match.arg(mating, c("random", "monogamy"))

  if (pairing == "separate") Waples.correction <- NULL

  if (!is.null(Waples.correction)) {
    if (!is.character(Waples.correction) || length(Waples.correction) != 1 ||
        !Waples.correction %in% c("nChromosomes", "genomeLength")) {
      stop(error(
        "  'Waples.correction' can only be either 'nChromosomes' or",
        "'genomeLength', or NULL for no correction.\n"
      ))
    }
    if (!is.numeric(Waples.correction.value) ||
        length(Waples.correction.value) != 1 ||
        is.na(Waples.correction.value) ||
        Waples.correction.value <= 0) {
      stop(error(
        "  'Waples.correction.value' should be a single positive number: the",
        "number of chromosomes for 'nChromosomes', or the genome length in cM",
        "for 'genomeLength'.\n"
      ))
    }
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
  rsq_sample <- Samp.Size <- pc.rsq_drift <- Mean_rsq <- Mean_rsq <- Samp.Size <- pc.rsq_drift <- rsq_sample <- NULL

  # Each call runs in its own directory under tempdir(). NeEstimator reads and
  # writes files with fixed names, so concurrent calls - forked parallel runs
  # share the parent's tempdir() - would otherwise overwrite each other's
  # input and output files.
  run.dir <- tempfile("LDNe_")
  dir.create(run.dir)
  on.exit(unlink(run.dir, recursive = TRUE), add = TRUE)

  # resolve outpath before the working directory changes, so that a relative
  # path such as '.' means the caller's working directory, as documented
  outpath <- normalizePath(outpath, winslash = "/", mustWork = FALSE)

  xx <- gl2genepop(x,
    outfile = "dummy.gen", outpath = run.dir,
    verbose = if (verbose >= 3) verbose else 0
  )

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

  # NeEstimator expects 0 for random mating and 1 for monogamy
  info[9] <- ifelse(mating == "random", 0, 1)

  con <- file(file.path(run.dir, "infodummy"), "w")
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
      file = file.path(run.dir, "ChrMap"), 
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
  
  con <- file(file.path(run.dir, "option"), "w")
  writeLines(option, con)
  close(con)

  sysname <- unname(Sys.info()["sysname"])
  if (sysname == "Windows") {
    prog <- "Ne2-1.exe"
    cmd <- "Ne2-1.exe i:infodummy o:option"
  } else if (sysname == "Linux") {
    prog <- "Ne2-1L"
    cmd <- "./Ne2-1L i:infodummy o:option"
  } else if (sysname == "Darwin") {
    prog <- "Ne2-1M"
    cmd <- "./Ne2-1M i:infodummy o:option"
  } else {
    stop(error(
      "  gl.LDNe does not know which NeEstimator executable to use on",
      sysname, ". The supported systems are Windows, Linux and Darwin",
      "(macOS).\n"
    ))
  }

  # check if file program can be found
  if (!file.exists(file.path(neest.path, prog))) {
    stop(error(
      "  Cannot find",
      prog,
      "in the specified folder given by neest.path:",
      neest.path,
      "\n"
    ))
  }
  file.copy(file.path(neest.path, prog),
    to = run.dir,
    overwrite = TRUE
  )

  # change into the run directory (run it there)
  old.path <- getwd()
  setwd(run.dir)
  on.exit(setwd(old.path), add = TRUE, after = FALSE)
  status <- system(cmd, ignore.stdout = verbose < 3)
  # a failed run must not be read from a file left by something else
  if (status != 0 || !file.exists(outfile)) {
    stop(error(
      "  NeEstimator (", prog, ") did not produce its output file",
      outfile, "and exited with status", status,
      ". Run with verbose = 3 to see its messages.\n"
    ))
  }
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

  # NeEstimator prints no jackknife line for such a population, so a
  # placeholder is inserted at its position. append() is used because
  # CI[1:(r - 1)] with r == 1 is c(1, 0), which selects the first element of
  # the next population instead of selecting nothing.
  if (length(ind_threshold) > 0) {
    for (r in unname(ind_threshold)) {
      CI_low_JackKnife <- append(CI_low_JackKnife, NA, after = r - 1)
      CI_high_JackKnife <- append(CI_high_JackKnife, NA, after = r - 1)
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

  saved <- file.copy(outfile, file.path(outpath, outfile), overwrite = TRUE)
  setwd(old.path)
  if (!saved && verbose >= 1) {
    cat(warn(
      "  Warning: could not copy the NeEstimator output to",
      file.path(outpath, outfile), "\n"
    ))
  }
  
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
    nNe <- naiveNe(tmp = run.dir, matingsys = mating)
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
    # lapply over seq_along drops the names, which callers index by
    names(pop_list) <- pops
  }
  
  # PLOTS. The plot is built whenever it is shown or saved: plot.file alone
  # used to reach utils.plot.save() with no plot object in scope.
  build.plot <- plot.out || !is.null(plot.file)

  if (build.plot) {
    # printing plots and reports assigning colors to populations

    # accept a palette function as well as a vector, and take one colour per
    # population; anything longer used to fail on the assignment below
    if (is.function(plot_colors_pop)) {
      plot_colors_pop <- plot_colors_pop(length(pops))
    }
    if (length(plot_colors_pop) < length(pops)) {
      stop(error(
        "  plot_colors_pop must supply at least one colour per population:",
        length(pops), "populations,", length(plot_colors_pop),
        "colours given.\n"
      ))
    }
    plot_colors_pop <- plot_colors_pop[seq_len(length(pops))]

    pop_list_plot <- lapply(pop_list, function(x) {
      stats::setNames(data.frame(t(x[, -1])), x[, 1])
    })

    pop_list_plot <- lapply(1:length(pops), function(i) {
      pop_temp <- pop_list_plot[[i]]
      pop_temp$pop <- pops[i]
      return(pop_temp)
    })

    pop_list_plot <- as.data.frame(rbindlist(pop_list_plot))
    pop_list_plot$pop <- factor(pop_list_plot$pop)#,levels= pop_list_plot$pop)
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

  if (verbose >= 2) {
    print(pop_list, row.names = FALSE)
  }

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
