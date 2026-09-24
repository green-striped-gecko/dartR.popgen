#' @name gl.sfs
#' @title Creates a site frequency spectrum based on a dartR or genlight object
#' @family demographic history
#' @description
#' Creates a folded or unfolded site frequency spectrum (SFS), for all
#' individuals together or as a multidimensional (joint) SFS with one
#' dimension per population.
#' @details
#' Only loci scored in every individual are used. A locus with missing calls
#' cannot be placed in a frequency class of the full sample: counted against
#' the full sample size it lands in a lower class, and a locus fixed in all
#' called individuals would appear polymorphic. Loci with missing calls are
#' therefore excluded and their number is reported (verbose >= 1); to keep
#' more loci, filter individuals with low call rate or impute
#' (\code{gl.impute}) before running the function.\cr
#' A folded multidimensional SFS is folded on the minor allele of all
#' populations combined (the joint minor allele frequency SFS of
#' fastsimcoal2): for each locus, if the alternative allele is the more
#' frequent one over all individuals, the counts of every population are
#' taken for the reference allele.\cr
#' For a multidimensional SFS, minbinsize sets to zero the cells whose total
#' count over all populations is below minbinsize (for minbinsize = 1, only
#' the cell of loci monomorphic in all populations); the array keeps its
#' full dimensions, so polymorphisms private to one population are kept.
#' @param x dartR/genlight object with SNP data [required].
#' @param minbinsize remove bins from the left of the sfs. For example to remove
#'  singletons (alleles only occurring once among all individuals) set
#'  minbinsize to 2. If set to zero, also monomorphic (d0) loci are returned.
#'  For a multidimensional sfs see details [default 0].
#' @param folded if set to TRUE (default) a folded sfs (minor allele frequency
#'  sfs) is returned. If set to FALSE then an unfolded (derived allele frequency
#'   sfs) is returned. It is assumed that 0 is homozygote for the reference and
#'   2 is homozygote for the derived allele. So you need to make sure your
#'   coding is correct [default TRUE].
#' @param singlepop switch to force to create a one-dimensional sfs, even
#' though the genlight/dartR object contains more than one population
#' [default FALSE].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension); only a one-dimensional sfs is plotted [default NULL]
#' @param plot.dir Directory in which to save files [default as specified by
#' the global working directory or tempdir()]
#' @param plot.theme User specified theme [default theme_dartR()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return returns a site frequency spectrum, either a one dimensional vector
#' (only a single population in the dartR/genlight object or singlepop=TRUE) or
#' an n-dimensional array (n is the number of populations in the genlight/dartR
#'  object). If the dartR/genlight object consists of several populations the
#'   multidimensional site frequency spectrum for each population is returned
#'    [=a multidimensional site frequency spectrum]. Be aware the
#'    multidimensional spectrum works only for a limited number of population
#'     and individuals [if too high the function stops with an error as the
#'     number of populations and individuals (and
#'     therefore dimensions) are too large]. To get a single sfs for a
#'     genlight/dartR object with multiple populations, you need to set
#'     singlepop to TRUE. The returned sfs can be used to analyse demographics,
#'      e.g. using fastsimcoal2.
#' @export
#' @examples
#' gl.sfs(bandicoot.gl, singlepop = TRUE)
#' gl.sfs(possums.gl[c(1:5, 31:33), ], minbinsize = 1)
#' @references Excoffier L., Dupanloup I., Huerta-Sanchez E., Sousa V. C. and
#'  Foll M. (2013) Robust demographic inference from genomic and SNP data. PLoS
#'  genetics 9(10)
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber & Carlo Pacioni --
#' Post to \url{https://groups.google.com/d/forum/dartr}


gl.sfs <- function(x,
                   minbinsize = 0,
                   folded = TRUE,
                   singlepop = FALSE,
                   plot.out = TRUE,
                   plot.file = NULL,
                   plot.dir = NULL,
                   plot.theme = theme_dartR(),
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

  if (datatype == "SilicoDArT") {
    stop(error("Fatal Error: Detected Presence/Absence (SilicoDArT) data. Please provide a SNP dataset\n"))
  }

  # only a single population....
  if (nPop(x) == 0) {
    if (verbose >= 2) {
      cat(
        warn(
          "  No population definition provided. I proceed, assuming your genlight/dartR object is a single population.\n"
        )
      )
    }
    pop(x) <- rep("A", nInd(x))
  }

  if (!singlepop & (prod(table(pop(x)) * 2 + 1)) > 2^30) {
    stop(
      error(
        "Cannot create a multidimensional sfs, due to too high dimensions. Reduce the number of populations/individuals or use singlepop=TRUE.\n"
      )
    )
  }

  # loci with missing calls cannot be placed in a class of the full sample
  mat <- as.matrix(x)
  miss <- colSums(is.na(mat)) > 0
  if (any(miss)) {
    if (verbose >= 1) {
      cat(
        warn(
          "  Warning:", sum(miss), "of", length(miss),
          "loci have missing calls and were excluded from the sfs. Filter",
          "(gl.filter.callrate) or impute (gl.impute) to keep more loci.\n"
        )
      )
    }
    mat <- mat[, !miss, drop = FALSE]
  }
  if (ncol(mat) == 0) {
    stop(error("Fatal Error: No loci without missing calls; cannot build the sfs.\n"))
  }

  # DO THE JOB
  if (nPop(x) == 1 | singlepop == TRUE) {
    mi <- nInd(x)
    if (!folded) {
      mi <- 2 * mi # double the number of slots...
    }
    cs <- colSums(mat)
    if (folded) {
      sfs0 <- table(mi - (abs(mi - cs)))
    } else {
      sfs0 <- table(cs)
    }
    sfsf <- rep(0, mi + 1)
    sfsf[as.numeric(names(sfs0)) + 1] <- sfs0
    names(sfsf) <- paste0("d", 0:mi)
    # delete minbinsize
    if (minbinsize > 0) {
      sfsf <- sfsf[-c(1:(minbinsize))]
    }
    sfs <- sfsf
    # multidimensional
  } else {
    pops <- pop(x)
    ni <- as.vector(table(pops))
    cs <- lapply(levels(pops), function(p) {
      colSums(mat[pops == p, , drop = FALSE])
    })
    if (folded) {
      # fold on the minor allele of all populations combined
      flip <- Reduce(`+`, cs) > sum(ni)
      cs <- lapply(seq_along(cs), function(i) {
        ifelse(flip, 2 * ni[i] - cs[[i]], cs[[i]])
      })
    }
    msfs0 <- do.call(table, lapply(seq_along(cs), function(i) {
      factor(cs[[i]], levels = 0:(2 * ni[i]))
    }))

    aa <- array(0, dim = table(pop(x)) * 2 + 1)
    dimnames(aa) <-
      sapply(dim(aa), function(x) {
        paste0("d", 0:(x - 1))
      }, simplify = F)
    aa[] <- as.numeric(msfs0)

    # delete minbinsize: cells whose total count is below minbinsize
    if (minbinsize > 0) {
      tot <- Reduce(`+`, lapply(seq_along(dim(aa)), function(d) {
        slice.index(aa, d) - 1
      }))
      aa[tot < minbinsize] <- 0
    }

    sfs <- aa
  }

  # PLOT
  gp <- NULL
  if (!is.array(sfs) && (plot.out || !is.null(plot.file))) {
    df <- data.frame(sfs = as.numeric(sfs),
                     names = as.numeric(substr(names(sfs), 2, 100)))
    gp <-
      ggplot(df, aes(x = names, y = sfs)) +
      geom_bar(stat = "identity") +
      xlab("bin") +
      ylab("Frequency") +
      plot.theme
  }
  if (plot.out) {
    if (!is.null(gp)) {
      print(gp)
    } else if (verbose >= 2) {
      cat(report(
        "  The sfs is multidimensional, therefore no plot is returned\n"
      ))
    }
  }

  # Optionally save the plot ---------------------

  if (!is.null(plot.file)) {
    if (!is.null(gp)) {
      tmp <- utils.plot.save(gp,
        dir = plot.dir,
        file = plot.file,
        verbose = verbose
      )
    } else if (verbose >= 1) {
      cat(warn("  Warning: the sfs is multidimensional, no plot is saved\n"))
    }
  }

  # FLAG SCRIPT END

  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(sfs)
}
