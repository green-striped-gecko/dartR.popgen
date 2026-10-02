#' Utility function to run Structure
#'
#' These functions were copied from package strataG, which is no longer on CRAN
#'  (maintained by Eric Archer)
#' @export
#' @author Bernd Gruber (bugs? Post to
#'  \url{https://groups.google.com/d/forum/dartr}); original implementation of
#'  Eric Archer \url{https://github.com/EricArcher/strataG}
#'
#' @param g a gtypes object [see \code{strataG}].
#' @param k.range vector of values to for \code{maxpop} in multiple runs. If set
#'   to \code{NULL}, a single STRUCTURE run is conducted with \code{maxpops}
#'   groups. If specified, do not also specify \code{maxpops}.
#' @param num.k.rep number of replicates for each value in \code{k.range}.
#' @param label label to use for input and output files (currently unused;
#'   runs are labelled \code{k<K>.r<replicate>}).
#' @param delete.files logical. Delete all files when STRUCTURE is finished?
#'   If FALSE, the files are kept in a time-stamped folder under
#'   \code{keep.dir}.
#' @param exec name of executable for STRUCTURE. Defaults to "structure".
#' @param burnin Number of burnin reps [default 10000].
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
#' or not to use strata information, one per individual in the order of
#' \code{ind.names} (or of the sorted ids when \code{ind.names} is NULL).
#' Only used when pop.prior = "usepopinfo" [default NULL].
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
#'  mean alphapriora x alphapriorb, and variance alphapriora x alphapriorb^2
#'  [default TRUE].
#' @param alphamax Maximum for uniform prior on alpha when
#' unifprioralpha = TRUE [default 20].
#' @param alphapriora Parameters of Gamma prior on alpha when
#' unifprioralpha = FALSE [default 0.05].
#' @param alphapriorb Parameters of Gamma prior on alpha when
#' unifprioralpha = FALSE [default 0.001].
#' @param ind.names Individual names in the order they should appear in the
#' results. When given, individuals are passed to STRUCTURE by their index
#' (STRUCTURE truncates labels at 11 characters and rejects spaces) and the
#' names are restored in \code{q.mat} and \code{prior.anc}; \code{q.mat} is
#' returned in this order [default NULL, ids passed as they are].
#' @param keep.dir Directory under which the STRUCTURE files are kept when
#' \code{delete.files = FALSE} [default tempdir()].
#' @param verbose Verbosity: 0, silent; 2, one line per run; 3 and above,
#' STRUCTURE's own output is shown [default 0].
#'
#' @return \describe{ \item{\code{structureRun}}{a list where each element is a
#' list with results from \code{utils.structure.read} and a vector of the filenames
#' used} \item{\code{structureWrite}}{a vector of the filenames used by
#' STRUCTURE} \item{\code{utils.structure.read}}{a list containing: \describe{
#' \item{\code{summary}}{new locus name, which is a combination of loci in
#' group} \item{\code{q.mat}}{data.frame of assignment probabilities for each
#' id} \item{\code{prior.anc}}{list of prior ancestry estimates for each
#' individual where population priors were used} \item{\code{files}}{vector of
#' input and output files used by STRUCTURE} \item{\code{label}}{label for the
#' run} } } }

utils.structure.run <- function(g,
                                k.range,
                                num.k.rep,
                                label,
                                delete.files = TRUE,
                                exec,
                                burnin,
                                numreps,
                                noadmix,
                                freqscorr,
                                randomize,
                                seed,
                                pop.prior,
                                locpriorinit,
                                maxlocprior,
                                gensback,
                                migrprior,
                                pfrompopflagonly,
                                popflag,
                                inferalpha,
                                alpha,
                                unifprioralpha,
                                alphamax,
                                alphapriora,
                                alphapriorb,
                                ind.names = NULL,
                                keep.dir = NULL,
                                verbose = 0) {
  ################################################################

  .alleles2integer <- function(g,
                               min.val = 0) {
    g$data %>%
      dplyr::group_by(.data$locus) %>%
      dplyr::mutate(allele = min.val - 1 + as.integer(factor(.data$allele))) %>%
      dplyr::ungroup()
  }

  ################################################################

  .stackedAlleles <- function(g,
                              alleles2integer = FALSE,
                              na.val = NULL,
                              ...) {
    x <- if (alleles2integer) {
      .alleles2integer(g, ...)
    } else {
      g$data
    }

    if (!is.null(na.val)) {
      x$allele[is.na(x$allele)] <- na.val
    }

    x %>%
      dplyr::arrange(.data$id, .data$locus) %>%
      dplyr::mutate(a = rep(1:g$ploidy, dplyr::n() / g$ploidy)) %>%
      tidyr::pivot_wider(names_from = "locus", values_from = "allele") %>%
      dplyr::rename(allele = "a") %>%
      dplyr::select("id", "stratum", "allele", dplyr::everything())
  }

  ####################################################

  structureWrite <- function(g,
                             label = NULL,
                             maxpops = 1:(dplyr::n_distinct(g$data$stratum)),
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
                             alphapriorb = 0.001) {
    if (!is.null(pop.prior)) {
      if (!pop.prior %in% c("locprior", "usepopinfo")) {
        stop(error("'pop.prior' must be 'locprior' or 'usepopinfo'."))
      }
    }

    # STRUCTURE reads labels up to the first space, so ids are sanitised
    # before popflag is matched to them
    ids <- gsub(" ", "_", unique(g$data$id))

    if (is.null(popflag)) {
      popflag <- rep(1, length(ids))
    }

    if (length(popflag) != length(ids)) {
      stop(error("  'popflag' should be the same length as the number of individuals in 'g'."))
    }
    if (!all(popflag %in% c(0, 1))) {
      stop(error("  All values in 'popflag' must be 0 or 1."))
    }

    if (is.null(names(popflag))) {
      names(popflag) <- ids
    } else {
      names(popflag) <- gsub(" ", "_", names(popflag))
    }

    in.file <- ifelse(is.null(label), "data", paste(label, "data", sep = "_"))
    out.file <- ifelse(is.null(label), "out", paste(label, "out", sep = "_"))
    main.file <- ifelse(is.null(label), "mainparams", paste(label, "mainparams", sep = "_"))
    extra.file <- ifelse(is.null(label), "extraparams", paste(label, "extraparams", sep = "_"))
    mat <- .stackedAlleles(g, alleles2integer = TRUE, na.val = -9) %>%
      dplyr::select(-"allele") %>%
      dplyr::mutate(
        id = gsub(" ", "_", .data$id),
        stratum = as.numeric(factor(.data$stratum)),
        popflag = popflag[.data$id]
      ) %>%
      dplyr::select("id", "stratum", "popflag", dplyr::everything()) %>%
      as.matrix()

    # the marker-name header must follow the column order of the matrix
    write(paste(colnames(mat)[-(1:3)], collapse = " "), file = in.file)

    for (i in 1:nrow(mat)) {
      write(paste(mat[i, ], collapse = " "),
        file = in.file,
        append = TRUE
      )
    }

    main.params <- c(
      paste("MAXPOPS", as.integer(maxpops)),
      paste("BURNIN", as.integer(burnin)),
      paste("NUMREPS", as.integer(numreps)),
      paste("INFILE", in.file),
      paste("OUTFILE", out.file),
      paste("NUMINDS", length(ids)),
      paste("NUMLOCI", length(unique(g$data$locus))),
      "MISSING -9", "LABEL 1", "POPDATA 1",
      "POPFLAG 1", "LOCDATA 0", "PHENOTYPE 0",
      "EXTRACOLS 0", "MARKERNAMES 1"
    )

    main.params <- paste("#define", main.params)

    write(main.params, file = main.file)

    extra.params <- c(
      paste("NOADMIX", as.integer(noadmix)),
      paste("FREQSCORR", as.integer(freqscorr)),
      paste("INFERALPHA", as.integer(inferalpha)),
      paste("ALPHA", as.numeric(alpha)),
      "FPRIORMEAN 0.01", "FPRIORSD 0.05", "LAMBDA 1.0",
      paste("UNIFPRIORALPHA", as.integer(unifprioralpha)),
      paste("ALPHAMAX", as.numeric(alphamax)),
      paste("ALPHAPRIORA", as.numeric(alphapriora)),
      paste("ALPHAPRIORB", as.numeric(alphapriorb)),
      "COMPUTEPROB 1",
      paste("ADMBURNIN", max(0, as.integer(burnin / 2))),
      "ALPHAPROPSD 0.025", "STARTATPOPINFO 0",
      paste("RANDOMIZE", as.integer(randomize)),
      paste("SEED", as.integer(seed)),
      "METROFREQ 10", "REPORTHITRATE 0"
    )

    if (!is.null(pop.prior)) {
      pop.prior <- tolower(pop.prior)
      prior.params <- if (pop.prior == "locprior") {
        c("LOCPRIOR 1", "LOCISPOP 1", paste(
          "LOCPRIORINIT",
          locpriorinit
        ), paste("MAXLOCPRIOR", maxlocprior))
      } else if (pop.prior == "usepopinfo") {
        c(
          "USEPOPINFO 1", paste("GENSBACK", trunc(gensback)),
          paste("MIGRPRIOR", migrprior), paste(
            "PFROMPOPFLAGONLY",
            as.integer(pfrompopflagonly)
          )
        )
      }
      extra.params <- c(extra.params, prior.params)
    }

    extra.params <- extra.params[!is.na(extra.params)]
    extra.params <- paste("#define", extra.params)

    write(extra.params, file = extra.file)

    invisible(list(
      files = c(
        data = in.file,
        mainparams = main.file,
        extraparams = extra.file,
        out = out.file
      ),
      pops = sort(unique(g$data$stratum))
    ))
  }

  ###########################################################
  # restore the individual names and order when ids were passed as an index
  .restoreIds <- function(result, ind.names) {
    if (is.null(ind.names) || is.null(result$q.mat)) {
      return(result)
    }
    idx <- as.integer(result$q.mat$id)
    result$q.mat$id <- ind.names[idx]
    result$q.mat <- result$q.mat[order(idx), , drop = FALSE]
    rownames(result$q.mat) <- NULL
    if (!is.null(result$prior.anc)) {
      anc.idx <- as.integer(names(result$prior.anc))
      names(result$prior.anc) <- ind.names[anc.idx]
      result$prior.anc <- result$prior.anc[order(anc.idx)]
    }
    result
  }

  ###########################################################

  exec <- normalizePath(path.expand(exec), mustWork = FALSE)
  if (!file.exists(exec)) {
    stop(error("  Cannot find the STRUCTURE executable:", exec, "\n"))
  }
  if (file.access(exec, mode = 1) != 0) {
    stop(error("  The STRUCTURE executable is not executable:", exec, "\n"))
  }

  # individuals are passed to STRUCTURE by index when their names are known:
  # STRUCTURE truncates labels at 11 characters and splits them at spaces
  if (!is.null(ind.names)) {
    if (!setequal(ind.names, unique(g$data$id))) {
      stop(error("  'ind.names' do not match the individual ids in 'g'.\n"))
    }
    id.index <- stats::setNames(seq_along(ind.names), ind.names)
    if (!is.null(popflag)) {
      if (is.null(names(popflag))) {
        if (length(popflag) != length(ind.names)) {
          stop(error("  'popflag' should be the same length as the number of individuals.\n"))
        }
        names(popflag) <- ind.names
      }
      names(popflag) <- as.character(id.index[names(popflag)])
    }
    g$data <- as.data.frame(g$data)
    g$data$id <- as.character(id.index[g$data$id])
  }

  # all files live in a per-call directory that is removed on exit; STRUCTURE
  # also writes seed.txt to the current directory, so the run happens there
  run.dir <- tempfile(pattern = "structureRun_")
  dir.create(run.dir)
  old.wd <- setwd(run.dir)
  on.exit(setwd(old.wd), add = TRUE)
  on.exit(unlink(run.dir, recursive = TRUE, force = TRUE), add = TRUE)

  if (is.null(k.range)) {
    k.range <- 1:(dplyr::n_distinct(g$data$stratum))
  }

  rep.df <- expand.grid(rep = 1:num.k.rep, k = k.range)
  rep.df$label <- paste0("k", rep.df$k, ".r", rep.df$rep)
  n.runs <- nrow(rep.df)

  run.result <- lapply(seq_len(n.runs), function(i) {
    run.label <- rep.df$label[i]
    if (verbose >= 2) {
      cat(report(
        "  Running STRUCTURE: K =", rep.df$k[i], ", replicate", rep.df$rep[i],
        paste0("(run ", i, " of ", n.runs, ")\n")
      ))
    }
    sw.out <- structureWrite(g,
      label = run.label,
      maxpops = rep.df$k[i],
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
    )

    files <- sw.out$files
    log.file <- paste(run.label, "log", sep = "_")
    args <- c(
      "-m", shQuote(files[["mainparams"]]),
      "-e", shQuote(files[["extraparams"]]),
      "-i", shQuote(files[["data"]]),
      "-o", shQuote(files[["out"]])
    )
    err.code <- system2(exec,
      args = args,
      stdout = if (verbose >= 3) "" else log.file,
      stderr = if (verbose >= 3) "" else log.file
    )
    if (err.code == 127) {
      stop(error(
        "  STRUCTURE could not be started (exit status 127). Command:\n ",
        paste(shQuote(exec), paste(args, collapse = " ")), "\n"
      ))
    } else if (err.code != 0) {
      tail.log <- if (file.exists(log.file)) {
        utils::tail(readLines(log.file, warn = FALSE), 10)
      } else {
        character()
      }
      stop(error(
        "  STRUCTURE exited with status", err.code, "for run", run.label,
        ".\n  Last lines of its output:\n ",
        paste(tail.log, collapse = "\n  "), "\n"
      ))
    }
    files["out"] <- paste(files["out"], "_f", sep = "")
    # shared with gl.read.structure
    result <- utils.structure.read(files[["out"]], sw.out$pops)
    result <- .restoreIds(result, ind.names)
    c(result, list(files = files, label = run.label))
  })

  names(run.result) <- rep.df$label

  if (delete.files) {
    run.result <- lapply(run.result, function(r) {
      r$files <- NULL
      r
    })
  } else {
    if (is.null(keep.dir)) {
      keep.dir <- tempdir()
    }
    keep.path <- file.path(
      keep.dir,
      paste0("structureRun_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    )
    dir.create(keep.path, recursive = TRUE, showWarnings = FALSE)
    file.copy(list.files(run.dir, full.names = TRUE), keep.path,
      overwrite = TRUE
    )
    run.result <- lapply(run.result, function(r) {
      r$files <- stats::setNames(
        file.path(keep.path, basename(r$files)),
        names(r$files)
      )
      r
    })
    if (verbose >= 2) {
      cat(report("  STRUCTURE files kept in", keep.path, "\n"))
    }
  }

  class(run.result) <- c("structure.result", class(run.result))
  run.result
}
