#' @name gl.run.assignpop
#' @title Runs assignPOP population assignment directly from a genlight object
#' @family linker
#'
#' @description
#' Provides a seamless interface between dartRverse genlight objects and the
#' \pkg{assignPOP} package, handling all data preparation and assignment steps
#' internally.
#'
#' \strong{Single-object mode} (\code{unknown.id = NULL}, \code{x.unknown =
#' NULL}): \code{x} contains individuals with known population assignments.
#' The returned object is structurally identical to the output of
#' \code{assignPOP::read.Structure()} and can be passed directly to
#' \code{assign.kfold()}, \code{assign.MC()}, or \code{assign.matrix()} for
#' cross-validation.
#'
#' \strong{Split mode} (\code{unknown.id} supplied): the simplest mode — a
#' single genlight object containing both reference and unknown individuals.
#' The unknown individuals are identified by name, extracted automatically,
#' and the remainder becomes the reference set. Reference populations with
#' fewer than \code{nmin} individuals are dropped. Population assignment is
#' performed by \code{assignPOP::assign.X()} internally.
#'
#' \strong{Two-object mode} (\code{x.unknown} supplied): \code{x} contains
#' the reference individuals (knowns) and \code{x.unknown} contains individuals
#' whose population of origin is to be determined. Reference populations with
#' fewer than \code{nmin} individuals are dropped. Population assignment is
#' performed by \code{assignPOP::assign.X()} internally.
#'
#' In all two-dataset modes the two objects are reduced to their common set of
#' loci before encoding.
#'
#' Internally, genotype dosage values (0, 1, 2) are one-hot encoded into allele
#' frequency pairs per locus, replicating the \code{structure_onehot()} step
#' performed by \code{read.Structure()}:
#' \itemize{
#'   \item dosage 0 (homozygous reference) \eqn{\rightarrow} \code{[1.0, 0.0]}
#'   \item dosage 1 (heterozygous)         \eqn{\rightarrow} \code{[0.5, 0.5]}
#'   \item dosage 2 (homozygous alternate)  \eqn{\rightarrow} \code{[0.0, 1.0]}
#'   \item \code{NA} (missing)             \eqn{\rightarrow} \code{[0.0, 0.0]}
#' }
#' Each locus contributes two columns to \code{DataMatrix}, named
#' \code{<locusname>_1} and \code{<locusname>_2}.
#'
#' @param x Genlight object. In single-object and split modes this is the only
#' object required. In two-object mode it contains the known reference
#' individuals [required].
#' @param unknown.id Character vector of individual names within \code{x} that
#' are to be treated as unknowns. These are extracted from \code{x}
#' automatically; the remainder of \code{x} becomes the reference set. Cannot
#' be used together with \code{x.unknown} [default NULL].
#' @param x.unknown Genlight object containing individuals whose population of
#' origin is unknown. Population labels in this object are ignored. Cannot be
#' used together with \code{unknown.id} [default NULL].
#' @param nmin Minimum number of individuals required in a reference population
#' for it to be retained. Populations with fewer individuals are dropped with a
#' warning. Applied only when unknowns are present [default 10].
#' @param outpath Directory where \code{assignPOP::assign.X()} writes its output
#' files (including \code{AssignmentResult.txt}). Only used in split and
#' two-object modes [default getwd()].
#' @param model Classification model passed to \code{assignPOP::assign.X()}.
#' One of \code{"svm"}, \code{"naiveBayes"}, \code{"randomForest"},
#' \code{"tree"}, or \code{"knn"} [default "svm"].
#' @param svm.kernel Kernel type for the SVM classifier, passed to
#' \code{assignPOP::assign.X()}. One of \code{"linear"}, \code{"radial"},
#' \code{"polynomial"}, or \code{"sigmoid"} [default "linear"].
#' @param svm.cost Cost parameter for the SVM classifier, passed to
#' \code{assignPOP::assign.X()} [default 1].
#' @param ntree Number of trees for the random forest classifier, passed to
#' \code{assignPOP::assign.X()} [default 50].
#' @param pca.PCs Number of principal components to retain for dimensionality
#' reduction inside \code{assign.X()}, or \code{"kaiser-guttman"} to select
#' automatically [default "kaiser-guttman"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#'
#' @author Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \donttest{
#' require("dartR.data")
#' if (isTRUE(getOption("dartR_fbm"))) possums.gl <- gl.gen2fbm(possums.gl)
#'
#' # Single-object mode: cross-validation on known individuals
#' dat <- gl.run.assignpop(possums.gl, verbose = 3)
#' # assignPOP::assign.kfold(dat, k.fold = 5, dir = tempdir())
#'
#' # Split mode: unknown extracted by name, assignment run automatically
#' # out <- gl.run.assignpop(gl, unknown.id = "AA011731", nmin = 10,
#' #                          dir = tempdir(), verbose = 3)
#' # out$results   # assignment result data frame
#' # out$train     # reference data object
#' # out$unknowns  # unknown data object
#' }
#'
#' @export
#' @return \strong{Single-object mode}: a named 3-element list identical in
#' structure to the output of \code{assignPOP::read.Structure()}:
#' \describe{
#'   \item{\code{DataMatrix}}{A data frame of one-hot encoded allele frequencies
#'     with a final \code{popNames_vector} column.}
#'   \item{\code{SampleID}}{Character vector of individual identifiers.}
#'   \item{\code{LocusName}}{Character vector of locus names.}
#' }
#' \strong{Split mode or two-object mode}: a named list with three elements:
#' \describe{
#'   \item{\code{results}}{A data frame containing the assignment result for
#'     each unknown individual: columns \code{Ind.ID}, \code{pred.pop}, and
#'     one column per reference population giving the posterior probability.
#'     Also written to \code{AssignmentResult.txt} in \code{dir}.}
#'   \item{\code{train}}{The assignPOP data object for the reference
#'     individuals, suitable for subsequent cross-validation with
#'     \code{assign.kfold()} or \code{assign.MC()}.}
#'   \item{\code{unknowns}}{The assignPOP data object for the unknown
#'     individuals.}
#' }

dir <- outpath

gl.run.assignpop <- function(x,
                             unknown.id = NULL,
                             x.unknown  = NULL,
                             nmin       = 10,
                             dir        = getwd(),
                             model      = "svm",
                             svm.kernel = "linear",
                             svm.cost   = 1,
                             ntree      = 50,
                             pca.PCs    = "kaiser-guttman",
                             verbose    = NULL) {

  # PRELIMINARIES -- checking ------------------------------------------------

  verbose  <- gl.check.verbosity(verbose)
  funname  <- match.call()[[1]]
  utils.flag.start(func = funname, build = "v.2023.3", verbose = verbose)
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # FUNCTION-SPECIFIC ERROR CHECKING -----------------------------------------

  if (!is.null(unknown.id) && !is.null(x.unknown)) {
    stop(
      "  Fatal Error: supply either unknown.id or x.unknown, not both.\n"
    )
  }

  if (is.null(pop(x))) {
    if (verbose >= 2) {
      cat(
        "  Population assignments not found in the genlight object.\n",
        "  All individuals assigned to a single population: 'Pop1'.\n"
      )
    }
    pop(x) <- rep("Pop1", nInd(x))
  }

  # SPLIT MODE: extract unknowns from x by individual name -------------------

  if (!is.null(unknown.id)) {

    missing_ids <- setdiff(unknown.id, indNames(x))
    if (length(missing_ids) > 0L) {
      stop(
        "  Fatal Error: the following unknown.id names were not found in x:\n",
        "  ", paste(missing_ids, collapse = ", "), "\n"
      )
    }

    if (verbose >= 2) {
      cat(
        "  Extracting", length(unknown.id),
        "unknown individual(s) from x:", paste(unknown.id, collapse = ", "), "\n"
      )
    }

    x.unknown <- x[indNames(x) %in% unknown.id, ]
    x         <- x[!indNames(x) %in% unknown.id, ]
  }

  # NMIN FILTERING: drop small reference populations -------------------------

  if (!is.null(x.unknown)) {

    pop_counts <- table(pop(x))
    pop_drop   <- names(pop_counts)[pop_counts < nmin]
    pop_keep   <- names(pop_counts)[pop_counts >= nmin]

    if (length(pop_drop) > 0L) {
      if (verbose >= 1) {
        cat(
          "  Dropping", length(pop_drop),
          "reference population(s) with fewer than", nmin, "individuals:\n",
          " ", paste(pop_drop, collapse = ", "), "\n"
        )
      }
      x      <- x[pop(x) %in% pop_keep, ]
      pop(x) <- droplevels(pop(x))
    }

    if (nInd(x) == 0L) {
      stop(
        "  Fatal Error: no reference individuals remain after nmin filtering.\n",
        "  Lower nmin or supply a larger reference dataset.\n"
      )
    }

    if (verbose >= 2) {
      cat(
        "  Reference set after nmin filtering:",
        nInd(x), "individuals,", nPop(x), "populations\n"
      )
    }
  }

  # LOCUS ALIGNMENT ----------------------------------------------------------

  if (!is.null(x.unknown)) {

    utils.check.datatype(x.unknown, accept = "SNP", verbose = 0)

    common_loci <- intersect(locNames(x), locNames(x.unknown))

    if (length(common_loci) == 0L) {
      stop(
        "  Fatal Error: no loci in common between reference and unknown sets.\n",
        "  Ensure both derive from the same SNP set.\n"
      )
    }

    n_dropped_x   <- nLoc(x)        - length(common_loci)
    n_dropped_unk <- nLoc(x.unknown) - length(common_loci)

    if (n_dropped_x > 0L && verbose >= 1) {
      cat(
        "  Warning:", n_dropped_x,
        "reference loci not in unknown set — dropped.\n"
      )
    }
    if (n_dropped_unk > 0L && verbose >= 1) {
      cat(
        "  Warning:", n_dropped_unk,
        "unknown loci not in reference set — dropped.\n"
      )
    }

    x         <- x[, locNames(x) %in% common_loci]
    x         <- x[, match(common_loci, locNames(x))]
    x.unknown <- x.unknown[, locNames(x.unknown) %in% common_loci]
    x.unknown <- x.unknown[, match(common_loci, locNames(x.unknown))]

    if (verbose >= 2) {
      cat("  Loci retained after alignment:", length(common_loci), "\n")
    }
  }

  # INTERNAL HELPER: build an assignPOP data object from a genlight object ---

  .build_assignpop_obj <- function(gl, pop_override = NULL, label = "x") {

    gmat   <- as.matrix(gl)
    n_ind  <- nInd(gl)
    n_loc  <- nLoc(gl)
    ind_id <- indNames(gl)
    loc_id <- locNames(gl)
    pop_id <- if (!is.null(pop_override)) {
      rep(pop_override, n_ind)
    } else {
      as.character(pop(gl))
    }

    if (verbose >= 2) {
      cat("  One-hot encoding allele frequencies for", label, "\n")
    }

    lookup1 <- c(1.0, 0.5, 0.0)
    lookup2 <- c(0.0, 0.5, 1.0)

    freq1         <- matrix(0.0, nrow = n_ind, ncol = n_loc)
    freq2         <- matrix(0.0, nrow = n_ind, ncol = n_loc)
    non_na        <- !is.na(gmat)
    freq1[non_na] <- lookup1[gmat[non_na] + 1L]
    freq2[non_na] <- lookup2[gmat[non_na] + 1L]

    dm_numeric <- matrix(0.0, nrow = n_ind, ncol = 2L * n_loc)
    dm_numeric[, seq(1L, 2L * n_loc, by = 2L)] <- freq1
    dm_numeric[, seq(2L, 2L * n_loc, by = 2L)] <- freq2
    colnames(dm_numeric) <- c(rbind(paste0(loc_id, "_1"),
                                    paste0(loc_id, "_2")))
    rownames(dm_numeric) <- ind_id

    dm                      <- as.data.frame(dm_numeric, stringsAsFactors = FALSE)
    dm[["popNames_vector"]] <- pop_id

    attr(dm, ".n_ind")    <- n_ind
    attr(dm, ".n_loc")    <- n_loc
    attr(dm, ".miss_pct") <- round(100 * sum(is.na(gmat)) / (n_ind * n_loc), 2)

    list(DataMatrix = dm, SampleID = ind_id, LocusName = loc_id)
  }

  # BUILD DATA OBJECTS -------------------------------------------------------

  if (verbose >= 2) {
    cat("  Building assignPOP data object for known individuals\n")
  }
  train_obj <- .build_assignpop_obj(x, label = "knowns")

  if (!is.null(x.unknown)) {
    if (verbose >= 2) {
      cat("  Building assignPOP data object for unknown individuals\n")
    }
    unk_obj <- .build_assignpop_obj(x.unknown, pop_override = "Unknown",
                                    label = "unknowns")
  }

  # RESULTS SUMMARY ----------------------------------------------------------

  if (verbose >= 3) {
    cat("\n  --- Known individuals (training) ---\n")
    cat("  Individuals  :", attr(train_obj$DataMatrix, ".n_ind"), "\n")
    cat("  Loci         :", attr(train_obj$DataMatrix, ".n_loc"), "\n")
    cat("  Populations  :", nPop(x), "\n")
    cat("  Missing data :", attr(train_obj$DataMatrix, ".miss_pct"), "%\n")

    if (!is.null(x.unknown)) {
      cat("\n  --- Unknown individuals ---\n")
      cat("  Individuals  :", attr(unk_obj$DataMatrix, ".n_ind"), "\n")
      cat("  Loci         :", attr(unk_obj$DataMatrix, ".n_loc"), "\n")
      cat("  Missing data :", attr(unk_obj$DataMatrix, ".miss_pct"), "%\n")
    }
  }

  # ASSEMBLE OUTPUT ----------------------------------------------------------

  if (is.null(x.unknown)) {

    # Single-object mode: return data object for cross-validation
    out <- train_obj

    if (verbose >= 1) {
      cat("\n  Next steps -- pass the returned object to assignPOP functions:\n")
      cat("\n    Optionally reduce loci by minor allele count:\n")
      cat("      dat <- assignPOP::reduce.allele(dat, min.mac = 5)\n")
      cat("\n    k-fold cross-validation assignment:\n")
      cat("      assignPOP::assign.kfold(dat, k.fold = 5, dir = tempdir())\n")
      cat("\n    Monte Carlo cross-validation assignment:\n")
      cat("      assignPOP::assign.MC(dat, train.inds = 0.7, dir = tempdir())\n")
      cat("\n    Compile and plot accuracy (after assign.kfold or assign.MC):\n")
      cat("      acc <- assignPOP::accuracy.kfold(dat, dir = tempdir())\n")
      cat("      assignPOP::accuracy.plot(acc)\n")
      cat("\n    Plot individual membership probabilities:\n")
      cat("      assignPOP::membership.plot(dir = tempdir())\n")
    }

  } else {

    # Two-dataset mode: run assign.X() internally
    if (!requireNamespace("assignPOP", quietly = TRUE)) {
      stop(
        "  Package 'assignPOP' is required but not installed.\n",
        "  Install it with: install.packages('assignPOP')\n"
      )
    }

    if (verbose >= 2) {
      cat("\n  Running assignPOP::assign.X()...\n")
    }

    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    # assign.X() requires a trailing slash
    dir_slash <- paste0(gsub("/*$", "", dir), "/")

    assignPOP::assign.X(
      x1         = train_obj,
      x2         = unk_obj,
      dir        = dir_slash,
      model      = model,
      svm.kernel = svm.kernel,
      svm.cost   = svm.cost,
      ntree      = ntree,
      pca.PCs    = pca.PCs,
      skipQ      = TRUE      # suppress interactive prompts
    )

    # Read results back and include in return value
    result_file <- file.path(dir, "AssignmentResult.txt")
    results     <- read.table(result_file, header = TRUE, sep = " ",
                              stringsAsFactors = FALSE, check.names = FALSE)

    if (verbose >= 3) {
      cat("\n  --- Assignment results ---\n")
      pops  <- names(results)[3:ncol(results)]
      for (i in seq_len(nrow(results))) {
        ind_id  <- results$Ind.ID[i]
        pred    <- results$pred.pop[i]
        posts   <- as.numeric(results[i, 3:ncol(results)])
        names(posts) <- pops
        top5    <- sort(posts, decreasing = TRUE)[1:min(5L, length(posts))]
        cat(sprintf("\n  Individual: %s  ->  %s\n", ind_id, pred))
        cat("  Top 5 posterior probabilities:\n")
        for (j in seq_along(top5)) {
          cat(sprintf("    %-30s %.4f\n", names(top5)[j], top5[j]))
        }
      }
    }

    out <- list(results = results, train = train_obj, unknowns = unk_obj)

    if (verbose >= 1) {
      cat("\n  Results written to:", file.path(dir, "AssignmentResult.txt"), "\n")
      cat("\n  To run cross-validation on the reference individuals:\n")
      cat("      assignPOP::assign.kfold(out$train, k.fold = 5, dir = tempdir())\n")
    }
  }

  # FLAG SCRIPT END ----------------------------------------------------------

  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN -------------------------------------------------------------------

  invisible(out)
}
