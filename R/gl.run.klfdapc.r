#' @name gl.run.klfdapc
#' @title Runs KLFDAPC population assignment directly from a genlight object
#' @family linker
#'
#' @description
#' Provides a seamless interface between dartRverse genlight objects and the
#' \pkg{KLFDAPC} package (Qin et al., 2022, Briefings in Bioinformatics,
#' doi:10.1093/bib/bbac202), which performs Kernel Local Fisher Discriminant
#' Analysis of Principal Components for population assignment.
#'
#' The function handles all data conversion and file management internally.
#' Temporary GDS (Genomic Data Structure) files are created as needed, used,
#' and deleted automatically on exit — the user never needs to manage files
#' on disk.
#'
#' \strong{Single-object mode} (\code{unknown.id = NULL},
#' \code{x.unknown = NULL}): \code{x} contains individuals with known
#' population assignments. The trained model object is returned and can be
#' passed to \code{KLFDAPC::predict.klfdapc()} for self-assignment or for
#' predicting new individuals outside this function.
#'
#' \strong{Split mode} (\code{unknown.id} supplied): the simplest mode — a
#' single genlight object containing both reference and unknown individuals.
#' The unknown individuals are identified by name, extracted automatically,
#' and the remainder becomes the reference set.
#'
#' \strong{Two-object mode} (\code{x.unknown} supplied): \code{x} contains
#' the reference individuals (knowns) used to train the model, and
#' \code{x.unknown} contains individuals whose population of origin is to be
#' determined. The two objects are reduced to their common set of loci. The
#' unknowns are projected into the training PC space using SNP loadings derived
#' from the training PCA (\code{SNPRelate::snpgdsPCASNPLoading()} and
#' \code{SNPRelate::snpgdsPCASampLoading()}), and \code{predict.klfdapc()} is
#' called internally. Both the trained model and the prediction results are
#' returned.
#'
#' \strong{KLFDAPC minimum data requirements:} the KLFDA algorithm requires
#' each reference population to have at least \code{max(knn, r + 1)}
#' individuals. With default parameters (\code{knn = 6}, \code{r = 3}) this
#' minimum is 6. Populations below this threshold are dropped automatically
#' when \code{nmin} is set, or the function stops with an informative message.
#' At least 2 reference populations must remain after any filtering.
#'
#' @param x Genlight object containing individuals with \strong{known}
#' population assignments, used as the reference (training) dataset [required].
#' @param unknown.id Character vector of individual names within \code{x} that
#' are to be treated as unknowns. These are extracted from \code{x}
#' automatically; the remainder of \code{x} becomes the reference set. Cannot
#' be used together with \code{x.unknown} [default NULL].
#' @param x.unknown Genlight object containing individuals whose population of
#' origin is unknown and is to be determined. Population labels in this object
#' are ignored. Cannot be used together with \code{unknown.id}. If both
#' \code{unknown.id} and \code{x.unknown} are \code{NULL}, single-object mode
#' is used and only the trained model is returned [default NULL].
#' @param nmin Minimum number of individuals required in a reference population
#' for it to be retained. Populations with fewer individuals are dropped with a
#' warning. Must be at least \code{max(knn, r + 1)} to satisfy KLFDAPC's
#' mathematical requirements. Applied in split and two-object modes
#' [default NULL, which sets the minimum to \code{max(knn, r + 1)}].
#' @param n.pc Number of principal components to retain from the SNP PCA and
#' pass to the kernel discriminant step. Must be >= \code{r}. A value between
#' 10 and 30 is typical [default 10].
#' @param r Number of discriminant dimensions to retain in the KLFDA
#' projection. Cannot exceed the number of populations minus one
#' [default 3].
#' @param kernel A \pkg{kernlab} kernel function object used to compute
#' pairwise similarities between individuals in PC space. The default
#' polynomial kernel of degree 1 (linear) matches the \code{KLFDAPC()}
#' default. Alternatives include \code{kernlab::rbfdot()} (Gaussian) and
#' \code{kernlab::polydot(degree = 2)} (quadratic)
#' [default \code{kernlab::polydot(degree = 1, scale = 1, offset = 1)}].
#' @param maf Minor allele frequency filter applied before PCA; SNPs with
#' MAF below this threshold are excluded. Set to \code{NaN} to skip
#' [default NaN].
#' @param missing.rate Locus missing-data rate filter applied before PCA;
#' SNPs with a missing rate above this threshold are excluded. Set to
#' \code{NaN} to skip [default NaN].
#' @param knn Number of nearest neighbours used to construct the affinity
#' matrix in the KLFDA step. Each reference population must contain at least
#' this many individuals [default 6].
#' @param reg Regularisation parameter added to the within-class scatter
#' matrix diagonal to ensure invertibility [default 0.001].
#' @param metric Type of transformation matrix returned by KLFDA. One of
#' \code{"weighted"}, \code{"orthonormalized"}, or \code{"plain"}
#' [default "weighted"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#'
#' @author Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \donttest{
#' require(dartR.data)
#' if (isTRUE(getOption("dartR_fbm"))) possums.gl <- gl.gen2fbm(possums.gl)
#' # Single-object mode: train model and self-assign
#' model <- gl.run.klfdapc(possums.gl, n.pc = 10, r = 3, verbose = 3)
#'
#' # Split mode: unknown extracted by name from a combined genlight object
#' # out <- gl.run.klfdapc(gl, unknown.id = "AA011731",
#' #                        n.pc = 10, r = 3, verbose = 3)
#' # out$predictions$posteriors.class1   # Bayesian assignment
#'
#' # Two-object mode: separate reference and unknown genlight objects
#' # out <- gl.run.klfdapc(gl.ref, x.unknown = gl.unknown,
#' #                        n.pc = 10, r = 3, verbose = 3)
#' }
#'
#' @export
#' @return \strong{Single-object mode} (\code{unknown.id = NULL} and
#' \code{x.unknown = NULL}): an object of class \code{"KLFDAPC"}, identical in
#' structure to the output of \code{KLFDAPC::KLFDAPC()}. It is a named list
#' with two elements:
#' \describe{
#'   \item{\code{KLFDAPC}}{A \code{"klfda"} object containing the trained
#'     discriminant model, including transformation matrix (\code{$T}),
#'     projected training data (\code{$Z}), class means and covariances,
#'     Naive Bayes classifier, and posterior probabilities.}
#'   \item{\code{PCA}}{The SNP PCA result object returned by
#'     \code{SNPRelate::snpgdsPCA()}, containing eigenvalues, eigenvectors,
#'     and per-individual PC scores.}
#' }
#' \strong{Split mode or two-object mode}: a named list with two elements:
#' \describe{
#'   \item{\code{model}}{The trained \code{"KLFDAPC"} object as above.}
#'   \item{\code{predictions}}{The output of
#'     \code{KLFDAPC::predict.klfdapc()} for the unknown individuals,
#'     containing \code{$posteriors.class1} (Bayesian assignment),
#'     \code{$posteriors1} (Bayesian posteriors), \code{$class} (Naive Bayes
#'     assignment), and \code{$posterior} (LDA-based posteriors).}
#' }

gl.run.klfdapc <- function(x,
                            unknown.id    = NULL,
                            x.unknown     = NULL,
                            nmin          = NULL,
                            n.pc          = 10,
                            r             = 3,
                            kernel        = kernlab::polydot(degree = 1,
                                                             scale  = 1,
                                                             offset = 1),
                            maf           = NaN,
                            missing.rate  = NaN,
                            knn           = 6,
                            reg           = 0.001,
                            metric        = c("weighted",
                                              "orthonormalized",
                                              "plain"),
                            verbose       = NULL) {

  # PRELIMINARIES -- checking ------------------------------------------------

  verbose  <- gl.check.verbosity(verbose)
  funname  <- match.call()[[1]]
  utils.flag.start(func = funname, build = "v.2023.3", verbose = verbose)
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # FUNCTION-SPECIFIC ERROR CHECKING -----------------------------------------

  if (!requireNamespace("SNPRelate", quietly = TRUE)) {
    stop(
      "  Package 'SNPRelate' is required but not installed.\n",
      "  Install it with: BiocManager::install('SNPRelate')\n"
    )
  }
  if (!requireNamespace("KLFDAPC", quietly = TRUE)) {
    stop(
      "  Package 'KLFDAPC' is required but not installed.\n",
      "  Install it with: remotes::install_github('xinghuq/KLFDAPC')\n"
    )
  }
  if (!requireNamespace("kernlab", quietly = TRUE)) {
    stop(
      "  Package 'kernlab' is required but not installed.\n",
      "  Install it with: install.packages('kernlab')\n"
    )
  }

  if (!is.null(unknown.id) && !is.null(x.unknown)) {
    stop(
      "  Fatal Error: supply either unknown.id or x.unknown, not both.\n"
    )
  }

  if (is.null(pop(x))) {
    if (verbose >= 2) {
      cat(
        "  Warning: no population assignments found in x.\n",
        "  Assigning all individuals to a single population: 'Pop1'.\n"
      )
    }
    pop(x) <- rep("Pop1", nInd(x))
  }

  # KLFDAPC HARD LIMITS ------------------------------------------------------
  # The KLFDA algorithm has two structural constraints that, if violated, cause
  # matrix operations to fail with cryptic errors:
  #
  # 1. knn constraint: lfda::getAffinityMatrix() constructs a k-NN affinity
  #    matrix within each class. Each class must have >= knn members or the
  #    k-NN graph is undefined.
  #
  # 2. Covariance constraint: cov(redZ[[i]]) computes an r x r covariance
  #    matrix from the r-dimensional KLFDA projections of class i. A covariance
  #    matrix from n observations has rank min(n-1, r), so solve() requires
  #    n >= r + 1 to avoid a singular matrix.
  #
  # Combined minimum per population: max(knn, r + 1)
  # Minimum number of populations: 2 (r >= 1 requires at least 2 classes)

  nmin_hard <- max(knn, r + 1L)

  if (is.null(nmin)) {
    nmin <- nmin_hard
  } else if (nmin < nmin_hard) {
    if (verbose >= 1) {
      cat(
        "  Warning: nmin =", nmin, "is below the KLFDAPC minimum of",
        nmin_hard, "(max(knn, r + 1) = max(", knn, ",", r + 1L, ")).\n",
        "  Resetting nmin to", nmin_hard, ".\n"
      )
    }
    nmin <- nmin_hard
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
          "reference population(s) with fewer than", nmin, "individuals\n",
          "  (KLFDAPC requires >= max(knn, r+1) =", nmin_hard,
          "per population):\n",
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
  }

  # Check minimum population count
  n_pop <- nPop(x)
  if (n_pop < 2L) {
    stop(
      "  Fatal Error: KLFDAPC requires at least 2 reference populations;\n",
      "  only ", n_pop, " remain after filtering.\n",
      "  Lower nmin or supply more populations.\n"
    )
  }

  # Check all remaining populations meet the per-population minimum
  pop_counts_final <- table(pop(x))
  pops_too_small   <- names(pop_counts_final)[pop_counts_final < nmin_hard]
  if (length(pops_too_small) > 0L) {
    stop(
      "  Fatal Error: the following reference populations have fewer than\n",
      "  max(knn, r+1) = ", nmin_hard, " individuals, which will cause\n",
      "  KLFDAPC matrix operations to fail:\n",
      "  ", paste(pops_too_small, collapse = ", "), "\n",
      "  Increase nmin or reduce knn / r.\n"
    )
  }

  if (verbose >= 2 && !is.null(x.unknown)) {
    cat(
      "  Reference set after filtering:",
      nInd(x), "individuals,", n_pop, "populations\n"
    )
  }

  # Adjust r and n.pc to stay within bounds ----------------------------------

  if (r > n_pop - 1L) {
    r <- n_pop - 1L
    if (verbose >= 2) {
      cat("  Warning: r cannot exceed nPop - 1. Resetting r to", r, "\n")
    }
  }

  if (n.pc < r) {
    n.pc <- r
    if (verbose >= 2) {
      cat(
        "  Warning: n.pc must be >= r. Resetting n.pc to", n.pc, "\n"
      )
    }
  }

  if (n.pc >= nInd(x)) {
    n.pc <- nInd(x) - 1L
    if (verbose >= 2) {
      cat(
        "  Warning: n.pc must be < nInd. Resetting n.pc to", n.pc, "\n"
      )
    }
  }

  metric <- match.arg(metric)

  # LOCUS ALIGNMENT ----------------------------------------------------------

  if (!is.null(x.unknown)) {

    utils.check.datatype(x.unknown, accept = "SNP", verbose = 0)

    common_loci <- intersect(locNames(x), locNames(x.unknown))

    if (length(common_loci) == 0L) {
      stop(
        "  Fatal Error: no loci in common between x and x.unknown.\n",
        "  Ensure both genlight objects were filtered from the same SNP set.\n"
      )
    }

    n_dropped_x   <- nLoc(x)         - length(common_loci)
    n_dropped_unk <- nLoc(x.unknown)  - length(common_loci)

    if (n_dropped_x > 0L && verbose >= 1) {
      cat(
        "  Warning:", n_dropped_x,
        "loci in x not found in x.unknown and will be dropped.\n"
      )
    }
    if (n_dropped_unk > 0L && verbose >= 1) {
      cat(
        "  Warning:", n_dropped_unk,
        "loci in x.unknown not found in x and will be dropped.\n"
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

  # INTERNAL HELPER: recode a genlight object to a GDS file -----------------

  .gl_to_gds <- function(gl, filepath) {

    gmat     <- as.matrix(gl)
    n_loc_gl <- nLoc(gl)

    # Recode dartR dosage (0=ref/ref, 1=het, 2=alt/alt, NA=missing)
    # to SNPRelate convention (2=AA, 1=AB, 0=BB, 3=missing)
    gmat_gds <- 2L - gmat
    storage.mode(gmat_gds) <- "integer"
    gmat_gds[is.na(gmat)]  <- 3L

    # Chromosome metadata
    has_chr <- FALSE
    if (!is.null(gl@other$loc.metrics)) {
      lm <- gl@other$loc.metrics
      if ("Chrom" %in% names(lm)) {
        chr     <- suppressWarnings(as.integer(lm$Chrom))
        chr[is.na(chr)] <- 0L
        has_chr <- TRUE
      } else if ("ChromName" %in% names(lm)) {
        chr     <- suppressWarnings(as.integer(lm$ChromName))
        chr[is.na(chr)] <- 0L
        has_chr <- TRUE
      }
    }
    if (!has_chr) {
      chr <- rep(1L, n_loc_gl)
    }

    # Position metadata
    if (!is.null(gl@other$loc.metrics) &&
        "ChromPos" %in% names(gl@other$loc.metrics)) {
      pos <- suppressWarnings(as.integer(gl@other$loc.metrics$ChromPos))
      pos[is.na(pos)] <- seq_len(n_loc_gl)[is.na(pos)]
    } else {
      pos <- seq_len(n_loc_gl)
    }

    SNPRelate::snpgdsCreateGeno(
      gds.fn         = filepath,
      genmat         = t(gmat_gds),
      sample.id      = indNames(gl),
      snp.id         = locNames(gl),
      snp.chromosome = chr,
      snp.position   = pos,
      snpfirstdim    = TRUE
    )

    invisible(gmat)
  }

  # DO THE JOB: TRAIN MODEL ON KNOWNS ----------------------------------------

  gmat  <- as.matrix(x)
  n_ind <- nInd(x)
  n_loc <- nLoc(x)

  if (verbose >= 2) {
    cat("  Writing temporary GDS file for known individuals\n")
  }

  tmp_gds_train <- tempfile(fileext = ".gds")
  on.exit({
    if (file.exists(tmp_gds_train)) file.remove(tmp_gds_train)
  }, add = TRUE)

  .gl_to_gds(x, tmp_gds_train)

  if (verbose >= 2) {
    cat("  Running KLFDAPC (PCA + kernel discriminant analysis) on knowns\n")
  }

  pop_labels <- as.character(pop(x))

  # KLFDAPC() opens the GDS, runs snpgdsPCA, builds the kernel, fits the
  # discriminant model, then closes the GDS.
  # autosome.only = FALSE because DArT SNPs are typically not annotated
  # by chromosome and have all been assigned to chr 1 above.
  model <- KLFDAPC::KLFDAPC(
    infile        = tmp_gds_train,
    y             = pop_labels,
    n.pc          = n.pc,
    autosome.only = FALSE,
    maf           = maf,
    missing.rate  = missing.rate,
    kernel        = kernel,
    r             = r,
    knn           = knn,
    reg           = reg,
    metric        = metric,
    verbose       = (verbose >= 5)
  )

  # DO THE JOB: PROJECT UNKNOWNS & PREDICT -----------------------------------

  if (!is.null(x.unknown)) {

    n_ind_unk <- nInd(x.unknown)
    gmat_unk  <- as.matrix(x.unknown)

    if (verbose >= 2) {
      cat("  Computing SNP loadings from training PCA\n")
    }

    # Reopen the training GDS to derive SNP loadings from the PCA result.
    # snpgdsPCASNPLoading() computes the per-SNP contribution to each PC,
    # allowing projection of new samples into the same PC space.
    gds_train    <- SNPRelate::snpgdsOpen(tmp_gds_train)
    snp_loadings <- SNPRelate::snpgdsPCASNPLoading(
      pcaobj  = model$PCA,
      gdsobj  = gds_train,
      verbose = (verbose >= 5)
    )
    SNPRelate::snpgdsClose(gds_train)

    if (verbose >= 2) {
      cat("  Writing temporary GDS file for unknown individuals\n")
    }

    tmp_gds_unk <- tempfile(fileext = ".gds")
    on.exit({
      if (file.exists(tmp_gds_unk)) file.remove(tmp_gds_unk)
    }, add = TRUE)

    .gl_to_gds(x.unknown, tmp_gds_unk)

    if (verbose >= 2) {
      cat("  Projecting unknown individuals into training PC space\n")
    }

    # snpgdsPCASampLoading() projects unknowns onto PC axes defined by the
    # training SNP loadings. It returns SCALED scores (unit-variance), whereas
    # the training eigenvectors stored in model$KLFDAPC$obj.trainData are RAW
    # (unit-length). We must invert the scaling before passing to
    # predict.klfdapc() so that kernel distances are computed in the same space.
    #   scaled   = raw * sqrt(ss / eigenval)
    #   raw      = scaled * sqrt(eigenval / ss),  where ss = (n-1) / TraceXTX
    gds_unk <- SNPRelate::snpgdsOpen(tmp_gds_unk)
    proj    <- SNPRelate::snpgdsPCASampLoading(
      loadobj = snp_loadings,
      gdsobj  = gds_unk,
      verbose = (verbose >= 5)
    )
    SNPRelate::snpgdsClose(gds_unk)

    n_train   <- length(model$PCA$sample.id)
    ss        <- (n_train - 1L) / model$PCA$TraceXTX
    inv_scale <- sqrt(snp_loadings$eigenval[seq_len(n.pc)] / ss)

    unknown_pcs <- proj$eigenvect[, seq_len(n.pc), drop = FALSE]
    unknown_pcs <- sweep(unknown_pcs, 2L, inv_scale, `*`)

    if (verbose >= 2) {
      cat("  Predicting population of origin for unknown individuals\n")
    }

    # KLFDAPC::predict.klfdapc() calls klaR::Mabayes(), which was moved to
    # the DA package in newer versions of klaR. Patch the function body at
    # runtime so that it uses DA::Mabayes() instead, without modifying any
    # installed package files.
    predict_fn <- KLFDAPC::predict.klfdapc
    if (requireNamespace("DA", quietly = TRUE) &&
        !exists("Mabayes", envir = asNamespace("klaR"), inherits = FALSE)) {
      body_lines   <- deparse(body(predict_fn))
      body_lines   <- gsub("klaR::Mabayes", "DA::Mabayes",
                           body_lines, fixed = TRUE)
      body(predict_fn) <- parse(text = paste(body_lines, collapse = "\n"))[[1L]]
      if (verbose >= 2) {
        cat("  Applied klaR/DA compatibility patch for predict.klfdapc\n")
      }
    }

    predictions <- predict_fn(model, newdata = unknown_pcs)
  }

  # RESULTS SUMMARY ----------------------------------------------------------

  if (verbose >= 3) {
    cat("\n  --- Known individuals (training) ---\n")
    cat("  Individuals  :", n_ind, "\n")
    cat("  Loci         :", n_loc, "\n")
    cat("  Populations  :", n_pop, "\n")
    cat("  PCs used     :", n.pc, "\n")
    cat("  Discriminant dimensions (r):", r, "\n")
    cat(
      "  Missing data :",
      round(100 * sum(is.na(gmat)) / (n_ind * n_loc), 2), "%\n"
    )
    cat("  Self-assignment confusion matrix:\n")
    tbl <- table(
      Predicted = model$KLFDAPC$posteriors.classZ,
      Actual    = pop_labels
    )
    print(tbl)

    if (!is.null(x.unknown)) {
      cat("\n  --- Unknown individuals ---\n")
      cat("  Individuals  :", n_ind_unk, "\n")
      cat(
        "  Missing data :",
        round(100 * sum(is.na(gmat_unk)) / (n_ind_unk * n_loc), 2), "%\n"
      )
      cat("  Bayesian assignment:\n")
      print(table(predictions$posteriors.class1))
    }
  }

  # ASSEMBLE OUTPUT & NEXT STEPS ---------------------------------------------

  if (is.null(x.unknown)) {

    out <- model

    if (verbose >= 1) {
      cat("\n  Next steps -- using the returned model object:\n")
      cat("\n    Self-assignment (training data):\n")
      cat("      pred <- KLFDAPC::predict.klfdapc(\n")
      cat("                model,\n")
      cat("                newdata = model$KLFDAPC$obj.trainData)\n")
      cat("\n    Access posterior probabilities:\n")
      cat("      pred$posteriors1   # Bayesian posteriors\n")
      cat("      pred$posterior     # LDA-based posteriors\n")
      cat("\n    Access training PC scores:\n")
      cat("      model$PCA$eigenvect[, 1:", n.pc, "]\n", sep = "")
    }

  } else {

    out <- list(model = model, predictions = predictions)

    if (verbose >= 1) {
      cat("\n  Next steps -- using the returned list:\n")
      cat("\n    Access Bayesian assigned populations:\n")
      cat("      out$predictions$posteriors.class1\n")
      cat("\n    Access Bayesian posterior probabilities:\n")
      cat("      out$predictions$posteriors1\n")
      cat("\n    Access LDA-based posterior probabilities:\n")
      cat("      out$predictions$posterior\n")
      cat("\n    Self-assignment accuracy on knowns:\n")
      cat("      KLFDAPC::predict.klfdapc(\n")
      cat("        out$model,\n")
      cat("        newdata = out$model$KLFDAPC$obj.trainData)\n")
    }
  }

  # FLAG SCRIPT END ----------------------------------------------------------

  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN -------------------------------------------------------------------

  invisible(out)
}
