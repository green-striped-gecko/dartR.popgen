#' Check a SNP panel for specified conservation monitoring parameters
#'
#' @description
#' Evaluates how well a SNP panel recreates specified parameters of
#' conservation concern compared to the full original dataset, or assesses
#' the panel's power for individual identification, parentage assignment,
#' population assignment, or F1 hybrid detection. Multiple parameters can be
#' evaluated simultaneously. All plots are combined on a single page with the
#' summary bar chart as the first panel in the same grid.
#'
#' @param x A \code{dartR} or \code{genlight} object containing the SNP panel
#'   (subset of loci).
#' @param xorig A \code{dartR} or \code{genlight} object containing the full
#'   original dataset for comparison. Default \code{NULL}. Required for all
#'   correlation/comparison parameters (\code{"Fst"}, \code{"He"},
#'   \code{"Ho"}, \code{"Fis"}, \code{"Nall"}, \code{"Ne"},
#'   \code{"Ho_ind"}, \code{"Fis_ind"}, \code{"relatedness"}). May be left
#'   \code{NULL} when only \code{"drift_resistance"}, \code{"id"},
#'   \code{"parentage"}, \code{"assignment"}, or \code{"hybridisation"} are
#'   requested, as these are evaluated on \code{x} alone.
#' @param parameter Character string or vector specifying the parameter(s) to
#'   evaluate. Options:
#'   \itemize{
#'     \item \code{"Fst"} — pairwise FST' (between-population).
#'     \item \code{"He"} — expected heterozygosity per population.
#'     \item \code{"Ho"} — observed heterozygosity per population.
#'     \item \code{"Fis"} — inbreeding coefficient per population.
#'     \item \code{"Nall"} or \code{"Ar"} — rarefaction-corrected allelic
#'       richness per population.
#'     \item \code{"Ne"} — effective population size per population.
#'     \item \code{"Ho_ind"} — individual observed heterozygosity.
#'       Correlation r² of per-individual Ho between the panel and full
#'       dataset, ignoring population membership.
#'     \item \code{"Fis_ind"} — individual within-population inbreeding
#'       coefficient. For each individual, \eqn{F_i = 1-H_{O,i}/H_{E,i}} is
#'       calculated from the full dataset and from the panel. Performance is
#'       the selected \code{corr.method} squared correlation between the two
#'       sets of individual estimates, consistent with the other
#'       correlation-based parameters.
#'     \item \code{"drift_resistance"} — absolute panel score from 0 to 1
#'       based on within-population minor allele frequencies. For locus
#'       \eqn{l} in population \eqn{k}, the score is
#'       \eqn{[2\min(p_{lk},1-p_{lk})]^2}. Scores are averaged across loci
#'       within each population and then equally across populations. A score
#'       of 1 means every panel locus has allele frequency 0.5 in every
#'       population; 0 means all panel loci are fixed.
#'     \item \code{"relatedness"} — pairwise genomic relatedness (VanRaden
#'       GRM), optionally rescaled to kinship via \code{metric}. The
#'       allele-frequency reference is controlled by \code{ref}: pooled
#'       across the dataset (\code{"global"}, population membership
#'       ignored) or computed within each population
#'       (\code{"by.pop"}, within-population pairs only).
#'     \item \code{"id"} — individual identification power.
#'     \item \code{"parentage"} — parentage assignment power.
#'     \item \code{"assignment"} — population assignment accuracy (LOO).
#'     \item \code{"hybridisation"} — F1 hybrid detection across all
#'       pairwise population crosses (requires at least 2 populations in
#'       \code{x}).
#'     \item \code{"all"} — all parameters except \code{"hybridisation"}.
#'     \item \code{"all-Ne"} — all except \code{"Ne"} and
#'       \code{"hybridisation"}.
#'   }
#' @param corr.method Character. Correlation method used as the performance
#'   indicator for correlation-based parameters. Options are
#'   \code{"spearman"} (default, preserving previous behaviour) and
#'   \code{"pearson"}. Applies to \code{"Fst"}, \code{"He"},
#'   \code{"Ho"}, \code{"Fis"}, \code{"Nall"}, \code{"Ne"},
#'   \code{"Ho_ind"}, \code{"Fis_ind"}, and \code{"relatedness"}.
#'   Accuracy/proportion
#'   metrics such as \code{"id"}, \code{"parentage"},
#'   \code{"assignment"}, and \code{"hybridisation"} are unchanged.
#' @param ref Character. Allele-frequency reference for
#'   \code{"relatedness"}. \code{"global"} (default) pools allele
#'   frequencies across the whole dataset and uses all pairs; between
#'   differentiated populations this produces strongly negative values
#'   that mainly reflect structure. \code{"by.pop"} computes allele
#'   frequencies within each population, references each population to
#'   its own base, and retains only within-population pairs — appropriate
#'   for within-population kinship. Populations with fewer than two
#'   individuals are skipped with a warning. Ignored by all other
#'   parameters.
#' @param metric Character. Reporting scale for \code{"relatedness"}.
#'   \code{"grm"} (default) returns the VanRaden genomic relationship
#'   matrix; \code{"kinship"} returns kinship = GRM / 2 (so an outbred
#'   self is ~0.5, parent--offspring/full-sib ~0.25). Because this is a
#'   constant rescaling applied to both panel and full dataset, it does
#'   not change the correlation-based performance value; it affects only
#'   the reported values, axis labels, and plot. Ignored by all other
#'   parameters.
#' @param inverse_dr Logical. If \code{TRUE} and \code{"drift_resistance"}
#'   is among the requested parameters, the performance is reported as
#'   \code{1 - drift_resistance} (drift sensitivity). A high value then
#'   indicates a panel whose loci are near fixation and therefore
#'   maximally responsive to allele-frequency change — desirable for
#'   genetic-erosion monitoring but not for long-lived stable panels.
#'   Default \code{FALSE}.
#' @param neest.path Character string. Path to NEstimator executable, required
#'   only for \code{"Ne"}.
#' @param error.rate Numeric. Per-allele genotyping error rate. Default
#'   \code{0.01}.
#' @param threshold Numeric. Confusion threshold for \code{"id"} and
#'   \code{"parentage"}. Default \code{0.001}.
#' @param n_sim_parents Integer. Simulated trios for \code{"parentage"}.
#'   Default \code{100}.
#' @param n_sim_hyb Integer. Simulated F1 individuals per pairwise
#'   population cross for \code{"hybridisation"}. Total F1 simulations
#'   are \code{n_sim_hyb * choose(nPop(x), 2)}. Default \code{100}.
#' @param n_cores Integer or \code{NULL}. Cores for parallel parentage
#'   simulation. \code{NULL} = auto. On Windows always sequential.
#' @param target Numeric (0--1) or \code{NULL}. Reference line in the summary
#'   bar chart. Default \code{0.9}.
#' @param plot.out Logical. Print combined plot. Default \code{TRUE}.
#' @param plot.file Character string. File name for saving the plot.
#' @param plot.dir Character string. Directory for saving the plot.
#' @param verbose Integer. \code{0} = silent, \code{1} = parameter names and
#'   performance scores only, \code{2} = full progress (default).
#'
#' @details
#' \strong{Performance metric (0--1):}
#' \itemize{
#'   \item Population and between-population parameters (\code{"Fst"},
#'     \code{"He"}, \code{"Ho"}, \code{"Fis"}, \code{"Nall"}, \code{"Ne"})
#'     — chosen \code{corr.method} r² between panel and full dataset
#'     values. Plots show both Spearman r² and Pearson R².
#'   \item \code{"Ho_ind"} — chosen \code{corr.method} r² of
#'     per-individual observed heterozygosity. Population membership is
#'     ignored.
#'   \item \code{"Fis_ind"} — agreement between individual inbreeding
#'     estimates from the panel and the full dataset. Within each population,
#'     expected heterozygosity is calculated as \eqn{H_{e,l}=2p_l(1-p_l)}.
#'     For each individual and its called loci,
#'     \eqn{F_i=1-H_{O,i}/H_{E,i}}. The returned data contain the actual
#'     \code{fis_orig} and \code{fis_panel} values for every individual.
#'     Performance is the chosen \code{corr.method} squared correlation
#'     between \code{fis_orig} and \code{fis_panel}, exactly as for the other
#'     correlation-based parameters. The summary also reports RMSE, MAE,
#'     bias, Pearson \eqn{R^2}, and Spearman \eqn{r^2}.
#'   \item \code{"drift_resistance"} — absolute panel score
#'     \eqn{\mathrm{mean}_{k}[\mathrm{mean}_{l}
#'     \{(2\min(p_{lk},1-p_{lk}))^2\}]}. This is not a correlation with
#'     \code{xorig}; it directly rewards loci with intermediate allele
#'     frequencies and strongly penalises rare or nearly fixed alleles.
#'   \item \code{"relatedness"} — chosen \code{corr.method} r² of all
#'     retained pairwise relatedness values (upper triangle), on the
#'     scale set by \code{metric} (GRM or kinship = GRM/2) and the
#'     reference set by \code{ref} (\code{"global"}: all pairs, pooled
#'     frequencies; \code{"by.pop"}: within-population pairs, local
#'     frequencies) (VanRaden 2008).
#'   \item \code{"id"} — fraction of individuals with MaxPmatch below
#'     \code{threshold}.
#'   \item \code{"parentage"} — fraction of simulated trios correctly
#'     assigned.
#'   \item \code{"assignment"} — fraction correctly assigned by LOO.
#'   \item \code{"hybridisation"} — correct F1-pair assignment rate,
#'     calculated as the proportion of simulated F1 hybrids assigned to the
#'     correct pairwise F1 class across all pairwise population crosses.
#' }
#'
#' \strong{Heterozygosity caching:} When any of \code{"He"}, \code{"Ho"},
#' \code{"Fis"} are requested, \code{gl.report.heterozygosity} is called
#' once per dataset and reused for all three.
#'
#' \strong{Relatedness reference and scale:} The VanRaden GRM is centred on
#' a reference allele-frequency vector, so a value of zero corresponds to
#' the average pair in that reference and roughly half of the off-diagonal
#' values are negative by construction. With \code{ref = "global"} the
#' reference is pooled across all populations; in structured data this makes
#' between-population pairs strongly negative, so the metric largely tracks
#' recovery of population structure. With \code{ref = "by.pop"} each
#' population is referenced to its own allele frequencies and only
#' within-population pairs are kept, giving within-population kinship; this
#' requires adequate within-population sample sizes. \code{metric = "kinship"}
#' rescales the matrix by 1/2 (numerator-relationship to kinship); as a
#' constant applied to both datasets it leaves the correlation performance
#' unchanged. The panel relatedness data frame columns are
#' \code{rel_orig} and \code{rel_panel}.
#'
#' @return
#' A named list with one element per parameter, each containing:
#' \code{name}, \code{performance}, \code{corr.method}, \code{data},
#' \code{summary}, \code{confusion}, and \code{plot}. \code{corr.method}
#' is the correlation method used for correlation-based parameters and
#' \code{NA} otherwise. For most parameters,
#' \code{summary} and \code{confusion} are \code{NULL}; for
#' \code{"drift_resistance"}, \code{summary} contains scores by population;
#' for \code{"Fis_ind"}, \code{summary} contains correlation and error
#' diagnostics; for \code{"hybridisation"}, \code{summary} gives
#' per-cross F1 detection
#' and correct-pair rates, and \code{confusion} gives the full true-cross by
#' assigned-class table. The \code{performance} value for
#' \code{"hybridisation"} is the correct F1-pair rate.
#'
#' @references
#' VanRaden, P.M. (2008). Efficient methods to compute genomic predictions.
#' \emph{Journal of Dairy Science}, 91(11), 4414--4423.
#'
#' @examples
#' panel <- gl.select.panel(possums.gl, method = "random", nl = 50)
#' res   <- gl.check.panel(panel, possums.gl, parameter = "all-Ne")
#' sapply(res, function(r) r$performance)
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_col geom_histogram geom_bar geom_vline geom_hline
#'   geom_abline geom_point geom_smooth geom_tile geom_text annotate labs scale_x_log10
#'   scale_fill_manual scale_fill_gradient scale_fill_gradient2 coord_cartesian coord_equal
#'   theme element_text
#' @importFrom patchwork wrap_plots

gl.check.panel <- function(x,
                           xorig       = NULL,
                           parameter   = "Fst",
                           corr.method = c("spearman", "pearson"),
                           ref         = c("global", "by.pop"),
                           metric      = c("grm", "kinship"),
                           inverse_dr  = FALSE,
                           neest.path  = NULL,
                           error.rate  = 0.01,
                           threshold   = 0.001,
                           n_sim_parents = 100,
                           n_sim_hyb   = 100,
                           n_cores     = NULL,
                           target      = 0.9,
                           plot.out    = TRUE,
                           plot.file   = NULL,
                           plot.dir    = NULL,
                           verbose     = NULL) {
  
  # ---------------------------------------------------------------
  # expand shortcuts and validate
  # ---------------------------------------------------------------
  all_params <- c("Fst","He","Ho","Fis","Nall","Ne",
                  "Ho_ind","Fis_ind","drift_resistance","relatedness",
                  "id","parentage","assignment")
  
  if (length(parameter) == 1 && parameter == "all")
    parameter <- all_params
  if (length(parameter) == 1 && parameter == "all-Ne")
    parameter <- setdiff(all_params, "Ne")
  
  parameter[parameter == "Ar"] <- "Nall"
  parameter <- unique(parameter)
  corr.method <- match.arg(tolower(corr.method), c("spearman", "pearson"))
  ref         <- match.arg(ref,    c("global", "by.pop"))
  metric      <- match.arg(metric, c("grm", "kinship"))
  correlation_params <- c("Fst", "He", "Ho", "Fis", "Nall", "Ne",
                          "Ho_ind", "Fis_ind", "relatedness")
  comparison_params  <- correlation_params
  
  valid_params <- c(all_params, "hybridisation")
  bad <- setdiff(parameter, valid_params)
  if (length(bad) > 0)
    stop(paste("Unknown parameter(s):", paste(bad, collapse=", "),
               "\nValid:", paste(valid_params, collapse=", ")))

  # xorig is needed only for correlation/comparison parameters.
  # drift_resistance and the power metrics run on x alone.
  if (is.null(xorig) && any(parameter %in% comparison_params)) {
    need <- intersect(parameter, comparison_params)
    stop(paste0("xorig is required for parameter(s): ",
                paste(need, collapse=", "),
                ".\nLeave xorig = NULL only when requesting ",
                "drift_resistance, id, parentage, assignment, or ",
                "hybridisation."))
  }

  if (error.rate  < 0 || error.rate  > 1) stop("error.rate must be in [0,1].")
  if (threshold   <= 0 || threshold  >= 1) stop("threshold must be in (0,1).")
  if (!is.numeric(n_sim_parents) || n_sim_parents < 1) stop("n_sim_parents must be a positive integer.")
  if (!is.numeric(n_sim_hyb) || n_sim_hyb < 1) stop("n_sim_hyb must be a positive integer.")
  if (!is.null(target) && (target <= 0 || target > 1)) stop("target must be in (0,1].")
  
  if ("hybridisation" %in% parameter && nPop(x) < 2)
    stop(paste0("'hybridisation' requires at least 2 populations. Found: ",
                nPop(x), "."))
  
  n_cores <- if (is.null(n_cores)) max(1L, parallel::detectCores()-1L) else
    max(1L, as.integer(n_cores))
  
  # resolve verbosity: 0=silent, 1=headers+scores, 2=full (default)
  verbose <- if (is.null(verbose)) 2L else as.integer(verbose)
  
  # ---------------------------------------------------------------
  # order populations
  # ---------------------------------------------------------------
  x     <- x[order(pop(x)),]
  if (!is.null(xorig)) xorig <- xorig[order(pop(xorig)),]
  
  # ---------------------------------------------------------------
  # shared helpers
  # ---------------------------------------------------------------
  make_trans <- function(e) {
    matrix(c(
      (1-e)^2,  2*e*(1-e),        e^2,
      e*(1-e),  (1-e)^2+e^2,      e*(1-e),
      e^2,      2*e*(1-e),        (1-e)^2
    ), nrow=3, byrow=TRUE)
  }
  
  r2_cor <- function(a, b, method = corr.method) {
    method <- match.arg(method, c("spearman", "pearson"))
    ok <- complete.cases(a, b)
    if (sum(ok) < 2L) return(NA_real_)
    as.numeric(suppressWarnings(
      cor(a[ok], b[ok], method = method)^2
    ))
  }
  
  r2_spearman <- function(a, b)
    r2_cor(a, b, method = "spearman")
  
  r2_pearson <- function(a, b)
    r2_cor(a, b, method = "pearson")
  
  r2_performance <- function(a, b)
    r2_cor(a, b, method = corr.method)
  
  corr_metric_label <- function(method = corr.method)
    if (method == "spearman") "Spearman r2" else "Pearson R2"
  
  r2_label <- function(sp, pe)
    sprintf("Spearman r\u00b2 = %.3f%s\nPearson R\u00b2 = %.3f%s",
            sp, if (corr.method == "spearman") " (perf)" else "",
            pe, if (corr.method == "pearson") " (perf)" else "")
  
  centre_theme <- function()
    theme(plot.title    = element_text(hjust=0.5),
          plot.subtitle = element_text(hjust=0.5))
  
  # formula = y ~ x suppresses the geom_smooth method message
  scatter_plot <- function(res, xcol, ycol, title, xlab, ylab) {
    sp <- r2_spearman(res[[xcol]], res[[ycol]])
    pe <- r2_pearson( res[[xcol]], res[[ycol]])
    ggplot(res, aes(x=.data[[xcol]], y=.data[[ycol]])) +
      geom_point(alpha=0.7) +
      geom_smooth(method="lm", formula=y~x, se=TRUE, colour="steelblue") +
      annotate("text", x=-Inf, y=Inf,
               label=r2_label(sp, pe),
               hjust=-0.1, vjust=1.3, size=3.5) +
      labs(title=title, x=xlab, y=ylab) +
      centre_theme()
  }
  
  # ---------------------------------------------------------------
  # VanRaden GRM with external reference allele frequencies
  # ---------------------------------------------------------------
  grm_vanraden <- function(gmat, p) {
    gmat_imp <- gmat
    na_loci  <- which(apply(gmat, 2, anyNA))
    for (l in na_loci) {
      na_idx              <- is.na(gmat[,l])
      gmat_imp[na_idx, l] <- 2*p[l]
    }
    Z     <- sweep(gmat_imp, 2, 2*p, "-")
    denom <- 2*sum(p*(1-p), na.rm=TRUE)
    (Z %*% t(Z)) / denom
  }
  
  # ---------------------------------------------------------------
  # cache heterozygosity — shared across He / Ho / Fis
  # ---------------------------------------------------------------
  het_cache <- if (any(c("He","Ho","Fis") %in% parameter)) {
    if (verbose >= 2) cat("Computing heterozygosity (shared for He, Ho, Fis)...\n")
    list(
      orig  = gl.report.heterozygosity(xorig, verbose=0),
      panel = gl.report.heterozygosity(x,     verbose=0)
    )
  } else NULL
  
  # ---------------------------------------------------------------
  # main loop
  # ---------------------------------------------------------------
  output <- vector("list", length(parameter))
  names(output) <- parameter
  
  for (param in parameter) {
    
    if (verbose >= 1) cat(sprintf("\n--- %s ---\n", param))
    
    gg           <- NULL
    res          <- NULL
    perf         <- NA_real_
    param_summary <- NULL
    confusion     <- NULL
    
    # ---------------------------------------------------------
    # Fst
    # ---------------------------------------------------------
    if (param == "Fst") {
      fs_orig  <- gl.report.fstat(xorig, verbose=0, plot.display=FALSE)
      fs_panel <- gl.report.fstat(x,     verbose=0, plot.display=FALSE)
      res <- data.frame(
        fst_orig  = as.numeric(fs_orig$Stat_matrices$Fstp[
          upper.tri(fs_orig$Stat_matrices$Fstp)]),
        fst_panel = as.numeric(fs_panel$Stat_matrices$Fstp[
          upper.tri(fs_panel$Stat_matrices$Fstp)])
      )
      res  <- res[complete.cases(res),]
      perf <- r2_performance(res$fst_orig, res$fst_panel)
      gg   <- scatter_plot(res, "fst_orig", "fst_panel",
                           "Fst\u2019",
                           "Fst\u2019 (original)", "Fst\u2019 (panel)")
    }
    
    # ---------------------------------------------------------
    # He
    # ---------------------------------------------------------
    if (param == "He") {
      res  <- data.frame(he_orig  = het_cache$orig$He,
                         he_panel = het_cache$panel$He)
      res  <- res[complete.cases(res),]
      perf <- r2_performance(res$he_orig, res$he_panel)
      gg   <- scatter_plot(res, "he_orig", "he_panel",
                           "He", "He (original)", "He (panel)")
    }
    
    # ---------------------------------------------------------
    # Ho
    # ---------------------------------------------------------
    if (param == "Ho") {
      res  <- data.frame(ho_orig  = het_cache$orig$Ho,
                         ho_panel = het_cache$panel$Ho)
      res  <- res[complete.cases(res),]
      perf <- r2_performance(res$ho_orig, res$ho_panel)
      gg   <- scatter_plot(res, "ho_orig", "ho_panel",
                           "Ho", "Ho (original)", "Ho (panel)")
    }
    
    # ---------------------------------------------------------
    # Fis
    # ---------------------------------------------------------
    if (param == "Fis") {
      res  <- data.frame(fis_orig  = het_cache$orig$FIS,
                         fis_panel = het_cache$panel$FIS)
      res  <- res[complete.cases(res),]
      perf <- r2_performance(res$fis_orig, res$fis_panel)
      gg   <- scatter_plot(res, "fis_orig", "fis_panel",
                           "Fis", "Fis (original)", "Fis (panel)")
    }
    
    # ---------------------------------------------------------
    # Nall
    # ---------------------------------------------------------
    if (param == "Nall") {
      ar_orig  <- gl.report.allelerich(xorig, verbose=0, plot.display=FALSE)
      ar_panel <- gl.report.allelerich(x,     verbose=0, plot.display=FALSE)
      res <- data.frame(
        ar_orig  = ar_orig$`Allelic Richness per population`$mean_corrected_richness,
        ar_panel = ar_panel$`Allelic Richness per population`$mean_corrected_richness
      )
      res  <- res[complete.cases(res),]
      perf <- r2_performance(res$ar_orig, res$ar_panel)
      gg   <- scatter_plot(res, "ar_orig", "ar_panel",
                           "Allelic richness",
                           "Ar (original)", "Ar (panel)")
    }
    
    # ---------------------------------------------------------
    # Ne
    # ---------------------------------------------------------
    if (param == "Ne") {
      if (is.null(neest.path)) stop("neest.path required for 'Ne'.")
      ne_orig  <- gl.LDNe(xorig, neest.path=neest.path, critical=0.05,
                          verbose=0, mating="random", plot.out = F)
      ne_panel <- gl.LDNe(x,     neest.path=neest.path, critical=0.05,
                          verbose=0, mating="random", singleton.rm=FALSE, plot.out = F)
      extract_ne <- function(nl)
        as.numeric(unlist(lapply(nl, function(z) z$`Frequency 1`[6])))
      res <- data.frame(ne_orig=extract_ne(ne_orig), ne_panel=extract_ne(ne_panel))
      res[res==Inf|res==-Inf] <- NA
      res  <- res[complete.cases(res),]
      perf <- r2_performance(res$ne_orig, res$ne_panel)
      gg   <- scatter_plot(res, "ne_orig", "ne_panel",
                           "Ne", "Ne (original)", "Ne (panel)")
    }
    
    # ---------------------------------------------------------
    # Ho_ind — individual observed heterozygosity
    # population membership ignored
    # ---------------------------------------------------------
    if (param == "Ho_ind") {
      ho_ind_orig  <- rowMeans(as.matrix(xorig)==1L, na.rm=TRUE)
      ho_ind_panel <- rowMeans(as.matrix(x)    ==1L, na.rm=TRUE)
      res <- data.frame(
        individual   = indNames(x),
        population   = as.character(pop(x)),
        ho_ind_orig  = ho_ind_orig,
        ho_ind_panel = ho_ind_panel
      )
      res  <- res[complete.cases(res[,c("ho_ind_orig","ho_ind_panel")]),]
      perf <- r2_performance(res$ho_ind_orig, res$ho_ind_panel)
      if (verbose >= 2)
        cat(sprintf("Individuals: %d  |  %s = %.4f\n",
                    nrow(res), corr_metric_label(), perf))
      gg   <- scatter_plot(res, "ho_ind_orig", "ho_ind_panel",
                           "Individual heterozygosity",
                           "Ho_ind (original)", "Ho_ind (panel)")
    }

    # ---------------------------------------------------------
    # Fis_ind — actual individual Fis from the panel versus the
    # full original dataset
    #
    # Within each population and for each individual's called loci:
    #   Ho_i  = number of heterozygous loci / number of called loci
    #   He_i  = sum[2 p_l (1-p_l)] / number of called loci
    #   Fis_i = 1 - Ho_i / He_i
    #
    # Performance is the selected corr.method squared correlation between
    # panel and full-data individual Fis, matching the implementation used
    # for Fst, He, Ho, Fis, Nall, Ne and relatedness.
    # ---------------------------------------------------------
    if (param == "Fis_ind") {
      ind_orig <- indNames(xorig)
      ind_pan  <- indNames(x)
      idx_pan_ind <- match(ind_orig, ind_pan)

      if (anyNA(idx_pan_ind)) {
        missing_ind <- ind_orig[is.na(idx_pan_ind)]
        stop(paste0(
          "Fis_ind: panel object is missing ", length(missing_ind),
          " individual(s) present in xorig: ",
          paste(head(missing_ind, 10), collapse = ", "),
          if (length(missing_ind) > 10) " ..." else ""
        ))
      }

      mat_orig <- as.matrix(xorig)
      mat_pan  <- as.matrix(x)[idx_pan_ind, , drop = FALSE]
      loci_all <- locNames(xorig)
      loci_pan <- locNames(x)
      idx_pan_loci <- match(loci_pan, loci_all)

      if (anyNA(idx_pan_loci)) {
        stop("Fis_ind: one or more panel loci are not present in xorig.")
      }

      pops_chr <- as.character(pop(xorig))
      upops    <- unique(pops_chr)
      n_ind    <- nrow(mat_orig)

      fis_orig <- rep(NA_real_, n_ind)
      fis_panel <- rep(NA_real_, n_ind)
      ho_orig <- rep(NA_real_, n_ind)
      ho_panel <- rep(NA_real_, n_ind)
      he_orig <- rep(NA_real_, n_ind)
      he_panel <- rep(NA_real_, n_ind)
      n_called_orig <- integer(n_ind)
      n_called_panel <- integer(n_ind)

      calc_ind_fis <- function(gmat, he_locus) {
        called <- !is.na(gmat)
        n_called <- rowSums(called)
        n_het <- rowSums(gmat == 1L, na.rm = TRUE)

        he_mat <- matrix(
          rep(he_locus, each = nrow(gmat)),
          nrow = nrow(gmat),
          ncol = ncol(gmat)
        )
        exp_het_sum <- rowSums(he_mat * called, na.rm = TRUE)

        Ho <- ifelse(n_called > 0L, n_het / n_called, NA_real_)
        He <- ifelse(n_called > 0L, exp_het_sum / n_called, NA_real_)
        Fis <- ifelse(is.finite(He) & He > 0,
                      1 - Ho / He,
                      NA_real_)

        list(Fis = Fis, Ho = Ho, He = He, n_called = n_called)
      }

      for (pp in upops) {
        idx <- which(pops_chr == pp)
        mo  <- mat_orig[idx, , drop = FALSE]
        mp  <- mat_pan[idx, , drop = FALSE]

        p_full <- colMeans(mo, na.rm = TRUE) / 2
        p_full[!is.finite(p_full)] <- NA_real_
        he_full <- 2 * p_full * (1 - p_full)
        he_pan  <- he_full[idx_pan_loci]

        fo <- calc_ind_fis(mo, he_full)
        fp <- calc_ind_fis(mp, he_pan)

        fis_orig[idx] <- fo$Fis
        fis_panel[idx] <- fp$Fis
        ho_orig[idx] <- fo$Ho
        ho_panel[idx] <- fp$Ho
        he_orig[idx] <- fo$He
        he_panel[idx] <- fp$He
        n_called_orig[idx] <- fo$n_called
        n_called_panel[idx] <- fp$n_called
      }

      res <- data.frame(
        individual = ind_orig,
        population = pops_chr,
        fis_orig = fis_orig,
        fis_panel = fis_panel,
        difference = fis_panel - fis_orig,
        ho_orig = ho_orig,
        ho_panel = ho_panel,
        he_orig = he_orig,
        he_panel = he_panel,
        n_called_orig = n_called_orig,
        n_called_panel = n_called_panel,
        stringsAsFactors = FALSE
      )
      res <- res[complete.cases(res[, c("fis_orig", "fis_panel")]), ]

      perf <- r2_performance(res$fis_orig, res$fis_panel)

      rmse <- if (nrow(res) > 0L)
        sqrt(mean((res$fis_panel - res$fis_orig)^2)) else NA_real_
      mae <- if (nrow(res) > 0L)
        mean(abs(res$fis_panel - res$fis_orig)) else NA_real_
      bias <- if (nrow(res) > 0L)
        mean(res$fis_panel - res$fis_orig) else NA_real_
      pearson_r2 <- r2_pearson(res$fis_orig, res$fis_panel)
      spearman_r2 <- r2_spearman(res$fis_orig, res$fis_panel)

      param_summary <- data.frame(
        n_individuals = nrow(res),
        n_loci_panel = nLoc(x),
        n_loci_original = nLoc(xorig),
        mean_fis_original = mean(res$fis_orig),
        mean_fis_panel = mean(res$fis_panel),
        corr_method = corr.method,
        performance = perf,
        rmse = rmse,
        mae = mae,
        bias_panel_minus_original = bias,
        pearson_R2 = pearson_r2,
        spearman_r2 = spearman_r2,
        stringsAsFactors = FALSE
      )

      if (verbose >= 2) {
        cat(sprintf(
          paste0(
            "Individuals: %d | mean Fis original: %.4f | ",
            "mean Fis panel: %.4f\n",
            "%s/performance: %.4f | RMSE: %.4f | MAE: %.4f | ",
            "bias: %.4f\n"
          ),
          nrow(res), mean(res$fis_orig), mean(res$fis_panel),
          corr_metric_label(), perf, rmse, mae, bias
        ))
      }

      gg <- scatter_plot(
        res,
        "fis_orig",
        "fis_panel",
        "Individual Fis",
        "Fis_ind (original)",
        "Fis_ind (panel)"
      )
    }

    # ---------------------------------------------------------
    # drift_resistance — absolute within-population MAF score
    #
    # For locus l in population k:
    #   score_lk = (2 * min(p_lk, 1 - p_lk))^2
    #
    # Loci are averaged within populations, then populations are
    # weighted equally. The score is 1 when every p_lk = 0.5 and
    # 0 when all loci are fixed within populations.
    # ---------------------------------------------------------
    if (param == "drift_resistance") {
      gmat     <- as.matrix(x)
      pops_chr <- as.character(pop(x))
      upops    <- unique(pops_chr)

      drift_rows <- lapply(upops, function(pp) {
        idx <- which(pops_chr == pp)
        gm  <- gmat[idx, , drop = FALSE]

        p_locus <- colMeans(gm, na.rm = TRUE) / 2
        p_locus[is.nan(p_locus)] <- NA_real_

        maf_locus   <- pmin(p_locus, 1 - p_locus)
        drift_score <- (2 * maf_locus)^2
        call_rate   <- colMeans(!is.na(gm))

        data.frame(
          population       = pp,
          locus            = locNames(x),
          allele_frequency = p_locus,
          maf              = maf_locus,
          call_rate        = call_rate,
          drift_score      = drift_score,
          stringsAsFactors = FALSE
        )
      })

      res <- do.call(rbind, drift_rows)
      rownames(res) <- NULL

      pop_summary <- do.call(
        rbind,
        lapply(upops, function(pp) {
          rr <- res[res$population == pp, , drop = FALSE]
          valid_score <- is.finite(rr$drift_score)
          data.frame(
            population       = pp,
            n_loci           = nrow(rr),
            n_loci_scored    = sum(valid_score),
            mean_maf         = if (any(valid_score))
              mean(rr$maf[valid_score]) else NA_real_,
            mean_drift_score = if (any(valid_score))
              mean(rr$drift_score[valid_score]) else NA_real_,
            stringsAsFactors = FALSE
          )
        })
      )

      valid_pop <- is.finite(pop_summary$mean_drift_score)
<<<<<<< HEAD
      raw_perf <- if (any(valid_pop))
        mean(pop_summary$mean_drift_score[valid_pop]) else NA_real_

      ## apply inversion if requested
      perf <- if (inverse_dr) 1 - raw_perf else raw_perf
      dr_label <- if (inverse_dr) "Drift sensitivity (1 \u2212 DR)"
                  else "Drift resistance"

      ## practical maximum: best achievable score from xorig
      ## DR: top-nl highest-MAF loci; inverse: top-nl lowest-MAF loci
=======
      perf <- if (any(valid_pop))
        mean(pop_summary$mean_drift_score[valid_pop]) else NA_real_

      ## practical maximum: score of the top-nLoc(x) MAF loci from xorig
>>>>>>> 6d2d167fb6f097250b1360fed5b015ada35b8f1d
      dr_max <- NA_real_
      if (!is.null(xorig)) {
        gmat_orig  <- as.matrix(xorig)
        pops_orig  <- as.character(pop(xorig))
        upops_orig <- unique(pops_orig)
        dr_all <- vapply(seq_len(ncol(gmat_orig)), function(l) {
          mean(vapply(upops_orig, function(pp) {
            p <- mean(gmat_orig[pops_orig == pp, l], na.rm = TRUE) / 2
            if (!is.finite(p)) return(NA_real_)
            (2 * min(p, 1 - p))^2
          }, numeric(1)), na.rm = TRUE)
        }, numeric(1))
<<<<<<< HEAD
        n_sel <- min(nLoc(x), length(dr_all))
        if (inverse_dr) {
          ## best for sensitivity = lowest DR loci
          bot_idx <- order(dr_all, decreasing = FALSE)[seq_len(n_sel)]
          dr_max  <- 1 - mean(dr_all[bot_idx], na.rm = TRUE)
        } else {
          top_idx <- order(dr_all, decreasing = TRUE)[seq_len(n_sel)]
          dr_max  <- mean(dr_all[top_idx], na.rm = TRUE)
        }
=======
        top_idx <- order(dr_all, decreasing = TRUE)[seq_len(min(nLoc(x), length(dr_all)))]
        dr_max  <- mean(dr_all[top_idx], na.rm = TRUE)
>>>>>>> 6d2d167fb6f097250b1360fed5b015ada35b8f1d
      }

      overall <- data.frame(
        population       = "Overall",
        n_loci           = nLoc(x),
        n_loci_scored    = sum(is.finite(res$drift_score)),
        mean_maf         = if (any(is.finite(res$maf)))
          mean(res$maf[is.finite(res$maf)]) else NA_real_,
        mean_drift_score = perf,
        practical_max    = dr_max,
<<<<<<< HEAD
        inverse_dr       = inverse_dr,
=======
>>>>>>> 6d2d167fb6f097250b1360fed5b015ada35b8f1d
        stringsAsFactors = FALSE
      )

      pop_summary$practical_max <- NA_real_
<<<<<<< HEAD
      pop_summary$inverse_dr    <- inverse_dr
      if (inverse_dr)
        pop_summary$mean_drift_score <- 1 - pop_summary$mean_drift_score
=======
>>>>>>> 6d2d167fb6f097250b1360fed5b015ada35b8f1d
      param_summary <- rbind(pop_summary, overall)

      if (verbose >= 2) {
        cat(sprintf(
          paste0(
<<<<<<< HEAD
            "%s score: %.4f | %d loci | %d populations",
            if (is.finite(dr_max)) " | practical max: %.4f" else "",
            if (inverse_dr) "\nScore = 1 - mean_pop(mean_locus[(2 * MAF)^2])\n"
            else "\nScore = mean_pop(mean_locus[(2 * MAF)^2])\n"
          ),
          dr_label, perf, nLoc(x), length(upops),
=======
            "Drift-resistance score: %.4f | %d loci | %d populations",
            if (is.finite(dr_max)) " | practical max: %.4f" else "",
            "\nScore = mean_pop(mean_locus[(2 * MAF)^2])\n"
          ),
          perf, nLoc(x), length(upops),
>>>>>>> 6d2d167fb6f097250b1360fed5b015ada35b8f1d
          if (is.finite(dr_max)) dr_max
        ))
      }

<<<<<<< HEAD
      plot_data <- pop_summary
      gg <- ggplot(
        plot_data,
=======
      gg <- ggplot(
        pop_summary,
>>>>>>> 6d2d167fb6f097250b1360fed5b015ada35b8f1d
        aes(x = population, y = mean_drift_score, fill = mean_drift_score)
      ) +
        geom_col(colour = "white", linewidth = 0.3) +
        scale_fill_gradient2(
<<<<<<< HEAD
          low = if (inverse_dr) "steelblue" else "salmon",
          mid = "goldenrod",
          high = if (inverse_dr) "salmon" else "steelblue",
=======
          low = "salmon", mid = "goldenrod", high = "steelblue",
>>>>>>> 6d2d167fb6f097250b1360fed5b015ada35b8f1d
          midpoint = 0.5, limits = c(0, 1), name = NULL
        ) +
        coord_cartesian(ylim = c(0, 1))

      if (is.finite(dr_max))
        gg <- gg +
        geom_hline(yintercept = dr_max, linetype = "dotted",
                   colour = "darkorange", linewidth = 0.6) +
        annotate("text", x = Inf, y = dr_max,
<<<<<<< HEAD
                 label = sprintf(" %s max=%.3f ",
                                 if (inverse_dr) "sens." else "DR", dr_max),
=======
                 label = sprintf(" max=%.3f ", dr_max),
>>>>>>> 6d2d167fb6f097250b1360fed5b015ada35b8f1d
                 hjust = 1, vjust = -0.4, size = 3, colour = "darkorange")

      gg <- gg +
        labs(
<<<<<<< HEAD
          title = dr_label,
          subtitle = sprintf(
            "Overall score = %.3f | %s%s",
            perf,
            if (inverse_dr) "1 - mean of within-population (2 x MAF)^2"
            else "mean of within-population (2 x MAF)^2",
            if (is.finite(dr_max))
              sprintf(" | practical max = %.3f (top-%d loci)", dr_max, nLoc(x))
            else ""
          ),
          x = "Population",
          y = if (inverse_dr) "Drift-sensitivity score (1-DR)" else "Drift-resistance score"
=======
          title = "Drift resistance",
          subtitle = sprintf(
            "Overall score = %.3f | mean of within-population (2 \u00d7 MAF)\u00b2%s",
            perf,
            if (is.finite(dr_max))
              sprintf(" | practical max = %.3f (top-%d MAF loci)", dr_max, nLoc(x))
            else ""
          ),
          x = "Population",
          y = "Drift-resistance score"
>>>>>>> 6d2d167fb6f097250b1360fed5b015ada35b8f1d
        ) +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)
        )
    }

    # ---------------------------------------------------------
    # relatedness — pairwise VanRaden GRM, optionally rescaled to
    # kinship (GRM / 2). Reference allele frequencies are either
    # pooled across the whole dataset (ref = "global") or computed
    # within each population (ref = "by.pop"). Under "by.pop" only
    # within-population pairs are retained and each population is
    # referenced to its own allele-frequency base.
    # ---------------------------------------------------------
    if (param == "relatedness") {
      scale_fac  <- if (metric == "kinship") 0.5 else 1
      rel_label  <- if (metric == "kinship") "Kinship" else "GRM"
      ref_label  <- if (ref == "by.pop") "within-population freq"
      else "full-dataset freq"

      # align panel individuals to xorig order so the upper-triangle
      # pairs from the two matrices correspond row-for-row
      mat_orig <- as.matrix(xorig)
      mat_pan  <- as.matrix(x)[indNames(xorig), , drop = FALSE]
      loc_pan  <- colnames(mat_pan)

      if (ref == "global") {
        if (verbose >= 2)
          cat("Computing global reference allele frequencies...\n")
        p_ref <- colMeans(mat_orig, na.rm = TRUE) / 2
        p_ref[is.nan(p_ref)] <- 0.5

        grm_orig    <- grm_vanraden(mat_orig, p_ref)            * scale_fac
        grm_panel   <- grm_vanraden(mat_pan,  p_ref[loc_pan])   * scale_fac
        pairs_orig  <- grm_orig[ upper.tri(grm_orig)]
        pairs_panel <- grm_panel[upper.tri(grm_panel)]

      } else {  # by.pop
        if (verbose >= 2)
          cat("Computing within-population reference allele frequencies...\n")
        pops_chr <- as.character(pop(xorig))
        upops    <- unique(pops_chr)
        pairs_orig <- pairs_panel <- numeric(0)
        skipped   <- character(0)

        for (pp in upops) {
          idx <- which(pops_chr == pp)
          if (length(idx) < 2) { skipped <- c(skipped, pp); next }
          mo   <- mat_orig[idx, , drop = FALSE]
          mp   <- mat_pan[ idx, , drop = FALSE]
          p_pp <- colMeans(mo, na.rm = TRUE) / 2
          p_pp[is.nan(p_pp)] <- 0.5
          g_o  <- grm_vanraden(mo, p_pp)            * scale_fac
          g_p  <- grm_vanraden(mp, p_pp[loc_pan])   * scale_fac
          pairs_orig  <- c(pairs_orig,  g_o[upper.tri(g_o)])
          pairs_panel <- c(pairs_panel, g_p[upper.tri(g_p)])
        }
        if (length(skipped) > 0)
          warning(sprintf(
            "relatedness (by.pop): %d population(s) with <2 individuals skipped: %s",
            length(skipped), paste(skipped, collapse = ", ")))
        if (length(pairs_orig) == 0)
          warning("relatedness (by.pop): no within-population pairs available.")
      }

      res  <- data.frame(rel_orig = pairs_orig, rel_panel = pairs_panel)
      res  <- res[complete.cases(res), ]
      perf <- r2_performance(res$rel_orig, res$rel_panel)

      if (verbose >= 2)
        cat(sprintf("Pairs: %d  |  %s (%s)  |  %s = %.4f\n",
                    nrow(res), rel_label, ref, corr_metric_label(), perf))

      plot_res <- if (nrow(res) > 10000L) res[sample(nrow(res), 10000L), ] else res
      sp <- r2_spearman(res$rel_orig, res$rel_panel)
      pe <- r2_pearson( res$rel_orig, res$rel_panel)

      gg <- ggplot(plot_res, aes(x = rel_orig, y = rel_panel)) +
        geom_point(alpha = 0.3, size = 0.8) +
        geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "steelblue") +
        annotate("text", x = -Inf, y = Inf,
                 label = r2_label(sp, pe),
                 hjust = -0.1, vjust = 1.3, size = 3.5) +
        labs(title    = sprintf("Relatedness (%s)", rel_label),
             subtitle = sprintf("%d pairs  |  %s, %s%s",
                                 nrow(res), ref, ref_label,
                                 if (nrow(res) > 10000L) "  |  10,000 pairs plotted" else ""),
             x = sprintf("%s (original)", rel_label),
             y = sprintf("%s (panel)",    rel_label)) +
        centre_theme()
    }
    
    # ---------------------------------------------------------
    # id
    # ---------------------------------------------------------
    if (param == "id") {
      trans <- make_trans(error.rate)
      gmat  <- as.matrix(x)
      n     <- nrow(gmat)
      L     <- ncol(gmat)
      
      if (verbose >= 2)
        cat(sprintf("Computing pairwise match probabilities (%d ind, %d loci)...\n", n, L))
      
      log_pm <- matrix(0, n, n)
      for (l in seq_len(L)) {
        g     <- gmat[,l]; valid <- !is.na(g)
        p_obs <- matrix(1/3, nrow=n, ncol=3)
        if (any(valid)) p_obs[valid,] <- trans[g[valid]+1,]
        pm_l <- p_obs %*% t(p_obs)
        if (any(!valid)) {
          na_idx         <- which(!valid)
          pm_l[na_idx,]  <- 1; pm_l[,na_idx] <- 1
        }
        log_pm <- log_pm + log(pmax(pm_l, .Machine$double.eps))
      }
      pm <- exp(log_pm); diag(pm) <- 0
      maxpm        <- apply(pm, 1, max)
      n_identified <- sum(maxpm < threshold)
      perf         <- n_identified / n
      
      if (verbose >= 2)
        cat(sprintf("Identified: %d / %d (%.1f%%)\n", n_identified, n, 100*perf))
      
      res <- data.frame(
        individual = indNames(x),
        population = as.character(pop(x)),
        maxpm      = maxpm,
        identified = maxpm < threshold
      )
      gg <- ggplot(res, aes(x=maxpm)) +
        geom_histogram(aes(fill=identified), bins=40,
                       colour="white", linewidth=0.2) +
        geom_vline(xintercept=threshold, linetype="dashed",
                   colour="firebrick", linewidth=0.8) +
        scale_x_log10() +
        scale_fill_manual(
          values=c("TRUE"="steelblue","FALSE"="salmon"),
          labels=c("TRUE"="Identified","FALSE"="Not identified"),
          name=NULL) +
        labs(title    = "Individual ID",
             subtitle = sprintf("%d / %d (%.1f%%) identified  (MaxPmatch < %g, e=%.3f)",
                                n_identified, n, 100*perf, threshold, error.rate),
             x="Max P(match) to nearest neighbour  [log10]", y="Count") +
        centre_theme()
    }
    
    # ---------------------------------------------------------
    # parentage — fork-parallel on Mac/Linux, sequential on Windows
    # ---------------------------------------------------------
    if (param == "parentage") {
      e     <- error.rate
      trans <- make_trans(e)
      gmat  <- as.matrix(x)
      n     <- nrow(gmat); L <- ncol(gmat)
      a1    <- c(0, 0.5, 1)
      T_0   <- outer(1-a1, 1-a1); T_2 <- outer(a1, a1); T_1 <- 1-T_0-T_2
      
      A0 <- A1 <- A2 <- vector("list", L)
      na_locus        <- vector("list", L)
      
      if (verbose >= 2)
        cat(sprintf("Precomputing A matrices (%d loci, %d ind)...\n", L, n))
      for (l in seq_len(L)) {
        g             <- gmat[,l]; valid <- !is.na(g)
        na_locus[[l]] <- which(!valid)
        p_obs         <- matrix(1/3, nrow=n, ncol=3)
        if (any(valid)) p_obs[valid,] <- trans[g[valid]+1,]
        A0[[l]] <- (p_obs %*% T_0) %*% t(p_obs)
        A1[[l]] <- (p_obs %*% T_1) %*% t(p_obs)
        A2[[l]] <- (p_obs %*% T_2) %*% t(p_obs)
      }
      
      pair_mat  <- which(upper.tri(matrix(0,n,n)), arr.ind=TRUE)
      cores_use <- if (.Platform$OS.type=="windows") {
        #if (n_cores>1) message("Parallel parentage disabled on Windows. Sequential.")
        1L
      } else min(n_cores, n_sim_parents)
      
      if (verbose >= 2)
        cat(sprintf("Simulating %d trios  (%s)...\n", n_sim_parents,
                    if (cores_use>1) sprintf("%d cores, fork",cores_use) else "sequential"))
      
      run_one_sim <- function(idx) {
        p_idx <- sample(n,2,replace=FALSE); p1 <- p_idx[1]; p2 <- p_idx[2]
        o_obs <- integer(L)
        for (l in seq_len(L)) {
          g1 <- gmat[p1,l]; g2 <- gmat[p2,l]
          ap1 <- if (is.na(g1)) 0.5 else a1[g1+1]
          ap2 <- if (is.na(g2)) 0.5 else a1[g2+1]
          ot  <- rbinom(1,1,ap1)+rbinom(1,1,ap2)
          o_obs[l] <- sample(0:2,1,prob=trans[ot+1,])
        }
        log_LL <- matrix(0,nrow=n,ncol=n)
        for (l in seq_len(L)) {
          p_o     <- trans[o_obs[l]+1,]
          contrib <- p_o[1]*A0[[l]]+p_o[2]*A1[[l]]+p_o[3]*A2[[l]]
          if (length(na_locus[[l]])>0) {
            idx_na           <- na_locus[[l]]
            contrib[idx_na,] <- 1; contrib[,idx_na] <- 1
          }
          log_LL <- log_LL+log(pmax(contrib,.Machine$double.eps))
        }
        log_LL    <- (log_LL+t(log_LL))/2
        L_true    <- log_LL[p1,p2]
        pair_LL   <- log_LL[pair_mat]
        true_pair <- pair_mat[,1]==min(p1,p2) & pair_mat[,2]==max(p1,p2)
        confusion <- mean(pair_LL[!true_pair]>=L_true)
        list(p1=p1, p2=p2, confusion=confusion, metric=1-confusion)
      }
      
      sim_results <- if (cores_use>1) {
        rp <- parallel::mclapply(seq_len(n_sim_parents), run_one_sim, mc.cores=cores_use)
        fl <- sapply(rp, inherits, "try-error")
        if (any(fl)) { warning(sum(fl)," sim(s) failed; rerunning."); rp[fl] <- lapply(which(fl),run_one_sim) }
        rp
      } else lapply(seq_len(n_sim_parents), run_one_sim)
      
      res <- data.frame(
        sim       = seq_len(n_sim_parents),
        p1        = sapply(sim_results,`[[`,"p1"),
        p2        = sapply(sim_results,`[[`,"p2"),
        confusion = sapply(sim_results,`[[`,"confusion"),
        metric    = sapply(sim_results,`[[`,"metric")
      )
      n_assigned <- sum(res$confusion<threshold)
      perf       <- n_assigned/n_sim_parents
      
      if (verbose >= 2)
        cat(sprintf("Assigned: %d / %d (%.1f%%)  (confusion < %g)\n",
                    n_assigned, n_sim_parents, 100*perf, threshold))
      
      gg <- ggplot(res, aes(x=metric)) +
        geom_histogram(aes(fill=confusion<threshold), bins=40,
                       colour="white", linewidth=0.2) +
        geom_vline(xintercept=1-threshold, linetype="dashed",
                   colour="firebrick", linewidth=0.8) +
        scale_fill_manual(
          values=c("TRUE"="steelblue","FALSE"="salmon"),
          labels=c("TRUE"="Assigned","FALSE"="Not assigned"),
          name=NULL) +
        labs(title    = "Parentage",
             subtitle = sprintf("%d / %d (%.1f%%) assigned  (confusion < %g, e=%.3f)",
                                n_assigned, n_sim_parents, 100*perf, threshold, e),
             x="1 \u2212 confusion fraction", y="Count") +
        centre_theme()
    }
    
    # ---------------------------------------------------------
    # assignment
    # ---------------------------------------------------------
    if (param == "assignment") {
      e     <- error.rate
      trans <- make_trans(e)
      gmat  <- as.matrix(x)
      n     <- nrow(gmat); L <- ncol(gmat)
      pop_f <- pop(x)
      empty_levels <- setdiff(levels(pop_f), as.character(unique(pop_f)))
      if (length(empty_levels) > 0) {
        warning(sprintf(
          "assignment: dropping %d population(s) with zero individuals: %s",
          length(empty_levels), paste(empty_levels, collapse=", ")))
        pop_f <- droplevels(pop_f)
      }
      pops  <- as.character(pop_f)
      upops <- levels(pop_f); K <- length(upops)
      
      if (K<2) stop("At least 2 populations required for 'assignment'.")
      
      af <- matrix(NA, nrow=K, ncol=L)
      for (k in seq_len(K)) {
        idx_k <- which(pops==upops[k])
        for (l in seq_len(L))
          af[k,l] <- pmax(pmin(mean(gmat[idx_k,l],na.rm=TRUE)/2,1-1e-6),1e-6)
      }
      
      ll_global <- array(NA, dim=c(L,3,K))
      for (k in seq_len(K)) for (l in seq_len(L)) {
        p <- af[k,l]; hwe <- c((1-p)^2,2*p*(1-p),p^2)
        ll_global[l,,k] <- as.numeric(hwe %*% trans)
      }
      
      af_loo <- matrix(NA, nrow=n, ncol=L)
      for (i in seq_len(n)) {
        k_i <- which(upops==pops[i])
        idx_loo <- setdiff(which(pops==pops[i]),i)
        for (l in seq_len(L))
          af_loo[i,l] <- if (length(idx_loo)==0) af[k_i,l] else {
            g_loo <- gmat[idx_loo,l]
            if (all(is.na(g_loo))) af[k_i,l] else
              pmax(pmin(mean(g_loo,na.rm=TRUE)/2,1-1e-6),1e-6)
          }
      }
      
      if (verbose >= 2) cat(sprintf("LOO assignment for %d individuals...\n", n))
      
      assigned_pop <- character(n)
      for (i in seq_len(n)) {
        k_i   <- which(upops==pops[i]); g_obs <- gmat[i,]; log_L <- numeric(K)
        for (l in seq_len(L)) {
          if (is.na(g_obs[l])) next
          ll_l <- ll_global[l,g_obs[l]+1,]
          p_loo <- af_loo[i,l]
          hwe_loo <- c((1-p_loo)^2,2*p_loo*(1-p_loo),p_loo^2)
          ll_l[k_i] <- as.numeric(hwe_loo %*% trans)[g_obs[l]+1]
          log_L     <- log_L+log(pmax(ll_l,.Machine$double.eps))
        }
        assigned_pop[i] <- upops[which.max(log_L)]
      }
      
      res <- data.frame(
        individual         = indNames(x),
        population         = pops,
        assigned_pop       = assigned_pop,
        correctly_assigned = pops==assigned_pop,
        stringsAsFactors   = FALSE
      )
      n_correct <- sum(res$correctly_assigned); perf <- n_correct/n
      
      if (verbose >= 2)
        cat(sprintf("Correctly assigned: %d / %d (%.1f%%)\n", n_correct, n, 100*perf))
      
      gg <- ggplot(res, aes(x=population, fill=correctly_assigned)) +
        geom_bar(colour="white", linewidth=0.2) +
        scale_fill_manual(
          values=c("TRUE"="steelblue","FALSE"="salmon"),
          labels=c("TRUE"="Correct","FALSE"="Misassigned"),
          name=NULL) +
        labs(title    = "Population assignment",
             subtitle = sprintf("%d / %d (%.1f%%) correctly assigned  (e=%.3f)",
                                n_correct, n, 100*perf, e),
             x="True population", y="Count") +
        theme(axis.text.x  = element_text(angle=45, hjust=1),
              plot.title   = element_text(hjust=0.5),
              plot.subtitle= element_text(hjust=0.5))
    }
    
    # ---------------------------------------------------------
    # hybridisation
    # ---------------------------------------------------------
    if (param == "hybridisation") {
      e     <- error.rate
      trans <- make_trans(e)
      gmat  <- as.matrix(x)
      L     <- ncol(gmat)
      pop_f <- as.factor(pop(x))
      empty_levels <- setdiff(levels(pop_f), as.character(unique(pop_f)))
      if (length(empty_levels) > 0) {
        warning(sprintf(
          "hybridisation: dropping %d population(s) with zero individuals: %s",
          length(empty_levels), paste(empty_levels, collapse=", ")))
        pop_f <- droplevels(pop_f)
      }
      pops  <- as.character(pop_f)
      upops <- levels(pop_f)
      K     <- length(upops)
      
      if (K < 2)
        stop("'hybridisation' requires at least 2 populations.")
      
      pop_idx <- setNames(lapply(upops, function(z) which(pops == z)), upops)
      
      pair_mat <- utils::combn(upops, 2, simplify = TRUE)
      pair_df <- data.frame(
        pop1       = pair_mat[1,],
        pop2       = pair_mat[2,],
        true_cross = paste(pair_mat[1,], pair_mat[2,], sep=" x "),
        f1_class   = paste0("F1: ", pair_mat[1,], " x ", pair_mat[2,]),
        stringsAsFactors = FALSE
      )
      n_pairs <- nrow(pair_df)
      
      # Allele frequencies for pure populations.
      af <- matrix(NA_real_, nrow=K, ncol=L,
                   dimnames=list(upops, locNames(x)))
      for (k in seq_len(K)) {
        idx_k <- pop_idx[[upops[k]]]
        af[k,] <- colMeans(gmat[idx_k,,drop=FALSE], na.rm=TRUE)/2
      }
      af[is.nan(af) | is.na(af)] <- 0.5
      af <- pmax(pmin(af, 1-1e-6), 1e-6)
      
      # Classifier contains all pure classes and all pairwise F1 classes.
      # This is essential for >2 populations because an F1 can be detected
      # as a hybrid but assigned to the wrong parental pair.
      class_labels <- c(upops, pair_df$f1_class)
      class_type   <- c(rep("pure", K), rep("F1", n_pairs))
      class_cross  <- c(rep(NA_character_, K), pair_df$true_cross)
      n_classes    <- length(class_labels)
      
      ll_class <- array(NA_real_, dim=c(L, 3, n_classes),
                        dimnames=list(locNames(x), as.character(0:2),
                                      class_labels))
      
      # Pure-population genotype likelihoods.
      for (k in seq_len(K)) {
        for (l in seq_len(L)) {
          p <- af[k,l]
          hwe <- c((1-p)^2, 2*p*(1-p), p^2)
          ll_class[l,,k] <- as.numeric(hwe %*% trans)
        }
      }
      
      # F1 genotype likelihoods for every pairwise population cross.
      for (pp in seq_len(n_pairs)) {
        i <- match(pair_df$pop1[pp], upops)
        j <- match(pair_df$pop2[pp], upops)
        class_i <- K + pp
        for (l in seq_len(L)) {
          p1 <- af[i,l]
          p2 <- af[j,l]
          f1 <- c((1-p1)*(1-p2),
                  p1*(1-p2) + (1-p1)*p2,
                  p1*p2)
          ll_class[l,,class_i] <- as.numeric(f1 %*% trans)
        }
      }
      
      a1_v <- c(0, 0.5, 1)
      sample_one <- function(idx) idx[sample.int(length(idx), 1)]
      
      sim_offspring <- function(par1, par2) {
        o <- integer(L)
        for (l in seq_len(L)) {
          g1 <- gmat[par1,l]
          g2 <- gmat[par2,l]
          ap1 <- if (is.na(g1)) 0.5 else a1_v[g1+1]
          ap2 <- if (is.na(g2)) 0.5 else a1_v[g2+1]
          ot  <- rbinom(1, 1, ap1) + rbinom(1, 1, ap2)
          o[l] <- sample(0:2, 1, prob=trans[ot+1,])
        }
        o
      }
      
      assign_class <- function(g_obs) {
        log_L <- numeric(n_classes)
        for (l in seq_len(L)) {
          if (is.na(g_obs[l])) next
          log_L <- log_L +
            log(pmax(ll_class[l, g_obs[l]+1, ], .Machine$double.eps))
        }
        which.max(log_L)
      }
      
      if (verbose >= 2)
        cat(sprintf("Simulating %d F1 offspring for each of %d population pairs (%d total F1 simulations)...\n",
                    n_sim_hyb, n_pairs, n_pairs*n_sim_hyb))
      
      records <- vector("list", n_pairs*n_sim_hyb)
      rec_i <- 0L
      
      # Simulate F1 hybrids for every pairwise population cross.
      for (pp in seq_len(n_pairs)) {
        idx_1 <- pop_idx[[pair_df$pop1[pp]]]
        idx_2 <- pop_idx[[pair_df$pop2[pp]]]
        
        for (s in seq_len(n_sim_hyb)) {
          ass_i <- assign_class(sim_offspring(sample_one(idx_1), sample_one(idx_2)))
          rec_i <- rec_i + 1L
          records[[rec_i]] <- data.frame(
            sim              = s,
            pop1             = pair_df$pop1[pp],
            pop2             = pair_df$pop2[pp],
            true_cross       = pair_df$true_cross[pp],
            true_class       = pair_df$f1_class[pp],
            assigned_class   = class_labels[ass_i],
            assigned_type    = class_type[ass_i],
            assigned_cross   = class_cross[ass_i],
            hybrid_detected  = class_type[ass_i] == "F1",
            correct_f1_pair  = class_labels[ass_i] == pair_df$f1_class[pp],
            stringsAsFactors = FALSE
          )
        }
      }
      
      res <- do.call(rbind, records)
      res$true_cross      <- factor(res$true_cross, levels=pair_df$true_cross)
      res$true_class      <- factor(res$true_class, levels=pair_df$f1_class)
      res$assigned_class  <- factor(res$assigned_class,
                                    levels=c(pair_df$f1_class, upops))
      res$assigned_type   <- factor(res$assigned_type, levels=c("F1", "pure"))
      
      # Performance requested here: proportion of all simulated F1s assigned
      # to the correct pairwise F1 class. The broader hybrid-detected rate is
      # kept only as an additional diagnostic.
      f1_detected_r <- mean(res$hybrid_detected)
      correct_pair_r <- mean(res$correct_f1_pair)
      perf <- correct_pair_r
      
      param_summary <- aggregate(
        cbind(hybrid_detected, correct_f1_pair) ~ true_cross,
        data = res,
        FUN = mean
      )
      names(param_summary)[names(param_summary) == "hybrid_detected"] <-
        "f1_detected_rate"
      names(param_summary)[names(param_summary) == "correct_f1_pair"] <-
        "correct_f1_pair_rate"
      param_summary$n_sim_hyb <- n_sim_hyb
      
      confusion <- as.data.frame(
        table(true_cross=res$true_cross,
              assigned_class=res$assigned_class),
        stringsAsFactors=FALSE
      )
      names(confusion)[names(confusion) == "Freq"] <- "n"
      totals <- aggregate(n ~ true_cross, data=confusion, sum)
      names(totals)[2] <- "total"
      confusion <- merge(confusion, totals, by="true_cross", all.x=TRUE)
      confusion$proportion <- confusion$n / confusion$total
      confusion$assigned_type <- ifelse(grepl("^F1: ", confusion$assigned_class),
                                        "F1", "pure")
      confusion$true_cross <- factor(confusion$true_cross,
                                     levels=pair_df$true_cross)
      confusion$assigned_class <- factor(confusion$assigned_class,
                                         levels=c(pair_df$f1_class, upops))
      confusion <- confusion[order(confusion$true_cross,
                                   confusion$assigned_class), ]
      
      if (verbose >= 2)
        cat(sprintf("F1 detected over all pairs: %.1f%%  |  correct F1 pair/performance: %.1f%%\n",
                    100*f1_detected_r, 100*correct_pair_r))
      
      gg <- ggplot(confusion,
                   aes(x=assigned_class, y=true_cross, fill=proportion)) +
        geom_tile(colour="white", linewidth=0.2) +
        scale_fill_gradient(low="white", high="steelblue", limits=c(0,1),
                            name="Proportion") +
        labs(title    = "Hybridisation",
             subtitle = sprintf(paste0("F1 detected: %.1f%%  |  correct F1 pair: %.1f%%  |  ",
                                      "%d pairs, %d F1 per pair"),
                                100*f1_detected_r, 100*correct_pair_r,
                                n_pairs, n_sim_hyb),
             x="Assigned class", y="True F1 cross") +
        theme(axis.text.x  = element_text(angle=45, hjust=1),
              plot.title   = element_text(hjust=0.5),
              plot.subtitle= element_text(hjust=0.5))
      
      if (nrow(confusion) <= 120) {
        gg <- gg + geom_text(aes(label=sprintf("%.2f", proportion)),
                             size=2.4)
      }
    }
    
    if (verbose >= 1) {
      if (param %in% correlation_params)
        cat(sprintf("Performance (%s, %s): %.4f\n",
                    param, corr_metric_label(), perf))
      else
        cat(sprintf("Performance (%s): %.4f\n", param, perf))
    }
    
    output[[param]] <- list(
      name        = param,
      performance = round(perf, 4),
      corr.method = if (param %in% correlation_params) corr.method else NA_character_,
      data        = res,
      summary     = param_summary,
      confusion   = confusion,
      plot        = gg
    )
  }
  
  # ---------------------------------------------------------------
  # combined output page
  # summary bar chart is first panel in the same grid — same size
  # as parameter plots rather than spanning the full width
  # ---------------------------------------------------------------
  if (plot.out) {
    
    perf_df <- data.frame(
      parameter   = factor(names(output), levels=names(output)),
      performance = sapply(output, `[[`, "performance"),
      stringsAsFactors = FALSE
    )
    
    ## relabel drift_resistance when inverted
    if (inverse_dr && "drift_resistance" %in% levels(perf_df$parameter))
      levels(perf_df$parameter)[levels(perf_df$parameter) == "drift_resistance"] <-
        "drift_sensitivity"
    
    summary_gg <- ggplot(perf_df, aes(x=parameter, y=performance,
                                      fill=performance)) +
      geom_col(colour="white", linewidth=0.3) +
      scale_fill_gradient2(low="salmon", mid="goldenrod", high="steelblue",
                           midpoint=0.7, limits=c(0,1), name=NULL) +
      coord_cartesian(ylim=c(0,1)) +
      labs(title    = "Performance summary",
           subtitle = sprintf("Mean: %.3f  |  %d loci  |  %d pops  |  corr: %s",
                              mean(perf_df$performance, na.rm=TRUE),
                              nLoc(x), nPop(x), corr.method),
           x=NULL, y="Performance (0\u20131)") +
      theme(legend.position = "none",
            axis.text.x     = element_text(angle=45, hjust=1),
            plot.title      = element_text(hjust=0.5),
            plot.subtitle   = element_text(hjust=0.5))
    
    if (!is.null(target))
      summary_gg <- summary_gg +
      geom_hline(yintercept=target, linetype="dashed",
                 colour="grey40", linewidth=0.7) +
      annotate("text", x=Inf, y=target,
               label=sprintf(" target=%.2f", target),
               hjust=1, vjust=-0.4, size=3, colour="grey40")
    
    param_plots <- Filter(Negate(is.null), lapply(output, function(r) r$plot))
    
    # summary is first cell in the same grid — same size as all other plots
    all_plots <- c(list(summary_gg), param_plots)
    print(patchwork::wrap_plots(all_plots, ncol=3L))
  }
  
  return(invisible(output))
}