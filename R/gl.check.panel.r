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
#'   original dataset for comparison. Not required for \code{"id"},
#'   \code{"parentage"}, \code{"assignment"}, or \code{"hybridisation"}.
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
#'     \item \code{"Ho_ind"} — individual observed heterozygosity. Spearman
#'       r² of per-individual Ho between panel and full dataset, ignoring
#'       population membership.
#'     \item \code{"relatedness"} — pairwise genomic relatedness (VanRaden
#'       GRM). Both GRMs computed using full-dataset allele frequencies as
#'       reference. Population membership ignored.
#'     \item \code{"id"} — individual identification power.
#'     \item \code{"parentage"} — parentage assignment power.
#'     \item \code{"assignment"} — population assignment accuracy (LOO).
#'     \item \code{"hybridisation"} — F1 hybrid detection (requires exactly
#'       2 populations in \code{x}).
#'     \item \code{"all"} — all parameters except \code{"hybridisation"}.
#'     \item \code{"all-Ne"} — all except \code{"Ne"} and
#'       \code{"hybridisation"}.
#'   }
#' @param neest.path Character string. Path to NEstimator executable, required
#'   only for \code{"Ne"}.
#' @param error.rate Numeric. Per-allele genotyping error rate. Default
#'   \code{0.01}.
#' @param threshold Numeric. Confusion threshold for \code{"id"} and
#'   \code{"parentage"}. Default \code{0.001}.
#' @param n_sim Integer. Simulated trios for \code{"parentage"}.
#'   Default \code{100}.
#' @param n_sim_hyb Integer. Simulated individuals per class for
#'   \code{"hybridisation"}. Default \code{100}.
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
#'     — Spearman r² between panel and full dataset values. Plots also show
#'     Pearson R².
#'   \item \code{"Ho_ind"} — Spearman r² of per-individual observed
#'     heterozygosity. Population membership ignored.
#'   \item \code{"relatedness"} — Spearman r² of all pairwise GRM values
#'     (upper triangle). Both GRMs use full-dataset allele frequencies as
#'     reference (VanRaden 2008). Population membership ignored.
#'   \item \code{"id"} — fraction of individuals with MaxPmatch below
#'     \code{threshold}.
#'   \item \code{"parentage"} — fraction of simulated trios correctly
#'     assigned.
#'   \item \code{"assignment"} — fraction correctly assigned by LOO.
#'   \item \code{"hybridisation"} — F1 detection rate.
#' }
#'
#' \strong{Heterozygosity caching:} When any of \code{"He"}, \code{"Ho"},
#' \code{"Fis"} are requested, \code{gl.report.heterozygosity} is called
#' once per dataset and reused for all three.
#'
#' \strong{Relatedness reference frequencies:} Both panel and full-dataset
#' GRMs are computed using allele frequencies from \code{xorig}, ensuring
#' both are on the same scale (VanRaden 2008).
#'
#' @return
#' A named list with one element per parameter, each containing:
#' \code{name}, \code{performance}, \code{data}, \code{plot}.
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
#' @importFrom ggplot2 ggplot aes geom_col geom_histogram geom_bar geom_vline
#'   geom_point geom_smooth annotate labs scale_x_log10 scale_fill_manual
#'   scale_fill_gradient2 coord_cartesian theme element_text
#' @importFrom patchwork wrap_plots

gl.check.panel <- function(x,
                           xorig,
                           parameter   = "Fst",
                           neest.path  = NULL,
                           error.rate  = 0.01,
                           threshold   = 0.001,
                           n_sim       = 100,
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
                  "Ho_ind","relatedness",
                  "id","parentage","assignment")
  
  if (length(parameter) == 1 && parameter == "all")
    parameter <- all_params
  if (length(parameter) == 1 && parameter == "all-Ne")
    parameter <- setdiff(all_params, "Ne")
  
  parameter[parameter == "Ar"] <- "Nall"
  parameter <- unique(parameter)
  
  valid_params <- c(all_params, "hybridisation")
  bad <- setdiff(parameter, valid_params)
  if (length(bad) > 0)
    stop(paste("Unknown parameter(s):", paste(bad, collapse=", "),
               "\nValid:", paste(valid_params, collapse=", ")))
  
  if (error.rate  < 0 || error.rate  > 1) stop("error.rate must be in [0,1].")
  if (threshold   <= 0 || threshold  >= 1) stop("threshold must be in (0,1).")
  if (!is.numeric(n_sim)     || n_sim     < 1) stop("n_sim must be a positive integer.")
  if (!is.numeric(n_sim_hyb) || n_sim_hyb < 1) stop("n_sim_hyb must be a positive integer.")
  if (!is.null(target) && (target <= 0 || target > 1)) stop("target must be in (0,1].")
  
  if ("hybridisation" %in% parameter && nPop(x) != 2)
    stop(paste0("'hybridisation' requires exactly 2 populations. Found: ",
                nPop(x), "."))
  
  n_cores <- if (is.null(n_cores)) max(1L, parallel::detectCores()-1L) else
    max(1L, as.integer(n_cores))
  
  # resolve verbosity: 0=silent, 1=headers+scores, 2=full (default)
  verbose <- if (is.null(verbose)) 2L else as.integer(verbose)
  
  # ---------------------------------------------------------------
  # order populations
  # ---------------------------------------------------------------
  x     <- x[order(pop(x)),]
  xorig <- xorig[order(pop(xorig)),]
  
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
  
  r2_spearman <- function(a, b)
    cor(a, b, use="complete.obs", method="spearman")^2
  
  r2_pearson <- function(a, b)
    cor(a, b, use="complete.obs", method="pearson")^2
  
  r2_label <- function(sp, pe)
    sprintf("Spearman r\u00b2 = %.3f\nPearson R\u00b2 = %.3f", sp, pe)
  
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
    
    gg   <- NULL
    res  <- NULL
    perf <- NA_real_
    
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
      perf <- r2_spearman(res$fst_orig, res$fst_panel)
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
      perf <- r2_spearman(res$he_orig, res$he_panel)
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
      perf <- r2_spearman(res$ho_orig, res$ho_panel)
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
      perf <- r2_spearman(res$fis_orig, res$fis_panel)
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
      perf <- r2_spearman(res$ar_orig, res$ar_panel)
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
                          verbose=0, mating="random")
      ne_panel <- gl.LDNe(x,     neest.path=neest.path, critical=0.05,
                          verbose=0, mating="random", singleton.rm=FALSE)
      extract_ne <- function(nl)
        as.numeric(unlist(lapply(nl, function(z) z$`Frequency 1`[6])))
      res <- data.frame(ne_orig=extract_ne(ne_orig), ne_panel=extract_ne(ne_panel))
      res[res==Inf|res==-Inf] <- NA
      res  <- res[complete.cases(res),]
      perf <- r2_spearman(res$ne_orig, res$ne_panel)
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
      perf <- r2_spearman(res$ho_ind_orig, res$ho_ind_panel)
      if (verbose >= 2)
        cat(sprintf("Individuals: %d  |  Spearman r2 = %.4f\n", nrow(res), perf))
      gg   <- scatter_plot(res, "ho_ind_orig", "ho_ind_panel",
                           "Individual heterozygosity",
                           "Ho_ind (original)", "Ho_ind (panel)")
    }
    
    # ---------------------------------------------------------
    # relatedness — pairwise VanRaden GRM
    # both GRMs use full-dataset allele frequencies
    # population membership ignored
    # ---------------------------------------------------------
    if (param == "relatedness") {
      if (verbose >= 2) cat("Computing reference allele frequencies from full dataset...\n")
      p_ref <- colMeans(as.matrix(xorig), na.rm=TRUE)/2
      p_ref[is.nan(p_ref)] <- 0.5
      
      if (verbose >= 2)
        cat(sprintf("Computing GRM: full dataset (%d ind, %d loci)...\n",
                    nInd(xorig), nLoc(xorig)))
      grm_orig  <- grm_vanraden(as.matrix(xorig), p_ref)
      
      if (verbose >= 2)
        cat(sprintf("Computing GRM: panel (%d ind, %d loci)...\n",
                    nInd(x), nLoc(x)))
      grm_panel <- grm_vanraden(as.matrix(x), p_ref[locNames(x)])
      
      pairs_orig  <- grm_orig[ upper.tri(grm_orig)]
      pairs_panel <- grm_panel[upper.tri(grm_panel)]
      
      res  <- data.frame(grm_orig=pairs_orig, grm_panel=pairs_panel)
      res  <- res[complete.cases(res),]
      perf <- r2_spearman(res$grm_orig, res$grm_panel)
      
      if (verbose >= 2)
        cat(sprintf("Pairs: %d  |  Spearman r2 = %.4f\n", nrow(res), perf))
      
      plot_res <- if (nrow(res)>10000L) res[sample(nrow(res),10000L),] else res
      sp <- r2_spearman(res$grm_orig, res$grm_panel)
      pe <- r2_pearson( res$grm_orig, res$grm_panel)
      
      gg <- ggplot(plot_res, aes(x=grm_orig, y=grm_panel)) +
        geom_point(alpha=0.3, size=0.8) +
        geom_smooth(method="lm", formula=y~x, se=TRUE, colour="steelblue") +
        annotate("text", x=-Inf, y=Inf,
                 label=r2_label(sp, pe),
                 hjust=-0.1, vjust=1.3, size=3.5) +
        labs(title    = "Relatedness (GRM)",
             subtitle = sprintf("%d pairs  |  ref freq from full dataset%s",
                                nrow(res),
                                if (nrow(res)>10000L) "  |  10,000 pairs plotted" else ""),
             x="GRM (original)", y="GRM (panel)") +
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
        if (n_cores>1) message("Parallel parentage disabled on Windows. Sequential.")
        1L
      } else min(n_cores, n_sim)
      
      if (verbose >= 2)
        cat(sprintf("Simulating %d trios  (%s)...\n", n_sim,
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
        rp <- parallel::mclapply(seq_len(n_sim), run_one_sim, mc.cores=cores_use)
        fl <- sapply(rp, inherits, "try-error")
        if (any(fl)) { warning(sum(fl)," sim(s) failed; rerunning."); rp[fl] <- lapply(which(fl),run_one_sim) }
        rp
      } else lapply(seq_len(n_sim), run_one_sim)
      
      res <- data.frame(
        sim       = seq_len(n_sim),
        p1        = sapply(sim_results,`[[`,"p1"),
        p2        = sapply(sim_results,`[[`,"p2"),
        confusion = sapply(sim_results,`[[`,"confusion"),
        metric    = sapply(sim_results,`[[`,"metric")
      )
      n_assigned <- sum(res$confusion<threshold)
      perf       <- n_assigned/n_sim
      
      if (verbose >= 2)
        cat(sprintf("Assigned: %d / %d (%.1f%%)  (confusion < %g)\n",
                    n_assigned, n_sim, 100*perf, threshold))
      
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
                                n_assigned, n_sim, 100*perf, threshold, e),
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
      pops  <- as.character(pop(x))
      upops <- levels(pop(x)); K <- length(upops)
      
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
      pops  <- as.character(pop(x)); upops <- levels(pop(x))
      idx1  <- which(pops==upops[1]); idx2 <- which(pops==upops[2])
      
      p1 <- pmax(pmin(colMeans(gmat[idx1,,drop=FALSE],na.rm=TRUE)/2,1-1e-6),1e-6)
      p2 <- pmax(pmin(colMeans(gmat[idx2,,drop=FALSE],na.rm=TRUE)/2,1-1e-6),1e-6)
      
      ll_class <- array(NA, dim=c(L,3,3))
      for (l in seq_len(L)) {
        hwe1 <- c((1-p1[l])^2,2*p1[l]*(1-p1[l]),p1[l]^2)
        hwe2 <- c((1-p2[l])^2,2*p2[l]*(1-p2[l]),p2[l]^2)
        f1   <- c((1-p1[l])*(1-p2[l]),p1[l]*(1-p2[l])+(1-p1[l])*p2[l],p1[l]*p2[l])
        ll_class[l,,1] <- as.numeric(hwe1 %*% trans)
        ll_class[l,,2] <- as.numeric(hwe2 %*% trans)
        ll_class[l,,3] <- as.numeric(f1   %*% trans)
      }
      
      a1_v <- c(0,0.5,1)
      
      sim_offspring <- function(par1, par2) {
        o <- integer(L)
        for (l in seq_len(L)) {
          g1 <- gmat[par1,l]; g2 <- gmat[par2,l]
          ap1 <- if(is.na(g1)) 0.5 else a1_v[g1+1]
          ap2 <- if(is.na(g2)) 0.5 else a1_v[g2+1]
          ot  <- rbinom(1,1,ap1)+rbinom(1,1,ap2)
          o[l] <- sample(0:2,1,prob=trans[ot+1,])
        }
        o
      }
      
      assign_class <- function(g_obs) {
        log_L <- numeric(3)
        for (l in seq_len(L)) {
          if (is.na(g_obs[l])) next
          log_L <- log_L+log(pmax(ll_class[l,g_obs[l]+1,],.Machine$double.eps))
        }
        which.max(log_L)
      }
      
      class_labels <- c(upops[1],upops[2],"F1")
      
      if (verbose >= 2)
        cat(sprintf("Simulating %d F1, %d pure %s, %d pure %s...\n",
                    n_sim_hyb, n_sim_hyb, upops[1], n_sim_hyb, upops[2]))
      
      true_class     <- rep(c("F1",upops[1],upops[2]), each=n_sim_hyb)
      assigned_class <- character(3*n_sim_hyb)
      
      for (s in seq_len(n_sim_hyb))
        assigned_class[s] <-
        class_labels[assign_class(sim_offspring(sample(idx1,1),sample(idx2,1)))]
      for (s in seq_len(n_sim_hyb))
        assigned_class[n_sim_hyb+s] <-
        class_labels[assign_class(sim_offspring(sample(idx1,1),sample(idx1,1)))]
      for (s in seq_len(n_sim_hyb))
        assigned_class[2*n_sim_hyb+s] <-
        class_labels[assign_class(sim_offspring(sample(idx2,1),sample(idx2,1)))]
      
      res <- data.frame(
        true_class=true_class, assigned_class=assigned_class,
        correct=true_class==assigned_class, stringsAsFactors=FALSE
      )
      hyb_r   <- mean(res$correct[res$true_class=="F1"])
      pure1_r <- mean(res$correct[res$true_class==upops[1]])
      pure2_r <- mean(res$correct[res$true_class==upops[2]])
      perf    <- hyb_r
      
      if (verbose >= 2)
        cat(sprintf("F1: %.1f%%  |  pure %s: %.1f%%  |  pure %s: %.1f%%\n",
                    100*hyb_r, upops[1], 100*pure1_r, upops[2], 100*pure2_r))
      
      plot_df <- as.data.frame(
        table(true_class=res$true_class, assigned_class=res$assigned_class),
        stringsAsFactors=FALSE
      )
      plot_df$true_class <- factor(plot_df$true_class,
                                   levels=c(upops[1],upops[2],"F1"))
      fill_cols <- setNames(c("steelblue","salmon","goldenrod"),
                            c(upops[1],upops[2],"F1"))
      
      gg <- ggplot(plot_df, aes(x=true_class,y=Freq,fill=assigned_class)) +
        geom_bar(stat="identity", position="fill",
                 colour="white", linewidth=0.2) +
        scale_fill_manual(values=fill_cols, name="Assigned to") +
        labs(title    = "Hybridisation",
             subtitle = sprintf("F1: %.1f%%  |  pure %s: %.1f%%  |  pure %s: %.1f%%",
                                100*hyb_r, upops[1], 100*pure1_r, upops[2], 100*pure2_r),
             x="True class", y="Proportion assigned") +
        centre_theme()
    }
    
    if (verbose >= 1) cat(sprintf("Performance (%s): %.4f\n", param, perf))
    
    output[[param]] <- list(
      name        = param,
      performance = round(perf, 4),
      data        = res,
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
    
    summary_gg <- ggplot(perf_df, aes(x=parameter, y=performance,
                                      fill=performance)) +
      geom_col(colour="white", linewidth=0.3) +
      scale_fill_gradient2(low="salmon", mid="goldenrod", high="steelblue",
                           midpoint=0.7, limits=c(0,1), name=NULL) +
      coord_cartesian(ylim=c(0,1)) +
      labs(title    = "Performance summary",
           subtitle = sprintf("Mean: %.3f  |  %d loci  |  %d pops",
                              mean(perf_df$performance, na.rm=TRUE),
                              nLoc(x), nPop(x)),
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