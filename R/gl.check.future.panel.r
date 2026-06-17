#' Evaluate SNP panel performance over future generations
#'
#' @description
#' Evaluates how well a SNP panel maintains its performance across multiple
#' generations, either by simulating genetic drift internally (Wright-Fisher
#' model) or by accepting user-supplied genlight objects representing the
#' population at each future generation.
#'
#' @param x A \code{dartR} or \code{genlight} object. The SNP panel (subset
#'   of loci from \code{xorig}).
#' @param xorig A \code{dartR} or \code{genlight} object. The full original
#'   dataset. All loci in \code{x} must be present in \code{xorig}.
#' @param parameter Character vector. Parameters to evaluate. Same options as
#'   \code{\link{gl.check.panel}}: \code{"Fst"}, \code{"He"}, \code{"Ho"},
#'   \code{"Fis"}, \code{"Nall"}, \code{"Ho_ind"}, \code{"relatedness"},
#'   \code{"id"}, \code{"parentage"}, \code{"assignment"},
#'   \code{"hybridisation"}, \code{"all"}, \code{"all-Ne"}.
#'   For \code{"hybridisation"}, all pairwise population crosses are
#'   evaluated when two or more populations are present.
#'   Default \code{c("Fst", "He")}.
#' @param corr.method Character. Correlation method passed to
#'   \code{gl.check.panel} for correlation-based performance metrics.
#'   Options are \code{"spearman"} (default) and \code{"pearson"}.
#' @param neest.path Character string. Path to NEstimator executable,
#'   required only when \code{"Ne"} is included in \code{parameter}.
#' @param type Character. \code{"drift"} (default) simulates genetic drift
#'   internally. \code{"user"} accepts user-supplied genlights via
#'   \code{user.gl}.
#' @param user.gl List of \code{dartR} or \code{genlight} objects. Required
#'   when \code{type = "user"}. Element \code{[[i]]} is the full population
#'   dataset at generation \code{i}. Panel loci are extracted automatically.
#' @param n_gen Integer. Number of generations to simulate
#'   (\code{type = "drift"} only). For \code{type = "user"}, \code{n_gen}
#'   is set to \code{length(user.gl)}. If \code{n_check} requests
#'   generations beyond \code{n_gen}, \code{n_gen} is extended automatically.
#'   Default \code{10}.
#' @param n_sim Integer. Number of replicate drift simulations
#'   (\code{type = "drift"} only). Default \code{5}.
#' @param n_check Integer vector or \code{NULL}. Generations at which to
#'   evaluate panel performance. Generation 0 (baseline) is always included.
#'   \code{NULL} (default) = every generation. For \code{type = "drift"},
#'   \code{n_gen} is extended automatically to cover all values in
#'   \code{n_check}. For \code{type = "user"}, values exceeding
#'   \code{length(user.gl)} are silently dropped.
#' @param error.rate Numeric. Per-allele error rate. Default \code{0.01}.
#' @param threshold Numeric. Confusion threshold for \code{"id"} and
#'   \code{"parentage"}. Default \code{0.001}.
#' @param n_sim_parent Integer. Simulated trios for \code{"parentage"}.
#'   Default \code{50}.
#' @param n_sim_hyb Integer. Simulated F1 individuals per pairwise
#'   population cross for \code{"hybridisation"}. Default \code{100}.
#' @param n_cores Integer or \code{NULL}. Cores for parallel parentage
#'   simulation. Default \code{NULL} = auto.
#' @param target Numeric (0--1) or \code{NULL}. Reference line in the plot.
#'   Default \code{0.9}.
#' @param plot.out Logical. Print the performance-over-time plot.
#'   Default \code{TRUE}.
#' @param plot.file Character string. File name for saving the plot.
#' @param plot.dir Character string. Directory for saving the plot.
#' @param verbose Integer. \code{0} = silent, \code{2} = full (default).
#'
#' @details
#' \strong{Drift model:} Allele frequencies tracked as K x L matrix.
#' Each generation: \eqn{p' = \text{Binomial}(2N, p) / 2N} per population.
#' No genlight objects created during drift — materialised only at check
#' generations as proper \code{"dartR"} objects.
#'
#' \strong{User mode:} User supplies full-dataset genlights. Panel loci
#' extracted from each element and \code{gl.check.panel} run at each check
#' generation. For correlation-based parameters, the chosen
#' \code{corr.method} is used at every generation.
#'
#' @return
#' A list: \code{$performance} (data frame), \code{$summary} (per-generation
#' summaries), \code{$final_panels} (panel genlights at last check
#' generation), \code{$plot}.
#'
#' @examples
#' panel <- gl.select.panel(possums.gl, method = "random", nl = 50)
#'
#' # check at generations 5, 10, 20 — n_gen extended automatically to 20
#' fut <- gl.check.future.panel(panel, possums.gl,
#'                               parameter = c("Fst", "He"),
#'                               n_gen  = 10,
#'                               n_sim  = 5,
#'                               n_check = c(5, 10, 20))
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point
#'   geom_hline annotate labs coord_cartesian theme element_text
#'   scale_colour_brewer scale_fill_brewer

gl.check.future.panel <- function(x,
                                  xorig,
                                  parameter    = c("Fst", "He"),
                                  corr.method  = c("spearman", "pearson"),
                                  neest.path   = NULL,
                                  type         = "drift",
                                  user.gl      = NULL,
                                  n_gen        = 10,
                                  n_sim        = 5,
                                  n_check      = NULL,
                                  error.rate   = 0.01,
                                  threshold    = 0.001,
                                  n_sim_parent = 50,
                                  n_sim_hyb    = 100,
                                  n_cores      = NULL,
                                  target       = 0.9,
                                  plot.out     = TRUE,
                                  plot.file    = NULL,
                                  plot.dir     = NULL,
                                  verbose      = NULL) {
  
  # ---------------------------------------------------------------
  # expand shortcuts
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
  corr.method <- match.arg(tolower(corr.method), c("spearman", "pearson"))
  
  # ---------------------------------------------------------------
  # validate
  # ---------------------------------------------------------------
  valid_params <- c(all_params, "hybridisation")
  bad <- setdiff(parameter, valid_params)
  if (length(bad) > 0)
    stop(paste("Unknown parameter(s):", paste(bad, collapse=", ")))
  if ("Ne" %in% parameter && is.null(neest.path))
    stop("neest.path is required when parameter includes 'Ne'.")
  
  if (!type %in% c("drift","user"))
    stop("type must be 'drift' or 'user'.")
  
  if (type == "user") {
    if (is.null(user.gl))
      stop("user.gl must be supplied when type = 'user'.")
    if (!is.list(user.gl) || length(user.gl) == 0)
      stop("user.gl must be a non-empty list of genlight objects.")
    if (!all(sapply(user.gl, inherits, "genlight")))
      stop("All elements of user.gl must be genlight or dartR objects.")
    missing_check <- sapply(user.gl, function(gl)
      length(setdiff(locNames(x), locNames(gl))))
    if (any(missing_check > 0))
      stop(paste0("Panel loci missing from user.gl element(s): ",
                  paste(which(missing_check > 0), collapse=", ")))
    n_gen <- length(user.gl)
    n_sim <- 1L
  }
  
  if (type == "drift") {
    if (n_gen < 1) stop("n_gen must be >= 1.")
    if (n_sim < 1) stop("n_sim must be >= 1.")
  }
  
  if (!is.null(target) && (target <= 0 || target > 1))
    stop("target must be in (0, 1].")
  
  panel_loci   <- locNames(x)
  missing_loci <- setdiff(panel_loci, locNames(xorig))
  if (length(missing_loci) > 0)
    stop(paste0(length(missing_loci), " panel loci not found in xorig."))
  
  if ("hybridisation" %in% parameter && nPop(x) < 2)
    stop(paste0("'hybridisation' requires at least 2 populations. Found: ",
                nPop(x), "."))
  
  verbose <- if (is.null(verbose)) 2L else as.integer(verbose)
  
  # ---------------------------------------------------------------
  # resolve check generations
  # for drift: n_gen is extended if n_check requests beyond it
  # for user:  values beyond length(user.gl) are silently dropped
  # ---------------------------------------------------------------
  if (is.null(n_check)) {
    check_gens <- seq_len(n_gen)
  } else {
    valid_checks <- as.integer(n_check[n_check >= 1])
    if (length(valid_checks) == 0)
      stop("n_check must contain at least one positive integer.")
    
    if (type == "drift") {
      # extend n_gen to cover all requested check generations
      n_gen      <- max(n_gen, max(valid_checks))
      check_gens <- sort(unique(valid_checks))
    } else {
      # user mode: filter to available genlights
      valid_checks <- valid_checks[valid_checks <= n_gen]
      if (length(valid_checks) == 0)
        stop(sprintf(
          "All n_check values exceed the number of user.gl genlights (%d).",
          n_gen))
      check_gens <- sort(unique(valid_checks))
    }
  }
  last_gen <- max(check_gens)
  
  # ---------------------------------------------------------------
  # precompute fixed quantities
  # ---------------------------------------------------------------
  xorig     <- xorig[order(pop(xorig)),]
  x         <- x[order(pop(x)),]
  pop_names <- levels(pop(xorig))
  K         <- length(pop_names)
  L_orig    <- nLoc(xorig)
  loc_orig  <- locNames(xorig)
  panel_idx <- which(loc_orig %in% panel_loci)
  
  cat(sprintf(
    "gl.check.future.panel [%s]: %d panel loci | %d full loci | %d pops | %d gen | %d sim | corr=%s | check: %s\n",
    type, length(panel_loci), L_orig, K, n_gen, n_sim, corr.method,
    if (length(check_gens) <= 10) paste(check_gens, collapse=", ")
    else sprintf("%d generations", length(check_gens))
  ))
  
  # ---------------------------------------------------------------
  # drift-specific setup
  # ---------------------------------------------------------------
  if (type == "drift") {
    
    n_per_pop <- setNames(
      sapply(pop_names, function(p) sum(as.character(pop(xorig))==p)),
      pop_names
    )
    
    af_init <- matrix(NA_real_, nrow=K, ncol=L_orig,
                      dimnames=list(pop_names, loc_orig))
    for (k in seq_len(K)) {
      idx_k        <- which(as.character(pop(xorig))==pop_names[k])
      af_k         <- colMeans(as.matrix(xorig[idx_k,]), na.rm=TRUE)/2
      af_k[is.nan(af_k)] <- NA_real_
      af_init[k,]  <- af_k
    }
    
    sim_drift <- function(af_mat) {
      for (k in seq_len(K)) {
        valid <- !is.na(af_mat[k,])
        if (any(valid)) {
          new_counts       <- rbinom(sum(valid), size=2L*n_per_pop[k],
                                     prob=af_mat[k,valid])
          af_mat[k,valid]  <- new_counts/(2L*n_per_pop[k])
        }
      }
      af_mat
    }
    
    af_to_genlight <- function(af_mat, loci_idx=NULL) {
      loci_use  <- if (is.null(loci_idx)) seq_len(L_orig) else loci_idx
      L_use     <- length(loci_use)
      n_total   <- sum(n_per_pop)
      full_mat  <- matrix(NA_integer_, nrow=n_total, ncol=L_use)
      pop_vec   <- rep(pop_names, n_per_pop)
      row_start <- 1L
      
      for (k in seq_len(K)) {
        n_k     <- n_per_pop[k]
        alf_k   <- af_mat[k,loci_use]
        valid   <- !is.na(alf_k)
        row_end <- row_start+n_k-1L
        if (any(valid)) {
          p_rep <- rep(alf_k[valid], each=n_k)
          draws <- matrix(rbinom(n_k*sum(valid), size=2L, prob=p_rep), nrow=n_k)
          full_mat[row_start:row_end, which(valid)] <- draws
        }
        row_start <- row_end+1L
      }
      
      ind_names <- paste0(rep(pop_names, n_per_pop), "_",
                          unlist(lapply(n_per_pop, seq_len)))
      template  <- if (is.null(loci_idx)) xorig else xorig[,loci_idx]
      
      gl <- new("dartR", gen=full_mat, fbm=NULL, ploidy=2L,
                ind.names=ind_names, loc.names=locNames(template),
                loc.all=template@loc.all, position=position(template),
                pop=factor(pop_vec, levels=pop_names))
      gl@other <- template@other
      if (!is.null(gl@other$loc.metrics))
        rownames(gl@other$loc.metrics) <- locNames(gl)
      gl@other$ind.metrics <- data.frame(id=ind_names, pop=pop_vec,
                                         stringsAsFactors=FALSE)
      gl
    }
  }
  
  # ---------------------------------------------------------------
  # check_perf: run gl.check.panel silently
  # ---------------------------------------------------------------
  check_perf <- function(x_panel, x_full) {
    result <- tryCatch(
      gl.check.panel(
        x            = x_panel,
        xorig        = x_full,
        parameter    = parameter,
        corr.method  = corr.method,
        error.rate   = error.rate,
        threshold    = threshold,
        n_sim        = n_sim_parent,
        n_sim_hyb    = n_sim_hyb,
        n_cores      = n_cores,
        target       = target,
        neest.path   = neest.path,
        plot.out     = FALSE,
        verbose      = 0
      ),
      error = function(e) {
        message("check_perf error: ", conditionMessage(e)); NULL
      }
    )
    if (is.null(result))
      return(setNames(rep(NA_real_, length(parameter)), parameter))
    vapply(result, function(r) {
      p <- r$performance
      if (is.null(p) || is.na(p)) NA_real_ else as.numeric(p)
    }, numeric(1))
  }
  
  # ---------------------------------------------------------------
  # generation 0 — baseline
  # ---------------------------------------------------------------
  cat("Generation 0: baseline...\n")
  base_perf <- check_perf(x, xorig)
  
  records      <- list()
  final_panels <- vector("list", n_sim)
  names(final_panels) <- paste0("rep_", seq_len(n_sim))
  
  for (p in parameter)
    records[[length(records)+1L]] <- data.frame(
      generation=0L, replicate=0L, parameter=p,
      performance=base_perf[p], stringsAsFactors=FALSE)
  
  # ---------------------------------------------------------------
  # main loop
  # ---------------------------------------------------------------
  if (type == "drift") {
    
    for (rep in seq_len(n_sim)) {
      cat(sprintf("Replicate %d / %d\n", rep, n_sim))
      af_curr <- af_init
      
      for (gen in seq_len(n_gen)) {
        af_curr <- sim_drift(af_curr)
        
        if (gen %in% check_gens) {
          cat(sprintf("\r  gen %3d / %d  [evaluating]   ", gen, n_gen))
          flush.console()
          xorig_sim <- af_to_genlight(af_curr)
          x_sim     <- af_to_genlight(af_curr, loci_idx=panel_idx)
          perf      <- check_perf(x_sim, xorig_sim)
          
          for (p in parameter)
            records[[length(records)+1L]] <- data.frame(
              generation=gen, replicate=rep, parameter=p,
              performance=perf[p], stringsAsFactors=FALSE)
          
          if (gen == last_gen)
            final_panels[[paste0("rep_",rep)]] <- x_sim
          
        } else {
          cat(sprintf("\r  gen %3d / %d               ", gen, n_gen))
          flush.console()
        }
      }
      cat("\n")
    }
    
  } else {  # user
    
    cat("Replicate 1 / 1  (user-supplied genlights)\n")
    
    for (gen in seq_len(n_gen)) {
      if (gen %in% check_gens) {
        cat(sprintf("\r  gen %3d / %d  [evaluating]   ", gen, n_gen))
        flush.console()
        xorig_sim <- user.gl[[gen]][order(pop(user.gl[[gen]])),]
        x_sim     <- gl.keep.loc(xorig_sim, loc.list=panel_loci, verbose=0)
        perf      <- check_perf(x_sim, xorig_sim)
        
        for (p in parameter)
          records[[length(records)+1L]] <- data.frame(
            generation=gen, replicate=1L, parameter=p,
            performance=perf[p], stringsAsFactors=FALSE)
        
        if (gen == last_gen)
          final_panels[["rep_1"]] <- x_sim
        
      } else {
        cat(sprintf("\r  gen %3d / %d               ", gen, n_gen))
        flush.console()
      }
    }
    cat("\n")
  }
  
  # ---------------------------------------------------------------
  # assemble results
  # ---------------------------------------------------------------
  df <- do.call(rbind, records)
  
  # ---------------------------------------------------------------
  # per-generation summary
  # ---------------------------------------------------------------
  all_gens <- c(0L, check_gens)
  
  summary_rows <- lapply(all_gens, function(gen) {
    lapply(parameter, function(p) {
      sub <- df[df$generation==gen & df$parameter==p,]
      pv  <- sub$performance[!is.na(sub$performance)]
      if (length(pv)==0) {
        data.frame(generation=gen, parameter=p, mean=NA, sd=NA,
                   min=NA, max=NA, q25=NA, q75=NA, stringsAsFactors=FALSE)
      } else {
        data.frame(generation=gen, parameter=p,
                   mean = mean(pv),
                   sd   = if (length(pv)>1) sd(pv) else 0,
                   min  = min(pv), max=max(pv),
                   q25  = as.numeric(quantile(pv,0.25)),
                   q75  = as.numeric(quantile(pv,0.75)),
                   stringsAsFactors=FALSE)
      }
    })
  })
  
  summary_df <- do.call(rbind, unlist(summary_rows, recursive=FALSE))
  summary_df <- summary_df[order(summary_df$parameter, summary_df$generation),]
  
  # ---------------------------------------------------------------
  # plot
  # ---------------------------------------------------------------
  gg <- NULL
  
  if (plot.out) {
    
    sub_title <- if (type=="drift")
      sprintf("%d panel loci  |  %d populations  |  %d replicates  |  corr=%s  |  ribbons: IQR (dark), range (light)",
              length(panel_loci), K, n_sim, corr.method)
    else
      sprintf("%d panel loci  |  %d populations  |  corr=%s  |  user-supplied simulation",
              length(panel_loci), K, corr.method)
    
    gg <- ggplot(summary_df, aes(x=generation, colour=parameter,
                                 fill=parameter)) +
      geom_ribbon(aes(ymin=min, ymax=max), alpha=0.10, colour=NA) +
      geom_ribbon(aes(ymin=q25, ymax=q75), alpha=0.25, colour=NA) +
      geom_line( aes(y=mean), linewidth=1.0) +
      geom_point(aes(y=mean), size=2.0) +
      coord_cartesian(ylim=c(0,1)) +
      scale_colour_brewer(palette="Set2", name="Parameter") +
      scale_fill_brewer(  palette="Set2", name="Parameter") +
      labs(title    = sprintf("Panel performance over %d generations [%s]",
                              n_gen, type),
           subtitle = sub_title,
           x="Generation", y="Performance (0\u20131)") +
      theme(plot.title    = element_text(hjust=0.5),
            plot.subtitle = element_text(hjust=0.5))
    
    if (!is.null(target))
      gg <- gg +
      geom_hline(yintercept=target, linetype="dashed",
                 colour="grey30", linewidth=0.7) +
      annotate("text", x=Inf, y=target,
               label=sprintf(" target=%.2f", target),
               hjust=1, vjust=-0.4, size=3.2, colour="grey30")
    
    print(gg)
  }
  
  cat(sprintf(
    "\nDone. %d records | %d parameters | %d check generations | %d replicate(s).\n",
    nrow(df), length(parameter), length(check_gens), n_sim))
  
  return(invisible(list(
    performance  = df,
    summary      = summary_df,
    corr.method  = corr.method,
    neest.path   = neest.path,
    final_panels = final_panels,
    plot         = gg
  )))
}