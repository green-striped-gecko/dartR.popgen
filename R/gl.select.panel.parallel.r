#' Run multiple simulated-annealing chains in parallel with synchronized
#' live plotting (burst-sync), returning the single best panel found.
#'
#' Runs \code{n_chains} independent SA optimisations of the same objective
#' as \code{\link{gl.select.panel.combined}}, each on its own core. Chains
#' advance in synchronized bursts: every \code{burst} iterations all chains
#' pause, the master collects their progress, redraws a combined convergence
#' plot showing every chain, then signals all chains to resume. This lets
#' you watch whether the chains are converging to a common performance level
#' (useful when you want confidence the panel has plateaued).
#'
#' Interruptible: pressing Escape stops all workers cleanly and returns the
#' best panel found across all chains up to the last completed burst.
#'
#' @param x A genlight object.
#' @param nl Panel size (number of loci).
#' @param parameter Character vector of parameters to optimise (see
#'   \code{gl.select.panel.combined}).
#' @param n_chains Number of parallel SA chains (recommended 3-5).
#' @param n_cores Cores to use. Default = \code{n_chains}.
#' @param n_iter Total SA iterations per chain.
#' @param burst Iterations each chain runs between synchronization/plot
#'   points. Smaller = smoother live plot but more overhead. Default 10.
#' @param weights,init.method,start_temp,cooling,restart_tol,stop.criterion
#'   As in \code{gl.select.panel.combined}. Each chain is seeded from a
#'   different random panel (init.method applied per chain with a distinct
#'   seed) unless \code{include.loci} forces some loci.
#' @param corr.method,ref,metric,inverse_dr,sample,include.loci,error.rate,threshold,n_sim_parents,n_sim_hyb,neest.path
#'   As in \code{gl.select.panel.combined}, passed to each chain's evaluator.
#' @param plot.out Logical. Draw the combined live plot at each burst.
#' @param verbose 0/1/2.
#'
#' @return A genlight (the single best panel across all chains) with the same
#'   \code{@other$sa_*} metadata slots as \code{gl.select.panel.combined},
#'   plus \code{@other$sa_all_chains} (per-chain best performance and history).
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_hline annotate labs
#'   scale_colour_brewer coord_cartesian theme element_text
#' @importFrom parallel makeCluster stopCluster clusterEvalQ clusterExport
#'   clusterApplyLB

gl.select.panel.parallel <- function(x,
                                     nl             = 10,
                                     parameter      = c("Fst", "He"),
                                     n_chains       = 4,
                                     n_cores        = NULL,
                                     n_iter         = 200,
                                     burst          = 10,
                                     neest.path     = NULL,
                                     corr.method    = c("spearman", "pearson"),
                                     ref            = c("global", "by.pop"),
                                     metric         = c("grm", "kinship"),
                                     inverse_dr     = FALSE,
                                     sample         = 1,
                                     include.loci   = NULL,
                                     weights        = NULL,
                                     init.method    = "random",
                                     start_temp     = 0.1,
                                     cooling        = NULL,
                                     restart_tol    = NULL,
                                     stop.criterion = NULL,
                                     error.rate     = 0.01,
                                     threshold      = 0.001,
                                     n_sim_parents  = 50,
                                     n_sim_hyb      = 100,
                                     plot.out       = TRUE,
                                     verbose        = 1) {

  corr.method <- match.arg(corr.method)
  ref         <- match.arg(ref)
  metric      <- match.arg(metric)

  if (is.null(n_cores)) n_cores <- n_chains
  if (is.null(cooling)) cooling <- (0.001 / start_temp)^(1 / n_iter)
  n_swap_max <- max(2L, round(nl * 0.3))

  # ---- console safety guard (restore sinks on any exit) ----
  .sink_depth_on_entry <- sink.number()
  on.exit({
    while (sink.number() > .sink_depth_on_entry) sink(type = "output")
  }, add = TRUE)

  # ============================================================
  # worker-side chain state + step function
  # This is shipped to each worker. Each worker holds ONE chain in
  # its global env (.chain) and advances it `burst` iterations per call.
  # ============================================================
  .worker_init <- function(chain_id, x, nl, parameter, weights, init.method,
                           include.loci, start_temp, cooling, n_swap_max,
                           n_iter, restart_tol, corr.method, ref, metric,
                           inverse_dr, sample, error.rate, threshold,
                           n_sim_parents, n_sim_hyb, neest.path) {
    set.seed(chain_id * 7919 + 1)   # distinct stream per chain

    w <- if (is.null(weights)) rep(1, length(parameter)) else weights
    names(w) <- parameter
    w_norm <- w / sum(w)
    has_dr <- "drift_resistance" %in% parameter

    eval_panel <- function(loci) {
      panel <- gl.keep.loc(x, loci, verbose = 0)
      tmp <- tempfile(); con <- file(tmp, open = "w")
      sink(con, type = "output")
      on.exit({
        if (sink.number() > 0) sink(type = "output")
        if (isOpen(con)) close(con)
        unlink(tmp)
      }, add = TRUE)
      res <- tryCatch(gl.check.panel(
        x = panel, xorig = x, parameter = parameter,
        corr.method = corr.method, ref = ref, metric = metric,
        inverse_dr = inverse_dr, sample = sample, error.rate = error.rate,
        threshold = threshold, n_sim_parents = n_sim_parents,
        n_sim_hyb = n_sim_hyb, plot.out = FALSE, neest.path = neest.path,
        verbose = 0), error = function(e) NULL)
      if (is.null(res)) return(list(composite = 0, dr = NA_real_))
      perfs <- vapply(res, function(r) {
        p <- r$performance
        if (is.null(p) || length(p) == 0 || is.na(p)) 0 else as.numeric(p)
      }, numeric(1))
      comp <- sum(w_norm[names(perfs)] * perfs[names(w_norm)])
      dr <- if (has_dr && "drift_resistance" %in% names(perfs))
        perfs[["drift_resistance"]] else NA_real_
      list(composite = comp, dr = dr)
    }

    protected <- if (is.null(include.loci)) character(0) else include.loci
    init_panel <- gl.select.panel(x, method = init.method, nl = nl,
                                  include.loci = include.loci,
                                  plot.out = FALSE, verbose = 0)
    curr_loci <- locNames(init_panel)
    ev        <- eval_panel(curr_loci)

    .chain <<- list(
      id = chain_id, x = x, nl = nl, all_loci = locNames(x),
      protected = protected, eval_panel = eval_panel,
      curr_loci = curr_loci, curr_perf = ev$composite,
      best_loci = curr_loci, best_perf = ev$composite,
      T_sa = start_temp, start_temp = start_temp, cooling = cooling,
      n_swap_max = n_swap_max, n_iter = n_iter, restart_tol = restart_tol,
      iter = 0L,
      hist_iter = 0L, hist_best = ev$composite, hist_curr = ev$composite
    )
    invisible(TRUE)
  }

  .worker_step <- function(n_steps) {
    ch <- .chain
    for (s in seq_len(n_steps)) {
      if (ch$iter >= ch$n_iter) break
      ch$iter <- ch$iter + 1L

      # restart-from-best
      if (!is.null(ch$restart_tol) &&
          ch$curr_perf < ch$best_perf - ch$restart_tol) {
        ch$curr_loci <- ch$best_loci
        ch$curr_perf <- ch$best_perf
      }

      n_swap <- max(1L, round(ch$n_swap_max * ch$T_sa / ch$start_temp))
      swappable <- setdiff(ch$curr_loci, ch$protected)
      n_swap_eff <- min(n_swap, length(swappable))
      if (n_swap_eff >= 1L) {
        out_locs <- sample(swappable, n_swap_eff)
        in_locs  <- sample(setdiff(ch$all_loci, ch$curr_loci), n_swap_eff)
        new_loci <- c(setdiff(ch$curr_loci, out_locs), in_locs)
        ev <- ch$eval_panel(new_loci)
        delta <- ev$composite - ch$curr_perf
        if (delta > 0 || runif(1) < exp(delta / ch$T_sa)) {
          ch$curr_loci <- new_loci
          ch$curr_perf <- ev$composite
        }
        if (ch$curr_perf > ch$best_perf) {
          ch$best_perf <- ch$curr_perf
          ch$best_loci <- ch$curr_loci
        }
      }
      ch$T_sa <- ch$T_sa * ch$cooling
      ch$hist_iter <- c(ch$hist_iter, ch$iter)
      ch$hist_best <- c(ch$hist_best, ch$best_perf)
      ch$hist_curr <- c(ch$hist_curr, ch$curr_perf)
    }
    .chain <<- ch
    list(id = ch$id, iter = ch$iter, best_perf = ch$best_perf,
         best_loci = ch$best_loci,
         hist_iter = ch$hist_iter, hist_best = ch$hist_best,
         hist_curr = ch$hist_curr, done = ch$iter >= ch$n_iter)
  }

  # ============================================================
  # master: set up cluster, init chains, run bursts, plot, collect
  # ============================================================
  cat(sprintf("Starting %d parallel SA chains on %d cores (%d iter, burst=%d)\n",
              n_chains, n_cores, n_iter, burst))

  cl <- parallel::makeCluster(min(n_cores, n_chains))
  # guarantee the cluster is torn down no matter how we exit (incl. Escape)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(cl, {
    suppressMessages(library(dartR.base))
    suppressMessages(library(dartR.sim))
  })
  # ship custom function versions if present in the master's global env
  fns <- character(0)
  for (fn in c("gl.check.panel", "gl.select.panel"))
    if (exists(fn, envir = globalenv())) fns <- c(fns, fn)
  if (length(fns)) parallel::clusterExport(cl, fns, envir = globalenv())
  parallel::clusterExport(cl, c(".worker_init", ".worker_step"),
                          envir = environment())

  # init one chain per worker.
  # Share all init args via a single exported list (avoids passing the
  # S4 genlight through clusterApply's positional arg machinery, which
  # fails with 'S4 class is not subsettable').
  .init_args <- list(
    x = x, nl = nl, parameter = parameter, weights = weights,
    init.method = init.method, include.loci = include.loci,
    start_temp = start_temp, cooling = cooling, n_swap_max = n_swap_max,
    n_iter = n_iter, restart_tol = restart_tol, corr.method = corr.method,
    ref = ref, metric = metric, inverse_dr = inverse_dr, sample = sample,
    error.rate = error.rate, threshold = threshold,
    n_sim_parents = n_sim_parents, n_sim_hyb = n_sim_hyb,
    neest.path = neest.path)
  parallel::clusterExport(cl, ".init_args", envir = environment())

  chain_ids <- seq_len(n_chains)
  parallel::clusterApply(cl, chain_ids, function(cid) {
    do.call(.worker_init, c(list(chain_id = cid), .init_args))
  })

  n_bursts <- ceiling(n_iter / burst)
  latest   <- vector("list", n_chains)
  stop_reason <- "completed"

  # ---- burst loop, wrapped so Escape returns best-so-far ----
  tryCatch({
    for (b in seq_len(n_bursts)) {
      steps <- min(burst, n_iter - (b - 1L) * burst)
      # advance every chain one burst. Export the step count as a global
      # so we don't pass extra named args through clusterApply.
      parallel::clusterExport(cl, "steps", envir = environment())
      latest <- parallel::clusterApply(cl, seq_len(n_chains),
                                       function(i) .worker_step(steps))

      # combined progress plot
      if (plot.out) .plot_chains(latest, n_iter, stop.criterion, parameter)

      # console summary
      if (verbose >= 1) {
        bests <- sapply(latest, function(z) z$best_perf)
        cat(sprintf("burst %d/%d | iter %d | chain best: %s | spread=%.4f\n",
                    b, n_bursts, latest[[1]]$iter,
                    paste(sprintf("%.3f", bests), collapse=" "),
                    max(bests) - min(bests)))
      }

      # optional early stop: all chains cleared the criterion
      if (!is.null(stop.criterion) &&
          all(sapply(latest, function(z) z$best_perf) >= stop.criterion)) {
        stop_reason <- "stop_criterion"
        cat(sprintf("All chains reached stop.criterion (%.3f). Stopping.\n",
                    stop.criterion))
        break
      }
    }
  }, interrupt = function(e) {
    stop_reason <<- "interrupted"
    cat("\nInterrupted - returning best panel found so far across all chains.\n")
  })

  # ---- pick global best across chains ----
  if (all(sapply(latest, is.null)))
    stop("No chain produced a result.")
  best_perfs <- sapply(latest, function(z) if (is.null(z)) -Inf else z$best_perf)
  win <- which.max(best_perfs)
  best_loci <- latest[[win]]$best_loci

  cat(sprintf("\nBest chain: #%d | performance = %.4f | reason = %s\n",
              win, best_perfs[win], stop_reason))

  # final combined plot
  if (plot.out) .plot_chains(latest, n_iter, stop.criterion, parameter,
                             final = TRUE, win = win)

  # build the winning panel genlight with metadata
  idx <- which(locNames(x) %in% best_loci)
  panel <- x[, idx]
  panel@other$sa_best_performance <- best_perfs[win]
  panel@other$sa_stop_reason      <- stop_reason
  panel@other$sa_all_chains <- lapply(latest, function(z) {
    if (is.null(z)) return(NULL)
    list(id = z$id, best_perf = z$best_perf,
         history = data.frame(iteration = z$hist_iter,
                              best = z$hist_best, current = z$hist_curr))
  })
  panel@other$sa_parameters <- list(
    parameter = parameter, weights = weights, n_chains = n_chains,
    n_iter = n_iter, corr.method = corr.method, ref = ref, metric = metric,
    inverse_dr = inverse_dr, sample = sample, include.loci = include.loci)

  cat(sprintf("Returning best panel: %d loci | best=%.4f | %s\n",
              nl, best_perfs[win], stop_reason))
  panel
}


# ---- combined chain plot (master side) ----
.plot_chains <- function(latest, n_iter, stop.criterion, parameter,
                         final = FALSE, win = NULL) {
  df <- do.call(rbind, lapply(latest, function(z) {
    if (is.null(z)) return(NULL)
    data.frame(chain = factor(z$id),
               iteration = z$hist_iter,
               best = z$hist_best,
               current = z$hist_curr)
  }))
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))

  ttl <- if (final) "Parallel SA - final (all chains)"
         else "Parallel SA - live (all chains)"
  sub <- sprintf("Params: %s | %d chains%s",
                 paste(parameter, collapse=", "),
                 length(unique(df$chain)),
                 if (!is.null(win)) sprintf(" | winner = chain %d", win) else "")

  gg <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = best,
                                         colour = chain)) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::geom_line(ggplot2::aes(y = current), linewidth = 0.3, alpha = 0.4) +
    ggplot2::scale_colour_brewer(palette = "Set2", name = "Chain") +
    ggplot2::coord_cartesian(xlim = c(0, n_iter), ylim = c(0, 1)) +
    ggplot2::labs(title = ttl, subtitle = sub,
                  x = "Iteration", y = "Composite performance") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5))

  if (!is.null(stop.criterion))
    gg <- gg + ggplot2::geom_hline(yintercept = stop.criterion,
                                   linetype = "dashed", colour = "firebrick")

  print(gg)
  invisible(gg)
}
