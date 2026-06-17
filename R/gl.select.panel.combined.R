#' Select a SNP panel optimised for specified conservation monitoring parameters
#'
#' @description
#' Selects a panel of \code{nl} loci optimised simultaneously across a
#' user-specified set of conservation monitoring parameters using simulated
#' annealing (SA) with adaptive swapping and optional restart-from-best.
#' Press Escape to stop — the best panel found is returned cleanly.
#'
#' @param x A \code{dartR} or \code{genlight} object. Full dataset from which
#'   the panel is selected and against which performance is evaluated. If
#'   \code{"hybridisation"} is in \code{parameter}, \code{x} must have
#'   exactly 2 populations.
#' @param nl Integer. Number of loci to select. Default \code{10}.
#' @param parameter Character vector. Parameters to optimise. Subset of
#'   \code{"Fst"}, \code{"He"}, \code{"Ho"}, \code{"Fis"}, \code{"Nall"},
#'   \code{"id"}, \code{"parentage"}, \code{"assignment"},
#'   \code{"hybridisation"}. \code{"Ne"} not supported. Default
#'   \code{c("Fst", "He")}.
#' @param weights Numeric vector, one per element of \code{parameter}.
#'   \code{NULL} (default) = equal weights.
#' @param init.method Initialisation method from \code{\link{gl.select.panel}}.
#'   Default \code{"random"}.
#' @param n_iter Integer. Maximum SA iterations. Default \code{100}.
#' @param start_temp Numeric. Starting temperature. Default \code{0.1}.
#' @param cooling Numeric (0,1) or \code{NULL}. Geometric cooling rate per
#'   iteration. If \code{NULL} (default), computed automatically:
#'   \deqn{\text{cooling} = \left(\frac{0.001}{\text{start\_temp}}\right)^{1/\text{n\_iter}}}
#'   Set explicitly to override.
#' @param restart_tol Numeric or \code{NULL}. Restart-from-best tolerance.
#'   If the current composite performance falls more than \code{restart_tol}
#'   below the best performance seen so far, the current panel is reset to
#'   the best panel before the next iteration. This prevents the search from
#'   wasting iterations in low-performance regions. \code{NULL} (default)
#'   disables restarting. A value of \code{0.1} is a reasonable starting
#'   point.
#' @param stop.criterion Numeric (0,1) or \code{NULL}. Auto-stop when
#'   composite performance reaches this value. Default \code{NULL}.
#' @param error.rate Numeric. Per-allele error rate for \code{"id"},
#'   \code{"parentage"}, \code{"assignment"}, \code{"hybridisation"}.
#'   Default \code{0.01}.
#' @param threshold Numeric. Confusion threshold for \code{"id"} and
#'   \code{"parentage"}. Default \code{0.001}.
#' @param n_sim_parent Integer. Number of simulated trios for
#'   \code{"parentage"}. Default \code{50}.
#' @param n_sim_hyb Integer. Simulated individuals per class for
#'   \code{"hybridisation"}. Default \code{100}.
#' @param plot.out Logical. Live convergence plot refreshed every
#'   \code{plot.every} iterations. Default \code{TRUE}.
#' @param plot.every Integer. Plot refresh interval in iterations.
#'   Default \code{10}.
#' @param plot.file Character string. File name for saving the final plot.
#' @param plot.dir Character string. Directory for saving the final plot.
#' @param verbose Integer controlling verbosity.
#'
#' @details
#' \strong{Adaptive swapping:}
#' The number of loci swapped per iteration scales linearly with temperature:
#' \deqn{n\_swap(t) = \max(1,\; \text{round}(\lfloor 0.3 \cdot nl \rceil \cdot T(t) / \text{start\_temp}))}
#' Up to 30\% of the panel is replaced early (broad exploration); by the end
#' only 1 locus is swapped (fine-tuning). Following Ingber (1989).
#'
#' \strong{Restart from best:}
#' When \code{restart_tol} is set, after each accept/reject step, if
#' \code{curr\_perf < best\_perf - restart\_tol} the current panel is reset
#' to the best panel found so far. This anchors exploration near the best
#' solution without disabling the temperature-controlled acceptance of new
#' proposals. Restarts are counted and reported.
#'
#' \strong{Auto-cooling:} When \code{cooling = NULL} the rate is chosen
#' automatically to match \code{n_iter}.
#'
#' \strong{Composite performance:}
#' \deqn{\text{composite} = \frac{\sum_p w_p \cdot \text{performance}_p}{\sum_p w_p}}
#'
#' @return
#' A \code{dartR} or \code{genlight} object with \code{nl} selected loci.
#' Stored in \code{x@other}:
#' \itemize{
#'   \item \code{sa_best_performance}
#'   \item \code{sa_history} — data frame: \code{iteration},
#'     \code{current_performance}, \code{best_performance}, \code{n_swap},
#'     \code{restarted}
#'   \item \code{sa_parameters} — list: \code{parameter}, \code{weights},
#'     \code{cooling}, \code{n_swap_max}, \code{restart_tol}
#'   \item \code{sa_stop_reason} — \code{"completed"},
#'     \code{"stop_criterion"}, or \code{"interrupted"}
#' }
#'
#' @references
#' Ingber, L. (1989). Very fast simulated re-annealing.
#' \emph{Mathematical and Computer Modelling}, 12(8), 967--973.
#'
#' Kirkpatrick, S., Gelatt, C. D., & Vecchi, M. P. (1983). Optimization by
#' simulated annealing. \emph{Science}, 220(4598), 671--680.
#'
#' @examples
#' # Basic run
#' panel <- gl.select.panel.combined(possums.gl, nl = 50,
#'                                    parameter = c("Fst", "He"))
#'
#' # With restart-from-best
#' panel <- gl.select.panel.combined(possums.gl, nl = 50,
#'                                    parameter    = c("Fst", "He", "Fis"),
#'                                    weights      = c(2, 1, 1),
#'                                    n_iter       = 200,
#'                                    restart_tol  = 0.1)
#'
#' panel@other$sa_best_performance
#' panel@other$sa_stop_reason
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_hline annotate labs
#'   scale_colour_manual scale_linewidth_manual coord_cartesian

gl.select.panel.combined <- function(x,
                                     nl             = 10,
                                     parameter      = c("Fst", "He"),
                                     weights        = NULL,
                                     init.method    = "random",
                                     n_iter         = 100,
                                     start_temp     = 0.1,
                                     cooling        = NULL,
                                     restart_tol    = NULL,
                                     stop.criterion = NULL,
                                     error.rate     = 0.01,
                                     threshold      = 0.001,
                                     n_sim_parent   = 50,
                                     n_sim_hyb      = 100,
                                     plot.out       = TRUE,
                                     plot.every     = 10,
                                     plot.file      = NULL,
                                     plot.dir       = NULL,
                                     verbose        = NULL) {
  
  # ---------------------------------------------------------------
  # input validation
  # ---------------------------------------------------------------
  valid_params <- c("Fst","He","Ho","Fis","Nall",
                    "Ho_ind","relatedness",
                    "id","parentage","assignment","hybridisation","Ne")
  
  #if ("Ne" %in% parameter)
  #  stop("'Ne' is not supported (too slow for repeated evaluation).")
  bad <- setdiff(parameter, valid_params)
  if (length(bad) > 0)
    stop(paste("Unknown parameter(s):", paste(bad, collapse=", "),
               "\nValid options:", paste(valid_params, collapse=", ")))
  if (length(parameter) == 0)  stop("At least one parameter must be specified.")
  if (nl > nLoc(x))            stop("nl exceeds the number of loci in x.")
  if (nl < 1)                  stop("nl must be at least 1.")
  if (n_iter < 1)              stop("n_iter must be at least 1.")
  if (start_temp <= 0)         stop("start_temp must be positive.")
  if (!is.null(restart_tol) && (restart_tol <= 0 || restart_tol >= 1))
    stop("restart_tol must be in (0, 1).")
  if (!is.null(stop.criterion) &&
      (stop.criterion <= 0 || stop.criterion > 1))
    stop("stop.criterion must be in (0, 1].")
  if ("hybridisation" %in% parameter && nPop(x) != 2)
    stop(paste0("'hybridisation' requires exactly 2 populations. ",
                "Found: ", nPop(x), ". Subset using gl.keep.pop()."))
  
  # auto-cooling: decay start_temp to ~0.001 over n_iter iterations
  if (is.null(cooling)) {
    cooling <- (0.001 / start_temp)^(1 / n_iter)
    cat(sprintf("Auto cooling rate: %.4f  (T: %.3f -> ~0.001 over %d iter)\n",
                cooling, start_temp, n_iter))
  }
  if (cooling <= 0 || cooling >= 1)
    stop("cooling must be between 0 and 1.")
  
  # adaptive swap: max swaps = 30% of panel
  n_swap_max <- max(2L, round(nl * 0.3))
  cat(sprintf("Adaptive swapping: %d loci/iter at start -> 1 at end\n",
              n_swap_max))
  
  if (!is.null(restart_tol))
    cat(sprintf("Restart from best when curr < best - %.3f\n", restart_tol))
  
  if (is.null(weights)) weights <- rep(1, length(parameter))
  if (length(weights) != length(parameter))
    stop(paste0("weights must be the same length as parameter (",
                length(parameter), ")."))
  if (any(weights < 0))  stop("All weights must be non-negative.")
  if (sum(weights) == 0) stop("At least one weight must be positive.")
  
  names(weights) <- parameter
  w_norm         <- weights / sum(weights)
  x              <- x[order(pop(x)),]
  
  # ---------------------------------------------------------------
  # evaluate_panel — all output suppressed via sink + finally
  # ---------------------------------------------------------------
  evaluate_panel <- function(loci) {
    panel  <- gl.keep.loc(x, loci, verbose = 0)
    tmp    <- tempfile()
    con    <- file(tmp, open = "w")
    result <- NULL
    sink(con, type = "output")
    tryCatch({
      result <- gl.check.panel(
        x          = panel,
        xorig      = x,
        parameter  = parameter,
        error.rate = error.rate,
        threshold  = threshold,
        n_sim      = n_sim_parent,
        n_sim_hyb  = n_sim_hyb,
        plot.out   = FALSE,
        neest.path = "d:/programms/NEestimator/",
        verbose    = 0
      )
    }, error = function(e) { },
    finally = {
      sink(type = "output")
      close(con)
      unlink(tmp)
    })
    if (is.null(result)) return(0)
    perfs <- vapply(result, function(r) {
      p <- r$performance
      if (is.null(p) || length(p) == 0 || is.na(p)) 0 else as.numeric(p)
    }, numeric(1))
    sum(w_norm[names(perfs)] * perfs[names(w_norm)])
  }
  
  # ---------------------------------------------------------------
  # convergence plot builder
  # ---------------------------------------------------------------
  build_conv_plot <- function(history, iter_done, best_perf, reason = "") {
    curr_hist <- history[seq_len(iter_done + 1), ]
    curr_hist <- curr_hist[!is.na(curr_hist$current_performance), ]
    if (nrow(curr_hist) == 0) return(invisible(NULL))
    
    plot_df <- rbind(
      data.frame(iteration   = curr_hist$iteration,
                 performance = curr_hist$current_performance,
                 type        = "Current"),
      data.frame(iteration   = curr_hist$iteration,
                 performance = curr_hist$best_performance,
                 type        = "Best")
    )
    
    ttl <- if (nchar(reason) > 0)
      sprintf("SA — %d loci  [%s | iter %d | best=%.4f]",
              nl, reason, iter_done, best_perf)
    else
      sprintf("SA — %d loci  [iter %d/%d | best=%.4f]",
              nl, iter_done, n_iter, best_perf)
    
    sub <- sprintf(
      "Params: %s  |  Weights: %s  |  cooling=%.4f  |  swap: %d->1%s",
      paste(parameter, collapse=", "),
      paste(names(w_norm), round(w_norm, 2), sep="=", collapse=", "),
      cooling, n_swap_max,
      if (!is.null(restart_tol))
        sprintf("  |  restart_tol=%.2f", restart_tol) else ""
    )
    
    gg <- ggplot(plot_df, aes(x = iteration, y = performance,
                              colour = type, linewidth = type)) +
      geom_line() +
      scale_colour_manual(
        values = c("Current" = "grey60", "Best" = "steelblue"),
        name   = NULL) +
      scale_linewidth_manual(
        values = c("Current" = 0.5, "Best" = 1.2),
        name   = NULL) +
      coord_cartesian(ylim = c(0, 1), xlim = c(0, n_iter))
    
    if (!is.null(stop.criterion))
      gg <- gg +
      geom_hline(yintercept = stop.criterion, linetype = "dashed",
                 colour = "firebrick", linewidth = 0.7) +
      annotate("text", x = 0, y = stop.criterion,
               label = sprintf(" stop=%.2f", stop.criterion),
               hjust = 0, vjust = -0.4, size = 3, colour = "firebrick")
    
    gg + labs(title = ttl, subtitle = sub,
              x = "Iteration", y = "Composite performance")
  }
  
  # ---------------------------------------------------------------
  # build_return — full result with SA metadata
  # ---------------------------------------------------------------
  build_return <- function(loci, perf, hist, iter_d, reason) {
    xx <- gl.keep.loc(x, loci, verbose = 0)
    xx@other$sa_best_performance <- perf
    xx@other$sa_history          <- hist[seq_len(iter_d + 1), ]
    xx@other$sa_parameters       <- list(parameter   = parameter,
                                         weights     = w_norm,
                                         cooling     = cooling,
                                         n_swap_max  = n_swap_max,
                                         restart_tol = restart_tol)
    xx@other$sa_stop_reason      <- reason
    cat(sprintf("Returning panel: %d loci | best=%.4f | %s\n",
                nl, perf, reason))
    xx
  }
  
  # ---------------------------------------------------------------
  # minimal_return — fallback inside interrupt handler
  # ---------------------------------------------------------------
  minimal_return <- function(loci, perf, reason) {
    idx <- which(locNames(x) %in% loci)
    xx  <- x[, idx]
    xx@other$sa_best_performance <- perf
    xx@other$sa_stop_reason      <- reason
    xx
  }
  
  # ---------------------------------------------------------------
  # initialise
  # ---------------------------------------------------------------
  cat(sprintf("Initialising with method = '%s'...\n", init.method))
  
  init_panel  <- gl.select.panel(x, method = init.method, nl = nl,
                                 plot.out = FALSE, verbose = 0)
  curr_loci   <- locNames(init_panel)
  curr_perf   <- evaluate_panel(curr_loci)
  best_loci   <- curr_loci
  best_perf   <- curr_perf
  all_loci    <- locNames(x)
  stop_reason <- "completed"
  iter_done   <- 0
  T_sa        <- start_temp
  n_restarts  <- 0L
  
  history <- data.frame(
    iteration           = 0:n_iter,
    current_performance = c(curr_perf, rep(NA_real_,    n_iter)),
    best_performance    = c(best_perf, rep(NA_real_,    n_iter)),
    n_swap              = c(n_swap_max, rep(NA_integer_, n_iter)),
    restarted           = c(FALSE,     rep(NA,          n_iter))
  )
  
  cat(sprintf("Initial composite performance: %.4f\n", curr_perf))
  cat(sprintf("Running up to %d iterations. Press Escape to stop early.\n",
              n_iter))
  if (!is.null(stop.criterion))
    cat(sprintf("Auto-stop when performance >= %.3f\n", stop.criterion))
  cat("\n")
  
  # ---------------------------------------------------------------
  # single outer tryCatch — interrupt returns panel cleanly
  # ---------------------------------------------------------------
  tryCatch(
    expr = {
      
      for (iter in seq_len(n_iter)) {
        
        # restart from best if current has drifted too far below best
        restarted <- FALSE
        if (!is.null(restart_tol) &&
            curr_perf < best_perf - restart_tol) {
          curr_loci <- best_loci
          curr_perf <- best_perf
          restarted  <- TRUE
          n_restarts <- n_restarts + 1L
          cat(sprintf("\n  iter %d: restart #%d (curr was %.4f < best %.4f - %.3f)\n",
                      iter, n_restarts, curr_perf, best_perf, restart_tol))
        }
        
        # adaptive swap size: scales linearly with temperature
        n_swap   <- max(1L, round(n_swap_max * T_sa / start_temp))
        
        # propose swap
        out_locs <- sample(curr_loci, n_swap, replace = FALSE)
        in_locs  <- sample(setdiff(all_loci, curr_loci), n_swap,
                           replace = FALSE)
        new_loci <- c(setdiff(curr_loci, out_locs), in_locs)
        
        # evaluate
        new_perf <- evaluate_panel(new_loci)
        delta    <- new_perf - curr_perf
        
        # accept / reject
        if (delta > 0 || runif(1) < exp(delta / T_sa)) {
          curr_loci <- new_loci
          curr_perf <- new_perf
        }
        
        # track best
        if (curr_perf > best_perf) {
          best_perf <- curr_perf
          best_loci <- curr_loci
          cat(sprintf("\n  iter %d: new best = %.4f  (swap=%d)\n",
                      iter, best_perf, n_swap))
        }
        
        history$current_performance[iter + 1] <- curr_perf
        history$best_performance[iter + 1]    <- best_perf
        history$n_swap[iter + 1]              <- n_swap
        history$restarted[iter + 1]           <- restarted
        T_sa      <- T_sa * cooling
        iter_done <- iter
        
        # progress line
        cat(sprintf(
          "\r  iter %4d/%d  T=%.5f  swap=%d  curr=%.4f  best=%.4f  restarts=%d  ",
          iter, n_iter, T_sa, n_swap, curr_perf, best_perf, n_restarts
        ))
        flush.console()
        
        # live plot
        if (plot.out && iter %% plot.every == 0) {
          cat("\n")
          gg_live <- build_conv_plot(history, iter, best_perf)
          if (!is.null(gg_live)) print(gg_live)
        }
        
        # stop criterion
        if (!is.null(stop.criterion) && best_perf >= stop.criterion) {
          cat(sprintf("\nStop criterion (%.3f) reached at iteration %d.\n",
                      stop.criterion, iter))
          stop_reason <- "stop_criterion"
          break
        }
      }
      
      if (stop_reason == "completed")
        cat(sprintf("\n\nSA complete after %d iterations.\n", n_iter))
      cat(sprintf("Best: %.4f  |  Total restarts: %d\n",
                  best_perf, n_restarts))
      
      if (plot.out) {
        gg_f <- build_conv_plot(history, iter_done, best_perf, stop_reason)
        if (!is.null(gg_f)) print(gg_f)
      }
      
      build_return(best_loci, best_perf, history, iter_done, stop_reason)
    },
    
    interrupt = function(e) {
      
      tryCatch(
        cat(sprintf("\n\nEscape: iter=%d  best=%.4f  restarts=%d\n",
                    iter_done, best_perf, n_restarts)),
        interrupt = function(e2) invisible(NULL)
      )
      
      stop_reason <<- "interrupted"
      
      tryCatch(
        build_return(best_loci, best_perf, history, iter_done, "interrupted"),
        
        interrupt = function(e2) {
          tryCatch(
            minimal_return(best_loci, best_perf, "interrupted"),
            
            interrupt = function(e3) {
              tryCatch(
                { idx <- which(locNames(x) %in% best_loci); x[, idx] },
                interrupt = function(e4) x,
                error     = function(e4) x
              )
            },
            error = function(e3) x
          )
        },
        error = function(e2) x
      )
    }
  )
}