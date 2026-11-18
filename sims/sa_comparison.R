## ------------------------------------------------------------------
## Step 2: SA optimisation comparison
##
## Does gl.select.panel.combined improve over rule-based selection?
## For each parameter, compare:
##   a) best rule-based method (from step 1)
##   b) SA with best method as init
##   c) SA with random init
##
## Reuses simulated datasets from method_comparison.R
##
## Requires: method_comparison.R sourced (for .run_one_method,
##           simulate_datasets, extract_method_results)
##           sim_single_pop_panel.R / sim_metapop_panel.R sourced
## ------------------------------------------------------------------

library(dartR.base)
library(dartR.sim)
library(ggplot2)

## ===================== FIND BEST METHODS =====================

#' From step 1 results, find the best method per parameter
#'
#' @param perf_df  data.frame from extract_method_results (step 1)
#' @return data.frame: parameter, best_method, mean_perf
find_best_methods <- function(perf_df) {
  ## mean performance per method × parameter (averaged over nl, reps, pop sizes)
  agg <- aggregate(performance ~ method + parameter, perf_df,
                   function(x) mean(x[!is.na(x)]))
  ## best method per parameter
  best <- do.call(rbind, lapply(split(agg, agg$parameter), function(sub) {
    sub <- sub[!is.na(sub$performance), ]
    if (nrow(sub) == 0) return(NULL)
    sub[which.max(sub$performance), ]
  }))
  names(best)[3] <- "mean_perf"
  rownames(best) <- NULL
  best
}


## ===================== DESIGN =====================

#' Build design for SA optimisation comparison
#'
#' For each parameter, creates rows for:
#'   - rule-based best (no SA, just the best method from step 1)
#'   - SA with best method as init
#'   - SA with random init
#'
#' @param best_methods data.frame from find_best_methods
#' @param scenario     "single" or "metapop"
#' @param nl           vector of panel sizes
#' @param n_reps       replicates
#' @param n_iter       SA iterations
#' @param n_target     census sizes (single pop, vector)
#' @param n_per_pop    per-pop sizes (metapop, vector)
#' @param n_pops       number of populations (metapop)
#' @param n_loci       loci in simulated dataset
#' @param n_founders   founders
#' @param n_burnin     burn-in generations
#' @return data.frame
build_design_sa <- function(best_methods,
                            scenario   = "single",
                            nl         = c(20, 50),
                            n_reps     = 3,
                            n_iter     = 200,
                            n_target   = 100,
                            n_per_pop  = 30,
                            n_pops     = 10,
                            n_loci     = 1000,
                            n_founders = 20,
                            n_burnin   = 4) {

  params <- best_methods$parameter

  ## three selection approaches per parameter
  approaches <- do.call(rbind, lapply(params, function(p) {
    bm <- best_methods$best_method[best_methods$parameter == p]
    data.frame(
      parameter   = p,
      approach    = c("rule_best", "SA_best_init", "SA_random_init"),
      init_method = c(bm, bm, "random"),
      use_sa      = c(FALSE, TRUE, TRUE),
      stringsAsFactors = FALSE
    )
  }))

  ## cross with pop sizes, nl, reps
  if (scenario == "single") {
    grid <- expand.grid(
      approach_id = seq_len(nrow(approaches)),
      nl          = nl,
      n_target    = n_target,
      rep         = seq_len(n_reps),
      stringsAsFactors = FALSE
    )
    grid <- cbind(grid, approaches[grid$approach_id, ])
    grid$approach_id <- NULL
    grid <- grid[order(grid$parameter, grid$n_target, grid$rep, grid$approach, grid$nl), ]
    grid$sim_id <- sprintf("S_%d_r%d", grid$n_target, grid$rep)
  } else {
    grid <- expand.grid(
      approach_id = seq_len(nrow(approaches)),
      nl          = nl,
      n_per_pop   = n_per_pop,
      rep         = seq_len(n_reps),
      stringsAsFactors = FALSE
    )
    grid <- cbind(grid, approaches[grid$approach_id, ])
    grid$approach_id <- NULL
    grid <- grid[order(grid$parameter, grid$n_per_pop, grid$rep, grid$approach, grid$nl), ]
    grid$sim_id <- sprintf("S_%d_r%d", grid$n_per_pop, grid$rep)
    grid$n_pops <- n_pops
  }

  grid$run_id     <- sprintf("SA%03d", seq_len(nrow(grid)))
  grid$scenario   <- scenario
  grid$n_loci     <- n_loci
  grid$n_founders <- n_founders
  grid$n_burnin   <- n_burnin
  grid$n_iter     <- n_iter
  grid$perf       <- NA_real_
  grid$done       <- FALSE

  rownames(grid) <- grid$run_id
  grid
}


## ===================== RUN SA COMPARISON =====================

#' Process one SA comparison row
#' @param d      one-row data.frame
#' @param xorig  genlight (gen_0)
#' @return list with run_id, performance, n_loci_selected
.run_one_sa <- function(d, xorig) {
  param <- d$parameter

  if (!d$use_sa) {
    ## rule-based only: select + check
    panel <- tryCatch(
      gl.select.panel(xorig, method = d$init_method, nl = d$nl, verbose = 0),
      error = function(e) NULL
    )
  } else {
    ## SA optimisation
    panel <- tryCatch(
      gl.select.panel.combined(
        x           = xorig,
        nl          = d$nl,
        parameter   = param,
        init.method = d$init_method,
        n_iter      = d$n_iter,
        plot.out    = FALSE,
        verbose     = 0
      ),
      error = function(e) NULL
    )
  }

  if (is.null(panel))
    return(list(run_id = d$run_id, performance = NA_real_,
                n_loci_selected = NA_integer_))

  ## check performance for this specific parameter
  power_params <- c("id", "parentage", "assignment", "hybridisation",
                    "drift_resistance")
  needs_xorig  <- !(param %in% power_params)

  check <- tryCatch(
    gl.check.panel(
      x         = panel,
      xorig     = if (needs_xorig) xorig else NULL,
      parameter = param,
      plot.out  = FALSE,
      verbose   = 0
    ),
    error = function(e) NULL
  )

  perf <- if (!is.null(check) && !is.null(check[[param]]$performance))
    as.numeric(check[[param]]$performance) else NA_real_

  list(run_id          = d$run_id,
       performance     = perf,
       n_loci_selected = nLoc(panel))
}


#' Run the SA comparison for all rows in the design
#'
#' @param design       data.frame from build_design_sa
#' @param sims         named list of simulations (reuse from step 1)
#' @param n_cores      parallel cores (1 = sequential)
#' @param source_files scripts to source on workers
#' @param verbose      0/1/2
#' @return list: $design (updated), $results
run_sa_comparison <- function(design, sims, n_cores = 1,
                              source_files = NULL, verbose = 1) {
  todo <- which(!design$done)
  if (length(todo) == 0) {
    cat("All rows already done.\n")
    return(list(design = design, results = list()))
  }

  work <- lapply(todo, function(i) {
    d   <- design[i, ]
    sim <- sims[[d$sim_id]]
    list(d = d, xorig = sim$sim_gls[["gen_0"]])
  })

  if (n_cores > 1) {
    cat(sprintf("Running %d SA tasks on %d cores...\n", length(work), n_cores))
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterEvalQ(cl, {
      suppressMessages(library(dartR.base))
      suppressMessages(library(dartR.sim))
    })

    ## export custom functions if sourced in main session
    fns_to_export <- ".run_one_sa"
    for (fn in c("gl.check.panel", "gl.select.panel",
                 "gl.select.panel.combined")) {
      if (exists(fn, envir = globalenv()))
        fns_to_export <- c(fns_to_export, fn)
    }
    parallel::clusterExport(cl, fns_to_export, envir = globalenv())

    if (!is.null(source_files)) {
      parallel::clusterExport(cl, "source_files", envir = environment())
      parallel::clusterEvalQ(cl, for (f in source_files) source(f))
    }

    out_list <- parallel::parLapply(cl, work, function(w) {
      .run_one_sa(w$d, w$xorig)
    })
  } else {
    out_list <- lapply(seq_along(work), function(j) {
      w <- work[[j]]
      if (verbose >= 1)
        cat(sprintf("[%s] %s | %s | nl=%d | iter=%s\n",
                    w$d$run_id, w$d$parameter, w$d$approach,
                    w$d$nl, ifelse(w$d$use_sa, w$d$n_iter, "—")))
      .run_one_sa(w$d, w$xorig)
    })
  }

  results <- list()
  for (k in seq_along(out_list)) {
    r   <- out_list[[k]]
    idx <- todo[k]
    results[[r$run_id]] <- r
    design$perf[idx]    <- r$performance
    design$done[idx]    <- TRUE
  }

  cat(sprintf("Done: %d tasks completed.\n", length(out_list)))
  list(design = design, results = results)
}


## ===================== PLOTTING =====================

#' Grouped barplot: rule-based vs SA for each parameter
#'
#' @param design  completed SA design (with perf column filled)
#' @param title   plot title
#' @return ggplot
plot_sa_comparison <- function(design, title = "SA optimisation: does it help?") {
  ## detect population size column
  pop_col <- if ("n_target" %in% names(design)) "n_target" else
    if ("n_per_pop" %in% names(design)) "n_per_pop" else NULL

  grp <- c("parameter", "approach", "nl")
  if (!is.null(pop_col)) grp <- c(grp, pop_col)

  agg <- aggregate(as.formula(paste("perf ~", paste(grp, collapse = " + "))),
                   design, function(x) {
                     x <- x[!is.na(x)]
                     if (length(x) == 0) return(c(m = NA, lo = NA, hi = NA))
                     c(m = mean(x), lo = min(x), hi = max(x))
                   })
  agg <- do.call(data.frame, agg)
  names(agg)[(ncol(agg)-2):ncol(agg)] <- c("m", "lo", "hi")
  agg <- agg[!is.na(agg$m), ]

  ## nice labels
  agg$approach <- factor(agg$approach,
                         levels = c("rule_best", "SA_best_init", "SA_random_init"),
                         labels = c("Rule-based\n(best method)", "SA\n(best init)", "SA\n(random init)"))
  agg$nl_label <- paste0("nl=", agg$nl)
  if (!is.null(pop_col)) agg$N_label <- paste0("N=", agg[[pop_col]])

  if (!is.null(pop_col)) {
    facet_form <- as.formula("nl_label + N_label ~ parameter")
  } else {
    facet_form <- as.formula("nl_label ~ parameter")
  }

  ggplot(agg, aes(x = approach, y = m, fill = approach)) +
    geom_col(colour = "grey30", linewidth = 0.2, alpha = 0.85, width = 0.7) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, linewidth = 0.4) +
    facet_grid(facet_form, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = c("Rule-based\n(best method)" = "grey70",
                                 "SA\n(best init)"           = "steelblue",
                                 "SA\n(random init)"         = "darkorange")) +
    labs(title = title,
         subtitle = "Bar = mean; whiskers = min-max across replicates",
         x = NULL, y = "Performance (0-1)") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none",
          strip.text = element_text(size = 9),
          plot.title = element_text(face = "bold"))
}


## ===================== EXAMPLE USAGE =====================

if (FALSE) {

  ## Assumes step 1 already run:
  ##   perf_sp, perf_mp from method_comparison.R
  ##   sims_sp, sims_mp from simulate_datasets()

  ## -------------------------------------------------------
  ## A. Single population
  ## -------------------------------------------------------

  ## find best methods from step 1
  best_sp <- find_best_methods(perf_sp)
  cat("Best methods (single pop):\n"); print(best_sp)

  ## build SA design
  design_sa_sp <- build_design_sa(
    best_methods = best_sp,
    scenario     = "single",
    nl           = c(20, 50),
    n_reps       = 3,
    n_iter       = 200,
    n_target     = c(100, 50, 20),
    n_loci       = 1000,
    n_founders   = 20,
    n_burnin     = 4
  )
  cat(sprintf("SA design (single): %d rows\n", nrow(design_sa_sp)))

  ## reuse simulated datasets from step 1
  out_sa_sp <- run_sa_comparison(design_sa_sp, sims_sp, n_cores = 10, verbose = 1)
  design_sa_sp <- out_sa_sp$design

  ## plot
  print(plot_sa_comparison(design_sa_sp, "Single pop: rule-based vs SA"))

  ## summary: improvement over rule-based
  cat("\n--- SA improvement summary (single pop) ---\n")
  imp <- aggregate(perf ~ parameter + approach + nl, design_sa_sp, mean, na.rm = TRUE)
  imp_wide <- reshape(imp, idvar = c("parameter", "nl"),
                      timevar = "approach", direction = "wide")
  num_cols <- sapply(imp_wide, is.numeric)
  imp_wide[num_cols] <- round(imp_wide[num_cols], 3)
  print(imp_wide)


  ## -------------------------------------------------------
  ## B. Multi-population
  ## -------------------------------------------------------

  best_mp <- find_best_methods(perf_mp)
  cat("\nBest methods (metapop):\n"); print(best_mp)

  design_sa_mp <- build_design_sa(
    best_methods = best_mp,
    scenario     = "metapop",
    nl           = c(20, 50),
    n_reps       = 3,
    n_iter       = 200,
    n_per_pop    = c(30, 50),
    n_pops       = 10,
    n_loci       = 1000,
    n_founders   = 20,
    n_burnin     = 4
  )
  cat(sprintf("SA design (metapop): %d rows\n", nrow(design_sa_mp)))

  out_sa_mp <- run_sa_comparison(design_sa_mp, sims_mp, n_cores = 1, verbose = 1)
  design_sa_mp <- out_sa_mp$design

  print(plot_sa_comparison(design_sa_mp, "Metapop: rule-based vs SA"))

  ## summary
  cat("\n--- SA improvement summary (metapop) ---\n")
  imp <- aggregate(perf ~ parameter + approach + nl, design_sa_mp, mean, na.rm = TRUE)
  imp_wide <- reshape(imp, idvar = c("parameter", "nl"),
                      timevar = "approach", direction = "wide")
  num_cols <- sapply(imp_wide, is.numeric)
  imp_wide[num_cols] <- round(imp_wide[num_cols], 3)
  print(imp_wide)
}
