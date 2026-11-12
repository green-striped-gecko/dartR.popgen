## ------------------------------------------------------------------
## Method comparison: which gl.select.panel method is best for each
## panel parameter?
##
## Workflow:
##   1. Simulate datasets (single-pop and metapop) — once per replicate
##   2. For each dataset × method × panel size: select panel, check perf
##   3. Aggregate into method × parameter performance table
##   4. Barplot: method vs performance per parameter
##
## Requires: sim_single_pop_panel.R and sim_metapop_panel.R sourced first
## ------------------------------------------------------------------

library(dartR.base)
library(dartR.sim)
library(ggplot2)

## ===================== DESIGN =====================

#' Build design for method comparison
#'
#' @param scenario   "single" or "metapop"
#' @param methods    character vector of gl.select.panel methods to test
#' @param parameters character vector of parameters to evaluate
#' @param nl         vector of panel sizes
#' @param n_reps     number of replicate datasets
#' @param n_target   census size (single pop)
#' @param n_per_pop  per-pop size (metapop)
#' @param n_pops     number of populations (metapop)
#' @param n_loci     loci in simulated dataset
#' @param n_founders founders
#' @param n_burnin   burn-in generations
#' @return data.frame with run_id, method, nl, rep, sim_id, etc.
build_design_methods <- function(scenario   = "single",
                                 methods    = c("pic", "random"),
                                 parameters = c("id", "drift_resistance"),
                                 nl         = c(20, 50),
                                 n_reps     = 3,
                                 n_target   = 100,
                                 n_per_pop  = 30,
                                 n_pops     = 10,
                                 n_loci     = 1000,
                                 n_founders = 20,
                                 n_burnin   = 4) {

  if (scenario == "single") {
    grid <- expand.grid(method   = methods,
                        nl       = nl,
                        n_target = n_target,
                        rep      = seq_len(n_reps),
                        stringsAsFactors = FALSE)
    grid <- grid[order(grid$n_target, grid$rep, grid$method, grid$nl), ]
    ## sim_id groups by n_target × rep (same population shared across methods)
    grid$sim_id <- sprintf("S_%d_r%d", grid$n_target, grid$rep)
  } else {
    grid <- expand.grid(method    = methods,
                        nl        = nl,
                        n_per_pop = n_per_pop,
                        rep       = seq_len(n_reps),
                        stringsAsFactors = FALSE)
    grid <- grid[order(grid$n_per_pop, grid$rep, grid$method, grid$nl), ]
    grid$sim_id <- sprintf("S_%d_r%d", grid$n_per_pop, grid$rep)
    grid$n_pops <- n_pops
  }

  grid$run_id     <- sprintf("C%03d", seq_len(nrow(grid)))
  grid$scenario   <- scenario
  grid$n_loci     <- n_loci
  grid$n_founders <- n_founders
  grid$n_burnin   <- n_burnin

  ## store parameter list as a single string for the design
  grid$parameters <- paste(parameters, collapse = ",")
  grid$done       <- FALSE

  rownames(grid) <- grid$run_id
  grid
}


## ===================== SIMULATION =====================

#' Simulate shared datasets (one per replicate)
#'
#' @param design  design data.frame from build_design_methods
#' @param verbose 0/1/2
#' @return named list of simulation results, keyed by sim_id
simulate_datasets <- function(design, verbose = 1) {
  sims     <- list()
  scenario <- design$scenario[1]
  sim_ids  <- unique(design$sim_id)

  for (sid in sim_ids) {
    if (sid %in% names(sims)) next
    d <- design[design$sim_id == sid, ][1, ]

    if (verbose >= 1)
      cat(sprintf("\nSimulating dataset %s (%s)...\n", sid, scenario))

    if (scenario == "single") {
      sims[[sid]] <- sim_single_pop(
        n_target   = d$n_target,
        n_loci     = d$n_loci,
        n_founders = d$n_founders,
        n_burnin   = d$n_burnin,
        n_gen      = 0,       # no forward sim needed, just the dataset
        verbose    = verbose
      )
    } else {
      sims[[sid]] <- sim_metapop(
        n_pops     = d$n_pops,
        n_per_pop  = d$n_per_pop,
        n_loci     = d$n_loci,
        n_founders = d$n_founders,
        n_burnin   = d$n_burnin,
        n_gen      = 0,
        verbose    = verbose
      )
    }
  }
  sims
}


## ===================== METHOD TESTING =====================

#' Process a single design row: select panel + check performance
#' (standalone function so it can be shipped to parallel workers)
#' @param d       one-row data.frame from the design
#' @param xorig   genlight (gen_0 from the simulation)
#' @return list with performance, n_loci_selected, error
.run_one_method <- function(d, xorig) {
  params <- strsplit(d$parameters, ",")[[1]]

  ## select panel
  panel <- tryCatch(
    gl.select.panel(xorig, method = d$method, nl = d$nl, verbose = 0),
    error = function(e) NULL
  )

  if (is.null(panel)) {
    return(list(
      run_id          = d$run_id,
      performance     = setNames(rep(NA_real_, length(params)), params),
      n_loci_selected = NA_integer_,
      error           = TRUE
    ))
  }

  ## determine which params need xorig
  power_params <- c("id", "parentage", "assignment", "hybridisation",
                    "drift_resistance")
  needs_xorig  <- setdiff(params, power_params)

  check_result <- tryCatch(
    gl.check.panel(
      x         = panel,
      xorig     = if (length(needs_xorig) > 0) xorig else NULL,
      parameter = params,
      plot.out  = FALSE,
      verbose   = 0
    ),
    error = function(e) NULL
  )

  if (is.null(check_result)) {
    perf_vec <- setNames(rep(NA_real_, length(params)), params)
  } else {
    perf_vec <- vapply(check_result, function(r) {
      p <- r$performance
      if (is.null(p) || length(p) == 0 || is.na(p)) NA_real_ else as.numeric(p)
    }, numeric(1))
  }

  list(
    run_id          = d$run_id,
    performance     = perf_vec,
    n_loci_selected = nLoc(panel),
    error           = is.null(check_result)
  )
}


#' Run panel selection + check for every row in the design
#'
#' @param design  design data.frame
#' @param sims    named list of simulated datasets (from simulate_datasets)
#' @param n_cores number of parallel cores (1 = sequential; >1 uses
#'                parallel::parLapply with a PSOCK cluster, works on Windows)
#' @param verbose 0/1/2 (verbose output only in sequential mode)
#' @return list: $design (updated), $results (named list with per-row results)
run_method_comparison <- function(design, sims, n_cores = 1, verbose = 1) {

  todo <- which(!design$done)
  if (length(todo) == 0) {
    cat("All rows already done.\n")
    return(list(design = design, results = list()))
  }

  ## build work list: each element is list(d = row, xorig = genlight)
  work <- lapply(todo, function(i) {
    d   <- design[i, ]
    sim <- sims[[d$sim_id]]
    list(d = d, xorig = sim$sim_gls[["gen_0"]])
  })

  if (n_cores > 1) {
    ## ---- parallel (PSOCK, Windows-safe) ----
    cat(sprintf("Running %d tasks on %d cores (PSOCK)...\n",
                length(work), n_cores))

    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    ## load dartR on each worker
    parallel::clusterEvalQ(cl, {
      suppressMessages(library(dartR.base))
      suppressMessages(library(dartR.sim))
    })

    ## export the worker function (defined at script top level = global env)
    parallel::clusterExport(cl, ".run_one_method", envir = globalenv())

    out_list <- parallel::parLapply(cl, work, function(w) {
      .run_one_method(w$d, w$xorig)
    })

  } else {
    ## ---- sequential ----
    out_list <- lapply(seq_along(work), function(j) {
      w <- work[[j]]
      if (verbose >= 1)
        cat(sprintf("[%s] method=%s, nl=%d, rep=%d\n",
                    w$d$run_id, w$d$method, w$d$nl, w$d$rep))
      .run_one_method(w$d, w$xorig)
    })
  }

  ## collect results
  results <- list()
  for (k in seq_along(out_list)) {
    r   <- out_list[[k]]
    idx <- todo[k]
    results[[r$run_id]] <- r
    design$done[idx]    <- TRUE
  }

  cat(sprintf("Done: %d tasks completed.\n", length(out_list)))
  list(design = design, results = results)
}


## ===================== EXTRACTION =====================

#' Extract method comparison results into a long data.frame
#'
#' @param design  design data.frame
#' @param results named list from run_method_comparison
#' @return data.frame: run_id, method, nl, rep, parameter, performance
extract_method_results <- function(design, results) {
  rows <- list()

  for (i in seq_len(nrow(design))) {
    d <- design[i, ]
    if (!d$done) next
    r <- results[[d$run_id]]
    if (is.null(r)) next

    params <- names(r$performance)
    for (p in params) {
      row_data <- data.frame(
        run_id    = d$run_id,
        scenario  = d$scenario,
        method    = d$method,
        nl        = d$nl,
        rep       = d$rep,
        parameter = p,
        performance = r$performance[p],
        n_loci_selected = r$n_loci_selected,
        error     = r$error,
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1]] <- row_data
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


## ===================== PLOTTING =====================

#' Barplot: method vs performance per parameter
#'
#' @param perf_df  data.frame from extract_method_results
#' @param title    plot title
#' @return ggplot object
plot_method_comparison <- function(perf_df, title = "Method comparison") {
  ## detect population size column
  pop_col <- if ("n_target" %in% names(perf_df)) "n_target" else
    if ("n_per_pop" %in% names(perf_df)) "n_per_pop" else NULL

  ## aggregate: mean + range per method × parameter × nl × pop_size
  grp <- c("method", "parameter", "nl")
  if (!is.null(pop_col)) grp <- c(grp, pop_col)

  agg <- aggregate(as.formula(paste("performance ~", paste(grp, collapse = " + "))),
                   perf_df, function(x) {
                     x <- x[!is.na(x)]
                     if (length(x) == 0) return(c(m = NA, lo = NA, hi = NA))
                     c(m = mean(x), lo = min(x), hi = max(x))
                   })
  agg <- do.call(data.frame, agg)
  names(agg)[(ncol(agg)-2):ncol(agg)] <- c("m", "lo", "hi")

  ## drop methods that failed everywhere for a parameter
  agg <- agg[!is.na(agg$m), ]

  agg$nl_label <- paste0("nl=", agg$nl)
  if (!is.null(pop_col))
    agg$N_label <- paste0("N=", agg[[pop_col]])

  ## build facet formula
  if (!is.null(pop_col)) {
    facet_form <- as.formula("nl_label + N_label ~ parameter")
  } else {
    facet_form <- as.formula("nl_label ~ parameter")
  }

  ggplot(agg, aes(x = method, y = m, fill = method)) +
    geom_col(colour = "grey30", linewidth = 0.2, alpha = 0.8) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.3, linewidth = 0.4) +
    facet_grid(facet_form, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = title,
         subtitle = "Bar = mean across replicates; whiskers = min-max range",
         x = NULL, y = "Performance (0-1)") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          strip.text = element_text(size = 9),
          plot.title = element_text(face = "bold"))
}


## ===================== EXAMPLE USAGE =====================

if (FALSE) {  # set to TRUE to run

  ## source the simulation frameworks first
  # source("sim_single_pop_panel.R")
  # source("sim_metapop_panel.R")

  ## -------------------------------------------------------
  ## A. Single-population methods
  ## -------------------------------------------------------
  design_sp <- build_design_methods(
    scenario   = "single",
    methods    = c("pic", "picdart", "random", "hafall", "stratified"),
    parameters = c("id", "parentage", "relatedness",
                   "Ho_ind", "Fis_ind", "drift_resistance"),
    nl         = c(20, 50, 100),
    n_target   = c(100, 50, 20),   # sweep over Ne
    n_reps     = 20,
    n_loci     = 1000,
    n_founders = 20,
    n_burnin   = 4
  )
  cat(sprintf("Single-pop design: %d rows (%d simulations)\n",
              nrow(design_sp), length(unique(design_sp$sim_id))))

  ## simulate shared datasets
  sims_sp <- simulate_datasets(design_sp, verbose = 1)

  ## run comparison (set n_cores > 1 for parallel on Windows)
  out_sp <- run_method_comparison(design_sp, sims_sp, n_cores = 10, verbose = 1)
  design_sp <- out_sp$design
  res_sp    <- out_sp$results

  ## extract and plot
  perf_sp <- extract_method_results(design_sp, res_sp)
  print(plot_method_comparison(perf_sp, "Single-population: method comparison"))

  ## summary table
  cat("\n--- Single-pop summary (mean performance) ---\n")
  tab <- reshape(
    aggregate(performance ~ method + parameter, perf_sp, mean, na.rm = TRUE),
    idvar = "method", timevar = "parameter", direction = "wide")
  num_cols <- sapply(tab, is.numeric)
  tab[num_cols] <- round(tab[num_cols], 3)
  print(tab)


  ## -------------------------------------------------------
  ## B. Multi-population methods
  ## -------------------------------------------------------
  design_mp <- build_design_methods(
    scenario   = "metapop",
    methods    = c("pic", "picdart", "random", "dapc", "pahigh",
                   "hafpop", "hafall", "stratified", "monopop"),
    parameters = c("Fst", "He", "Ho", "Nall", "Fis",
                   "assignment", "hybridisation", "drift_resistance"),
    nl         = c(20, 50),
    n_per_pop  = c(30),#, 50),   # sweep over per-pop size
    n_reps     = 5,
    n_pops     = 10,
    n_loci     = 1000,
    n_founders = 20,
    n_burnin   = 4
  )
  cat(sprintf("Metapop design: %d rows (%d simulations)\n",
              nrow(design_mp), length(unique(design_mp$sim_id))))

  ## simulate shared datasets
  sims_mp <- simulate_datasets(design_mp, verbose = 1)

  ## run comparison
  out_mp <- run_method_comparison(design_mp, sims_mp, n_cores = 10, verbose = 1)
  design_mp <- out_mp$design
  res_mp    <- out_mp$results

  ## extract and plot
  perf_mp <- extract_method_results(design_mp, res_mp)
  print(plot_method_comparison(perf_mp, "Multi-population: method comparison"))

  ## summary table
  cat("\n--- Metapop summary (mean performance) ---\n")
  tab <- reshape(
    aggregate(performance ~ method + parameter, perf_mp, mean, na.rm = TRUE),
    idvar = "method", timevar = "parameter", direction = "wide")
  num_cols <- sapply(tab, is.numeric)
  tab[num_cols] <- round(tab[num_cols], 3)
  print(tab)


  ## -------------------------------------------------------
  ## C. Combined overview heatmap
  ## -------------------------------------------------------
  #perf_all <- rbind(perf_sp, perf_mp)
  #perf_all <- perf_sp
  #perf_all <- perf_mp
  
  ## heatmap: method x parameter, fill = mean performance
  heat_df <- aggregate(performance ~ method + parameter, perf_all,
                       mean, na.rm = TRUE)
  heat_df$performance[is.nan(heat_df$performance)] <- NA

  print(
    ggplot(heat_df, aes(parameter, method, fill = performance)) +
      geom_tile(colour = "white", linewidth = 0.8) +
      geom_text(aes(label = ifelse(is.na(performance), "—",
                                    sprintf("%.2f", performance))),
                size = 3) +
      scale_fill_viridis_c(option = "D", limits = c(0, 1),
                           na.value = "grey85", name = "Performance") +
      labs(title = "Method × Parameter performance overview",
           subtitle = "Mean across panel sizes and replicates; — = method not applicable",
           x = NULL, y = NULL) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1),
            plot.title = element_text(face = "bold"))
  )
}
