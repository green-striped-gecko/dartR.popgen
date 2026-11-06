## ------------------------------------------------------------------
## Step 3: Does drift_resistance optimisation improve panel durability?
##
## For each parameter, compare two SA-optimised panels:
##   a) parameter only:  weights = c(1, 0) for c(param, drift_resistance)
##   b) parameter + DR:  weights = c(1, 1) for c(param, drift_resistance)
##
## Then track performance over 20 generations using
## gl.check.future.panel(type = "user").
##
## Requires: sim_single_pop_panel.R, sim_metapop_panel.R,
##           and the custom panel functions sourced.
## ------------------------------------------------------------------

library(dartR.base)
library(dartR.sim)
library(ggplot2)

## ===================== DESIGN =====================

#' Build design for durability comparison (with vs without drift resistance)
#'
#' @param parameters  vector of parameters to test (one SA per parameter)
#' @param best_inits  named vector: parameter -> best init method
#'                    (from step 1/2, e.g. c(id="pic", Fst="hafpop"))
#' @param scenario    "single" or "metapop"
#' @param nl          vector of panel sizes
#' @param n_reps      replicates
#' @param n_iter      SA iterations
#' @param n_target    census sizes (single, vector)
#' @param n_per_pop   per-pop sizes (metapop, vector)
#' @param n_pops      populations (metapop)
#' @param n_loci      loci
#' @param n_founders  founders
#' @param n_burnin    burn-in generations
#' @param n_gen       forward generations to track
#' @return data.frame
build_design_durability <- function(parameters,
                                    best_inits,
                                    scenario   = "single",
                                    nl         = c(20, 50),
                                    n_reps     = 3,
                                    n_iter     = 200,
                                    n_target   = 100,
                                    n_per_pop  = 30,
                                    n_pops     = 10,
                                    n_loci     = 1000,
                                    n_founders = 20,
                                    n_burnin   = 4,
                                    n_gen      = 20) {

  ## two approaches per parameter: without DR (w=0) and with DR (w=1)
  approaches <- do.call(rbind, lapply(parameters, function(p) {
    init <- if (p %in% names(best_inits)) best_inits[p] else "random"
    data.frame(
      parameter   = p,
      approach    = c("param_only", "param_plus_DR"),
      dr_weight   = c(0, 1),
      init_method = init,
      stringsAsFactors = FALSE
    )
  }))

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
    grid$sim_id <- sprintf("D_%d_r%d", grid$n_target, grid$rep)
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
    grid$sim_id <- sprintf("D_%d_r%d", grid$n_per_pop, grid$rep)
    grid$n_pops <- n_pops
  }

  grid$run_id     <- sprintf("DR%03d", seq_len(nrow(grid)))
  grid$scenario   <- scenario
  grid$n_loci     <- n_loci
  grid$n_founders <- n_founders
  grid$n_burnin   <- n_burnin
  grid$n_gen      <- n_gen
  grid$n_iter     <- n_iter
  grid$done       <- FALSE

  grid <- grid[order(grid$parameter, grid$run_id), ]
  rownames(grid) <- grid$run_id
  grid
}


## ===================== SIMULATE WITH FORWARD GENERATIONS =====================

#' Simulate datasets with forward generations for durability testing
#'
#' @param design  design data.frame
#' @param verbose 0/1/2
#' @return named list of simulations keyed by sim_id
simulate_durability_datasets <- function(design, verbose = 1) {
  sims     <- list()
  scenario <- design$scenario[1]
  sim_ids  <- unique(design$sim_id)

  for (sid in sim_ids) {
    if (sid %in% names(sims)) next
    d <- design[design$sim_id == sid, ][1, ]

    if (verbose >= 1)
      cat(sprintf("\nSimulating %s (n_gen=%d)...\n", sid, d$n_gen))

    if (scenario == "single") {
      sims[[sid]] <- sim_single_pop(
        n_target   = d$n_target,
        n_loci     = d$n_loci,
        n_founders = d$n_founders,
        n_burnin   = d$n_burnin,
        n_gen      = d$n_gen,
        verbose    = verbose
      )
    } else {
      sims[[sid]] <- sim_metapop(
        n_pops     = d$n_pops,
        n_per_pop  = d$n_per_pop,
        n_loci     = d$n_loci,
        n_founders = d$n_founders,
        n_burnin   = d$n_burnin,
        n_gen      = d$n_gen,
        verbose    = verbose
      )
    }
  }
  sims
}


## ===================== RUN DURABILITY COMPARISON =====================

#' Process one durability row: SA panel selection + future check
#' @param d      one-row data.frame
#' @param sim    full simulation result (with sim_gls)
#' @return list with run_id, panel_perf_t0, future_performance (data.frame)
.run_one_durability <- function(d, sim) {
  xorig     <- sim$sim_gls[["gen_0"]]
  param     <- d$parameter
  sa_params <- c(param, "drift_resistance")
  sa_weights <- c(1, d$dr_weight)

  ## select panel via SA (even with dr_weight=0, SA optimises the param)
  panel <- tryCatch(
    gl.select.panel.combined(
      x           = xorig,
      nl          = d$nl,
      parameter   = sa_params,
      weights     = sa_weights,
      init.method = d$init_method,
      n_iter      = d$n_iter,
      plot.out    = FALSE,
      verbose     = 0
    ),
    error = function(e) NULL
  )

  if (is.null(panel))
    return(list(run_id = d$run_id, panel_perf_t0 = NA_real_,
                future_performance = NULL))

  ## check t0 performance for the target parameter
  power_params <- c("id", "parentage", "assignment", "hybridisation",
                    "drift_resistance")
  needs_xorig  <- !(param %in% power_params)

  t0_check <- tryCatch(
    gl.check.panel(
      x         = panel,
      xorig     = if (needs_xorig) xorig else NULL,
      parameter = c(param, "drift_resistance"),
      plot.out  = FALSE,
      verbose   = 0
    ),
    error = function(e) NULL
  )

  t0_perf <- if (!is.null(t0_check) && !is.null(t0_check[[param]]))
    as.numeric(t0_check[[param]]$performance) else NA_real_
  t0_dr   <- if (!is.null(t0_check) && !is.null(t0_check[["drift_resistance"]]))
    as.numeric(t0_check[["drift_resistance"]]$performance) else NA_real_

  ## run future panel check
  dur <- tryCatch(
    gl.check.future.panel(
      x         = panel,
      xorig     = xorig,
      parameter = c(param, "drift_resistance"),
      type      = "user",
      user.gl   = sim$sim_gls[-1],
      plot.out  = FALSE,
      verbose   = 0
    ),
    error = function(e) NULL
  )

  future_perf <- if (!is.null(dur) && !is.null(dur$performance))
    dur$performance else NULL

  list(run_id             = d$run_id,
       panel_perf_t0      = t0_perf,
       panel_dr_t0        = t0_dr,
       future_performance = future_perf)
}


#' Run durability comparison for all design rows
#'
#' @param design       data.frame from build_design_durability
#' @param sims         named list from simulate_durability_datasets
#' @param n_cores      parallel cores
#' @param verbose      0/1/2
#' @return list: $design (updated), $results
run_durability_comparison <- function(design, sims, n_cores = 1, verbose = 1) {
  todo <- which(!design$done)
  if (length(todo) == 0) {
    cat("All rows already done.\n")
    return(list(design = design, results = list()))
  }

  work <- lapply(todo, function(i) {
    d   <- design[i, ]
    sim <- sims[[d$sim_id]]
    list(d = d, sim = sim)
  })

  if (n_cores > 1) {
    cat(sprintf("Running %d durability tasks on %d cores...\n",
                length(work), n_cores))
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterEvalQ(cl, {
      suppressMessages(library(dartR.base))
      suppressMessages(library(dartR.sim))
    })

    fns_to_export <- ".run_one_durability"
    for (fn in c("gl.check.panel", "gl.select.panel",
                 "gl.select.panel.combined", "gl.check.future.panel")) {
      if (exists(fn, envir = globalenv()))
        fns_to_export <- c(fns_to_export, fn)
    }
    parallel::clusterExport(cl, fns_to_export, envir = globalenv())

    out_list <- parallel::parLapply(cl, work, function(w) {
      .run_one_durability(w$d, w$sim)
    })
  } else {
    out_list <- lapply(seq_along(work), function(j) {
      w <- work[[j]]
      if (verbose >= 1)
        cat(sprintf("[%s] %s | %s | nl=%d | dr_w=%d\n",
                    w$d$run_id, w$d$parameter, w$d$approach,
                    w$d$nl, w$d$dr_weight))
      .run_one_durability(w$d, w$sim)
    })
  }

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

#' Extract durability results into a long data.frame
#'
#' @param design   design data.frame
#' @param results  named list from run_durability_comparison
#' @return data.frame: run_id, parameter, approach, dr_weight, nl,
#'         generation, perf_parameter, perf_dr
extract_durability <- function(design, results) {
  rows <- list()
  drop_cols <- c("done", "init_method", "n_iter", "n_loci",
                 "n_founders", "n_burnin")

  for (i in seq_len(nrow(design))) {
    d <- design[i, ]
    if (!d$done) next
    r <- results[[d$run_id]]
    if (is.null(r) || is.null(r$future_performance)) next

    fp   <- r$future_performance
    param <- d$parameter

    ## extract the two tracked metrics from the future performance df
    fp_param <- fp[fp$parameter == param, , drop = FALSE]
    fp_dr    <- fp[fp$parameter == "drift_resistance", , drop = FALSE]

    if (nrow(fp_param) == 0) next

    out <- data.frame(
      run_id    = d$run_id,
      parameter = d$parameter,
      approach  = d$approach,
      dr_weight = d$dr_weight,
      nl        = d$nl,
      rep       = d$rep,
      scenario  = d$scenario,
      stringsAsFactors = FALSE
    )

    ## add pop size column
    if ("n_target" %in% names(d)) out$n_target <- d$n_target
    if ("n_per_pop" %in% names(d)) out$n_per_pop <- d$n_per_pop

    out <- out[rep(1, nrow(fp_param)), , drop = FALSE]
    out$generation      <- fp_param$generation
    out$perf_parameter  <- fp_param$performance

    ## match drift_resistance by generation
    out$perf_dr <- fp_dr$performance[match(out$generation, fp_dr$generation)]

    rows[[length(rows) + 1]] <- out
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}


## ===================== PLOTTING =====================

#' Plot durability: parameter performance over generations,
#' with vs without drift resistance optimisation
#'
#' @param dur_df  data.frame from extract_durability
#' @param title   plot title
#' @return ggplot
plot_durability <- function(dur_df, title = "Panel durability: with vs without drift resistance") {

  pop_col <- if ("n_target" %in% names(dur_df)) "n_target" else
    if ("n_per_pop" %in% names(dur_df)) "n_per_pop" else NULL

  ## aggregate across replicates
  grp <- c("parameter", "approach", "nl", "generation")
  if (!is.null(pop_col)) grp <- c(grp, pop_col)

  agg <- aggregate(
    as.formula(paste("cbind(perf_parameter, perf_dr) ~", paste(grp, collapse = " + "))),
    dur_df, function(x) mean(x, na.rm = TRUE)
  )

  agg$approach_label <- factor(agg$approach,
    levels = c("param_only", "param_plus_DR"),
    labels = c("Without DR", "With DR"))
  agg$nl_label <- paste0("nl=", agg$nl)
  if (!is.null(pop_col)) agg$N_label <- paste0("N=", agg[[pop_col]])

  ## parameter performance plot
  if (!is.null(pop_col)) {
    facet_form <- as.formula("nl_label + N_label ~ parameter")
  } else {
    facet_form <- as.formula("nl_label ~ parameter")
  }

  p1 <- ggplot(agg, aes(generation, perf_parameter,
                         colour = approach_label, linetype = approach_label)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0.9, linetype = 3, colour = "grey50") +
    facet_grid(facet_form) +
    scale_colour_manual(values = c("Without DR" = "steelblue",
                                   "With DR"    = "darkorange")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = title,
         subtitle = "Solid/dashed = with/without drift resistance in SA weights",
         x = "Generation", y = "Parameter performance",
         colour = "SA strategy", linetype = "SA strategy") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))

  ## drift resistance comparison plot
  p2 <- ggplot(agg, aes(generation, perf_dr,
                         colour = approach_label, linetype = approach_label)) +
    geom_line(linewidth = 0.8) +
    facet_grid(facet_form) +
    scale_colour_manual(values = c("Without DR" = "steelblue",
                                   "With DR"    = "darkorange")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Drift resistance over generations",
         subtitle = "Higher = more buffer remaining; DR-optimised panels should stay higher longer",
         x = "Generation", y = "Drift resistance",
         colour = "SA strategy", linetype = "SA strategy") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))

  list(parameter_plot = p1, dr_plot = p2)
}


## ===================== EXAMPLE USAGE =====================

if (FALSE) {

  ## source the frameworks first
  # source("sim_single_pop_panel.R")
  # source("sim_metapop_panel.R")

  ## -------------------------------------------------------
  ## A. Single population
  ## -------------------------------------------------------

  ## best init methods from step 1 (adjust to your results)
  best_inits_sp <- c(id = "picdart", parentage = "picdart",
                     relatedness = "picdart")

  ## build design
  design_dur_sp <- build_design_durability(
    parameters = c("id", "parentage", "relatedness"),
    best_inits = best_inits_sp,
    scenario   = "single",
    nl         = c(20, 50),
    n_reps     = 3,
    n_iter     = 200,
    n_target   = c(100, 50, 20),
    n_loci     = 1000,
    n_founders = 20,
    n_burnin   = 4,
    n_gen      = 20
  )
  cat(sprintf("Durability design (single): %d rows\n", nrow(design_dur_sp)))

  ## simulate (needs forward generations)
  sims_dur_sp <- simulate_durability_datasets(design_dur_sp, verbose = 1)

  ## run comparison
  out_dur_sp <- run_durability_comparison(design_dur_sp, sims_dur_sp,
                                          n_cores = 10, verbose = 1)
  design_dur_sp <- out_dur_sp$design
  res_dur_sp    <- out_dur_sp$results

  ## extract and plot
  dur_df_sp <- extract_durability(design_dur_sp, res_dur_sp)
  plots_sp  <- plot_durability(dur_df_sp, "Single pop: DR optimisation effect")
  print(plots_sp$parameter_plot)
  print(plots_sp$dr_plot)


  ## -------------------------------------------------------
  ## B. Multi-population
  ## -------------------------------------------------------

  best_inits_mp <- c(Fst = "hafpop", He = "pic",
                     assignment = "dapc", hybridisation = "dapc")

  design_dur_mp <- build_design_durability(
    parameters = c("Fst", "He", "assignment", "hybridisation"),
    best_inits = best_inits_mp,
    scenario   = "metapop",
    nl         = c(20, 50),
    n_reps     = 3,
    n_iter     = 200,
    n_per_pop  = 30,
    n_pops     = 10,
    n_loci     = 1000,
    n_founders = 20,
    n_burnin   = 4,
    n_gen      = 20
  )

  sims_dur_mp <- simulate_durability_datasets(design_dur_mp, verbose = 1)
  out_dur_mp  <- run_durability_comparison(design_dur_mp, sims_dur_mp,
                                           n_cores = 1, verbose = 1)
  design_dur_mp <- out_dur_mp$design
  res_dur_mp    <- out_dur_mp$results

  dur_df_mp <- extract_durability(design_dur_mp, res_dur_mp)
  plots_mp  <- plot_durability(dur_df_mp, "Metapop: DR optimisation effect")
  print(plots_mp$parameter_plot)
  print(plots_mp$dr_plot)
}
