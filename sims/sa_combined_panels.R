## ------------------------------------------------------------------
## Step 4: Combined multi-objective panels
##
## Can parameters be combined into a single panel?
##   - Compatible groups (all PIC or all diagnostic) combine easily
##   - Diverging objectives (PIC + diagnostic) need larger panels
##   - The SA resolves conflicts by finding compromise loci
##
## For each "recipe" (parameter combination) x panel size:
##   gl.select.panel.combined -> gl.check.panel on ALL recipe params
##
## Key figure: performance per parameter vs panel size, per recipe.
## Shows the minimum panel size for acceptable multi-objective performance.
##
## Requires: sim_single_pop_panel.R, sim_metapop_panel.R sourced
## ------------------------------------------------------------------

library(dartR.base)
library(dartR.sim)
library(ggplot2)

## ===================== RECIPE DEFINITIONS =====================

#' Define panel recipes: named lists of parameter combinations
#' Each recipe has: parameters, init method, scenario, description
panel_recipes <- list(

  ## compatible single-pop (all PIC-friendly)
  single_id_par = list(
    parameters  = c("id", "parentage", "drift_resistance"),
    init_method = "pic",
    scenario    = "single",
    description = "ID + Parentage (compatible, PIC)"
  ),

  single_all = list(
    parameters  = c("id", "parentage", "relatedness", "drift_resistance"),
    init_method = "pic",
    scenario    = "single",
    description = "All single-pop (compatible, PIC)"
  ),

  ## compatible multi-pop: differentiation group
  meta_diff = list(
    parameters  = c("Fst", "assignment", "hybridisation", "drift_resistance"),
    init_method = "dapc",
    scenario    = "metapop",
    description = "Differentiation group (compatible, DAPC)"
  ),

  ## compatible multi-pop: diversity group
  meta_diversity = list(
    parameters  = c("He", "Ho", "Nall", "drift_resistance"),
    init_method = "pic",
    scenario    = "metapop",
    description = "Diversity monitoring (compatible, PIC)"
  ),

  ## DIVERGING: single-pop power + multi-pop differentiation
  cross_id_assign = list(
    parameters  = c("id", "parentage", "assignment", "drift_resistance"),
    init_method = "hafpop",
    scenario    = "metapop",
    description = "ID + Assignment (diverging: PIC vs diagnostic)"
  ),

  ## DIVERGING: full management panel
  management = list(
    parameters  = c("id", "parentage", "assignment", "hybridisation",
                    "drift_resistance"),
    init_method = "hafpop",
    scenario    = "metapop",
    description = "Full management (diverging: PIC vs diagnostic)"
  ),

  ## universal: everything
  universal = list(
    parameters  = c("id", "parentage", "assignment", "He", "Fst",
                    "drift_resistance"),
    init_method = "hafpop",
    scenario    = "metapop",
    description = "Universal (maximum conflict)"
  )
)


## ===================== DESIGN =====================

#' Build design for combined panel comparison
#'
#' @param recipes    named list of recipe definitions
#' @param nl         vector of panel sizes (the key sweep axis)
#' @param n_reps     replicates
#' @param n_iter     SA iterations
#' @param n_target   single-pop census size
#' @param n_per_pop  metapop per-pop size
#' @param n_pops     metapop populations
#' @param n_loci     total loci
#' @param n_founders founders
#' @param n_burnin   burn-in
#' @return data.frame
build_design_combined <- function(recipes    = panel_recipes,
                                  nl         = c(10, 20, 30, 50, 75, 100),
                                  n_reps     = 3,
                                  n_iter     = 200,
                                  n_target   = 100,
                                  n_per_pop  = 30,
                                  n_pops     = 10,
                                  n_loci     = 1000,
                                  n_founders = 20,
                                  n_burnin   = 4) {

  rows <- list()
  for (rname in names(recipes)) {
    rec <- recipes[[rname]]
    for (nl_val in nl) {
      for (r in seq_len(n_reps)) {
        row <- data.frame(
          recipe      = rname,
          scenario    = rec$scenario,
          description = rec$description,
          parameters  = paste(rec$parameters, collapse = ","),
          init_method = rec$init_method,
          nl          = nl_val,
          rep         = r,
          stringsAsFactors = FALSE
        )
        if (rec$scenario == "single") {
          row$n_target  <- n_target
          row$sim_id    <- sprintf("CS_%d_r%d", n_target, r)
        } else {
          row$n_per_pop <- n_per_pop
          row$n_pops    <- n_pops
          row$sim_id    <- sprintf("CM_%d_r%d", n_per_pop, r)
        }
        rows[[length(rows) + 1]] <- row
      }
    }
  }

  grid <- do.call(rbind, rows)
  grid$run_id     <- sprintf("CP%03d", seq_len(nrow(grid)))
  grid$n_loci     <- n_loci
  grid$n_founders <- n_founders
  grid$n_burnin   <- n_burnin
  grid$n_iter     <- n_iter
  grid$done       <- FALSE

  grid <- grid[order(grid$recipe, grid$nl, grid$rep), ]
  rownames(grid) <- grid$run_id
  grid
}


## ===================== SIMULATE =====================

#' Simulate shared datasets for combined panel tests
#'
#' @param design  design data.frame
#' @param verbose 0/1/2
#' @return named list keyed by sim_id
simulate_combined_datasets <- function(design, verbose = 1) {
  sims    <- list()
  sim_ids <- unique(design$sim_id)

  for (sid in sim_ids) {
    if (sid %in% names(sims)) next
    d <- design[design$sim_id == sid, ][1, ]

    if (verbose >= 1) cat(sprintf("\nSimulating %s...\n", sid))

    if (d$scenario == "single") {
      sims[[sid]] <- sim_single_pop(
        n_target   = d$n_target,
        n_loci     = d$n_loci,
        n_founders = d$n_founders,
        n_burnin   = d$n_burnin,
        n_gen      = 0,
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


## ===================== RUN =====================

#' Process one combined-panel row
.run_one_combined <- function(d, xorig) {
  params <- strsplit(d$parameters, ",")[[1]]

  ## SA with equal weights on all parameters
  panel <- tryCatch(
    gl.select.panel.combined(
      x           = xorig,
      nl          = d$nl,
      parameter   = params,
      init.method = d$init_method,
      n_iter      = d$n_iter,
      plot.out    = FALSE,
      verbose     = 0
    ),
    error = function(e) NULL
  )

  if (is.null(panel))
    return(list(run_id = d$run_id,
                performance = setNames(rep(NA_real_, length(params)), params),
                composite   = NA_real_))

  ## check all parameters
  power_params <- c("id", "parentage", "assignment", "hybridisation",
                    "drift_resistance")
  needs_xorig  <- any(!params %in% power_params)

  check <- tryCatch(
    gl.check.panel(
      x         = panel,
      xorig     = if (needs_xorig) xorig else NULL,
      parameter = params,
      plot.out  = FALSE,
      verbose   = 0
    ),
    error = function(e) NULL
  )

  if (is.null(check)) {
    perf_vec <- setNames(rep(NA_real_, length(params)), params)
  } else {
    perf_vec <- vapply(check, function(r) {
      p <- r$performance
      if (is.null(p) || length(p) == 0 || is.na(p)) NA_real_ else as.numeric(p)
    }, numeric(1))
  }

  composite <- mean(perf_vec[!is.na(perf_vec)])

  list(run_id      = d$run_id,
       performance = perf_vec,
       composite   = composite)
}


#' Run combined panel tests
#'
#' @param design       data.frame from build_design_combined
#' @param sims         named list of simulations
#' @param n_cores      parallel cores
#' @param verbose      0/1/2
#' @return list: $design, $results
run_combined_comparison <- function(design, sims, n_cores = 1, verbose = 1) {
  todo <- which(!design$done)
  if (length(todo) == 0) {
    cat("All done.\n")
    return(list(design = design, results = list()))
  }

  work <- lapply(todo, function(i) {
    d   <- design[i, ]
    sim <- sims[[d$sim_id]]
    list(d = d, xorig = sim$sim_gls[["gen_0"]])
  })

  if (n_cores > 1) {
    cat(sprintf("Running %d combined tasks on %d cores...\n",
                length(work), n_cores))
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterEvalQ(cl, {
      suppressMessages(library(dartR.base))
      suppressMessages(library(dartR.sim))
    })
    fns_to_export <- ".run_one_combined"
    for (fn in c("gl.check.panel", "gl.select.panel",
                 "gl.select.panel.combined")) {
      if (exists(fn, envir = globalenv()))
        fns_to_export <- c(fns_to_export, fn)
    }
    parallel::clusterExport(cl, fns_to_export, envir = globalenv())

    out_list <- parallel::parLapply(cl, work, function(w) {
      .run_one_combined(w$d, w$xorig)
    })
  } else {
    out_list <- lapply(seq_along(work), function(j) {
      w <- work[[j]]
      if (verbose >= 1)
        cat(sprintf("[%s] %s | nl=%d | rep=%d\n",
                    w$d$run_id, w$d$recipe, w$d$nl, w$d$rep))
      .run_one_combined(w$d, w$xorig)
    })
  }

  results <- list()
  for (k in seq_along(out_list)) {
    r   <- out_list[[k]]
    idx <- todo[k]
    results[[r$run_id]] <- r
    design$done[idx]    <- TRUE
  }

  cat(sprintf("Done: %d tasks.\n", length(out_list)))
  list(design = design, results = results)
}


## ===================== EXTRACTION =====================

#' Extract combined panel results into long data.frame
extract_combined_results <- function(design, results) {
  rows <- list()
  for (i in seq_len(nrow(design))) {
    d <- design[i, ]
    if (!d$done) next
    r <- results[[d$run_id]]
    if (is.null(r)) next

    for (p in names(r$performance)) {
      row <- data.frame(
        run_id      = d$run_id,
        recipe      = d$recipe,
        description = d$description,
        scenario    = d$scenario,
        nl          = d$nl,
        rep         = d$rep,
        parameter   = p,
        performance = r$performance[p],
        composite   = r$composite,
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1]] <- row
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


## ===================== PLOTTING =====================

#' Performance vs panel size per recipe, faceted by parameter
#'
#' Shows the minimum panel size where all parameters meet threshold
plot_combined_vs_size <- function(comb_df,
                                  title = "Combined panels: performance vs panel size",
                                  threshold = 0.9) {
  agg <- aggregate(performance ~ recipe + description + nl + parameter,
                   comb_df, function(x) {
                     x <- x[!is.na(x)]
                     if (length(x) == 0) return(c(m = NA, lo = NA, hi = NA))
                     c(m = mean(x), lo = min(x), hi = max(x))
                   })
  agg <- do.call(data.frame, agg)
  names(agg)[5:7] <- c("m", "lo", "hi")
  agg <- agg[!is.na(agg$m), ]

  ggplot(agg, aes(nl, m, colour = recipe, fill = recipe)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.1, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = threshold, linetype = 3, colour = "grey40") +
    facet_wrap(~ parameter, scales = "free_y") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title    = title,
         subtitle = sprintf("Dotted = %.0f%% threshold; lines = recipes (parameter combos)",
                            threshold * 100),
         x = "Panel size (number of loci)", y = "Performance") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}


#' Heatmap: recipe x parameter at a given panel size
plot_combined_heatmap <- function(comb_df, nl_show = 50,
                                  title = NULL) {
  sub <- comb_df[comb_df$nl == nl_show, ]
  if (nrow(sub) == 0) stop("No data for nl = ", nl_show)

  agg <- aggregate(performance ~ recipe + description + parameter, sub,
                   mean, na.rm = TRUE)

  if (is.null(title)) title <- sprintf("Combined panel performance (nl = %d)", nl_show)

  ggplot(agg, aes(parameter, recipe, fill = performance)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = ifelse(is.na(performance), "\u2014",
                                  sprintf("%.2f", performance))),
              size = 3.2) +
    scale_fill_viridis_c(option = "D", limits = c(0, 1),
                         na.value = "grey85", name = "Perf") +
    scale_y_discrete(labels = function(x) {
      descs <- setNames(agg$description, agg$recipe)
      ifelse(x %in% names(descs), descs[x], x)
    }) +
    labs(title = title,
         subtitle = "Grey = parameter not in recipe",
         x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          plot.title = element_text(face = "bold"))
}


#' Minimum panel size to reach threshold on ALL parameters in a recipe
plot_min_panel_size <- function(comb_df, threshold = 0.9,
                                title = "Minimum panel size for all-pass") {
  ## for each recipe × nl × rep: does every parameter meet threshold?
  agg <- aggregate(performance ~ recipe + description + nl + rep,
                   comb_df, function(x) all(!is.na(x) & x >= threshold))
  names(agg)[5] <- "all_pass"

  ## fraction of reps that pass, per recipe × nl
  pass_rate <- aggregate(all_pass ~ recipe + description + nl, agg, mean)

  ggplot(pass_rate, aes(nl, all_pass, colour = recipe)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
    annotate("text", x = max(pass_rate$nl), y = 0.52,
             label = "50% pass rate", hjust = 1, size = 3, colour = "grey50") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = title,
         subtitle = sprintf("Fraction of replicates where ALL parameters >= %.0f%%",
                            threshold * 100),
         x = "Panel size (number of loci)",
         y = "Pass rate") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}


## ===================== EXAMPLE USAGE =====================

if (FALSE) {

  # source("sim_single_pop_panel.R")
  # source("sim_metapop_panel.R")

  ## build design — sweep panel sizes to find the crossover
  design_cp <- build_design_combined(
    recipes    = panel_recipes,
    nl         = c(10, 20, 30, 50, 75, 100),
    n_reps     = 3,
    n_iter     = 200,
    n_target   = 100,
    n_per_pop  = 30,
    n_pops     = 10,
    n_loci     = 1000,
    n_founders = 20,
    n_burnin   = 4
  )
  cat(sprintf("Combined design: %d rows\n", nrow(design_cp)))

  ## simulate
  sims_cp <- simulate_combined_datasets(design_cp, verbose = 1)

  ## run
  out_cp <- run_combined_comparison(design_cp, sims_cp, n_cores = 1, verbose = 1)
  design_cp <- out_cp$design
  res_cp    <- out_cp$results

  ## extract
  comb_df <- extract_combined_results(design_cp, res_cp)

  ## --- plots ---

  ## 1. Performance vs panel size, per parameter
  print(plot_combined_vs_size(comb_df))

  ## 2. Heatmap at nl=50
  print(plot_combined_heatmap(comb_df, nl_show = 50))

  ## 3. Minimum panel size for all-pass
  print(plot_min_panel_size(comb_df))

  ## 4. Summary table
  cat("\n--- Mean performance by recipe x panel size ---\n")
  tab <- aggregate(performance ~ recipe + nl, comb_df, mean, na.rm = TRUE)
  tab_wide <- reshape(tab, idvar = "recipe", timevar = "nl", direction = "wide")
  num_cols <- sapply(tab_wide, is.numeric)
  tab_wide[num_cols] <- round(tab_wide[num_cols], 3)
  print(tab_wide)
}
