## ------------------------------------------------------------------
## SNP panel durability simulation framework
##
## sim_single_pop()    — runs one WF simulation, returns genlight list + Ne
## build_design()      — creates a design data.frame from parameter grid
## run_simulations()   — runs all rows of the design, stores in results list
## run_panel_checks()  — runs gl.check.future.panel on each row
## extract_performance() — pulls performance into a long data.frame
##
## Workflow:
##   design  <- build_design(...)
##   out     <- run_simulations(design)
##   design  <- out$design
##   res     <- out$results
##   out2    <- run_panel_checks(design, res, ...)
##   design  <- out2$design
##   res     <- out2$results
##   perf_df <- extract_performance(design, res)
## ------------------------------------------------------------------

library(dartR.base)
library(dartR.sim)
library(ggplot2)

## ===================== CORE SIMULATION =====================

#' True Wright-Fisher breeding (multinomial resampling)
#' @param x       genlight of current generation
#' @param n_target constant population size
#' @param pop_name population label
#' @return list: $gl (genlight), $Vk, $Ne_Vk
breed_generation <- function(x, n_target, pop_name) {
  n    <- nInd(x)
  gmat <- as.matrix(x)
  L    <- ncol(gmat)

  dad_idx <- sample(n, n_target, replace = TRUE)
  mum_idx <- sample(n, n_target, replace = TRUE)

  k     <- tabulate(c(dad_idx, mum_idx), nbins = n)
  k_bar <- mean(k)
  Vk    <- var(k)
  Ne_Vk <- (n * k_bar - 1) / (k_bar - 1 + Vk / k_bar)

  a_dad <- matrix(rbinom(n_target * L, 1L, gmat[dad_idx, ] / 2), n_target, L)
  a_mum <- matrix(rbinom(n_target * L, 1L, gmat[mum_idx, ] / 2), n_target, L)
  offspring_mat <- a_dad + a_mum
  offspring_mat[is.na(gmat[dad_idx, ]) | is.na(gmat[mum_idx, ])] <- NA_integer_

  colnames(offspring_mat) <- colnames(gmat)
  gl_new <- new("genlight", gen = offspring_mat, ploidy = 2L,
                ind.names = paste0(pop_name, "_", seq_len(n_target)),
                loc.names = colnames(gmat), loc.all = x@loc.all)
  pop(gl_new) <- factor(rep(pop_name, n_target))

  list(gl = gl_new, Vk = Vk, Ne_Vk = Ne_Vk)
}


#' Ne from He decay (consecutive-generation harmonic mean)
estimate_Ne_He <- function(sim_gls) {
  n_gen <- length(sim_gls) - 1
  He_vec <- vapply(sim_gls, function(gl) {
    af <- colMeans(as.matrix(gl), na.rm = TRUE) / 2
    mean(2 * af * (1 - af), na.rm = TRUE)
  }, numeric(1))

  r <- He_vec[-1] / He_vec[-(n_gen + 1)]
  ne <- 1 / (2 * (1 - r))
  valid <- is.finite(ne) & ne > 0
  if (sum(valid) == 0) return(NA_real_)
  1 / mean(1 / ne[valid])
}


#' Run a single-population WF simulation
#'
#' @param n_target   census / population size
#' @param n_loci     number of loci
#' @param n_founders number of unrelated founders
#' @param n_burnin   burn-in generations (builds pedigree structure)
#' @param n_gen      forward generations to track
#' @param pop_name   population label
#' @param verbose    0 = silent, 1 = summary, 2 = per-generation
#' @return list: $sim_gls, $Ne_Vk, $Ne_He, $Vk_mean, $He_trajectory
sim_single_pop <- function(n_target,
                           n_loci     = 1000,
                           n_founders = 20,
                           n_burnin   = 4,
                           n_gen      = 20,
                           pop_name   = "pop1",
                           verbose    = 1) {

  ## 1. Founders
  founders <- gl.sim.Neconst(n_founders, nlocs = n_loci, verbose = 0)
  pop(founders) <- factor(rep(pop_name, n_founders))
  founders <- gl.compliance.check(founders, verbose = 0)

  ## 2. Burn-in
  x <- founders
  res_b <- breed_generation(x, n_target, pop_name)
  x <- gl.compliance.check(res_b$gl, verbose = 0)

  for (b in seq_len(n_burnin)) {
    if (verbose >= 2) cat(sprintf("  burn-in %d\n", b))
    res_b <- breed_generation(x, n_target, pop_name)
    x <- gl.compliance.check(res_b$gl, verbose = 0)
  }

  ## 3. Forward simulation
  sim_gls   <- vector("list", n_gen + 1)
  sim_gls[[1]] <- x
  names(sim_gls)[1] <- "gen_0"
  ne_vk_vec <- vk_vec <- numeric(n_gen)

  for (g in seq_len(n_gen)) {
    if (verbose >= 2) cat(sprintf("  gen %d\n", g))
    res_g <- breed_generation(x, n_target, pop_name)
    x <- gl.compliance.check(res_g$gl, verbose = 0)
    sim_gls[[g + 1]] <- x
    names(sim_gls)[g + 1] <- paste0("gen_", g)
    ne_vk_vec[g] <- res_g$Ne_Vk
    vk_vec[g]    <- res_g$Vk
  }

  ## 4. Ne estimates
  Ne_Vk <- 1 / mean(1 / ne_vk_vec[is.finite(ne_vk_vec) & ne_vk_vec > 0])
  Ne_He <- estimate_Ne_He(sim_gls)

  ## He trajectory
  He_traj <- vapply(sim_gls, function(gl) {
    af <- colMeans(as.matrix(gl), na.rm = TRUE) / 2
    mean(2 * af * (1 - af), na.rm = TRUE)
  }, numeric(1))

  if (verbose >= 1)
    cat(sprintf("  N=%d | Ne_Vk=%.1f | Ne_He=%.1f | Vk=%.2f\n",
                n_target, Ne_Vk, Ne_He, mean(vk_vec)))

  list(sim_gls       = sim_gls,
       Ne_Vk         = Ne_Vk,
       Ne_He         = Ne_He,
       Vk_mean       = mean(vk_vec),
       ne_vk_per_gen = ne_vk_vec,
       He_trajectory = He_traj)
}


## ===================== DESIGN + EXECUTION =====================

#' Build a design data.frame from parameter grid
#'
#' @param n_targets  vector of census sizes
#' @param n_loci     vector of locus counts (can be >1 for panel-size sweep)
#' @param n_reps     number of replicates
#' @param n_founders founders per simulation
#' @param n_burnin   burn-in generations
#' @param n_gen      forward generations
#' @return data.frame with run_id (zero-padded), parameters, result placeholders
build_design <- function(n_targets,
                         n_loci     = 1000,
                         n_reps     = 3,
                         n_founders = 20,
                         n_burnin   = 4,
                         n_gen      = 20) {

  grid <- expand.grid(n_target  = n_targets,
                      n_loci    = n_loci,
                      rep       = seq_len(n_reps),
                      stringsAsFactors = FALSE)
  grid <- grid[order(grid$n_target, grid$n_loci, grid$rep), ]

  grid$run_id     <- sprintf("%03d", seq_len(nrow(grid)))
  grid$n_founders <- n_founders
  grid$n_burnin   <- n_burnin
  grid$n_gen      <- n_gen
  grid$Ne_Vk      <- NA_real_
  grid$Ne_He      <- NA_real_
  grid$Vk_mean    <- NA_real_
  grid$sim_done   <- FALSE
  grid$panel_done <- FALSE

  rownames(grid) <- grid$run_id
  grid[, c("run_id", "n_target", "n_loci", "rep",
           "n_founders", "n_burnin", "n_gen",
           "Ne_Vk", "Ne_He", "Vk_mean",
           "sim_done", "panel_done")]
}


#' Run all simulations in the design
#'
#' @param design  data.frame from build_design()
#' @param verbose 0/1/2
#' @return list: $design (updated with Ne), $results (named list)
run_simulations <- function(design, verbose = 1) {
  results <- list()

  for (i in seq_len(nrow(design))) {
    d <- design[i, ]
    if (d$sim_done) { cat(sprintf("[%s] already done, skipping\n", d$run_id)); next }

    if (verbose >= 1)
      cat(sprintf("\n[%s] N=%d, L=%d, rep=%d\n",
                  d$run_id, d$n_target, d$n_loci, d$rep))

    res <- sim_single_pop(
      n_target   = d$n_target,
      n_loci     = d$n_loci,
      n_founders = d$n_founders,
      n_burnin   = d$n_burnin,
      n_gen      = d$n_gen,
      verbose    = verbose
    )

    results[[d$run_id]] <- res
    design$Ne_Vk[i]    <- res$Ne_Vk
    design$Ne_He[i]    <- res$Ne_He
    design$Vk_mean[i]  <- res$Vk_mean
    design$sim_done[i] <- TRUE
  }

  list(design = design, results = results)
}


#' Run gl.check.future.panel on each completed simulation
#'
#' @param design     data.frame (sim_done rows will be checked)
#' @param results    named list from run_simulations
#' @param nl         panel size (number of loci to select)
#' @param method     gl.select.panel method, OR "combined" to use
#'                   gl.select.panel.combined (SA optimiser)
#' @param parameter  parameters to evaluate (also used as optimisation
#'                   targets when method = "combined")
#' @param weights    numeric vector of weights for combined (NULL = equal)
#' @param n_iter     SA iterations for combined (default 100)
#' @param init.method seed method for combined (default "random")
#' @param start_temp  SA starting temperature (default 0.1)
#' @param cooling     SA cooling rate (NULL = auto)
#' @param restart_tol SA restart tolerance (NULL = no restarts)
#' @param ...        extra args passed to gl.check.future.panel
#'                   (e.g. ref, metric, target, corr.method)
#' @return list: $design (updated), $results (panel_check added)
run_panel_checks <- function(design, results,
                             nl          = 50,
                             method      = "random",
                             parameter   = c("id", "parentage",
                                             "relatedness", "Ho_ind",
                                             "Fis_ind",
                                             "drift_resistance"),
                             weights     = NULL,
                             n_iter      = 100,
                             init.method = "random",
                             start_temp  = 0.1,
                             cooling     = NULL,
                             restart_tol = NULL,
                             ...) {

  for (i in seq_len(nrow(design))) {
    d <- design[i, ]
    if (!d$sim_done)  { warning(sprintf("[%s] sim not done, skip", d$run_id)); next }
    if (d$panel_done) { cat(sprintf("[%s] panel already done, skip\n", d$run_id)); next }

    r      <- results[[d$run_id]]
    xorig  <- r$sim_gls[["gen_0"]]

    ## population size label (works for both single and metapop designs)
    n_label <- if ("n_target" %in% names(d)) d$n_target
               else if ("n_per_pop" %in% names(d)) d$n_per_pop
               else "?"

    if (method == "combined") {
      cat(sprintf("\n[%s] SA panel: N=%s, L=%d, rep=%d, nl=%d, iter=%d\n",
                  d$run_id, n_label, d$n_loci, d$rep, nl, n_iter))

      panel <- gl.select.panel.combined(
        x           = xorig,
        nl          = nl,
        parameter   = parameter,
        weights     = weights,
        init.method = init.method,
        n_iter      = n_iter,
        start_temp  = start_temp,
        cooling     = cooling,
        restart_tol = restart_tol,
        verbose     = 0,
        plot.out    = FALSE
      )
    } else {
      cat(sprintf("\n[%s] panel check: N=%s, L=%d, rep=%d, nl=%d, method=%s\n",
                  d$run_id, n_label, d$n_loci, d$rep, nl, method))

      panel <- gl.select.panel(xorig, method = method, nl = nl, verbose = 0)
    }

    ## ... passes only gl.check.future.panel args (ref, metric, target, etc.)
    dur <- gl.check.future.panel(
      x         = panel,
      xorig     = xorig,
      parameter = parameter,
      type      = "user",
      user.gl   = r$sim_gls[-1],
      plot.out  = FALSE,
      verbose   = 0,
      ...
    )

    results[[d$run_id]]$panel_check  <- dur
    results[[d$run_id]]$panel_nl     <- nl
    results[[d$run_id]]$panel_method <- method
    design$panel_done[i] <- TRUE
  }

  list(design = design, results = results)
}


## ===================== RESULTS EXTRACTION =====================

#' Extract panel performance into a long data.frame for plotting
#'
#' @param design   design data.frame
#' @param results  results list (with $panel_check)
#' @return data.frame: run_id, n_target, n_loci, rep, Ne_Vk, Ne_He,
#'         nl, generation, parameter, performance
extract_performance <- function(design, results) {
  rows <- list()
  ## columns to carry from design (drop housekeeping flags)
  drop_cols <- c("sim_done", "panel_done")

  for (i in seq_len(nrow(design))) {
    d <- design[i, ]
    if (!d$panel_done) next
    r   <- results[[d$run_id]]
    dur <- r$panel_check
    if (is.null(dur) || is.null(dur$performance)) next

    perf <- dur$performance
    ## design columns (one row, recycled across performance rows)
    dcols <- d[, setdiff(names(d), drop_cols), drop = FALSE]
    dcols$nl <- r$panel_nl
    rownames(dcols) <- NULL

    rows[[length(rows) + 1]] <- cbind(
      dcols[rep(1, nrow(perf)), , drop = FALSE],
      perf[, c("generation", "parameter", "performance")],
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


## ===================== EXAMPLE USAGE =====================

if (FALSE) {  # set to TRUE to run

  ## 1. Build design
  design <- build_design(
    n_targets  = c(100, 50, 20),
    n_loci     = 1000,
    n_reps     = 3,
    n_founders = 20,
    n_burnin   = 4,
    n_gen      = 20
  )
  print(design)

  ## 2. Run simulations
  out <- run_simulations(design, verbose = 1)
  design <- out$design
  res    <- out$results

  ## quick Ne check
  print(design[, c("run_id", "n_target", "rep", "Ne_Vk", "Ne_He", "Vk_mean")])

  ## 3. Run panel checks — compare random vs SA-optimised
  ##    (reuse the same simulations, just different panel selection)

  ## 3a. Random panel
  out_rnd <- run_panel_checks(
    design, res,
    nl        = 50,
    method    = "random",
    parameter = c("id", "parentage", "relatedness", "Ho_ind",
                  "Fis_ind", "drift_resistance"),
    ref       = "by.pop",
    metric    = "kinship",
    target    = 0.9
  )
  design_rnd <- out_rnd$design
  res_rnd    <- out_rnd$results
  perf_rnd   <- extract_performance(design_rnd, res_rnd)
  perf_rnd$selection <- "random"

  ## 3b. SA-optimised panel (reset panel_done first)
  design$panel_done <- FALSE
  out_sa <- run_panel_checks(
    design, res,
    nl          = 50,
    method      = "combined",
    parameter   = c("id", "parentage", "relatedness", "Ho_ind",
                    "Fis_ind", "drift_resistance"),
    n_iter      = 100,
    init.method = "random",
    ref         = "by.pop",
    metric      = "kinship",
    target      = 0.9
  )
  design_sa <- out_sa$design
  res_sa    <- out_sa$results
  perf_sa   <- extract_performance(design_sa, res_sa)
  perf_sa$selection <- "combined (SA)"

  ## 4. Compare
  perf_all <- rbind(perf_rnd, perf_sa)

  perf_agg <- aggregate(performance ~ n_target + generation + parameter + selection,
                        perf_all, function(x)
                          c(m = mean(x), lo = quantile(x, 0.1), hi = quantile(x, 0.9)))
  perf_agg <- do.call(data.frame, perf_agg)
  names(perf_agg)[5:7] <- c("m", "lo", "hi")
  perf_agg$label <- factor(paste0("N=", perf_agg$n_target))

  print(
    ggplot(perf_agg, aes(generation, m, colour = label, fill = label,
                         linetype = selection)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.1, colour = NA) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 0.9, linetype = 3) +
      facet_wrap(~ parameter, scales = "free_y") +
      labs(title = "Panel durability: random vs SA-optimised",
           x = "Generation", y = "Performance",
           colour = "Census N", fill = "Census N",
           linetype = "Selection") +
      theme_bw()
  )

  ## Ne diagnostic plot
  print(
    ggplot(design, aes(n_target, Ne_Vk)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
      geom_point(aes(colour = "Vk"), size = 3) +
      geom_point(aes(y = Ne_He, colour = "He"), size = 3, shape = 17) +
      coord_equal(xlim = c(0, max(design$n_target) * 1.1),
                  ylim = c(0, max(design$n_target) * 1.1)) +
      labs(title = "Census N vs Ne", x = "Census N", y = "Ne",
           colour = "Method") +
      theme_bw()
  )

  ## He decay plot
  he_df <- do.call(rbind, lapply(seq_len(nrow(design)), function(i) {
    r <- res[[design$run_id[i]]]
    data.frame(run_id = design$run_id[i], n_target = design$n_target[i],
               rep = design$rep[i], gen = 0:design$n_gen[i],
               He = r$He_trajectory)
  }))
  he_df <- do.call(rbind, lapply(
    split(he_df, he_df$run_id),
    function(s) { s$He_rel <- s$He / s$He[s$gen == 0]; s }
  ))
  he_agg <- aggregate(He_rel ~ n_target + gen, he_df, function(x)
    c(m = mean(x), lo = quantile(x, 0.1), hi = quantile(x, 0.9)))
  he_agg <- do.call(data.frame, he_agg)
  names(he_agg)[3:5] <- c("m", "lo", "hi")
  he_agg$label <- factor(paste0("N=", he_agg$n_target))

  print(
    ggplot(he_agg, aes(gen, m, colour = label, fill = label)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
      geom_line(linewidth = 0.8) +
      labs(title = "He retention", x = "Generation",
           y = expression(He/He[0]), colour = "N", fill = "N") +
      theme_bw()
  )
}
