## ------------------------------------------------------------------
## Metapopulation panel durability simulation framework
##
## sim_metapop()        — one metapop WF simulation with emigration
## build_design_metapop() — design data.frame from parameter grid
## run_simulations_metapop() — runs all rows, stores in results list
##
## Uses run_panel_checks() and extract_performance() from the
## single-pop framework (source that script first or combine).
##
## Connectivity structure (default):
##   Cluster 1 (A-D):  chain, high gene flow
##   Cluster 2 (E-G):  chain, moderate gene flow
##   Bridge    (D-E):  weak link between clusters
##   Isolated  (H-J):  minimal or no gene flow
##
## Workflow:
##   design  <- build_design_metapop(...)
##   out     <- run_simulations_metapop(design)
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

## ===================== EMIGRATION MATRIX =====================

#' Build a symmetric integer emigration matrix
#'
#' @param n_pops    number of populations
#' @param pop_names character vector of population names
#' @param cl1_m     migrants per gen within cluster 1 (A-B-C-D chain)
#' @param cl2_m     migrants per gen within cluster 2 (E-F-G chain)
#' @param bridge_m  migrants per gen on the D-E bridge
#' @param iso_m     migrants between H-I (0 = fully isolated)
#' @return integer matrix [i,j] = migrants FROM pop j TO pop i
build_emi_matrix <- function(n_pops, pop_names,
                             cl1_m    = 3,
                             cl2_m    = 2,
                             bridge_m = 1,
                             iso_m    = 0) {
  emi <- matrix(0L, nrow = n_pops, ncol = n_pops,
                dimnames = list(pop_names, pop_names))
  link <- function(a, b, m) {
    emi[a, b] <<- as.integer(m)
    emi[b, a] <<- as.integer(m)
  }
  if (n_pops >= 4) { link("A","B",cl1_m); link("B","C",cl1_m); link("C","D",cl1_m) }
  if (n_pops >= 7) { link("E","F",cl2_m); link("F","G",cl2_m) }
  if (n_pops >= 5) link("D","E",bridge_m)
  if (n_pops >= 9 && iso_m > 0) link("H","I",iso_m)
  emi
}


## ===================== CORE SIMULATION =====================

#' True WF breeding within one population (multinomial resampling)
#' @param gl_pop   genlight of one population (post-emigration)
#' @param n_target constant population size
#' @param popname  population label
#' @param pop_levels all population levels (to preserve factor)
#' @return list: $gl, $Vk, $Ne_Vk
breed_pop_wf <- function(gl_pop, n_target, popname, pop_levels) {
  n <- nInd(gl_pop)
  if (n < 2) {
    warning(sprintf("Pop %s has %d individual(s) - cloning.", popname, n))
    gl_out <- gl_pop[sample(n, n_target, replace = TRUE), ]
    pop(gl_out) <- factor(rep(popname, n_target), levels = pop_levels)
    return(list(gl = gl_out, Vk = NA, Ne_Vk = NA))
  }

  gmat <- as.matrix(gl_pop)
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
                ind.names = paste0(popname, "_", seq_len(n_target)),
                loc.names = colnames(gmat), loc.all = gl_pop@loc.all)
  pop(gl_new) <- factor(rep(popname, n_target), levels = pop_levels)

  list(gl = gl_new, Vk = Vk, Ne_Vk = Ne_Vk)
}


#' Ne from He decay per population (harmonic mean of consecutive ratios)
estimate_Ne_He_metapop <- function(sim_gls, pop_names) {
  n_gen <- length(sim_gls) - 1
  ne_by_pop <- setNames(numeric(length(pop_names)), pop_names)

  for (pp in pop_names) {
    He_vec <- vapply(sim_gls, function(gl) {
      gl_p <- gl[pop(gl) == pp, ]
      af <- colMeans(as.matrix(gl_p), na.rm = TRUE) / 2
      mean(2 * af * (1 - af), na.rm = TRUE)
    }, numeric(1))

    r <- He_vec[-1] / He_vec[-(n_gen + 1)]
    ne <- 1 / (2 * (1 - r))
    valid <- is.finite(ne) & ne > 0
    ne_by_pop[pp] <- if (sum(valid) > 0) 1 / mean(1 / ne[valid]) else NA_real_
  }
  ne_by_pop
}


#' Run a metapopulation WF simulation with emigration
#'
#' @param n_pops     number of populations
#' @param n_per_pop  constant population size per deme
#' @param n_loci     number of loci
#' @param n_founders number of unrelated founders (total, shared pool)
#' @param n_burnin   burn-in generations (with emigration, builds structure)
#' @param n_gen      forward generations to track
#' @param emi_table  emigration matrix (from build_emi_matrix)
#' @param pop_names  population labels
#' @param verbose    0 = silent, 1 = summary, 2 = per-generation
#' @return list: $sim_gls, $Ne_Vk_by_pop, $Ne_He_by_pop, $Vk_by_pop,
#'               $emi_table, $pop_names
sim_metapop <- function(n_pops     = 10,
                        n_per_pop  = 30,
                        n_loci     = 1000,
                        n_founders = 20,
                        n_burnin   = 4,
                        n_gen      = 20,
                        emi_table  = NULL,
                        pop_names  = LETTERS[1:n_pops],
                        verbose    = 1) {

  if (is.null(emi_table))
    emi_table <- build_emi_matrix(n_pops, pop_names)

  ## 1. Founders - shared pool, split into populations
  ##    All pops descend from the same allele frequency base;
  ##    differentiation builds during burn-in via drift + emigration.
  n_total <- n_pops * n_founders
  founders_all <- gl.sim.Neconst(n_total, nlocs = n_loci, verbose = 0)
  pop(founders_all) <- factor(rep(pop_names, each = n_founders),
                              levels = pop_names)
  x <- gl.compliance.check(founders_all, verbose = 0)

  ## 2. Burn-in: emigrate + breed, builds differentiation
  for (b in seq_len(n_burnin + 1)) {
    if (verbose >= 2) cat(sprintf("  burn-in %d\n", b))

    ## emigration (skip for first expansion round)
    if (b > 1) x <- gl.sim.emigration(x, emi.table = emi_table)

    ## breed within each pop to target size
    pop_gls <- list()
    for (pp in pop_names) {
      gl_p <- x[pop(x) == pp, ]
      if (nInd(gl_p) == 0) {
        warning(sprintf("Burn-in %d: pop %s extinct.", b, pp)); next
      }
      pop_gls[[pp]] <- breed_pop_wf(gl_p, n_per_pop, pp, pop_names)$gl
    }
    x <- do.call(rbind, pop_gls)
    x <- gl.compliance.check(x, verbose = 0)
  }

  ## 3. Forward simulation
  sim_gls   <- vector("list", n_gen + 1)
  sim_gls[[1]] <- x
  names(sim_gls)[1] <- "gen_0"

  ## per-pop Ne tracking
  ne_vk_mat <- matrix(NA_real_, nrow = n_gen, ncol = n_pops,
                      dimnames = list(NULL, pop_names))
  vk_mat    <- ne_vk_mat

  for (g in seq_len(n_gen)) {
    if (verbose >= 2) cat(sprintf("  gen %d\n", g))

    ## emigration
    x <- gl.sim.emigration(x, emi.table = emi_table)

    ## breed within each pop
    pop_gls <- list()
    for (pp in pop_names) {
      gl_p <- x[pop(x) == pp, ]
      if (nInd(gl_p) == 0) {
        warning(sprintf("Gen %d: pop %s extinct.", g, pp)); next
      }
      res_p <- breed_pop_wf(gl_p, n_per_pop, pp, pop_names)
      pop_gls[[pp]] <- res_p$gl
      ne_vk_mat[g, pp] <- res_p$Ne_Vk
      vk_mat[g, pp]    <- res_p$Vk
    }
    x <- do.call(rbind, pop_gls)
    x <- gl.compliance.check(x, verbose = 0)
    sim_gls[[g + 1]] <- x
    names(sim_gls)[g + 1] <- paste0("gen_", g)
  }

  ## 4. Ne estimates
  ## Vk-based: harmonic mean across generations, per pop
  Ne_Vk_by_pop <- apply(ne_vk_mat, 2, function(col) {
    v <- col[is.finite(col) & col > 0]
    if (length(v) == 0) NA_real_ else 1 / mean(1 / v)
  })

  ## He-decay-based: per pop
  Ne_He_by_pop <- estimate_Ne_He_metapop(sim_gls, pop_names)

  Vk_by_pop <- colMeans(vk_mat, na.rm = TRUE)

  if (verbose >= 1) {
    cat(sprintf("  %d pops x %d ind | %d loci | %d gen\n",
                n_pops, n_per_pop, n_loci, n_gen))
    for (pp in pop_names)
      cat(sprintf("    %s: Ne_Vk=%.1f  Ne_He=%.1f  Vk=%.2f\n",
                  pp, Ne_Vk_by_pop[pp], Ne_He_by_pop[pp], Vk_by_pop[pp]))
  }

  list(sim_gls       = sim_gls,
       Ne_Vk_by_pop  = Ne_Vk_by_pop,
       Ne_He_by_pop  = Ne_He_by_pop,
       Vk_by_pop     = Vk_by_pop,
       emi_table     = emi_table,
       pop_names     = pop_names)
}


## ===================== DESIGN + EXECUTION =====================

#' Build a design data.frame for metapopulation simulations
#'
#' @param n_pops     number of populations
#' @param n_per_pop  vector of per-pop sizes to sweep
#' @param n_loci     vector of locus counts
#' @param n_reps     replicates per combination
#' @param n_founders founders (total, shared pool)
#' @param n_burnin   burn-in generations
#' @param n_gen      forward generations
#' @param cl1_m      cluster 1 migration rate(s) to sweep
#' @param cl2_m      cluster 2 migration rate
#' @param bridge_m   bridge migration rate
#' @param iso_m      isolated migration rate
#' @return data.frame with run_id, parameters, result placeholders
build_design_metapop <- function(n_pops     = 10,
                                 n_per_pop  = 30,
                                 n_loci     = 1000,
                                 n_reps     = 3,
                                 n_founders = 20,
                                 n_burnin   = 4,
                                 n_gen      = 20,
                                 cl1_m      = 3,
                                 cl2_m      = 2,
                                 bridge_m   = 1,
                                 iso_m      = 0) {

  grid <- expand.grid(n_per_pop = n_per_pop,
                      n_loci    = n_loci,
                      cl1_m     = cl1_m,
                      rep       = seq_len(n_reps),
                      stringsAsFactors = FALSE)
  grid <- grid[order(grid$n_per_pop, grid$n_loci, grid$cl1_m, grid$rep), ]

  grid$run_id     <- sprintf("M%03d", seq_len(nrow(grid)))
  grid$n_pops     <- n_pops
  grid$n_founders <- n_founders
  grid$n_burnin   <- n_burnin
  grid$n_gen      <- n_gen
  grid$cl2_m      <- cl2_m
  grid$bridge_m   <- bridge_m
  grid$iso_m      <- iso_m
  grid$Ne_Vk_mean <- NA_real_
  grid$Ne_He_mean <- NA_real_
  grid$Vk_mean    <- NA_real_
  grid$sim_done   <- FALSE
  grid$panel_done <- FALSE

  rownames(grid) <- grid$run_id
  grid[, c("run_id", "n_pops", "n_per_pop", "n_loci", "cl1_m", "rep",
           "n_founders", "n_burnin", "n_gen",
           "cl2_m", "bridge_m", "iso_m",
           "Ne_Vk_mean", "Ne_He_mean", "Vk_mean",
           "sim_done", "panel_done")]
}


#' Run all metapopulation simulations in the design
#'
#' @param design  data.frame from build_design_metapop()
#' @param verbose 0/1/2
#' @return list: $design (updated), $results (named list)
run_simulations_metapop <- function(design, verbose = 1) {
  results <- list()

  for (i in seq_len(nrow(design))) {
    d <- design[i, ]
    if (d$sim_done) {
      cat(sprintf("[%s] already done, skipping\n", d$run_id)); next
    }

    pop_names <- LETTERS[1:d$n_pops]
    emi_table <- build_emi_matrix(d$n_pops, pop_names,
                                  cl1_m = d$cl1_m, cl2_m = d$cl2_m,
                                  bridge_m = d$bridge_m, iso_m = d$iso_m)

    if (verbose >= 1)
      cat(sprintf("\n[%s] pops=%d, N/pop=%d, L=%d, cl1=%d, rep=%d\n",
                  d$run_id, d$n_pops, d$n_per_pop, d$n_loci, d$cl1_m, d$rep))

    res <- sim_metapop(
      n_pops     = d$n_pops,
      n_per_pop  = d$n_per_pop,
      n_loci     = d$n_loci,
      n_founders = d$n_founders,
      n_burnin   = d$n_burnin,
      n_gen      = d$n_gen,
      emi_table  = emi_table,
      pop_names  = pop_names,
      verbose    = verbose
    )

    results[[d$run_id]] <- res
    design$Ne_Vk_mean[i] <- mean(res$Ne_Vk_by_pop, na.rm = TRUE)
    design$Ne_He_mean[i] <- mean(res$Ne_He_by_pop, na.rm = TRUE)
    design$Vk_mean[i]    <- mean(res$Vk_by_pop, na.rm = TRUE)
    design$sim_done[i]   <- TRUE
  }

  list(design = design, results = results)
}


## ===================== EXAMPLE USAGE =====================

if (FALSE) {  # set to TRUE to run

  ## 1. Build design
  design_mp <- build_design_metapop(
    n_pops    = 10,
    n_per_pop = 30,
    n_loci    = 1000,
    n_reps    = 3,
    n_founders = 20,
    n_burnin  = 4,
    n_gen     = 20,
    cl1_m     = 3,
    cl2_m     = 2,
    bridge_m  = 1,
    iso_m     = 0
  )
  print(design_mp)

  ## 2. Run simulations
  out_mp <- run_simulations_metapop(design_mp, verbose = 1)
  design_mp <- out_mp$design
  res_mp    <- out_mp$results

  ## Ne summary
  print(design_mp[, c("run_id", "n_per_pop", "cl1_m", "rep",
                       "Ne_Vk_mean", "Ne_He_mean")])

  ## 3. Run panel checks - multi-population parameters
  ##    (requires run_panel_checks from sim_single_pop_panel.R)
  out2_mp <- run_panel_checks(
    design_mp, res_mp,
    nl        = 50,
    method    = "combined",
    parameter = c("Fst", "He", "assignment", "hybridisation",
                  "drift_resistance"),
    n_iter    = 200,
    target    = 0.9
  )
  design_mp <- out2_mp$design
  res_mp    <- out2_mp$results

  ## 4. Extract and plot
  perf_mp <- extract_performance(design_mp, res_mp)

  perf_agg <- aggregate(performance ~ n_per_pop + generation + parameter,
                        perf_mp, function(x)
                          c(m = mean(x), lo = quantile(x, 0.1),
                            hi = quantile(x, 0.9)))
  perf_agg <- do.call(data.frame, perf_agg)
  names(perf_agg)[4:6] <- c("m", "lo", "hi")

  print(
    ggplot(perf_agg, aes(generation, m, colour = parameter)) +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = parameter),
                  alpha = 0.1, colour = NA) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 0.9, linetype = 3) +
      labs(title = "Metapop panel durability",
           x = "Generation", y = "Performance") +
      theme_bw()
  )

  ## Fst trajectory diagnostic
  fst_rows <- list()
  for (rid in names(res_mp)) {
    r <- res_mp[[rid]]
    for (g in seq(0, design_mp$n_gen[1], by = 5)) {
      gl_g <- r$sim_gls[[paste0("gen_", g)]]
      if (is.null(gl_g)) next
      
      fst <- gl.report.fstat(gl_g,  verbose = 0, plot.display = F)$Stat_matrices$Fstp
      fst_rows[[length(fst_rows) + 1]] <- data.frame(
        run_id = rid, gen = g, Fst = mean(fst[lower.tri(fst)]))
    }
  }
  fst_df <- do.call(rbind, fst_rows)
  print(
    ggplot(fst_df, aes(gen, Fst)) +
      geom_line(aes(group = run_id), alpha = 0.5) +
      geom_smooth(se = FALSE) +
      labs(title = "Fst over generations",
           x = "Generation", y = "Fst") +
      theme_bw()
  )
}




