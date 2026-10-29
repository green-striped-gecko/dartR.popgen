#' @name gl.opt.ind.all
# HEADER INFORMATION -------------------
#' @title Finds the order of individuals that captures the most alleles
#' @family population genetics
#'
#' @description
#' Ranks individuals so that each one added captures as many new alleles as
#' possible (greedy allele-capture order). Use it to decide how many
#' individuals, and from which populations, are needed to keep a set
#' proportion of the alleles in a sample. Typical uses are choosing founders
#' for a captive breeding program, picking individuals for translocation or
#' population mixing, or setting priorities for sequencing or banking.
#'
#' @param x Name of the genlight object containing SNP data [required].
#' @param first How to choose the first individual. 'optimal' starts with the
#' individual carrying the most alleles and runs once; 'random' starts with a
#' random individual and is repeated \code{bootstrap} times; 'all' runs once
#' for each individual as the starting individual [default 'optimal'].
#' @param bootstrap Number of runs when \code{first = 'random'}. Ignored for
#' 'optimal' and 'all' [default 100].
#' @param percent Proportion of alleles to retain (between 0 and 1). Used to
#' count how many individuals (in total and per population) are needed
#' [default 0.95].
#' @param ncores Number of cores used for 'random' and 'all'. Values above 1
#' use \code{future::multisession}, which works on Windows, Linux and macOS
#' [default 1].
#' @param plot.display Specify if the plot is to be displayed in the graphics
#' window [default TRUE].
#' @param plot.theme User specified theme [default theme_dartR()].
#' @param plot.colors Vector with two colors: the first for lines and borders,
#' the second for fills [default gl.colors(2)].
#' @param plot.colors.pop A color palette for populations or a vector with as
#' many colors as there are populations [default gl.colors("dis")].
#' @param plot.dir Directory to save the plot RDS file [default as specified
#' by the global working directory or tempdir()].
#' @param plot.file Filename (minus extension) for the RDS plot file
#' [required for plot save].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' \strong{How alleles are counted}
#'
#' At each locus an individual carries one distinct allele if homozygous
#' (0 or 2) and two if heterozygous (1). The function picks the individual
#' that adds the most alleles not yet captured, records them as captured and
#' repeats. A heterozygote captures both alleles at a locus, so the locus is
#' then complete. A homozygote captures one allele; afterwards only
#' individuals that carry the other allele can add anything at that locus.
#' The search stops when all alleles are captured. Individuals that would add
#' no new alleles are therefore not listed in the output.
#'
#' Missing genotypes (NA) are treated as unknown and add nothing. Individuals
#' with many missing values will look less valuable than they may be, so
#' filter for call rate first (e.g. \code{gl.filter.callrate(x, method = "ind")}).
#'
#' \strong{Choice of the first individual}
#'
#' The greedy order depends on where it starts. 'optimal' gives one ranking.
#' 'random' and 'all' show how stable the result is when the start changes,
#' and give a spread for the number of individuals needed per population.
#' When several individuals add the same number of alleles, one is picked at
#' random, so use \code{set.seed()} for reproducible results. Results are
#' reproducible for the same seed whatever the value of \code{ncores}.
#'
#' \strong{Individuals needed}
#'
#' The individuals needed to reach \code{percent} are the first k individuals
#' in the order, where k is the first position at which the cumulative
#' proportion of alleles reaches \code{percent}. This includes the individual
#' that crosses the threshold.
#'
#' \strong{Plot}
#'
#' The upper panel shows the cumulative proportion of alleles captured as
#' individuals are added (one curve per run for 'random' and 'all', with the
#' median in bold). The dashed line marks \code{percent}. The lower panel
#' shows the number of individuals needed from each population: bars for
#' 'optimal' (with the population size shown in light grey behind), and
#' boxplots across runs for 'random' and 'all'.
#'
#' If a plot.file is given, the ggplot arising from this function is saved as
#' an "RDS" binary file using saveRDS(); it can be reloaded with readRDS(). If
#' a plot directory (plot.dir) is specified, the ggplot binary is saved to
#' that directory; otherwise to the tempdir().
#'
#' @return A list (returned invisibly when \code{verbose} is 0) with:
#' \itemize{
#' \item \code{table} data.frame with one row per individual and run: run,
#' start (id of the first individual), order (position in the order),
#' ind.index (row number in x), id, pop, nall.added (new alleles added),
#' nall.cum (cumulative alleles), prop.cum (cumulative proportion of all
#' alleles) and needed (TRUE if the individual is among those needed to reach
#' \code{percent}).
#' \item \code{summary} data.frame with one row per population: n.ind
#' (individuals in x), and the number of individuals needed to reach
#' \code{percent} as mean_n, median_n, sd, min_n and max_n across runs, plus
#' prop (share of the total needed). For 'optimal' mean, median, min and max
#' are identical and sd is NA.
#' \item \code{by.run} matrix with one row per run and one column per
#' population: number of individuals needed to reach \code{percent}.
#' \item \code{n.needed} vector with the total number of individuals needed
#' in each run.
#' \item \code{plot} the ggplot (patchwork) object.
#' \item \code{settings} list of the input settings and data dimensions
#' (first, bootstrap, percent, ncores, n.runs, n.ind, n.loc, n.pop,
#' total.alleles).
#' }
#'
#' @author Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \donttest{
#' test <- possums.gl[1:60, 1:200]
#' # one run, starting with the individual with the most alleles
#' res <- gl.opt.ind.all(test, first = "optimal", percent = 0.95)
#' res$summary
#' head(res$table)
#' # 20 runs, each starting with a random individual
#' set.seed(42)
#' res2 <- gl.opt.ind.all(test, first = "random", bootstrap = 20)
#' res2$summary
#' # save the table
#' # write.csv(res2$table, file = "opt.ind.csv", row.names = FALSE)
#' }
#' # small example that runs quickly
#' res <- gl.opt.ind.all(possums.gl[1:20, 1:50], plot.display = FALSE)
#'
#' @seealso \code{\link{gl.filter.callrate}}
#' @importFrom future plan multisession sequential availableCores
#' @importFrom furrr future_map furrr_options
#' @importFrom stats median sd aggregate
#' @export

gl.opt.ind.all <- function(x,
                           first = "optimal",
                           bootstrap = 100,
                           percent = 0.95,
                           ncores = 1,
                           plot.display = TRUE,
                           plot.theme = theme_dartR(),
                           plot.colors = NULL,
                           plot.colors.pop = NULL,
                           plot.dir = NULL,
                           plot.file = NULL,
                           verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Elastic Elapid",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  first <- tolower(first)
  if (length(first) != 1 || !(first %in% c("optimal", "random", "all"))) {
    stop(error("  Fatal Error: parameter first must be one of 'optimal', 'random' or 'all'.\n"))
  }
  
  if (!is.numeric(percent) || length(percent) != 1 ||
      percent <= 0 || percent > 1) {
    stop(error("  Fatal Error: parameter percent must be a single value larger than 0 and not larger than 1.\n"))
  }
  
  if (!is.numeric(ncores) || length(ncores) != 1 || ncores < 1) {
    stop(error("  Fatal Error: parameter ncores must be a single number of 1 or more.\n"))
  }
  ncores <- as.integer(ncores)
  maxcores <- future::availableCores()
  if (ncores > maxcores) {
    if (verbose >= 2) {
      cat(warn("  Warning: ncores set to", ncores, "but only", maxcores,
               "cores available. ncores set to", maxcores, "\n"))
    }
    ncores <- as.integer(maxcores)
  }
  
  if (first == "random") {
    if (is.null(bootstrap) || !is.numeric(bootstrap) ||
        length(bootstrap) != 1 || bootstrap < 1) {
      stop(error("  Fatal Error: parameter bootstrap must be a single number of 1 or more when first = 'random'.\n"))
    }
    bootstrap <- as.integer(bootstrap)
  } else {
    if (!is.null(bootstrap) && !missing(bootstrap) && verbose >= 2) {
      cat(warn("  Warning: parameter bootstrap is ignored when first =", first, "\n"))
    }
    bootstrap <- NA
  }
  
  if (first == "optimal" && ncores > 1 && verbose >= 2) {
    cat(warn("  Warning: first = 'optimal' runs once, ncores is not used.\n"))
  }
  
  if (nInd(x) < 2) {
    stop(error("  Fatal Error: at least two individuals are needed.\n"))
  }
  
  if (is.null(pop(x)) || length(pop(x)) != nInd(x)) {
    if (verbose >= 2) {
      cat(warn("  Warning: no population assignments found, all individuals assigned to pop1.\n"))
    }
    pop(x) <- rep("pop1", nInd(x))
  }
  
  if (verbose == 0) {
    plot.display <- FALSE
  }
  
  # SET COLOURS
  if (is.null(plot.colors)) {
    plot.colors <- gl.colors(2, verbose = 0)
  } else if (length(plot.colors) > 2) {
    if (verbose >= 2) {
      cat(warn("  Warning: More than 2 colors specified, only the first 2 are used\n"))
    }
    plot.colors <- plot.colors[1:2]
  } else if (length(plot.colors) == 1) {
    plot.colors <- rep(plot.colors, 2)
  }
  
  if (is.null(plot.colors.pop)) {
    plot.colors.pop <- gl.colors("dis", verbose = 0)
  }
  if (is.function(plot.colors.pop)) {
    plot.colors.pop <- plot.colors.pop(nPop(x))
  }
  if (length(plot.colors.pop) < nPop(x)) {
    if (verbose >= 2) {
      cat(warn("  Warning: fewer colors than populations in plot.colors.pop, default colors used.\n"))
    }
    plot.colors.pop <- gl.colors("dis", verbose = 0)(nPop(x))
  }
  plot.colors.pop <- plot.colors.pop[1:nPop(x)]
  names(plot.colors.pop) <- popNames(x)
  
  if (verbose >= 3 && sum(is.na(as.matrix(x))) > 0) {
    cat(warn("  Warning: data contain missing values. Missing genotypes add no alleles, consider filtering for call rate first.\n"))
  }
  
  # DO THE JOB
  
  xx <- as.matrix(x)
  dimnames(xx) <- NULL
  ids <- indNames(x)
  pops <- factor(pop(x), levels = popNames(x))
  
  # starting individuals for each run (NA = individual with most alleles)
  starts <- switch(first,
                   optimal = NA_integer_,
                   random = sample.int(nInd(x), bootstrap, replace = TRUE),
                   all = seq_len(nInd(x)))
  nruns <- length(starts)
  
  if (verbose >= 2) {
    cat(report("  Finding allele-capture order for", nInd(x), "individuals and",
               nLoc(x), "loci;", nruns, ifelse(nruns == 1, "run", "runs"),
               "(first =", paste0("'", first, "')"), "\n"))
  }
  
  if (nruns == 1) {
    runs <- list(utils.opt.ind.order(xx, start = starts[1]))
  } else {
    # keep the user's plan and restore it on exit
    oplan <- future::plan()
    on.exit(future::plan(oplan), add = TRUE)
    if (ncores > 1) {
      future::plan(future::multisession, workers = ncores)
    } else {
      future::plan(future::sequential)
    }
    # Pass the helper to the workers as an argument and detach it (and the
    # wrapper) from the package namespace and this function's frame. Workers
    # then need neither dartR.popgen (important with devtools::load_all() or
    # an older installed version) nor a copy of the genlight object.
    ord.fun <- utils.opt.ind.order
    environment(ord.fun) <- baseenv()
    worker <- function(s, f, m) f(m, start = s)
    environment(worker) <- baseenv()
    runs <- furrr::future_map(starts,
                              worker,
                              f = ord.fun,
                              m = xx,
                              .options = furrr::furrr_options(seed = TRUE,
                                                              packages = character(0)))
  }
  
  # build the output table
  total.alleles <- sum(runs[[1]]$nalltot)
  tabs <- lapply(seq_len(nruns), function(r) {
    o <- runs[[r]]$orderind
    na <- runs[[r]]$nalltot
    cum <- cumsum(na)
    prop <- cum / sum(na)
    k <- which(prop >= percent - 1e-12)[1]
    data.frame(run = r,
               start = ids[o[1]],
               order = seq_along(o),
               ind.index = o,
               id = ids[o],
               pop = pops[o],
               nall.added = na,
               nall.cum = cum,
               prop.cum = round(prop, 4),
               needed = seq_along(o) <= k,
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, tabs)
  rownames(tab) <- NULL
  
  # individuals needed per population and run (populations not needed count 0)
  by.run <- do.call(rbind, lapply(tabs, function(d) {
    table(factor(d$pop[d$needed], levels = popNames(x)))
  }))
  by.run <- matrix(by.run, nrow = nruns,
                   dimnames = list(paste0("run", seq_len(nruns)), popNames(x)))
  n.needed <- rowSums(by.run)
  names(n.needed) <- rownames(by.run)
  
  summ <- data.frame(pop = popNames(x),
                     n.ind = as.integer(table(pops)),
                     mean_n = round(colMeans(by.run), 2),
                     median_n = apply(by.run, 2, stats::median),
                     sd = if (nruns > 1) round(apply(by.run, 2, stats::sd), 2) else NA_real_,
                     min_n = apply(by.run, 2, min),
                     max_n = apply(by.run, 2, max),
                     prop = round(colMeans(by.run) / sum(colMeans(by.run)), 3),
                     row.names = NULL,
                     stringsAsFactors = FALSE)
  
  # PRINTOUT
  if (verbose >= 3) {
    cat("  Total number of alleles:", total.alleles, "\n")
    if (nruns == 1) {
      cat("  Individuals needed to retain", percent * 100, "% of alleles:",
          n.needed, "of", nInd(x), "\n")
    } else {
      cat("  Individuals needed to retain", percent * 100, "% of alleles: median",
          stats::median(n.needed), "(range", min(n.needed), "-", max(n.needed),
          ") of", nInd(x), "\n")
    }
    cat("  Individuals needed per population:\n")
    print(summ)
  }
  
  # PLOT
  run <- order <- prop.cum <- pop <- n <- n.ind <- NULL # for R CMD check
  
  plot.title <- paste0("Allele capture (first = '", first, "', ",
                       nruns, ifelse(nruns == 1, " run", " runs"), ")")
  
  if (nruns == 1) {
    k <- n.needed[1]
    p1 <- ggplot(tab, aes(x = order, y = prop.cum)) +
      geom_line(color = "grey60") +
      geom_point(aes(color = pop), size = 2) +
      geom_hline(yintercept = percent, linetype = "dashed",
                 color = plot.colors[1]) +
      geom_vline(xintercept = k, linetype = "dashed", color = plot.colors[1]) +
      scale_color_manual(values = plot.colors.pop, drop = FALSE) +
      labs(x = "Number of individuals", y = "Proportion of alleles captured",
           color = "Population", title = plot.title,
           subtitle = paste0(k, " of ", nInd(x), " individuals needed to retain ",
                             percent * 100, "% of ", total.alleles, " alleles")) +
      plot.theme
    
    dfbar <- data.frame(pop = factor(summ$pop, levels = popNames(x)),
                        n = summ$median_n, n.ind = summ$n.ind)
    p2 <- ggplot(dfbar, aes(x = pop)) +
      geom_col(aes(y = n.ind), fill = "grey90", color = "grey70") +
      geom_col(aes(y = n, fill = pop), color = plot.colors[1]) +
      scale_fill_manual(values = plot.colors.pop, drop = FALSE) +
      labs(x = "Population", y = "Individuals needed",
           subtitle = "Needed (colored) out of available (grey)") +
      plot.theme +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    medcurve <- stats::aggregate(prop.cum ~ order, data = tab, FUN = stats::median)
    p1 <- ggplot(tab, aes(x = order, y = prop.cum, group = run)) +
      geom_line(color = plot.colors[2], alpha = max(0.05, min(0.5, 10 / nruns))) +
      geom_line(data = medcurve, aes(x = order, y = prop.cum),
                inherit.aes = FALSE, color = plot.colors[1], linewidth = 1) +
      geom_hline(yintercept = percent, linetype = "dashed",
                 color = plot.colors[1]) +
      labs(x = "Number of individuals", y = "Proportion of alleles captured",
           title = plot.title,
           subtitle = paste0("Individuals needed to retain ", percent * 100,
                             "% of ", total.alleles, " alleles: median ",
                             stats::median(n.needed), " (range ", min(n.needed),
                             "-", max(n.needed), ") of ", nInd(x))) +
      plot.theme
    
    dfbox <- data.frame(pop = factor(rep(popNames(x), each = nruns),
                                     levels = popNames(x)),
                        n = as.vector(by.run))
    p2 <- ggplot(dfbox, aes(x = pop, y = n, fill = pop)) +
      geom_boxplot(color = plot.colors[1], outlier.shape = NA) +
      geom_jitter(width = 0.15, height = 0.1, size = 0.8, alpha = 0.4) +
      scale_fill_manual(values = plot.colors.pop, drop = FALSE) +
      labs(x = "Population", y = "Individuals needed",
           subtitle = paste0("Individuals needed per population to retain ",
                             percent * 100, "% of alleles")) +
      plot.theme +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  p3 <- (p1 / p2) + plot_layout(heights = c(3, 2))
  
  if (plot.display) {
    suppressWarnings(print(p3))
  }
  
  # Optionally save the plot
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p3,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }
  
  out <- list(table = tab,
              summary = summ,
              by.run = by.run,
              n.needed = n.needed,
              plot = p3,
              settings = list(first = first,
                              bootstrap = bootstrap,
                              percent = percent,
                              ncores = ncores,
                              n.runs = nruns,
                              n.ind = nInd(x),
                              n.loc = nLoc(x),
                              n.pop = nPop(x),
                              total.alleles = total.alleles))
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  if (verbose == 0) {
    return(invisible(out))
  }
  return(out)
}

# Internal helper: greedy allele-capture order for one run.
# xx: genotype matrix (individuals x loci, coded 0/1/2/NA).
# start: row index of the first individual, or NA to start with the
#        individual carrying the most alleles.
# Kept outside gl.opt.ind.all so that parallel workers only receive the
# matrix, not the genlight object.
utils.opt.ind.order <- function(xx, start = NA) {
  n <- nrow(xx)
  orderind <- integer(0)
  nalltot <- numeric(0)
  
  for (i in seq_len(n)) {
    # alleles still to be captured: 1 for homozygotes, 2 for heterozygotes
    nall_ind <- rowSums((xx %% 2) + 1, na.rm = TRUE)
    if (i == 1 && !is.na(start)) {
      nextind <- as.integer(start)
    } else {
      best <- which(nall_ind == max(nall_ind))
      nextind <- if (length(best) > 1) best[sample.int(length(best), 1)] else best
    }
    if (i > 1 && nall_ind[nextind] == 0) break
    
    g <- xx[nextind, ]
    xx[nextind, ] <- NA
    # heterozygote: both alleles captured, locus done
    xx[, which(g == 1)] <- NA
    # homozygote: same genotype adds nothing; heterozygotes now carry only
    # the other allele (recoded as the other homozygote)
    hom <- which(g == 0 | g == 2)
    if (length(hom) > 0) {
      sub <- xx[, hom, drop = FALSE]
      gm <- matrix(g[hom], nrow = n, ncol = length(hom), byrow = TRUE)
      sub[which(sub == gm)] <- NA
      het <- which(sub == 1)
      sub[het] <- 2 - gm[het]
      xx[, hom] <- sub
    }
    
    orderind <- c(orderind, nextind)
    nalltot <- c(nalltot, nall_ind[nextind])
    if (all(is.na(xx))) break
  }
  list(orderind = orderind, nalltot = nalltot)
}
