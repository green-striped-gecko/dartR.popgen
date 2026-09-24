#' @name gl.TajimasD
#' @title Calculation of Tajima's D
#' @family demographic history
#' @description
#' This function calculates Tajima's D for each population in the
#' dartR/genlight object, with p-values from the beta distribution, the
#' standard normal distribution and, if \code{rep} is given, from
#' simulations of a neutral population with unlinked SNPs.
#' @details
#' D is calculated from the allele counts at each locus: pi (the average
#' number of pairwise differences) and Watterson's estimator use the number
#' of called sequences at each locus, so loci with missing calls are used
#' with their own sample size. The constants in the variance of D use the
#' mean sample size (reported as N, the mean number of genotyped individuals).
#' \cr
#' Pval.normal and Pval.beta follow Tajima (1989), which assumes that all
#' segregating sites lie on one non-recombining locus. SNPs from DArT data
#' are mostly unlinked, and for unlinked sites D varies much less under
#' neutrality, so these two p-values are too large (conservative) for such
#' data. sim_pval is computed against a neutral null of unlinked SNPs: for
#' each replicate, the derived-allele count of each of the S segregating sites
#' is drawn from the neutral distribution P(k) proportional to 1/k
#' (k = 1, ..., 2N - 1), and D is computed from the simulated sites.
#' sim_pval is the proportion of replicates with |D| at least as large as
#' the observed |D|.\cr
#' DArT loci are selected because they are polymorphic, which shifts the
#' site frequency spectrum towards intermediate frequencies and D upwards;
#' interpret D with this ascertainment in mind.
#' @param x Name of the genlight object containing the SNP data [required].
#' @param ms.path Not needed any more: the neutral simulation is done in R.
#' Accepted for compatibility and ignored [default NULL].
#' @param simulation.out Directory in which to save the simulated Tajima's D
#' values (one file per population), if rep is given [default NULL].
#' @param rep Number of simulated replicates for sim_pval. If NULL, no
#' simulation is run [default NULL].
#' @param seeds Seed for the random number generator of the simulation (the
#' first value is used) [default NULL, not set].
#' @param cleanup Not needed any more (no temporary files are written);
#' accepted for compatibility [default TRUE].
#' @param plot.dir Directory in which to save files [default as specified by
#' the global working directory or tempdir()].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot.file Name for the RDS binary file to save (base name only,
#'  exclude extension) [default NULL].
#' @param plot_theme Theme of the plot [default theme_dartR()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return
#' A data frame with, for each population, pi, the number of segregating
#' sites S, D, Pval.normal, Pval.beta, N (mean number of genotyped
#' individuals), theta_per_site (pi divided by the number of loci, not by
#' the number of base pairs) and, if rep is given, sim_pval. A plot of D by
#' population and, if rep is given, the simulated distributions are shown.
#' @export
#' @importFrom stringr str_split_fixed
#' @importFrom terra split
#' @references
#' \itemize{
#' \item Tajima, F. (1989). Statistical method for testing the neutral
#' mutation hypothesis by DNA polymorphism. Genetics, 123(3), 585-595.
#' \item Fu, Y. X. (1995). Statistical properties of segregating sites.
#' Theoretical Population Biology, 48(2), 172-197.
#' \item Paradis, E. (2010). pegas: an R package for population genetics with
#' an integrated-modular approach. Bioinformatics, 26(3), 419-420.
#' }
#' @author Author(s): Renee Catullo. Custodian: Ching Ching Lau -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' Tajima <- gl.TajimasD(x=bandicoot.gl)
#' # with a neutral null of unlinked SNPs
#' \donttest{
#' Tajima <- gl.TajimasD(x=bandicoot.gl, rep=100, seeds=1)
#' }
#' @importFrom stats D pbeta pnorm sd

gl.TajimasD <- function(x,
                        ms.path = NULL,
                        simulation.out = NULL,
                        rep = NULL,
                        seeds = NULL,
                        cleanup = TRUE,
                        plot.dir = NULL,
                        plot.out = TRUE,
                        plot.file = NULL,
                        plot_theme = NULL,
                        verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  if (datatype == "SilicoDArT") {
    stop(error("Fatal Error: Detected Presence/Absence (SilicoDArT) data. Please provide a SNP dataset\n"))
  }
  if (!is.null(rep) && (!is.numeric(rep) || length(rep) != 1 || rep < 1)) {
    stop(error("Fatal Error: rep must be a positive number of replicates\n"))
  }
  if (!is.null(ms.path) && verbose >= 1) {
    cat(warn("  Warning: ms is no longer needed; ms.path is ignored and the neutral simulation is done in R\n"))
  }
  if (!is.null(ms.path) && is.null(rep)) {
    stop(error("Fatal Error: rep (number of simulated replicates) is required for the simulation\n"))
  }

  if (is.null(plot_theme)) {
    plot_theme <- theme_dartR()
  }

  population <- D  <- meane <- NULL

  # constants of Tajima's D for n sequences
  tajima_const <- function(n) {
    tmp <- 1:(n - 1)
    a1 <- sum(1 / tmp)
    a2 <- sum(1 / tmp^2)
    b1 <- (n + 1) / (3 * (n - 1))
    b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
    c1 <- b1 - 1 / a1
    c2 <- b2 - (n + 2) / (a1 * n) + a2 / a1^2
    list(a1 = a1, a2 = a2, e1 = c1 / a1, e2 = c2 / (a1^2 + a2))
  }

  # get Tajima's D (code adopted from Paradis, E. (2010) -- pegas package --
  # written by Renee Catullo)

  get_tajima_D <- function(x) {
    # allele counts and called individuals for every locus in every population
    allele_freqs <- utils.get.allele.freq(x, verbose = 0)

    #split each population
    allele_freqs_by_pop <- split(allele_freqs, allele_freqs$popn)

    get_tajima_D_for_one_pop <- function(af) {
      # per-locus number of called sequences and exact frequencies
      n_i <- af$nobs * 2
      p1 <- af$sum / n_i
      p2 <- 1 - p1
      h <- (n_i / (n_i - 1)) * (1 - (p1^2 + p2^2))
      pi <- sum(h[n_i > 1], na.rm = T)

      #number of segregating sites, ignoring missing data
      seg <- !is.na(p1) & p1 > 0 & p1 < 1
      S <- sum(seg)
      if (S == 0) {
        warning("No segregating sites")
      }

      n <- mean(n_i)
      k <- tajima_const(n)

      # Watterson's estimator with the sample size of each site
      a1_i <- vapply(n_i[seg], function(m) sum(1 / (1:(m - 1))), numeric(1))
      thetaW <- sum(1 / a1_i)

      # calculate D and do beta testing
      D <- (pi - thetaW) / sqrt(k$e1 * S + k$e2 * S * (S - 1))
      Dmin <- (2 / n - 1 / k$a1) / sqrt(k$e2)
      Dmax <- ((n / (2 * (n - 1))) - 1 / k$a1) / sqrt(k$e2)
      tmp1 <- 1 + Dmin * Dmax
      tmp2 <- Dmax - Dmin
      a <- -tmp1 * Dmax / tmp2
      b <- tmp1 * Dmin / tmp2
      p <- pbeta((D - Dmin) / tmp2, b, a)
      p <- ifelse(p < 0.5, 2 * p, 2 * (1 - p))

      data.frame(
        pi = pi,
        S = S,
        D = D,
        Pval.normal = 2 * pnorm(-abs(D)),
        Pval.beta = p,
        N = n / 2
      )
    }

    output <- do.call("rbind",
                      lapply(allele_freqs_by_pop, get_tajima_D_for_one_pop))
    data.frame(population = rownames(output),
               output,
               row.names = NULL)
  }

  tmp_tajD <- get_tajima_D(x)
  tmp_tajD$theta_per_site <- tmp_tajD$pi / x@n.loc

  # NEUTRAL NULL OF UNLINKED SNPs when rep is provided
  # the derived count k of an unlinked segregating site in n sequences has
  # P(k) proportional to 1/k under the standard neutral model
  simulate_D <- function(n, S, rep) {
    k <- tajima_const(n)
    kk <- 1:(n - 1)
    cnt <- matrix(sample(kk, S * rep, replace = TRUE, prob = 1 / kk),
                  nrow = S)
    p <- cnt / n
    pi <- colSums((n / (n - 1)) * 2 * p * (1 - p))
    (pi - S / k$a1) / sqrt(k$e1 * S + k$e2 * S * (S - 1))
  }

  if (!is.null(rep)) {
    if (!is.null(seeds)) {
      # use the seed without changing the user's random number stream
      if (exists(".Random.seed", envir = globalenv())) {
        old.seed <- get(".Random.seed", envir = globalenv())
        on.exit(assign(".Random.seed", old.seed, envir = globalenv()), add = TRUE)
      }
      set.seed(seeds[1])
    }
    tmp_tajD$sim_pval <- NA
    sim_sum <- data.frame(row.names = seq_len(rep))

    for (p in tmp_tajD$population) {
      row <- which(tmp_tajD$population == p)
      S <- tmp_tajD$S[row]
      nsam <- round(2 * tmp_tajD$N[row])
      if (S == 0 || nsam < 2) {
        sim_sum[[p]] <- NA_real_
        next
      }
      simD <- simulate_D(nsam, S, rep)
      sim_sum[[p]] <- simD
      tmp_tajD$sim_pval[row] <- mean(abs(simD) >= abs(tmp_tajD$D[row]))

      # if choose to output the simulation
      if (!is.null(simulation.out)) {
        write.table(
          simD,
          file.path(simulation.out, paste0("Sim_TajimasD_", p, ".txt")),
          row.names = F,
          col.names = F,
          quote = F
        )
      }
    }

    # plot Tajima's D distribution
    plot.list <- list()
    for (p in tmp_tajD$population) {
      plot.list[[p]] <- ggplot(sim_sum, aes(x = .data[[p]])) + geom_histogram(bins = 20) +
        geom_vline(xintercept = tmp_tajD[which(tmp_tajD$population == p), 'D'], colour =
                     "red") +
        plot_theme +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
        labs(subtitle = p)
    }

    sim_plot <- plot.list %>% purrr::map(function(x) {
      ggplot2::ggplot_gtable(ggplot2::ggplot_build(x))
    })

    maxWidth <- do.call(grid::unit.pmax,
                        purrr::map(sim_plot, function(x)
                          x$widths[2:3]))

    for (i in 1:length(sim_plot)) {
      sim_plot[[i]]$widths[2:3] <- maxWidth
    }

    sim_plot$bottom <- "Simulated Tajima's D"
    sim_plot$left <- "Count"
  }

  p1 <-
    ggplot(tmp_tajD, aes(x = population, y = D, fill = population)) +
    geom_bar(position = "dodge", stat = "identity", color = "black") +
    plot_theme +
    theme(axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +
    labs(fill = "Population", x = "Populations") +
    ggtitle("Tajima's D by Population")

  if (plot.out & exists("sim_plot")) {
    print(p1)
    do.call(gridExtra::grid.arrange, sim_plot)
  } else if (plot.out) {
    print(p1)
  }

  if (!is.null(plot.file) & exists("plot.list")) {
    tmp <- utils.plot.save(
      plot.list,
      dir = plot.dir,
      file = paste0(plot.file, "_distribution"),
      verbose = verbose
    )
    tmp2 <- utils.plot.save(
      p1,
      dir = plot.dir,
      file = paste0(plot.file, "_TajimasD"),
      verbose = verbose
    )
  } else if (!is.null(plot.file) & exists("plot.list") == F) {
    tmp <- utils.plot.save(
      p1,
      dir = plot.dir,
      file = paste0(plot.file, "_TajimasD"),
      verbose = verbose
    )
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(invisible(tmp_tajD))

}
