#' @name gl.ld.haplotype
#' @title Visualize patterns of linkage disequilibrium and identification of
#' haplotypes
#' @description
#' This function plots a Linkage disequilibrium (LD) heatmap, where the colour
#' shading indicates the strength of LD. Chromosome positions (Mbp) are shown on
#'  the horizontal axis, and haplotypes appear as triangles and delimited by
#'  dark yellow vertical lines. Numbers identifying each haplotype are shown in
#'  the upper part of the plot.
#'
#'  The heatmap also shows heterozygosity for each SNP.
#'
#'  When \code{haplo_id = TRUE}, the function identifies haplotypes as runs of
#'  adjacent SNPs whose pairwise LD is at least \code{ld_threshold_haplo} and
#'  that contain at least \code{min_snps} SNPs. With the default
#'  \code{haplo_id = FALSE} no haplotypes are identified and the returned table
#'  is empty.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param pop_name Name of the population to analyse. If NULL all the
#' populations are analysed [default NULL].
#' @param chrom_name Name of the chromosome to analyse. If NULL all the
#' chromosomes are analysed [default NULL].
#' @param ld_max_pairwise Maximum distance in number of base pairs at which LD
#' should be calculated [default 10000000].
#' @param maf Minor allele frequency (by population) threshold to filter out
#' loci. If a value > 1 is provided it will be interpreted as MAC (i.e. the
#'  minimum number of times an allele needs to be observed) [default 0.05].
#' @param ld_stat The LD measure to be calculated: "LLR", "OR", "Q", "Covar",
#'   "D.prime", "R.squared", and "R". See \code{\link[snpStats]{ld}}
#'    (package snpStats) for details [default "R.squared"].
#' @param ind.limit Minimum number of individuals that a population must have
#' to be analysed. Populations with fewer individuals are skipped, with a
#' warning; the function stops if every population is skipped
#' [default 10].
#' @param haplo_id Whether to identify haplotypes [default FALSE].
#' @param min_snps Minimum number of SNPs that a haplotype must contain to be
#' called [default 10].
#' @param ld_threshold_haplo Minimum LD between adjacent SNPs to call a
#' haplotype [default 0.5].
#' @param plot_het Whether to plot heterozygosity [default TRUE].
#' @param snp_pos Whether to plot SNP positions. The SNP position track is
#' drawn only when no haplotypes are identified [default TRUE].
#' @param target.snp1 Vector of position(s) of target SNP(s) in base pairs;
#' the closest SNP to each position is highlighted in the SNP position track
#' [default NULL].
#' @param target.snp2 Vector of position(s) of target SNP(s) in base pairs
#' [default NULL].
#' @param target.snp3 Vector of position(s) of target SNP(s) in base pairs
#' [default NULL].
#' @param col.all Color of line indicating position for all SNPs
#' [default "black"].
#' @param col.target1 Color of line indicating position for target.snp1
#' [default "green"].
#' @param col.target2 Color of line indicating position for target.snp2
#'  [default "blue"].
#' @param col.target3 Color of line indicating position for target.snp3
#'  [default "red"].
#' @param coordinates A vector of two elements with the start and end
#' coordinates in base pairs to which restrict the
#' analysis e.g. c(1,1000000) [default NULL].
#' @param color_haplo Color palette for haplotype plot. Options are: "magma",
#' "inferno", "plasma", "viridis", "cividis", "rocket", "mako" and "turbo"
#'  [default "viridis"].
#' @param color_het Color for heterozygosity [default "deeppink"].
#' @param plot.out Specify if heatmap plot is to be produced [default TRUE].
#' @param plot.dir Directory in which to save files [default tempdir(), or the
#' directory set with gl.set.wd()].
#' @param plot.save Whether to save each plot in pdf format, as
#' <population>_<chromosome>.pdf in plot.dir [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' The information for SNP's position should be stored in the genlight accessor
#'   "@@position" and the SNP's chromosome name in the accessor "@@chromosome"
#'   (see examples). The function will then calculate LD within each chromosome.
#'   Chromosomes with fewer than four SNPs after filtering are skipped.
#'
#' The output of the function includes a table with the haplotypes
#'  that were identified and their location.
#'
#'  Colors of the heatmap (\code{color_haplo}) are based on the function
#'    \code{\link[viridis]{scale_fill_viridis}} from  package \code{viridis}.
#'    Other color palettes options are "magma", "inferno", "plasma", "viridis",
#'     "cividis", "rocket", "mako" and "turbo".
#' @return A data frame with one row per haplotype identified: population,
#' chromosome, haplotype number, start and end (base pairs), the corresponding
#' x coordinates in the LD plot, midpoints and a label in Mbp. Returned
#' invisibly; empty when \code{haplo_id = FALSE} or no haplotype is found.
#' @family ld functions
#' @examples
#' require("dartR.data")
#' x <- platypus.gl
#' x <- gl.filter.callrate(x, threshold = 1)
#' # only the first 15 individuals because of speed during tests
#' x <- gl.keep.pop(x, pop.list = "TENTERFIELD")[1:15, ]
#' x$chromosome <- as.factor(x$other$loc.metrics$Chrom_Platypus_Chrom_NCBIv1)
#' x$position <- x$other$loc.metrics$ChromPos_Platypus_Chrom_NCBIv1
#' ld_res <- gl.ld.haplotype(x,
#'   chrom_name = "NC_041728.1_chromosome_1",
#'   ld_max_pairwise = 10000000
#' )
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @export

gl.ld.haplotype <- function(x,
                            pop_name = NULL,
                            chrom_name = NULL,
                            ld_max_pairwise = 10000000,
                            maf = 0.05,
                            ld_stat = "R.squared",
                            ind.limit = 10,
                            haplo_id = FALSE,
                            min_snps = 10,
                            ld_threshold_haplo = 0.5,
                            plot_het = TRUE,
                            snp_pos = TRUE,
                            target.snp1 = NULL,
                            target.snp2 = NULL,
                            target.snp3 = NULL,
                            col.all = "black",
                            col.target1 = "green",
                            col.target2 = "blue",
                            col.target3 = "red",
                            coordinates = NULL,
                            color_haplo = "viridis",
                            color_het = "deeppink",
                            plot.out = TRUE,
                            plot.save = FALSE,
                            plot.dir = NULL,
                            verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(
    func = funname,
    verbose = verbose
  )

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING

  # check if packages are installed
  for (pkg in c("snpStats", "sp", "raster", "scales", "viridis")) {
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      stop(error(
        "Package", pkg,
        "needed for this function to work. Please install it.\n"
      ))
    }
  }

  # chromosome and position information is required
  if (length(x@chromosome) != nLoc(x) || length(x@position) != nLoc(x) ||
    all(is.na(x@position))) {
    stop(error(
      "  SNP chromosome names and positions must be stored in x@chromosome",
      "and x@position (see examples).\n"
    ))
  }

  # DO THE JOB

  # bind variables used in ggplot aesthetics
  group <- lat <- long <- het <- position <- colorL <- y <- NULL
  start_ld_plot <- end_ld_plot <- xintercept <- NULL

  if (!is.null(chrom_name)) {
    chrom_missing <- setdiff(chrom_name, as.character(unique(x@chromosome)))
    if (length(chrom_missing) > 0) {
      stop(error(
        "  Chromosome(s)", paste(chrom_missing, collapse = ", "),
        "not found in x@chromosome.\n"
      ))
    }
    chrom_tmp <- which(as.character(x@chromosome) %in% chrom_name)
    x <- gl.keep.loc(x, loc.list = locNames(x)[chrom_tmp], verbose = 0)
  }

  if (!is.null(pop_name)) {
    pop_missing <- setdiff(pop_name, popNames(x))
    if (length(pop_missing) > 0) {
      stop(error(
        "  Population(s)", paste(pop_missing, collapse = ", "),
        "not found in x.\n"
      ))
    }
    x <- gl.keep.pop(x, pop.list = pop_name, verbose = 0)
  }

  if (!is.null(coordinates)) {
    if (verbose >= 2) {
      cat(report(
        "  Restricting the analysis from", coordinates[1], "to",
        coordinates[2], "base pairs\n"
      ))
    }
    loc_in_range <- which(x@position >= coordinates[1] &
      x@position <= coordinates[2])
    if (length(loc_in_range) == 0) {
      stop(error(
        "  No loci found between", coordinates[1], "and", coordinates[2],
        "base pairs.\n"
      ))
    }
    x <- gl.keep.loc(x, loc.list = locNames(x)[loc_in_range], verbose = 0)
  }

  x_list <- seppop(x)

  haplo_table <- data.frame(
    population = character(), chromosome = character(),
    haplotype = integer(), start = numeric(), end = numeric(),
    start_ld_plot = numeric(), end_ld_plot = numeric(),
    midpoint = numeric(), midpoint_ld_plot = numeric(),
    labels = character(), stringsAsFactors = FALSE
  )

  chr_list <- as.character(unique(x@chromosome))

  # why each skipped population was skipped; reported after the loop so
  # callers that only see conditions (e.g. DartRShiny) learn of it
  skipped <- character(0)

  for (pop_n in seq_along(x_list)) {
    pop_ld <- x_list[[pop_n]]
    pop_label <- popNames(pop_ld)

    if (nInd(pop_ld) < ind.limit) {
      if (verbose >= 1) {
        cat(warn(
          "  Skipping population", pop_label,
          "from analysis because it has fewer than", ind.limit,
          "individuals.\n"
        ))
      }
      skipped <- c(skipped, paste0(
        pop_label, " (", nInd(pop_ld), " individuals, fewer than ind.limit = ",
        ind.limit, ")"
      ))
      next
    }

    if (verbose >= 2) {
      cat(report("  Calculating pairwise LD in population", pop_label, "\n"))
    }
    # ordering SNPs by chromosome and position; loc.metrics re-subset from
    # the original object so they stay in step with the genotypes
    loc_order <- order(pop_ld@chromosome, pop_ld@position)
    hold <- pop_ld
    pop_ld <- hold[, loc_order]
    pop_ld@other$loc.metrics <- hold@other$loc.metrics[loc_order, , drop = FALSE]
    pop_ld <- gl.recalc.metrics(pop_ld, verbose = 0)
    if (maf > 0) {
      pop_ld <- gl.filter.maf(pop_ld, threshold = maf, verbose = 0)
    }
    if (nLoc(pop_ld) < 4) {
      if (verbose >= 1) {
        cat(warn(
          "  Skipping population", pop_label,
          "because fewer than 4 SNPs remain after filtering.\n"
        ))
      }
      skipped <- c(skipped, paste0(
        pop_label, " (fewer than 4 SNPs after filtering)"
      ))
      next
    }

    # per-locus heterozygosity for the whole population, computed once and
    # subset per chromosome below
    het_pop <- as.numeric(as.matrix(utils.basic.stats(pop_ld)$Hs)[, 1])

    # PLINK intermediates: written to and read from the same per-call path,
    # independent of the working directory set with gl.set.wd()
    plink_prefix <- tempfile(pattern = "gl_plink_")
    gl2plink(pop_ld,
      outfile = basename(plink_prefix),
      outpath = dirname(plink_prefix),
      verbose = 0
    )

    # Read a pedfile as "SnpMatrix" object using a modified version of the
    # function read.pedfile from package snpStats
    snp_stats <-
      utils.read.ped(
        file = paste0(plink_prefix, ".ped"),
        snps = paste0(plink_prefix, ".map"),
        sep = " ",
        show_warnings = FALSE,
        na.strings = NA
      )
    unlink(paste0(plink_prefix, c(".ped", ".map")))

    ld_map <- snp_stats$map
    colnames(ld_map) <-
      c(
        "chr",
        "snp.name",
        "null",
        "loc_bp",
        "allele.1",
        "allele.2"
      )
    ld_map$chr <- as.character(pop_ld@chromosome)
    genotype <- snp_stats$genotypes
    colnames(genotype@.Data) <- ld_map$loc_bp

    for (chr_name in chr_list) {
      ld_loci <- which(ld_map$chr == chr_name)
      if (length(ld_loci) == 0) {
        next
      }
      ld_map_loci <- ld_map[ld_loci, , drop = FALSE]
      genotype_loci <- genotype[, ld_loci, drop = FALSE]
      het_loci <- het_pop[ld_loci]
      # removing loci that have the same location
      dupl_loci <- which(duplicated(ld_map_loci$loc_bp))
      if (length(dupl_loci) > 0) {
        ld_map_loci <- ld_map_loci[-dupl_loci, , drop = FALSE]
        genotype_loci <- genotype_loci[, -dupl_loci, drop = FALSE]
        het_loci <- het_loci[-dupl_loci]
      }
      n_snps_chr <- nrow(ld_map_loci)
      # the rotated LD raster needs at least four SNPs to form a polygon
      if (n_snps_chr < 4) {
        if (verbose >= 1) {
          cat(warn(
            "  Skipping chromosome", chr_name, "in population", pop_label,
            "because it has fewer than 4 SNPs after filtering.\n"
          ))
        }
        next
      }
      if (verbose >= 2) {
        cat(report("  Analysing chromosome", chr_name, "\n"))
      }
      loc_bp <- ld_map_loci$loc_bp

      # this is the mean distance between each snp which is used to determine
      # the depth at which LD analyses are performed
      mean_dis <- mean(diff(loc_bp))
      ld_depth_b <- ceiling((ld_max_pairwise / mean_dis)) - 1

      if (ld_depth_b < 5) {
        if (verbose >= 1) {
          cat(warn(
            "  The maximum distance at which LD should be calculated",
            "(ld_max_pairwise) is too short for chromosome", chr_name,
            ". Setting this distance to", round(mean_dis * 5, 0), "bp\n"
          ))
        }
        ld_depth_b <- 5
      }
      # snpStats::ld cannot look further than the number of SNPs minus one
      ld_depth_b <- min(ld_depth_b, n_snps_chr - 1)
      # function to calculate LD
      ld_snps <- snpStats::ld(genotype_loci,
        depth = ld_depth_b,
        stats = ld_stat
      )
      ld_dense <- as.matrix(ld_snps)
      dimnames(ld_dense) <- NULL
      # the sparse result stores 0 both for pairs beyond the depth band and
      # for computed pairs with no LD; cells are therefore selected by band
      # index, so LD of exactly 0 and negative LD are drawn
      band <- abs(col(ld_dense) - row(ld_dense))
      in_band <- which(band >= 1 & band <= ld_depth_b, arr.ind = TRUE)
      ld_columns_2 <- data.frame(
        Var1 = in_band[, "row"],
        Var2 = in_band[, "col"],
        Freq = ld_dense[in_band]
      )
      # remove cases where LD was not calculated
      ld_columns_2 <- ld_columns_2[complete.cases(ld_columns_2), ]
      raster_haplo <- raster::rasterFromXYZ(ld_columns_2)
      polygon_haplo <- suppressMessages(
        raster::rasterToPolygons(
          raster_haplo,
          fun = NULL,
          n = 4,
          na.rm = TRUE,
          digits = 12,
          dissolve = TRUE
        )
      )

      polygon_haplo <- sp::elide(polygon_haplo, rotate = 45)
      polygon_haplo$id <- rownames(as.data.frame(polygon_haplo))

      # this only has the coordinates
      polygon_haplo.pts <- fortify(polygon_haplo)
      # add the attributes back
      polygon_haplo.df <-
        merge(polygon_haplo.pts,
          as.data.frame(polygon_haplo),
          by = "id"
        )
      width_poly <- round(max(polygon_haplo.df$long), 0)
      height_poly <- max(polygon_haplo.df$lat)

      # x coordinate of each SNP along the rotated LD plot
      snp_x <- scales::rescale(seq_len(n_snps_chr),
        from = c(1, n_snps_chr),
        to = c(0, width_poly)
      )

      # SNP heterozygosity track for the SNPs of this chromosome
      snp_het_alone <- data.frame(
        position = snp_x,
        het = scales::rescale(het_loci, to = c(0, height_poly))
      )

      # identifying haplotypes: runs of adjacent SNPs whose pairwise LD is at
      # least ld_threshold_haplo, with at least min_snps SNPs
      haplo_blocks <- NULL
      if (haplo_id) {
        adjacent_ld <- ld_dense[cbind(1:(n_snps_chr - 1), 2:n_snps_chr)]
        in_ld <- !is.na(adjacent_ld) & adjacent_ld >= ld_threshold_haplo
        runs <- rle(in_ld)
        run_end <- cumsum(runs$lengths)
        run_start <- run_end - runs$lengths + 1
        haplo_blocks <- data.frame(
          start_idx = run_start[runs$values],
          end_idx = run_end[runs$values] + 1
        )
        haplo_blocks$n_snps <- haplo_blocks$end_idx - haplo_blocks$start_idx + 1
        haplo_blocks <- haplo_blocks[haplo_blocks$n_snps >= min_snps, ,
          drop = FALSE
        ]
        if (nrow(haplo_blocks) == 0) {
          haplo_blocks <- NULL
          if (verbose >= 2) {
            cat(warn(
              "  No haplotypes with at least", min_snps,
              "SNPs were found for chromosome", chr_name,
              ". Try using a lower threshold.\n"
            ))
          }
        }
      }

      if (!is.null(haplo_blocks)) {
        locations_temp_2 <- data.frame(
          start = loc_bp[haplo_blocks$start_idx],
          end = loc_bp[haplo_blocks$end_idx],
          start_ld_plot = snp_x[haplo_blocks$start_idx],
          end_ld_plot = snp_x[haplo_blocks$end_idx]
        )
        locations_temp_2$midpoint <-
          (locations_temp_2$start + locations_temp_2$end) / 2
        locations_temp_2$midpoint_ld_plot <-
          (locations_temp_2$start_ld_plot + locations_temp_2$end_ld_plot) / 2
        locations_temp_2$labels <-
          paste0(
            as.character(round(locations_temp_2$start / 1000000, 0)), "-",
            as.character(round(locations_temp_2$end / 1000000, 0))
          )

        # alternate background shading between consecutive haplotypes
        haplo_temp_a <-
          locations_temp_2[seq_len(nrow(locations_temp_2)) %% 2 == 1, ,
            drop = FALSE
          ]
        haplo_temp_b <-
          locations_temp_2[seq_len(nrow(locations_temp_2)) %% 2 == 0, ,
            drop = FALSE
          ]

        ticks_breaks <-
          c(
            locations_temp_2$start_ld_plot,
            locations_temp_2$end_ld_plot
          )
        ticks_lab <- c(locations_temp_2$start, locations_temp_2$end)
        ticks_joint <- data.frame(
          ticks_breaks = ticks_breaks[order(ticks_breaks)],
          ticks_lab = as.character(round(ticks_lab[order(ticks_lab)] / 1000000, 0)),
          stringsAsFactors = FALSE
        )
        ticks_joint <- ticks_joint[!duplicated(ticks_joint$ticks_lab), ]

        colors_plot <-
          c(
            "Heterozygosity" = color_het,
            "Haplotypes limits" = "lightgoldenrod3"
          )
        labels_haplo <- as.character(seq_len(nrow(locations_temp_2)))

        haplo_table_tmp <- data.frame(
          population = pop_label,
          chromosome = chr_name,
          haplotype = seq_len(nrow(locations_temp_2)),
          locations_temp_2,
          stringsAsFactors = FALSE
        )
        haplo_table <- rbind(haplo_table, haplo_table_tmp)

        y_min <- min(polygon_haplo.df$lat) - 30
        y_max <- max(polygon_haplo.df$lat) + 30

        p_temp <- ggplot() +
          geom_rect(
            data = haplo_temp_a,
            aes(
              xmin = start_ld_plot,
              xmax = end_ld_plot
            ),
            ymin = y_min,
            ymax = y_max,
            color = "cornsilk3",
            fill = "cornsilk3"
          ) +
          geom_rect(
            data = haplo_temp_b,
            aes(
              xmin = start_ld_plot,
              xmax = end_ld_plot
            ),
            ymin = y_min,
            ymax = y_max,
            color = "cornsilk4",
            fill = "cornsilk4"
          ) +
          geom_polygon(
            data = polygon_haplo.df,
            aes(long, lat,
              group = group,
              fill = polygon_haplo.df[, "Freq"]
            )
          ) +
          viridis::scale_fill_viridis(name = ld_stat, option = color_haplo) +
          geom_vline(
            data = data.frame(xintercept = ticks_breaks),
            aes(
              xintercept = xintercept,
              color = "Haplotypes limits"
            ),
            linewidth = 1
          ) +
          annotate(
            "text",
            x = locations_temp_2$midpoint_ld_plot,
            y = max(polygon_haplo.df$lat) + 13,
            label = labels_haplo,
            size = 3,
            color = "black"
          ) +
          annotate(
            "text",
            x = locations_temp_2$midpoint_ld_plot,
            y = min(polygon_haplo.df$lat) - 13,
            label = locations_temp_2$labels,
            size = 3,
            color = "black"
          ) +
          labs(
            x = "Chromosome location (Mbp)",
            y = "Het",
            title = paste(
              "Population", pop_label, "Chromosome", chr_name, "-",
              n_snps_chr, "SNPs"
            )
          ) +
          scale_x_continuous(
            breaks = ticks_joint$ticks_breaks,
            labels = ticks_joint$ticks_lab
          ) +
          scale_colour_manual(name = "", values = colors_plot) +
          theme_void() +
          theme(
            legend.position = "top",
            legend.text = element_text(size = 10),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_text(hjust = 0.5),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
          ) +
          coord_fixed(ratio = 1 / 1)

        if (plot_het) {
          p_temp <- p_temp +
            geom_line(
              data = snp_het_alone,
              aes(x = position, y = het, color = "Heterozygosity"),
              inherit.aes = FALSE,
              linewidth = 1 / 2,
              alpha = 1
            )
        }
      } else {
        if (haplo_id && verbose >= 2) {
          cat(warn(
            "  No haplotypes were identified for chromosome", chr_name, "\n"
          ))
        }

        colors_plot <- c("Heterozygosity" = color_het)

        p_temp <- ggplot() +
          geom_polygon(
            data = polygon_haplo.df,
            aes(long, lat,
              group = group,
              fill = polygon_haplo.df[, "Freq"]
            )
          ) +
          viridis::scale_fill_viridis(name = ld_stat, option = color_haplo) +
          labs(
            x = "Chromosome location (Mbp)",
            y = "Het",
            title = paste(
              "Population", pop_label, "Chromosome", chr_name, "-",
              n_snps_chr, "SNPs"
            )
          ) +
          theme_void() +
          theme(
            legend.position = "top",
            legend.text = element_text(size = 10),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
          ) +
          coord_fixed(ratio = 1 / 1)

        if (plot_het) {
          p_temp <- p_temp +
            geom_line(
              data = snp_het_alone,
              aes(x = position, y = het, color = "Heterozygosity"),
              inherit.aes = FALSE,
              linewidth = 1 / 2,
              alpha = 1
            ) +
            scale_colour_manual(name = "", values = colors_plot)
        }

        if (snp_pos) {
          # SNP position track: each SNP joins its index position (top) to
          # its physical position (bottom), both on the chromosome's scale
          snp <- data.frame(snp = loc_bp)
          snp$order <- seq_len(nrow(snp))
          snp$scale <- scales::rescale(snp$snp, to = c(1, n_snps_chr))
          snp1 <- snp
          snp1$y <- 2
          snp2 <- snp
          snp2$y <- 1
          snp2$order <- snp1$scale
          snp_fin <- rbind(snp1, snp2)

          axis_start <- if (!is.null(coordinates)) coordinates[1] else 1
          labels_tmp <- c(
            seq(axis_start, max(loc_bp), max(loc_bp) / 10),
            max(loc_bp)
          )
          labels_plot <- round(labels_tmp, -6) / 1000000
          breaks_plot <- scales::rescale(labels_plot * 1000000,
            from = range(loc_bp),
            to = c(1, n_snps_chr)
          )

          snp_fin$colorL <- "col.all"

          nearest_snp <- function(targets) {
            vapply(targets, function(target) {
              snp_fin$snp[which.min(abs(snp_fin$snp - target))]
            }, numeric(1))
          }
          if (!is.null(target.snp1)) {
            snp_fin[snp_fin$snp %in% nearest_snp(target.snp1), "colorL"] <-
              "col.target1"
          }
          if (!is.null(target.snp2)) {
            snp_fin[snp_fin$snp %in% nearest_snp(target.snp2), "colorL"] <-
              "col.target2"
          }
          if (!is.null(target.snp3)) {
            snp_fin[snp_fin$snp %in% nearest_snp(target.snp3), "colorL"] <-
              "col.target3"
          }

          p_pos <- ggplot(snp_fin, aes(
            x = order,
            y = y,
            group = snp,
            color = colorL
          )) +
            geom_line() +
            scale_color_manual(values = c(
              "col.all" = col.all,
              "col.target1" = col.target1,
              "col.target2" = col.target2,
              "col.target3" = col.target3
            )) +
            theme(
              panel.background = element_rect(fill = "transparent", colour = NA),
              plot.background = element_rect(fill = "transparent", colour = NA),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              plot.margin = unit(c(0, 0, 0, 0), "null"),
              panel.spacing = unit(c(0, 0, 0, 0), "null"),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              legend.position = "none"
            ) +
            scale_x_continuous(
              breaks = breaks_plot,
              labels = labels_plot
            ) +
            labs(x = "Chromosome location (Mbp)")

          layout <- c(
            area(t = 1, l = 1, b = 5, r = 1),
            area(t = 5.5, l = 1, b = 5.5, r = 1)
          )
          p_temp <- p_temp / p_pos +
            plot_layout(design = layout)
        }
      }

      # PRINTING OUTPUTS
      if (plot.out) {
        print(p_temp)
      }

      # Optionally save the plot
      if (isTRUE(plot.save)) {
        file.name <- file.path(plot.dir, paste0(pop_label, "_", chr_name, ".pdf"))
        ggsave(
          filename = file.name,
          plot = p_temp,
          width = 15,
          height = 5,
          units = "in",
          dpi = "retina",
          bg = "transparent",
          limitsize = FALSE
        )
        if (verbose >= 2) {
          cat(report("  Plot saved to", file.name, "\n"))
        }
      }
    }
  }

  rownames(haplo_table) <- NULL

  if (length(skipped) == length(x_list)) {
    stop(error(
      "  No population was analysed. Skipped:",
      paste(skipped, collapse = "; "), "\n"
    ))
  }
  if (length(skipped) > 0) {
    warning(paste(
      "Populations skipped:", paste(skipped, collapse = "; ")
    ), call. = FALSE)
  }

  if (verbose >= 3) {
    if (nrow(haplo_table) > 0) {
      print(haplo_table, row.names = FALSE)
    } else {
      cat(report("  No haplotypes in the results table\n"))
    }
  }

  # FLAG SCRIPT END

  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN

  return(invisible(haplo_table))
}
