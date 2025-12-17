#' @name gl.plot.snmf
#'
#' @title Plots ancestry coefficient from snmf
#'
#' @description
#' This function takes a Q matrix (output from
#'  \code{\link{gl.run.snmf}}) and plots the typical structure bar
#'   plot that visualize the Q matrix of a structure run.
#' @param snmf.result run object from \code{\link{gl.run.snmf}} [required].
#' @param border.ind The width of the border line between individuals
#' [default 0.25].
#' @param plot.K The number for K of the Q matrix that should be plotted. Needs to
#'  be within you simulated range of K's in your snmf run object [required].
#' @param plot.theme Theme for the plot. See Details for options
#' [default NULL].
#' @param color.clusters A color palette for clusters (K) or a list with
#' as many colors as there are clusters (K) [default NULL].
#' @param ind.name Whether to plot individual names [default TRUE].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot.file Name for the RDS binary file to save (base name only, exclude
#' extension) [default NULL].
#' @param plot.dir Directory in which to save files [default = working directory].
#' @param den Whether to include a dendrogram. It is necessary to include the 
#' original genlight object used in gl.run.structure in the parameter x 
#' [default FALSE].
#' @param inverse.den Flip dendrogram upside down [default TRUE].
#' @param x The original genlight object used in gl.run.structure description
#' [default NULL]. 
#' @param plot.colors.pop A color palette for population plots or a list with
#' as many colors as there are populations in the dataset 
#' [default gl.colors("dis")].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log ; 3, progress and results summary; 5, full report [default
#'   2, unless specified using gl.set.verbosity].
#'
#' @details The function outputs a barplot which is the typical output of
#'  snmf.
#'  Plots and table are saved to the working directory specified in plot.dir
#'  if plot.file is set.
#'
#' Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#'The Q matrices can be input to other R packages for plotting ancestry proportion, e.g. FSTruct
#'\url{https://github.com/MaikeMorrison/FSTruct}
#' @return Q-matrix
#'
#' @author Ching Ching Lau (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' # examples need LEA to be installed on the system (see above)
#' \dontrun{
#' m <- gl.run.snmf(x=bandicoot.gl, minK=1, 
#' maxK=5, rep=10)
#' Q <- gl.plot.snmf(snmf.result=m, plot.K = 3, ind.name=T)
#' gl.map.snmf(bandicoot.gl, qmat=Q)
#' # move population 4 (out of 5) 0.5 degrees to the right and populations 1
#' # 0.3 degree to the north of the map.
#' mp <- data.frame(lon=c(0,0,0,0.5,0), lat=c(-0.3,0,0,0,0))
#' gl.map.snmf(bandicoot.gl, qmat=Q, movepops=mp)
#' }
#' @export
#' @seealso \code{gl.run.snmf}, \code{gl.plot.snmf}
#' @references
#' \itemize{
#' \item Frichot E, Mathieu F, Trouillon T, Bouchard G, Francois O. (2014). Fast and Efficient 
#' Estimation of Individual Ancestry Coefficients. Genetics, 194(4): 973--983.
#' 
#' }

gl.plot.snmf <- function(snmf.result,
                         border.ind = 0.25,
                         plot.K,
                         plot.theme = NULL,
                         color.clusters = NULL,
                         ind.name = TRUE,
                         plot.out = TRUE,
                         plot.file = NULL,
                         plot.dir = NULL,
                         den = FALSE,
                         inverse.den = TRUE,
                         x = NULL,
                         plot.colors.pop = gl.colors("dis"),
                         verbose = 2) {
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(
    func  = funname,
    build = "Jody",
    verbose = verbose
  )
  
  if (is.null(plot.theme)) {
    plot.theme <- theme_dartR()
  }
  
  # Extract admixture analysis
  Q <- snmf.result$matrix[paste0("K", plot.K)][[1]]

  Q_long <- tidyr::pivot_longer(
    Q,
    cols = dplyr::starts_with("Pop_"),
    names_to = "K",
    values_to = "values"
  )
  
  Q_long$K <- as.factor(Q_long$K)
  
  if (is.null(color.clusters)) {
    color.clusters <- gl.colors(type = "structure")[1:max(plot.K)]
  }
  
  # Dendrogram 
  if (den) {
    
    if (is.null(x)) {
      stop("If 'den = TRUE', you must provide a genlight object in 'x'.")
    }
    
    # ---------- DISTANCE + HCLUST + DENDRO ----------
    res <- gl.dist.ind(x, method = "Manhattan", plot.display = FALSE, verbose = 0)
    
    reorderfun <- function(d, w) stats::reorder(d, w, agglo.FUN = mean)
    
    distr <- stats::dist(res)
    hcr   <- stats::hclust(distr)
    ddr   <- stats::as.dendrogram(hcr)
    ddr   <- reorderfun(ddr, TRUE)
    
    # ---------- DENDRO DATA ----------
    dd_data   <- ggdendro::dendro_data(ddr, type = "rectangle")
    segments  <- dd_data$segments   # segments of dendrogram
    labels_df <- dd_data$labels     # columns: x, y, label (tip labels)
    
    # ---------- DOMINANT K PER IND, KEEPING Pop ----------
    domK <- Q_long %>%
      dplyr::group_by(Label, Pop, K) %>%
      dplyr::summarise(mean_val = mean(values), .groups = "drop") %>%
      dplyr::group_by(Label, Pop) %>%        # one K per (Label, Pop)
      dplyr::slice_max(mean_val, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    
    lab2pop <- stats::setNames(domK$Pop, domK$Label)
    
    # ---------- ASSIGN Pop TO TERMINAL BRANCHES ----------
    segments$Pop_group <- NA_character_ 
    
    for (i in seq_len(nrow(labels_df))) {
      lab   <- labels_df$label[i]
      x_pos <- labels_df$x[i]
      
      this_pop <- lab2pop[[lab]]
      if (is.null(this_pop)) next   # skip labels not in domK
      
      # terminal segment for this leaf: yend == 0
      idx <- which(segments$x == x_pos & segments$yend == 0)
      if (length(idx) > 0) {
        segments$Pop_group[idx] <- this_pop
      }
    }
    
    # ---------- POP COLOUR PALETTE ----------
    pop_levels <- sort(unique(stats::na.omit(segments$Pop_group)))
    
    if (is(plot.colors.pop, "function")) {
      all_pops <- levels(pop(x))
      pop_cols_all <- plot.colors.pop(length(all_pops))
      names(pop_cols_all) <- all_pops
    } else {
      pop_cols_all <- plot.colors.pop
      # if (is.null(names(pop_cols_all))) {
        names(pop_cols_all) <- levels(pop(x))[seq_along(pop_cols_all)]
      # }
    }
    
    # subset/reorder to only pops present in the tree
    pop_cols <- pop_cols_all[pop_levels]
    
    # Split internal vs terminal segments
    seg_internal <- segments[is.na(segments$Pop_group), , drop = FALSE]
    seg_terminal <- segments[!is.na(segments$Pop_group), , drop = FALSE]
    
    # ---------- DENDROGRAM PLOT
    labels_df <- dd_data$labels
    
    # Split internal vs terminal segments
    seg_internal <- segments[is.na(segments$Pop_group), , drop = FALSE]
    seg_terminal <- segments[!is.na(segments$Pop_group), , drop = FALSE]
    
    # --- TIP LINES DATA (one per leaf) ---
    
    # leaves are segments with yend == 0 and with Pop_group
    tip_base <- seg_terminal[seg_terminal$yend == 0, c("x", "yend", "Pop_group")]
    tip_base <- tip_base[!duplicated(tip_base$x), , drop = FALSE]
    
    # height of the little vertical lines
    tip_height <- max(segments$y) * 0.05
    
    tip_lines <- transform(
      tip_base,
      y    = yend,                 # start at the leaf tip (usually 0)
      yend = yend - tip_height     # draw a short line downward
    )
    
    tip_lines$yend <- -0.5
    
    p_den <- ggplot() +
      geom_segment(
        data = seg_internal,
        aes(x = x, y = y, xend = xend, yend = yend),
        colour   = "black",
        linewidth = 0.4
      ) +
      geom_segment(
        data = seg_terminal,
        aes(x = x, y = y, xend = xend, yend = -0.5),
        linewidth = 0.4, colour = "black"
      ) +
      geom_segment(
        data = tip_lines,
        aes(x = x, xend = x, y = y, yend = yend, colour = Pop_group),
        linewidth = 1.2   # thicker so they stand out
      ) +
      scale_colour_manual(
        name   = "",
        values = pop_cols
      ) +
      scale_x_continuous(
        limits = c(0.5, nInd(x) + 0.5),
        expand = c(0, 0)
      ) +
      theme_dendro() +
      theme(
        legend.position = "bottom",
        legend.title    = element_text(size = 11),
        legend.text     = element_text(size = 10),
        legend.key.size = unit(0.8, "lines"),
        plot.margin = margin(t = 0, r = 5.5, b = 0, l = 5.5, unit = "pt")
      )
    
    if(inverse.den){
      p_den <- p_den +
      scale_y_reverse(expand = expansion(mult = c(0.02, 0.02))) 
    }
    
    # ---------- ORDER INDIVIDUALS TO MATCH DENDROGRAM ----------
    # Get label order from dendrogram
    rowInd    <- stats::order.dendrogram(ddr)
    lab_order <- indNames(x)[rowInd]   # vector of labels in dendrogram order
    
    rowInd_df <- data.frame(
      Label = lab_order,
      stringsAsFactors = FALSE
    )
    
    # Align Q_long with dendrogram order:
    Q_long <- Q_long %>%
      dplyr::mutate(
        Label = factor(Label, levels = lab_order),
        Order = as.integer(Label)   # <-- THIS is now dendrogram order
      ) %>%
      dplyr::arrange(Order)
    
    # For the barplot we don't want Pop facets (single panel),
    Q_long$Pop <- ""
    
    Q_long <- Q_long %>%
      dplyr::mutate(
        Label = factor(Label, levels = lab_order),    # ensure order matches tree
        id    = as.integer(Label)                     # 1..n in dendrogram order
      ) %>%
      dplyr::arrange(id)
  }
  
  # -----------------------------
  # STRUCTURE-LIKE BARPLOT
  # -----------------------------
  if (!den) {
    # keep the current order but make Label a factor explicitly
    Q_long$Label <- factor(Q_long$Label, levels = unique(Q_long$Label))
  }
  
  Q_long_plot <- Q_long %>%
    group_by(Pop) %>%
    arrange(Label, .by_group = TRUE) %>%
    mutate(
      pos_pop = row_number()  # 1..n within each Pop
    ) %>%
    ungroup()
  
  lab_df <- Q_long_plot %>%
    distinct(Pop, pos_pop, Label)
  
  p3 <- ggplot(Q_long, aes_(x = ~Order, y = ~values, fill = ~K)) +
    geom_col(aes(colour = K), linewidth = border.ind, width = 1) +
    facet_grid(
      ~ factor(Pop, levels = unique(Q_long$Pop)),
      scales = "free",
      space  = "free"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = color.clusters) +
    scale_colour_manual(values = color.clusters, guide = "none") +
    plot.theme +
    theme(
      panel.spacing = grid::unit(0, "lines"),
      panel.border  = element_rect(
        color    = "black",
        fill     = NA ,
        linewidth = 1
      ),
      strip.background = element_blank(),
      strip.text.x     = element_text(size = 12, angle = 90),
      axis.title.x     = element_blank(),
      axis.text.x      = element_text(
        size  = 8,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.title.y  = element_blank(),
      axis.text.y   = element_blank(),
      axis.ticks.y  = element_blank(),
      legend.position = "none"
    )
  
  if(den){
    p3 <- p3 +
      scale_x_continuous(
      limits = c(0.5, nInd(x) + 0.5),
      breaks = seq_len(nInd(x)),
      labels = unique(Q_long$Label),
      expand = c(0, 0)) +
        theme(
        plot.margin = margin(t = 0, r = 5.5, b = 0, l = 5.5, unit = "pt")
        )
  }else {
    lab_df <- lab_df[!duplicated(lab_df$Label),]
    lab_df$pos_pop <- 1:nrow(x)
    p3 <- p3 +  scale_x_continuous(
      expand = c(0, 0),
      breaks = lab_df$pos_pop,
      labels = lab_df$Label
    ) +
      theme(
        plot.margin = margin(t = 0, r = 5.5, b = 0, l = 5.5, unit = "pt")
      )

  }
  
  if (!ind.name) {
    p3 <- p3 + theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  
  # -----------------------------
  # Combine Dendrogram + Barplot
  # -----------------------------
  if (den) {
    p3 <-  p3 / p_den  + 
      plot_layout(heights = c(1, 3))
  }
  
  # -----------------------------
  # Print / Save
  # -----------------------------
  if (plot.out) {
    print(p3)
  }
  
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(
      p3,
      dir     = plot.dir,
      file    = plot.file,
      verbose = verbose
    )
  }
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  invisible(Q)
}
