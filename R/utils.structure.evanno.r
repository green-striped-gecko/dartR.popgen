#' @name utils.structure.evanno
#' @title Evanno statistics and plots from a STRUCTURE run object
#' @description
#' Computes the Evanno et al. (2005) statistics from a structure run object
#' (\code{\link{gl.run.structure}}) and builds their plots. Called by
#' \code{\link{gl.evanno}} and \code{\link{gl.run.structure}}.
#'
#' The code was adapted from package strataG (Eric Archer), which is no longer
#' on CRAN. See \code{\link{gl.evanno}} for the statistics returned.
#' @param sr Structure run object from \code{\link{gl.run.structure}}
#' [required].
#' @param plot Whether the combined plot is printed [default TRUE].
#' @param plot.theme Theme added to every plot; NULL leaves the ggplot2
#' default [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  brief progress messages; 3, progress and results summary; 5, full report
#'  [default 0].
#' @return A list with the statistics (element df) and a list of plots
#' (element plots: mean.ln.k, ln.pk, ln.ppk, delta.k when it can be computed,
#' and combined).
#' @author Author(s): Bernd Gruber, original implementation by Eric Archer
#' (\url{https://github.com/EricArcher/strataG}). Custodian: Bernd Gruber --
#' Post to \url{https://groups.google.com/d/forum/dartr}
#' @export

utils.structure.evanno <- function(sr,
                                   plot = TRUE,
                                   plot.theme = NULL,
                                   verbose = 0) {
  if (!is(sr, "structure.result")) {
    stop(error(
      "sr is not a structure.result object returned by gl.run.structure.\n"
    ))
  }
  k.tbl <- table(sapply(sr, function(x) x$summary["k"]))
  if (length(k.tbl) < 3) {
    stop(error(
      "The Evanno method needs at least three values of K;",
      length(k.tbl), "found.\n"
    ))
  }
  sr.smry <- t(sapply(sr, function(x) x$summary))
  ln.k <- tapply(sr.smry[, "est.ln.prob"], sr.smry[
    ,
    "k"
  ], mean)
  sd.ln.k <- tapply(sr.smry[, "est.ln.prob"], sr.smry[
    ,
    "k"
  ], stats::sd)

  k <- as.numeric(names(ln.k))
  n <- length(k)
  # LnP'(K) needs K - 1, and LnP''(K) and delta K need K - 1 and K + 1
  has.prev <- c(FALSE, diff(k) == 1)
  has.next <- c(diff(k) == 1, FALSE)
  both <- has.prev & has.next
  ln.pk <- ln.ppk <- delta.k <- rep(NA_real_, n)
  ln.pk[has.prev] <- ln.k[has.prev] - ln.k[which(has.prev) - 1]
  i <- which(both)
  ln.ppk[i] <- abs(ln.pk[i + 1] - ln.pk[i])
  delta.k[i] <- ln.ppk[i] / sd.ln.k[i]

  if (any(diff(k) != 1) && verbose >= 1) {
    # the smallest and largest K never have delta K; name the others
    no.nb <- k[!both & seq_len(n) != 1 & seq_len(n) != n]
    cat(warn(
      "  Warning: the values of K are not consecutive. LnP''(K) and delta K",
      "need K - 1 and K + 1, so they are NA for K =",
      paste(no.nb, collapse = ", "), "\n"
    ))
  }

  one.rep <- k[both & is.na(sd.ln.k)]
  if (length(one.rep) > 0 && verbose >= 1) {
    cat(warn(
      "  Warning: delta K needs at least two replicates per K; it is NA for",
      "K =", paste(one.rep, collapse = ", "), "\n"
    ))
  }

  # replicates with identical LnP(K) have sd = 0 and delta K is undefined
  zero.sd <- both & !is.na(sd.ln.k) & sd.ln.k == 0
  delta.k[zero.sd] <- NA
  if (any(zero.sd) && verbose >= 1) {
    cat(warn(
      "  Warning: all replicates agree on LnP(K) (sd = 0) for K =",
      paste(k[zero.sd], collapse = ", "), "so delta K is undefined (NA)",
      "there. Read mean LnP(K) instead.\n"
    ))
  }

  df <- data.frame(
    k = k,
    reps = as.numeric(table(sr.smry[, "k"])),
    mean.ln.k = as.numeric(ln.k),
    sd.ln.k = as.numeric(sd.ln.k),
    ln.pk = ln.pk,
    ln.ppk = ln.ppk,
    delta.k = delta.k
  )

  rownames(df) <- NULL
  df$sd.min <- df$mean.ln.k - df$sd.ln.k
  df$sd.max <- df$mean.ln.k + df$sd.ln.k
  plot.list <- list(
    mean.ln.k =
      ggplot2::ggplot(
        df,
        ggplot2::aes(
          x = .data$k,
          y = .data$mean.ln.k
        )
      ) +
        ggplot2::ylab("mean LnP(K)") +
        ggplot2::geom_segment(ggplot2::aes(
          x = .data$k,
          xend = .data$k,
          y = .data$sd.min,
          yend = .data$sd.max
        )),
    ln.pk = ggplot2::ggplot(
      df[!is.na(df$ln.pk), ],
      ggplot2::aes(
        x = .data$k,
        y = .data$ln.pk
      )
    ) +
      ggplot2::ylab("LnP'(K)"),
    ln.ppk = ggplot2::ggplot(
      df[!is.na(df$ln.ppk), ],
      ggplot2::aes(
        x = .data$k,
        y = .data$ln.ppk
      )
    ) +
      ggplot2::ylab("LnP''(K)")
  )

  if (!all(is.na(df$delta.k))) {
    plot.list$delta.k <- ggplot2::ggplot(
      df[!is.na(df$delta.k), ],
      ggplot2::aes(x = .data$k, y = .data$delta.k)
    ) +
      ggplot2::ylab(expression(Delta(K)))
  }
  for (i in 1:length(plot.list)) {
    plot.list[[i]] <- plot.list[[i]] + ggplot2::geom_line() +
      ggplot2::geom_point(
        fill = "white", shape = 21,
        size = 3
      ) + ggplot2::xlim(c(1, max(df$k)))
    if (!is.null(plot.theme)) {
      plot.list[[i]] <- plot.list[[i]] + plot.theme
    }
    plot.list[[i]] <- plot.list[[i]] +
      ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }

  # the four panels in two columns, K labelled once below
  plot.list$combined <- patchwork::wrap_plots(plot.list, ncol = 2) +
    patchwork::plot_annotation(
      caption = "K",
      theme = ggplot2::theme(
        plot.caption = ggplot2::element_text(hjust = 0.5)
      )
    )

  if (plot) {
    print(plot.list$combined)
  }
  df$sd.min <- df$sd.max <- NULL
  invisible(list(df = df, plots = plot.list))
}
