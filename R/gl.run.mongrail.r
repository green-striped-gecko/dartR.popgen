#' @name gl.run.mongrail
#' @title Classify hybrids with MONGRAIL from a genlight object
#' @family hybridisation
#'
#' @description
#' Runs MONGRAIL, a Bayesian full-likelihood method that classifies individuals
#' into one of six two-generation genealogical classes -- pure population A,
#' pure population B, F1, F2, backcross to A, and backcross to B -- from
#' multilocus SNP genotypes of two reference populations plus a set of test
#' individuals. The function writes the three MONGRAIL \code{.GT} input files
#' from a genlight object, invokes the \code{mongrail2} executable, and returns
#' the posterior probability of each class for every test individual together
#' with the most probable class.
#'
#' @param x Name of the genlight object [required].
#' @param pop.A Population name(s) in \code{pop(x)} forming reference population A
#' [required].
#' @param pop.B Population name(s) in \code{pop(x)} forming reference population B
#' [required].
#' @param test.pop Population name(s) in \code{pop(x)} to classify. NULL uses all
#' individuals not assigned to pop.A or pop.B [default NULL].
#' @param mongrail.path Directory holding the mongrail2 executable
#' [default getOption("mongrail.path")].
#' @param outpath Directory for the .GT input files and the output
#' [default tempdir()].
#' @param output.name Base name for the output file [default "mongrail"].
#' @param recomb.rate Recombination rate per base pair, passed as -r
#' [default 1e-8].
#' @param max.snps.per.chr Maximum SNPs retained per chromosome; chromosomes
#' with more are thinned to an evenly spaced subset. The MONGRAIL .GT hard limit
#' is 32, and exact enumeration is used at or below 16 [default 16].
#' @param prune.threshold Posterior probability threshold for pruning haplotype
#' pairs, passed as -p (0-1). NULL omits it [default NULL].
#' @param max.recomb Maximum recombinations for sparse enumeration, passed as -k.
#' Required when any chromosome holds more than 16 loci; NULL lets the function
#' set it when needed [default NULL].
#' @param phased If TRUE, heterozygotes are written phased (a|b); if FALSE they
#' are written unphased (a/b), which suits DArTseq data [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#'
#' @details
#' MONGRAIL models linkage and recombination along chromosomes, so loci are
#' grouped by their chromosome and position (from \code{x@chromosome} and
#' \code{x@position}). When those slots are empty -- as for an anonymous SNP
#' panel with no reference genome -- every locus is placed on its own
#' chromosome, i.e. treated as unlinked, which is the correct default in the
#' absence of map information. Because the method scales steeply with the number
#' of loci per chromosome, panels mapped to a genome are thinned to
#' \code{max.snps.per.chr} loci per chromosome.
#'
#' The .GT format tolerates unphased genotypes, so the two-allele DArTseq calls
#' are written directly (0/0, 0/1, 1/1) with missing genotypes as ./. and no
#' phasing step is required.
#'
#' MONGRAIL and its \code{mongrail2} executable are the work of the MONGRAIL
#' authors; this function only prepares its input and parses its output. See
#' \url{https://github.com/mongrail/mongrail2}.
#'
#' @author Custodian: Peter J. Unmack -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \dontrun{
#' res <- gl.run.mongrail(testset.gl, pop.A = "pop1", pop.B = "pop2",
#'                        mongrail.path = getOption("mongrail.path"))
#' head(res$assignments)
#' }
#'
#' @export
#' @return A list with \code{$assignments} (a data.frame of test individual,
#' most probable class, and the six posterior probabilities), \code{$files}
#' (the input and output file paths) and \code{$raw} (the parsed output table),
#' returned invisibly.

gl.run.mongrail <- function(x,
                            pop.A,
                            pop.B,
                            test.pop = NULL,
                            mongrail.path = NULL,
                            outpath = NULL,
                            output.name = "mongrail",
                            recomb.rate = 1e-8,
                            max.snps.per.chr = 16,
                            prune.threshold = NULL,
                            max.recomb = NULL,
                            phased = FALSE,
                            verbose = NULL) {

  # PRELIMINARIES -- checking ----------------
  funname <- match.call()[[1]]
  verbose <- gl.check.verbosity(verbose)
  outpath <- gl.check.wd(outpath, verbose = 0)
  utils.flag.start(func = funname, build = "v.2026.1", verbose = verbose)
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  if (is.null(mongrail.path)) {
    mongrail.path <- getOption("mongrail.path")
  }
  if (is.null(mongrail.path) || is.na(mongrail.path)) {
    stop(error(
      "Fatal error: mongrail.path is not set. Give the folder holding the",
      "mongrail2 executable, or register options(mongrail.path=...).\n"
    ))
  }
  os <- Sys.info()[["sysname"]]
  exe.name <- if (os == "Windows") "mongrail2.exe" else "mongrail2"
  exe <- file.path(path.expand(mongrail.path), exe.name)
  if (!file.exists(exe)) {
    stop(error("  MONGRAIL executable not found:", exe,
               "\n  Set mongrail.path to the folder that contains", exe.name,
               "\n"))
  }

  if (missing(pop.A) || missing(pop.B)) {
    stop(error("Fatal error: pop.A and pop.B are both required\n"))
  }
  all.pops <- as.character(unique(pop(x)))
  if (!all(pop.A %in% all.pops) || !all(pop.B %in% all.pops)) {
    stop(error("Fatal error: pop.A / pop.B not found in pop(x). Available:",
               paste(all.pops, collapse = ", "), "\n"))
  }

  # DO THE JOB --------------------------------
  # Partition individuals into the three MONGRAIL groups.
  pv <- as.character(pop(x))
  idx.A <- which(pv %in% pop.A)
  idx.B <- which(pv %in% pop.B)
  if (is.null(test.pop)) {
    idx.T <- which(!(pv %in% c(pop.A, pop.B)))
  } else {
    idx.T <- which(pv %in% test.pop)
  }
  if (length(idx.A) == 0 || length(idx.B) == 0) {
    stop(error("Fatal error: a reference population has no individuals\n"))
  }
  if (length(idx.T) == 0) {
    stop(error("Fatal error: no test individuals to classify\n"))
  }

  # Chromosome and position per locus; fall back to one locus per chromosome
  # (unlinked) when no map information is present.
  chrom <- as.character(x@chromosome)
  pos <- x@position
  if (length(chrom) != nLoc(x) || all(is.na(chrom))) {
    chrom <- paste0("L", seq_len(nLoc(x)))
  }
  if (length(pos) != nLoc(x) || all(is.na(pos))) {
    pos <- seq_len(nLoc(x))
  }
  pos[is.na(pos)] <- 0

  # Thin to max.snps.per.chr loci per chromosome (evenly spaced).
  keep <- unlist(lapply(split(seq_len(nLoc(x)), chrom), function(ix) {
    if (length(ix) <= max.snps.per.chr) return(ix)
    ix[round(seq(1, length(ix), length.out = max.snps.per.chr))]
  }), use.names = FALSE)
  keep <- sort(keep)
  if (verbose >= 2 && length(keep) < nLoc(x)) {
    cat(report("  Thinned from", nLoc(x), "to", length(keep),
               "loci to respect max.snps.per.chr =", max.snps.per.chr, "\n"))
  }
  chrom <- chrom[keep]
  pos <- pos[keep]

  # Dosage -> MONGRAIL genotype string.
  sep.het <- if (phased) "|" else "/"
  dose2gt <- function(d) {
    g <- rep(NA_character_, length(d))
    g[d == 0] <- paste0("0", sep.het, "0")
    g[d == 1] <- paste0("0", sep.het, "1")
    g[d == 2] <- paste0("1", sep.het, "1")
    g[is.na(d)] <- paste0(".", sep.het, ".")
    g
  }

  # Write one .GT file for a set of individuals. Rows are loci; each line is
  # "chrom:pos:" followed by space-separated genotypes in individual order.
  write.gt <- function(ind.idx, file) {
    m <- as.matrix(x[ind.idx, keep])          # individuals x loci (thinned)
    tag <- paste0(chrom, ":", pos, ":")
    con <- file(file, "w")
    on.exit(close(con), add = TRUE)
    for (j in seq_len(ncol(m))) {
      writeLines(paste(c(tag[j], dose2gt(m[, j])), collapse = " "), con)
    }
    invisible(file)
  }

  f.A <- file.path(outpath, paste0(output.name, "_popA.GT"))
  f.B <- file.path(outpath, paste0(output.name, "_popB.GT"))
  f.T <- file.path(outpath, paste0(output.name, "_hybrids.GT"))
  write.gt(idx.A, f.A)
  write.gt(idx.B, f.B)
  write.gt(idx.T, f.T)
  if (verbose >= 2) {
    cat(report("  Wrote .GT files:", length(idx.A), "in A,", length(idx.B),
               "in B,", length(idx.T), "test individuals\n"))
  }

  # Build the command line.
  args <- character(0)
  if (!phased) {
    # .GT with unphased genotypes is the default path (no -c / VCF).
  }
  args <- c(args, "-r", format(recomb.rate, scientific = TRUE))
  if (!is.null(prune.threshold)) args <- c(args, "-p", prune.threshold)
  max.per.chr <- max(table(chrom))
  if (max.per.chr > 16) {
    if (is.null(max.recomb)) max.recomb <- 2
    args <- c(args, "-k", max.recomb)
    if (verbose >= 2) {
      cat(report("  >16 loci on a chromosome; using sparse enumeration -k",
                 max.recomb, "\n"))
    }
  }
  if (verbose >= 2) args <- c(args, "-v")
  args <- c(args, shQuote(f.A), shQuote(f.B), shQuote(f.T))

  out.file <- file.path(outpath, paste0(output.name, ".out"))
  if (verbose >= 2) {
    cat(report("  Running:", exe.name, paste(args, collapse = " "), "\n"))
  }
  status <- system2(exe, args = args, stdout = out.file,
                    stderr = if (verbose >= 2) "" else FALSE)
  if (is.na(status) || status != 0) {
    stop(error(paste0("  MONGRAIL failed to run (exit status ", status,
                      "): ", exe, "\n")))
  }
  if (!file.exists(out.file) || length(readLines(out.file, warn = FALSE)) < 2) {
    stop(error("  MONGRAIL produced no output. Run with verbose >= 2 to see",
               "its console messages.\n"))
  }

  # Parse the output: header row then one row per test individual, with the six
  # class labels and their posterior probabilities.
  # MONGRAIL writes a #-commented banner and a #-commented header line
  # ("# indiv Ma .. Mf P(Ma) .. P(Mf)"), then one data row per test individual
  # (the indiv column is a 0-based index, in the order the hybrids were written).
  lines <- readLines(out.file, warn = FALSE)
  hdr.i <- grep("indiv", lines, fixed = TRUE)
  hdr.i <- hdr.i[grepl("P(Ma)", lines[hdr.i], fixed = TRUE)]
  if (length(hdr.i) == 0) {
    stop(error("  Could not find the MONGRAIL header line in", out.file, "\n"))
  }
  col.names <- strsplit(trimws(sub("^#", "", lines[hdr.i[1]])),
                        "[[:space:]]+")[[1]]
  data.lines <- lines[!grepl("^[[:space:]]*#", lines) &
                        nzchar(trimws(lines))]
  raw <- utils::read.table(text = data.lines, col.names = col.names,
                           stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(raw) != length(idx.T)) {
    cat(warn("  Warning: MONGRAIL returned", nrow(raw), "rows for",
             length(idx.T), "test individuals\n"))
  }
  prob.cols <- which(startsWith(colnames(raw), "P("))
  probs <- as.matrix(raw[, prob.cols, drop = FALSE])
  class.codes <- c("Ma", "Mb", "Mc", "Md", "Me", "Mf")
  class.labels <- c(Ma = "pure B", Mb = "backcross to A", Mc = "F1",
                    Md = "pure A", Me = "backcross to B", Mf = "F2")
  best <- class.codes[max.col(probs, ties.method = "first")]

  assignments <- data.frame(
    id = indNames(x)[idx.T],
    class = unname(class.labels[best]),
    class.code = best,
    stringsAsFactors = FALSE
  )
  colnames(probs) <- paste0("P.", class.codes)
  assignments <- cbind(assignments, probs)

  if (verbose >= 3) {
    tab <- table(assignments$class)
    cat(report("  MONGRAIL classifications:\n"))
    for (nm in names(tab)) cat(report(sprintf("    %-16s %d\n", nm, tab[[nm]])))
  }

  # FLAG SCRIPT END ---------------------------
  if (verbose >= 1) {
    cat(report("Completed:", as.character(funname), "\n"))
  }

  # RETURN
  invisible(list(
    assignments = assignments,
    files = list(popA = f.A, popB = f.B, hybrids = f.T, output = out.file),
    raw = raw
  ))
}
