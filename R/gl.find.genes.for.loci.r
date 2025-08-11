#' @name gl.find.genes.for.loci
#'
#' @title Map loci (SNPs) to the nearest gene feature from a GFF
#'
#' @description
#' Given a SNP genlight object and a GFF3 annotation file, find the closest gene
#' (or transcript, if requested) for each input locus. If a locus falls within a
#' gene interval, that gene is considered the closest with distance 0.
#'
#' @param x A SNP genlight object with mapped loci. Must contain per-locus
#'   x$chromosome and x$position. [required]
#' @param gff.file Path to a GFF3 file (either plain or with a .gz alongside).
#'   [required]
#' @param loci Character vector of locus names to map. Must match locNames(x).
#'   [required]
#' @param include_types Character vector of GFF types to treat as "gene"
#'   features. Defaults to c("gene","pseudogene").
#' @param fallback_to_mrna Logical. If no rows match include_types, use
#'   transcript features c("mRNA","transcript") as proxies. [default TRUE]
#' @param save2tmp Logical: save the result table to tempdir() (retrievable
#'   with gl.list.reports / gl.print.reports). [default FALSE]
#' @param verbose Verbosity: 0-5 (see gl.set.verbosity()). [default from
#'   gl.check.verbosity()]
#'
#' @details
#' The function parses common keys in the GFF attributes column (e.g., ID,
#' Name, gene, product, Parent) to provide informative gene labels.
#' Closeness is measured on the same sequence (chromosome/contig) as:
#' - 0 if the locus is within [gene_start, gene_end]
#' - otherwise, the minimum bp distance to the interval edges
#'
#' If multiple genes are exactly equally close, a deterministic tie-break is
#' applied: closest to gene midpoint, then shorter gene length, then lexicographic gene_id.
#'
#' @return A data.table with one row per input locus and columns:
#'   locus, chrom, pos, gene_start, gene_end, gene_type, gene_id, gene_name,
#'   gene_symbol, gene_product, gene_attributes, distance_bp, nearest_side.
#'   `distance_bp` is the absolute distance in bp; `nearest_side` is "inside",
#'   "left" (locus < gene_start), or "right" (locus > gene_end) in coordinate space.
#'
#' @examples
#' \dontrun{
#' res <- gl.find.genes.for.loci(
#'   x = testset.gl,
#'   gff.file = "species.gff3",
#'   loci = c("locus_12","locus_51","locus_89")
#' )
#' }
#'
#' @family annotation and mapping helpers
#'
#' @importFrom ape read.gff
#' @importFrom data.table as.data.table data.table setkey foverlaps :=
#' @importFrom stringr str_match str_detect
#'
#' @export
gl.find.genes.for.loci <- function(x,
                                   gff.file,
                                   loci,
                                   include_types = c("gene","pseudogene"),
                                   fallback_to_mrna = TRUE,
                                   save2tmp = FALSE,
                                   verbose = NULL) {
  
  # Avoid R CMD check NOTES for data.table NSE vars
  start <- end <- ID <- Name <- gene_k <- gene_id <- product <- Parent <-
    gene_sym <- type <- seqid <- gene_seqid <- gene_start <- gene_end <- NA
  chrom <- locus <- pos <- gene_type <- gene_name <- gene_symbol <-
    gene_product <- gene_attributes <- NA
  join_coord <- which_end <- gene_mid <- gene_len <- distance_bp <-
    nearest_side <- i.pos <- i.locus <- i.chrom <- NA
  
  # VERBOSITY ------------------------------------------------------------------
  verbose <- gl.check.verbosity(verbose)
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)
  
  # CHECKS ---------------------------------------------------------------------
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)
  
  for (pkg in c("ape", "data.table", "stringr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(error("Package ", pkg, " is required. Please install it.\n"))
    }
  }
  
  if (missing(gff.file) || !is.character(gff.file) || length(gff.file) != 1) {
    stop(error("Argument 'gff.file' must be a length-1 character path.\n"))
  }
  if (missing(loci) || !is.character(loci) || length(loci) < 1) {
    stop(error("Argument 'loci' must be a non-empty character vector.\n"))
  }
  
  if (is.null(x$chromosome) || is.null(x$position)) {
    stop(error("Input 'x' must contain per-locus 'chromosome' and 'position'.\n"))
  }
  if (length(x$chromosome) != nLoc(x) || length(x$position) != nLoc(x)) {
    stop(error("Lengths of 'x$chromosome' and 'x$position' must equal nLoc(x).\n"))
  }
  
  all_names <- locNames(x)
  in_x <- loci %in% all_names
  if (!all(in_x)) {
    missing_loci <- loci[!in_x]
    warning(warn(length(missing_loci), " locus name(s) not found in x: ",
                 paste(utils::head(missing_loci, 10), collapse = ", "),
                 if (length(missing_loci) > 10) " ..." else "", "\n"))
    loci <- loci[in_x]
  }
  if (length(loci) == 0L) {
    if (verbose >= 1) cat(report("No valid loci to map; returning empty table.\n"))
    return(data.table::data.table(
      locus=character(), chrom=character(), pos=integer(),
      gene_start=integer(), gene_end=integer(),
      gene_type=character(), gene_id=character(), gene_name=character(),
      gene_symbol=character(), gene_product=character(),
      gene_attributes=character(), distance_bp=integer(),
      nearest_side=character()
    ))
  }
  
  # LOAD GFF -------------------------------------------------------------------
  if (verbose >= 2) cat(report("Loading GFF and parsing attributes...\n"))
  gff_plain <- gff.file
  gff_gz    <- paste0(gff_plain, ".gz")
  
  has_plain <- file.exists(gff_plain)
  has_gz    <- file.exists(gff_gz)
  if (!has_plain && !has_gz) {
    stop(error("Cannot find '", gff_plain, "' or compressed '", gff_gz, "'.\n"))
  }
  
  gff <- if (has_plain) ape::read.gff(gff_plain) else ape::read.gff(gzfile(gff_gz))
  gff_dt <- data.table::as.data.table(gff)
  gff_dt[, ':='(start = as.integer(start), end = as.integer(end))]
  
  # ATTRIBUTE PARSERS ----------------------------------------------------------
  extract_attr <- function(x_attr, key) {
    m <- stringr::str_match(x_attr, paste0("(^|;)", key, "=([^;]+)"))
    m[, 3]
  }
  
  gff_dt[, ID        := extract_attr(attributes, "ID")]
  gff_dt[, Name      := extract_attr(attributes, "Name")]
  gff_dt[, gene_k    := extract_attr(attributes, "gene")]
  gff_dt[, gene_id   := data.table::fcoalesce(gene_k, ID)]
  gff_dt[, product   := extract_attr(attributes, "product")]
  gff_dt[, Parent    := extract_attr(attributes, "Parent")]
  gff_dt[, gene_sym  := extract_attr(attributes, "gene_symbol")]
  
  # CHOOSE GENE FEATURES -------------------------------------------------------
  gene_feats <- gff_dt[type %in% include_types]
  
  if (nrow(gene_feats) == 0L && isTRUE(fallback_to_mrna)) {
    if (verbose >= 2) cat(warn("No 'gene' features found. Falling back to transcripts (mRNA/transcript).\n"))
    gene_feats <- gff_dt[type %in% c("mRNA", "transcript")]
  }
  
  if (nrow(gene_feats) == 0L) {
    warning(warn("No suitable gene (or transcript) features found in the GFF.\n"))
    return(data.table::data.table(
      locus=character(), chrom=character(), pos=integer(),
      gene_start=integer(), gene_end=integer(),
      gene_type=character(), gene_id=character(), gene_name=character(),
      gene_symbol=character(), gene_product=character(),
      gene_attributes=character(), distance_bp=integer(),
      nearest_side=character()
    ))
  }
  
  # BUILD INTERVALS ------------------------------------------------------------
  gene_iv <- gene_feats[, list(
    gene_seqid = as.character(seqid),
    gene_start = as.integer(start),
    gene_end   = as.integer(end),
    gene_type  = as.character(type),
    gene_id    = data.table::fcoalesce(gene_id, Name, ID),
    gene_name  = data.table::fcoalesce(Name, gene_k, ID),
    gene_symbol= gene_sym,
    gene_product = product,
    gene_attributes = attributes
  )]
  gene_iv <- gene_iv[!is.na(gene_seqid) & !is.na(gene_start) & !is.na(gene_end)]
  data.table::setkey(gene_iv, gene_seqid, gene_start, gene_end)
  
  # LOCI TABLE (subset to requested loci) --------------------------------------
  loci_idx <- match(loci, all_names)
  loci_dt <- data.table::data.table(
    chrom = as.character(x$chromosome)[loci_idx],
    start = as.integer(x$position)[loci_idx],
    end   = as.integer(x$position)[loci_idx],
    pos   = as.integer(x$position)[loci_idx],
    locus = loci
  )
  loci_dt <- loci_dt[!is.na(chrom) & !is.na(start)]
  data.table::setkey(loci_dt, chrom, start, end)
  
  if (nrow(loci_dt) == 0L) {
    if (verbose >= 1) cat(report("No mappable loci (missing chrom/pos). Returning empty table.\n"))
    return(data.table::data.table(
      locus=character(), chrom=character(), pos=integer(),
      gene_start=integer(), gene_end=integer(),
      gene_type=character(), gene_id=character(), gene_name=character(),
      gene_symbol=character(), gene_product=character(),
      gene_attributes=character(), distance_bp=integer(),
      nearest_side=character()
    ))
  }
  
  # 1) OVERLAPS (distance = 0) -------------------------------------------------
  hits0 <- data.table::foverlaps(
    x = loci_dt, y = gene_iv,
    by.x = c("chrom","start","end"),
    by.y = c("gene_seqid","gene_start","gene_end"),
    nomatch = 0L
  )
  
  if (nrow(hits0)) {
    hits0[, gene_mid := (gene_start + gene_end)/2]
    hits0[, gene_len := gene_end - gene_start]
    # deterministic pick: closest to midpoint, then shorter gene, then gene_id
    best_overlap <- hits0[order(locus, abs(pos - gene_mid), gene_len, gene_id),
                          .SD[1L], by = locus]
    best_overlap[, `:=`(distance_bp = 0L, nearest_side = "inside")]
    best_overlap_out <- best_overlap[, list(
      locus, chrom, pos, gene_start, gene_end, gene_type, gene_id, gene_name,
      gene_symbol, gene_product, gene_attributes, distance_bp, nearest_side
    )]
  } else {
    best_overlap_out <- data.table::data.table(
      locus=character(), chrom=character(), pos=integer(),
      gene_start=integer(), gene_end=integer(),
      gene_type=character(), gene_id=character(), gene_name=character(),
      gene_symbol=character(), gene_product=character(),
      gene_attributes=character(), distance_bp=integer(),
      nearest_side=character()
    )
  }
  
  overlapped_loci <- unique(best_overlap_out$locus)
  non_overlap <- loci_dt[!(locus %chin% overlapped_loci)]
  
  # 2) NEAREST BY INTERVAL EDGES (non-overlapping loci) -----------------------
  # Build a table of gene boundary points (starts and ends)
  gene_pts <- rbind(
    gene_iv[, list(gene_seqid, join_coord = gene_start,
                gene_start, gene_end, gene_type, gene_id, gene_name,
                gene_symbol, gene_product, gene_attributes, which_end = "start")],
    gene_iv[, list(gene_seqid, join_coord = gene_end,
                gene_start, gene_end, gene_type, gene_id, gene_name,
                gene_symbol, gene_product, gene_attributes, which_end = "end")]
  )
  data.table::setkey(gene_pts, gene_seqid, join_coord)
  
  if (nrow(non_overlap)) {
    points <- non_overlap[, list(chrom, join_coord = pos, pos, locus)]
    data.table::setkey(points, chrom, join_coord)
    
    nearest_edge <- gene_pts[points, roll = "nearest"]
    
    non_out <- nearest_edge[, {
      d <- abs(i.pos - join_coord)
      side <- ifelse(i.pos < gene_start, "left", "right")
      data.table::data.table(
        locus       = i.locus,
        chrom       = i.chrom,
        pos         = i.pos,
        gene_start  = gene_start,
        gene_end    = gene_end,
        gene_type   = gene_type,
        gene_id     = gene_id,
        gene_name   = gene_name,
        gene_symbol = gene_symbol,
        gene_product= gene_product,
        gene_attributes = gene_attributes,
        distance_bp = as.integer(d),
        nearest_side= side
      )
    }]
  } else {
    non_out <- data.table::data.table(
      locus=character(), chrom=character(), pos=integer(),
      gene_start=integer(), gene_end=integer(),
      gene_type=character(), gene_id=character(), gene_name=character(),
      gene_symbol=character(), gene_product=character(),
      gene_attributes=character(), distance_bp=integer(),
      nearest_side=character()
    )
  }
  
  # COMBINE (one row per locus) -----------------------------------------------
  out <- rbind(best_overlap_out, non_out, fill = TRUE)
  
  # Safety: ensure each requested locus appears once (if present in x)
  # If a locus ended up duplicated somehow, keep the smallest distance.
  if (nrow(out)) {
    out <- out[order(locus, distance_bp, abs(pos - ((gene_start + gene_end)/2)),
                     (gene_end - gene_start), gene_id),
               .SD[1L], by = locus]
  }
  
  if (verbose >= 2) {
    cat(report("Assigned nearest gene for ", nrow(out), " locus/loci.\n"))
  }
  
  # SAVE (OPTIONAL) ------------------------------------------------------------
  if (isTRUE(save2tmp)) {
    fn <- tempfile(pattern = "dartR_table_loci2nearestgene_", fileext = ".rds")
    saveRDS(out, file = fn)
    if (verbose >= 2) {
      cat(report(" Saved table to ", fn, " (via saveRDS)\n"))
      cat(report(" Retrieve with gl.list.reports() / gl.print.reports()\n"))
    }
  }
  
  # END ------------------------------------------------------------------------
  return(out)
}
