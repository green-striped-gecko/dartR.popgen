#' @name gl.find.genes.for.loci
#'
#' @title Map loci (SNPs) to the nearest gene feature from a GFF
#'
#' @family annotation and mapping helpers
#'
#' @description
#' Given a SNP genlight object and a GFF3 annotation file, find the closest gene
#' (or transcript, if requested) for each input locus. If a locus falls within a
#' gene interval, that gene is considered the closest with distance 0.
#'
#' @details
#' The function parses common keys in the GFF attributes column (ID, Name,
#' gene, product, Parent, gene_symbol) to provide informative gene labels.
#' When a gene feature carries no product attribute, the product of its
#' child features (mRNA, then CDS or exon) is used, which is where NCBI and
#' Ensembl GFF3 files place it.
#'
#' Closeness is measured on the same sequence (chromosome/contig) as:
#' \itemize{
#' \item 0 if the locus is within [gene_start, gene_end]
#' \item otherwise, the minimum bp distance to the interval edges
#' }
#'
#' If multiple genes are exactly equally close, a deterministic tie-break is
#' applied: closest to gene midpoint, then shorter gene length, then
#' lexicographic gene_id. A locus on a sequence that has no gene feature is
#' returned with NA in the gene columns.
#'
#' The sequence names in x$chromosome must match the first column of the
#' GFF. Loci without a position, and sequence names absent from the GFF,
#' are reported at verbose >= 1.
#'
#' @param x A SNP genlight object with mapped loci. Must contain per-locus
#'   x$chromosome and x$position. [required]
#' @param gff.file Path to a GFF3 file, plain or gzip-compressed. A path
#'   without the .gz extension is also accepted when only the .gz file
#'   exists. [required]
#' @param loci Character vector of locus names to map, matching locNames(x).
#'   [default NULL, all loci with a chromosome and position]
#' @param include_types Character vector of GFF types to treat as "gene"
#'   features. [default c("gene","pseudogene")]
#' @param fallback_to_mrna Logical. If no rows match include_types, use
#'   transcript features c("mRNA","transcript") as proxies. [default TRUE]
#' @param save2tmp Logical: save the result table to tempdir() (retrievable
#'   with gl.list.reports / gl.print.reports). [default FALSE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'   brief progress messages; 3, progress and results summary; 5, full report
#'   [default 2, unless specified using gl.set.verbosity].
#'
#' @return A data.table with one row per locus that has a chromosome and
#'   position, and columns:
#'   locus, chrom, pos, gene_start, gene_end, gene_type, gene_id, gene_name,
#'   gene_symbol, gene_product, gene_attributes, distance_bp, nearest_side,
#'   gene_strand, relative_position.
#'   `distance_bp` is the absolute distance in bp; `nearest_side` is "inside",
#'   "left" (locus < gene_start) or "right" (locus > gene_end) in coordinate
#'   space; `relative_position` is "inside", "upstream" or "downstream" with
#'   respect to the gene strand (NA when the strand is unknown).
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # A small GFF3 with two genes on one sequence
#' gff <- tempfile(fileext = ".gff")
#' writeLines(c(
#'   "##gff-version 3",
#'   paste("chr1", "src", "gene", 100, 200, ".", "+", ".",
#'         "ID=gene-A;Name=A", sep = "\t"),
#'   paste("chr1", "src", "gene", 500, 800, ".", "-", ".",
#'         "ID=gene-B;Name=B", sep = "\t")
#' ), gff)
#' x <- testset.gl
#' x@chromosome <- factor(c("chr1", "chr1", rep(NA, nLoc(x) - 2)))
#' x@position <- c(150L, 300L, rep(NA_integer_, nLoc(x) - 2))
#' gl.find.genes.for.loci(x, gff.file = gff, loci = locNames(x)[1:2])
#'
#' @importFrom ape read.gff
#' @importFrom data.table as.data.table data.table setkey foverlaps :=
#' @importFrom stringr str_match str_detect
#'
#' @export
gl.find.genes.for.loci <- function(x,
                                   gff.file,
                                   loci = NULL,
                                   include_types = c("gene","pseudogene"),
                                   fallback_to_mrna = TRUE,
                                   save2tmp = FALSE,
                                   verbose = NULL) {

  # Avoid R CMD check NOTES for data.table NSE vars
  start <- end <- strand <- ID <- Name <- gene_k <- gene_id <- product <-
    Parent <- gene_sym <- type <- seqid <- gene_seqid <- gene_start <-
    gene_end <- NA
  chrom <- locus <- pos <- gene_type <- gene_name <- gene_symbol <-
    gene_product <- gene_attributes <- gene_strand <- relative_position <- NA
  join_coord <- which_end <- gene_mid <- gene_len <- distance_bp <-
    nearest_side <- parent_gene <- NA

  # VERBOSITY ------------------------------------------------------------------
  verbose <- gl.check.verbosity(verbose)
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECKS ---------------------------------------------------------------------
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  if (missing(gff.file) || !is.character(gff.file) || length(gff.file) != 1) {
    stop(error("Argument 'gff.file' must be a length-1 character path.\n"))
  }

  if (is.null(x$chromosome) || is.null(x$position)) {
    stop(error("Input 'x' must contain per-locus 'chromosome' and 'position'.\n"))
  }
  if (length(x$chromosome) != nLoc(x) || length(x$position) != nLoc(x)) {
    stop(error("Lengths of 'x$chromosome' and 'x$position' must equal nLoc(x).\n"))
  }

  all_names <- locNames(x)
  if (is.null(loci)) {
    loci <- all_names
  }
  if (!is.character(loci) || length(loci) < 1) {
    stop(error("Argument 'loci' must be NULL or a non-empty character vector.\n"))
  }

  empty_out <- function() {
    data.table::data.table(
      locus=character(), chrom=character(), pos=integer(),
      gene_start=integer(), gene_end=integer(),
      gene_type=character(), gene_id=character(), gene_name=character(),
      gene_symbol=character(), gene_product=character(),
      gene_attributes=character(), distance_bp=integer(),
      nearest_side=character(), gene_strand=character(),
      relative_position=character()
    )
  }

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
    return(empty_out())
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

  # ape::read.gff reads gzip through an R connection, so a .gz path passed
  # directly is handled by the first branch.
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
    return(empty_out())
  }

  # PRODUCT FROM CHILD FEATURES ------------------------------------------------
  # NCBI and Ensembl GFF3 put 'product' on mRNA, CDS or exon rows, not on the
  # gene row. Walk one or two Parent levels up from every row with a product
  # and give each gene the first product found among its descendants.
  if (any(is.na(gene_feats$product))) {
    parent_of <- gff_dt[!is.na(ID) & !is.na(Parent), list(ID, Parent)]
    parent_of <- parent_of[!duplicated(ID)]
    with_prod <- gff_dt[!is.na(product) & !is.na(Parent), list(Parent, product)]
    with_prod[, parent_gene := ifelse(Parent %in% gene_feats$ID, Parent,
                                      parent_of$Parent[match(Parent, parent_of$ID)])]
    with_prod <- with_prod[!is.na(parent_gene) & parent_gene %in% gene_feats$ID]
    with_prod <- with_prod[!duplicated(parent_gene)]
    gene_feats[is.na(product),
               product := with_prod$product[match(ID, with_prod$parent_gene)]]
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
    gene_attributes = attributes,
    gene_strand = as.character(strand)
  )]
  gene_iv <- gene_iv[!is.na(gene_seqid) & !is.na(gene_start) & !is.na(gene_end)]
  gene_iv[!(gene_strand %in% c("+", "-")), gene_strand := NA_character_]
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
  no_position <- loci_dt[is.na(chrom) | is.na(start), locus]
  loci_dt <- loci_dt[!is.na(chrom) & !is.na(start)]
  data.table::setkey(loci_dt, chrom, start, end)

  if (length(no_position) && verbose >= 1) {
    cat(warn(paste0("  ", length(no_position), " of ", length(loci),
                    " requested loci have no chromosome or position and are",
                    " not reported (e.g. ", no_position[1], ")\n")))
  }

  if (nrow(loci_dt) == 0L) {
    if (verbose >= 1) cat(report("No mappable loci (missing chrom/pos). Returning empty table.\n"))
    return(empty_out())
  }

  # Sequence names that never occur in the GFF give NA genes; say so.
  seq_missing <- setdiff(unique(loci_dt$chrom), unique(gene_iv$gene_seqid))
  if (length(seq_missing) && verbose >= 1) {
    n_affected <- sum(loci_dt$chrom %in% seq_missing)
    cat(warn(paste0("  ", length(seq_missing), " sequence name(s) of ",
                    n_affected, " loci do not occur in the GFF (e.g. ",
                    seq_missing[1], "); these loci get NA gene columns\n")))
    if (length(seq_missing) == length(unique(loci_dt$chrom))) {
      cat(warn(paste0("  No locus sequence name matches the GFF: check that",
                      " x$chromosome uses the same names as the GFF (e.g. ",
                      gene_iv$gene_seqid[1], ")\n")))
    }
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
      gene_symbol, gene_product, gene_attributes, distance_bp, nearest_side,
      gene_strand
    )]
  } else {
    best_overlap_out <- empty_out()[, list(
      locus, chrom, pos, gene_start, gene_end, gene_type, gene_id, gene_name,
      gene_symbol, gene_product, gene_attributes, distance_bp, nearest_side,
      gene_strand
    )]
  }

  overlapped_loci <- unique(best_overlap_out$locus)
  non_overlap <- loci_dt[!(locus %chin% overlapped_loci)]

  # 2) NEAREST BY INTERVAL EDGES (non-overlapping loci) -----------------------
  # Build a table of gene boundary points (starts and ends)
  gene_pts <- rbind(
    gene_iv[, list(gene_seqid, join_coord = gene_start,
                gene_start, gene_end, gene_type, gene_id, gene_name,
                gene_symbol, gene_product, gene_attributes, gene_strand,
                which_end = "start")],
    gene_iv[, list(gene_seqid, join_coord = gene_end,
                gene_start, gene_end, gene_type, gene_id, gene_name,
                gene_symbol, gene_product, gene_attributes, gene_strand,
                which_end = "end")]
  )
  data.table::setkey(gene_pts, gene_seqid, join_coord)

  if (nrow(non_overlap)) {
    points <- non_overlap[, list(chrom, join_coord = pos, pos, locus)]
    data.table::setkey(points, chrom, join_coord)

    # Rolling joins to the nearest gene edge at or before the locus
    # (roll = Inf) and at or after it (roll = -Inf), on the same sequence.
    # In the joined table the join columns keep the names of gene_pts
    # (gene_seqid, join_coord) but hold the values of points, and the other
    # columns of points (pos, locus) keep their own names; the distance is
    # therefore taken from the gene interval. Taking both sides lets the
    # documented tie-break decide between two genes at equal distance.
    candidates <- rbind(
      gene_pts[points, roll = Inf, nomatch = NA],
      gene_pts[points, roll = -Inf, nomatch = NA]
    )
    candidates <- candidates[!is.na(gene_start)]

    if (nrow(candidates)) {
      candidates[, distance_bp := as.integer(pmin(abs(pos - gene_start),
                                                  abs(pos - gene_end)))]
      candidates[, gene_mid := (gene_start + gene_end)/2]
      candidates[, gene_len := gene_end - gene_start]
      nearest <- candidates[order(locus, distance_bp, abs(pos - gene_mid),
                                  gene_len, gene_id), .SD[1L], by = locus]
      non_out <- nearest[, list(
        locus       = locus,
        chrom       = gene_seqid,
        pos         = pos,
        gene_start  = gene_start,
        gene_end    = gene_end,
        gene_type   = gene_type,
        gene_id     = gene_id,
        gene_name   = gene_name,
        gene_symbol = gene_symbol,
        gene_product= gene_product,
        gene_attributes = gene_attributes,
        distance_bp = distance_bp,
        nearest_side= ifelse(pos < gene_start, "left", "right"),
        gene_strand = gene_strand
      )]
    } else {
      non_out <- best_overlap_out[0L]
    }

    # Loci on a sequence without any gene feature: keep the row, NA gene.
    no_gene <- non_overlap[!(locus %chin% non_out$locus)]
    if (nrow(no_gene)) {
      non_out <- rbind(non_out, no_gene[, list(
        locus, chrom, pos,
        gene_start = NA_integer_, gene_end = NA_integer_,
        gene_type = NA_character_, gene_id = NA_character_,
        gene_name = NA_character_, gene_symbol = NA_character_,
        gene_product = NA_character_, gene_attributes = NA_character_,
        distance_bp = NA_integer_, nearest_side = NA_character_,
        gene_strand = NA_character_
      )])
    }
  } else {
    non_out <- best_overlap_out[0L]
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

  # Position relative to the gene strand: "left" of a plus-strand gene and
  # "right" of a minus-strand gene are both upstream.
  out[, relative_position := ifelse(
    is.na(nearest_side) | is.na(gene_strand), NA_character_,
    ifelse(nearest_side == "inside", "inside",
           ifelse((nearest_side == "left") == (gene_strand == "+"),
                  "upstream", "downstream")))]
  out[nearest_side %in% "inside", relative_position := "inside"]

  if (verbose >= 3) {
    cat(report(paste0(
      "  Assigned genes for ", nrow(out), " loci: ",
      sum(out$nearest_side %in% "inside"), " inside a gene, ",
      sum(out$nearest_side %in% c("left", "right")), " near a gene, ",
      sum(is.na(out$nearest_side)), " with no gene on their sequence\n")))
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

  # FLAG SCRIPT END ------------------------------------------------------------
  if (verbose > 0) {
    cat(report("Completed:", funname, "\n"))
  }

  # END ------------------------------------------------------------------------
  return(out)
}
