#' @name gl.find.loci.in.genes
#'
#' @title Find loci that fall within genes matching a pattern (from a GFF)
#'
#' @family annotation and mapping helpers
#'
#' @description
#' Given a genlight object with mapped loci (chromosome and position) and a
#' gene annotation file (GFF3, plain or gzip-compressed), this function
#' returns the names of the loci that fall inside genes whose annotation
#' matches a pattern (e.g. "MHC", "major histocompatibility").
#'
#' @details
#' The pattern is tested against the whole attributes column of every GFF
#' row. A gene or pseudogene is selected when its own row matches, or when
#' any of its descendant features matches (mRNA, lnc_RNA, CDS, exon and so
#' on, followed up through their Parent attribute, or linked by the gene=
#' key). NCBI and Ensembl GFF3 files put the product description on these
#' child rows, not on the gene row.
#'
#' Because the whole attribute string is searched, a short pattern can match
#' more than intended: "TAP" also matches WTAP, METAP1 or STAP1, and "ID"
#' matches every feature. Use word boundaries (e.g. "\\bTAP1\\b") and
#' "(?i)" for case-insensitive matching.
#'
#' Loci are placed by x$chromosome and x$position, which must be the
#' sequence name used in the first column of the GFF and a position in
#' genome coordinates. dartR sets x$position to the SNP position within the
#' tag when the data are read, so replace it with the genome position first
#' (e.g. the ChromPos locus metric plus SnpPosition). Loci without a
#' position, and sequence names absent from the GFF, are reported at
#' verbose >= 1.
#'
#' To learn which gene each locus falls in, pass the result to
#' \code{\link{gl.find.genes.for.loci}}.
#'
#' @param x Name of the genlight object containing SNP or SilicoDArT data,
#'   with per-locus x$chromosome and x$position [required].
#' @param gff.file Path to a GFF3 file, plain or gzip-compressed. A path
#'   without the .gz extension is also accepted when only the .gz file
#'   exists [required].
#' @param gene Regular expression identifying the target genes, e.g.
#'   "(?i)major histocompatibility|\\bMHC\\b" [required].
#' @param save2tmp Logical: save the table of loci and the genes they fall
#'   in to tempdir() (retrievable with gl.list.reports and
#'   gl.print.reports) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'   brief progress messages; 3, progress and results summary; 5, full report
#'   [default 2, unless specified using gl.set.verbosity].
#'
#' @return A character vector of the names of the loci that fall inside the
#'   matching genes, in genome order.
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # A small GFF3 with one MHC gene and one other gene on one sequence
#' gff <- tempfile(fileext = ".gff")
#' writeLines(c(
#'   "##gff-version 3",
#'   paste("chr1", "src", "gene", 100, 200, ".", "+", ".",
#'         "ID=gene-A;Name=A", sep = "\t"),
#'   paste("chr1", "src", "mRNA", 100, 200, ".", "+", ".",
#'         "ID=rna-A;Parent=gene-A;product=MHC class I antigen", sep = "\t"),
#'   paste("chr1", "src", "gene", 500, 800, ".", "-", ".",
#'         "ID=gene-B;Name=B;description=kinase", sep = "\t")
#' ), gff)
#' x <- testset.gl
#' x@chromosome <- factor(c("chr1", "chr1", rep(NA, nLoc(x) - 2)))
#' x@position <- c(150L, 600L, rep(NA_integer_, nLoc(x) - 2))
#' gl.find.loci.in.genes(x, gff.file = gff,
#'                       gene = "(?i)major histocompatibility|\\bMHC\\b")
#'
#' @importFrom ape read.gff
#' @importFrom data.table as.data.table data.table setkey foverlaps :=
#' @importFrom stringr str_match str_detect
#'
#' @export
gl.find.loci.in.genes <- function(x,
                                  gff.file,
                                  gene,
                                  save2tmp = FALSE,
                                  verbose = NULL) {

  # Avoid R CMD check NOTES for data.table NSE vars
  ID <- Name <- Parent <- gene_k <- chrom <- end <- locus <- seqid <-
    start <- type <- is_match <- gene_id <- gene_name <- gene_type <-
    i.start <- NA

  # SET VERBOSITY --------------------------------------------------------------
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START ----------------------------------------------------------
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK DATATYPE -------------------------------------------------------------
  datatype <- utils.check.datatype(x, accept = c("SNP", "SilicoDArT"),
                                   verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING -------------------------------------------
  if (missing(gff.file) || !is.character(gff.file) || length(gff.file) != 1 ||
      is.na(gff.file)) {
    stop(error("Argument 'gff.file' must be a length-1 character path.\n"))
  }
  if (missing(gene) || !is.character(gene) || length(gene) != 1 ||
      is.na(gene) || !nzchar(gene)) {
    stop(error("Argument 'gene' must be a single non-empty character pattern.\n"))
  }
  if (is.null(x$chromosome) || is.null(x$position)) {
    stop(error("Input 'x' must contain per-locus 'chromosome' and 'position'.\n"))
  }
  if (length(x$chromosome) != nLoc(x) || length(x$position) != nLoc(x)) {
    stop(error("Lengths of 'x$chromosome' and 'x$position' must equal nLoc(x).\n"))
  }

  gff_plain <- gff.file
  gff_gz    <- paste0(gff_plain, ".gz")
  has_plain <- file.exists(gff_plain)
  has_gz    <- file.exists(gff_gz)
  if (!has_plain && !has_gz) {
    stop(error(paste0("Cannot find '", gff_plain, "' or compressed '", gff_gz,
                      "'.\n")))
  }

  # Features treated as gene intervals
  gene_types <- c("gene", "pseudogene")

  extract_attr <- function(x_attr, key) {
    # Value of 'key' in a GFF attributes string "key1=val1;key2=val2"; NA if
    # the key is absent
    m <- stringr::str_match(x_attr, paste0("(^|;)", key, "=([^;]+)"))
    m[, 3]
  }

  # DO THE JOB -----------------------------------------------------------------
  # Load GFF annotation (gz or plain). ape::read.gff reads gzip through an R
  # connection, so a .gz path passed directly is handled by the first branch.
  if (verbose >= 2) cat(report("  Loading GFF file and parsing attributes\n"))
  gff_tab <- if (has_plain) {
    ape::read.gff(gff_plain)
  } else {
    ape::read.gff(gzfile(gff_gz))
  }
  gff_dt <- data.table::as.data.table(gff_tab)
  gff_dt[, `:=`(start = as.integer(start), end = as.integer(end))]

  gff_dt[, ID     := extract_attr(attributes, "ID")]
  gff_dt[, Parent := extract_attr(attributes, "Parent")]
  gff_dt[, Name   := extract_attr(attributes, "Name")]
  gff_dt[, gene_k := extract_attr(attributes, "gene")]

  # Rows whose attributes match the pattern, at any feature level
  gff_dt[, is_match := !is.na(attributes) &
           stringr::str_detect(attributes, gene)]
  hits_dt <- gff_dt[is_match == TRUE]

  # Follow Parent links from matching child rows (mRNA, lnc_RNA, CDS, exon)
  # up to their gene. A row can list several parents separated by commas.
  gene_ids_all <- unique(gff_dt[type %in% gene_types & !is.na(ID), ID])
  parent_map <- unique(gff_dt[!is.na(ID) & !is.na(Parent), list(ID, Parent)],
                       by = "ID")
  parent_of <- stats::setNames(parent_map$Parent, parent_map$ID)
  split_ids <- function(p) unlist(strsplit(p[!is.na(p)], ",", fixed = TRUE))

  matched_ids <- character(0)
  todo <- split_ids(hits_dt[!(type %in% gene_types), Parent])
  seen <- character(0)
  while (length(todo) > 0) {
    todo <- setdiff(unique(todo), seen)
    seen <- c(seen, todo)
    matched_ids <- c(matched_ids, todo[todo %in% gene_ids_all])
    todo <- split_ids(parent_of[todo[!(todo %in% gene_ids_all)]])
  }

  # GFFs without Parent links tie child rows to genes by the gene= key only
  matched_keys <- unique(stats::na.omit(hits_dt[!(type %in% gene_types),
                                                gene_k]))

  genes_dt <- gff_dt[type %in% gene_types & (
    is_match |
      (!is.na(ID) & ID %in% matched_ids) |
      (!is.na(gene_k) & gene_k %in% matched_keys) |
      (!is.na(Name) & Name %in% matched_keys)
  )]

  if (nrow(genes_dt) == 0 && verbose >= 1) {
    cat(warn("  Warning: no gene or pseudogene in the GFF matches the pattern",
             paste0("'", gene, "'."),
             "Check the pattern, or the attribute keys in the GFF.\n"))
  }

  # Gene interval table
  genes_iv <- genes_dt[, list(chrom = as.character(seqid),
                              start = as.integer(start),
                              end   = as.integer(end),
                              gene_id = ID,
                              gene_name = data.table::fcoalesce(gene_k, Name),
                              gene_type = as.character(type))]
  genes_iv <- genes_iv[!is.na(chrom) & !is.na(start) & !is.na(end)]
  data.table::setkey(genes_iv, chrom, start, end)

  # Loci table (1-bp intervals at 'position')
  loci_dt <- data.table::data.table(
    chrom = as.character(x$chromosome),
    start = as.integer(x$position),
    end   = as.integer(x$position),
    locus = locNames(x)
  )
  no_position <- loci_dt[is.na(chrom) | is.na(start), locus]
  loci_dt <- loci_dt[!is.na(chrom) & !is.na(start)]
  data.table::setkey(loci_dt, chrom, start, end)

  if (length(no_position) && verbose >= 1) {
    cat(warn(paste0("  ", length(no_position), " of ", nLoc(x),
                    " loci have no chromosome or position and are not",
                    " tested (e.g. ", no_position[1], ")\n")))
  }

  # Sequence names that never occur in the GFF cannot overlap any gene
  seq_missing <- setdiff(unique(loci_dt$chrom), unique(as.character(gff_dt$seqid)))
  if (length(seq_missing) && verbose >= 1) {
    n_affected <- sum(loci_dt$chrom %in% seq_missing)
    cat(warn(paste0("  ", length(seq_missing), " sequence name(s) of ",
                    n_affected, " loci do not occur in the GFF (e.g. ",
                    seq_missing[1], ")\n")))
    if (nrow(loci_dt) > 0 &&
        length(seq_missing) == length(unique(loci_dt$chrom))) {
      cat(warn(paste0("  No locus sequence name matches the GFF: check that",
                      " x$chromosome uses the names in the first column of",
                      " the GFF (e.g. ", unique(gff_dt$seqid)[1], ")\n")))
    }
  }

  # Overlap join: loci within any matching gene interval
  hits <- data.table::foverlaps(loci_dt, genes_iv, nomatch = 0L)
  loci_in_genes <- unique(hits$locus)

  if (verbose >= 3) {
    cat(report("  Genes matching the pattern:", nrow(genes_dt),
               paste0("(", sum(genes_dt$type == "pseudogene"),
                      " pseudogenes)\n")))
    cat(report("  Loci inside matching genes:", length(loci_in_genes), "\n"))
  }

  # SAVE (OPTIONAL) ------------------------------------------------------------
  if (isTRUE(save2tmp)) {
    out_tab <- hits[, list(locus, chrom, pos = i.start, gene_id, gene_name,
                           gene_type, gene_start = start, gene_end = end)]
    fn <- tempfile(pattern = "dartR_table_lociingenes_", fileext = ".rds")
    saveRDS(out_tab, file = fn)
    if (verbose >= 2) {
      cat(report("  Saved table to ", fn, " (via saveRDS)\n"))
      cat(report("  Retrieve with gl.list.reports() / gl.print.reports()\n"))
    }
  }

  # FLAG SCRIPT END ------------------------------------------------------------
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN ---------------------------------------------------------------------
  return(loci_in_genes)
}
