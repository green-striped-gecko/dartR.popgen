#' @name gl.find.loci.in.genes
#'
#' @title Find loci that fall within genes matching a pattern (from a GFF)
#'
#' @description
#' Given a genlight object with mapped loci (chromosome and position) and a
#' gene annotation file (GFF, plain or gz), this function returns the locus
#' names whose genomic positions overlap features of type "gene" whose
#' attributes match a user-supplied pattern (e.g., "MHC",
#' "major histocompatibility").
#'
#' @param x Name of the genlight object containing SNP data [required].
#' @param gff.file Path to a GFF3 file (plain or with a companion .gz)
#'   [required].
#' @param gene Character pattern to detect target genes. Interpreted as a
#'   regular expression, case-insensitive matching is recommended via `(?i)`
#'   [required].
#' @param save2tmp Logical: save intermediate tables to tempdir for
#'   retrieval with gl.list.reports and gl.print.reports [default FALSE].
#' @param verbose Verbosity: 0=silent/fatal; 1=begin/end; 2=progress;
#'   3=progress+summary; 5=full report [default 2 or as set by
#'   `gl.set.verbosity()`].
#'
#' @details
#' The function parses the GFF "attributes" column to extract common keys
#' (e.g., `Name`, `gene`, `product`) and flags any "gene" features whose
#' attributes match the supplied `gene` pattern. It then uses interval
#' overlap to identify input loci that fall inside those genes.
#'
#' Required fields in `x` for overlap are per-locus chromosome and base
#' position, accessible as `x$chromosome` and `x$position`, and locus names
#' via `locNames(x)`.
#'
#' @return A character vector of locus names overlapping the matching gene
#'   intervals (in genomic coordinates).
#'
#' @author Luis Mijangos (post to
#'   \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' \dontrun{
#' # Regex for case-insensitive MHC:
#' mhc_loci <- gl.find.loci.in.genes(
#'   x         = testset.gl,
#'   gff.file  = "species.gff3",
#'   gene      = "(?i)major histocompatibility|\\bMHC\\b"
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
gl.find.loci.in.genes <- function(x,
                        gff.file,
                        gene,
                        save2tmp = FALSE,
                        verbose = NULL) {
  Name <- cds_gene<- chrom<- end<- product<- seqid <-start <-type<- NA
  
  first_present <- function(candidates, where_names) {
    for (nm in candidates) if (nm %in% where_names) return(nm)
    NA_character_
  }
  
  extract_attr <- function(x, key) {
    # Extracts the value for 'key' in a GFF attributes string "key1=val1;key2=val2"
    # Returns NA if key not present
    m <- stringr::str_match(x, paste0("(^|;)", key, "=([^;]+)"))
    m[, 3]
  }
  
  # ---------------------------------------------------------------------------
  # 3) Load GFF annotation (gz or plain)
  # ---------------------------------------------------------------------------
  cat("LOADING GFF FILE...\n")
  
  gff_plain <- gff.file
  gff_gz    <- paste0(gff_plain, ".gz")
  
  if (file.exists(gff_plain)) {
    gff <- read.gff(gff_plain)
  } else if (file.exists(gff_gz)) {
    gff <- read.gff(gzfile(gff_gz))
  }
  
  # Convert to data.table and ensure numeric coords
  gff_dt <- as.data.table(gff)
  # Some GFFs load 'start'/'end' as character; coerce safely
  gff_dt[, `:=`(start = as.integer(start), end = as.integer(end))]
  
  # ---------------------------------------------------------------------------
  # 4) Parse attributes and identify MHC genes
  # ---------------------------------------------------------------------------
  # Extract useful keys for 'gene' and 'gene features'
  gff_dt[, Name    := extract_attr(attributes, "Name")]
  gff_dt[, gene    := extract_attr(attributes, "gene")]
  gff_dt[, product := extract_attr(attributes, "product")]
  
  # Candidate MHC pattern: "MHC" or "major histocompatibility" (case-insensitive)
  mhc_pat <- gene
  
  # CDS entries: gather gene names hinted as MHC by product/attributes
  cds_dt <- gff_dt[type == "CDS"]
  cds_dt[, cds_gene := extract_attr(attributes, "gene")]
  # Safe logicals (handle NAs explicitly)
  prod_is_mhc <- !is.na(cds_dt$product)    & str_detect(cds_dt$product, mhc_pat)
  attr_is_mhc <- !is.na(cds_dt$attributes) & str_detect(cds_dt$attributes, mhc_pat)
  
  mhc_gene_names <- unique(na.omit(c(
    cds_dt[prod_is_mhc, cds_gene],
    cds_dt[attr_is_mhc, cds_gene]
  )))
  
  # Gene features: keep any gene whose attributes match MHC pattern OR whose
  # name/gene appears among the MHC-like CDS gene names
  genes_dt <- gff_dt[type == "gene"]
  genes_dt[, gene := data.table::fcoalesce(gene, Name)]  # prefer 'gene'; fall back to 'Name'
  
  mhc_genes_dt <- genes_dt[
    !is.na(attributes) & (
      str_detect(attributes, mhc_pat) |
        (!is.na(gene) & gene %in% mhc_gene_names) |
        (!is.na(Name) & Name %in% mhc_gene_names)
    )
  ]
  
  # If nothing found via attributes, this may help diagnose naming mismatches
  if (nrow(mhc_genes_dt) == 0) {
    warning("No MHC genes detected via attributes. Consider widening 'mhc_pat' or checking attribute keys.")
  }
  
  # ---------------------------------------------------------------------------
  # 5) Interval overlap: locate BLAST-mapped loci that fall within MHC genes
  # ---------------------------------------------------------------------------
  # Build MHC interval table
  mhc_dt <- mhc_genes_dt[, list(chrom = as.character(seqid),
                             start = as.integer(start),
                             end   = as.integer(end))]
  # Guard against NAs
  mhc_dt <- mhc_dt[!is.na(chrom) & !is.na(start) & !is.na(end)]
  setkey(mhc_dt, chrom, start, end)
  
  # Build loci table (1-bp intervals at 'position')
  loci_dt <- data.table(
    chrom = as.character(x$chromosome),
    start = as.integer(x$position),
    end   = as.integer(x$position),
    locus = locNames(x)
  )
  loci_dt <- loci_dt[!is.na(chrom) & !is.na(start)]
  setkey(loci_dt, chrom, start, end)
  
  # Overlap join: loci within any MHC interval
  hits <- foverlaps(loci_dt, mhc_dt, nomatch = 0L)
  
  # Unique locus IDs that map inside MHC intervals
  mhc_genes_db <- unique(hits$locus)
  
  # ---------------------------------------------------------------------------
  # 6) Save results (optional) and print summary
  
  cat(sprintf("MHC genes detected: %d\n", nrow(mhc_genes_dt)))
  cat(sprintf("Loci overlapping MHC intervals: %d\n", length(mhc_genes_db)))
  
  # Return/print vector of loci (as in your original script)
  mhc_genes_db
}
