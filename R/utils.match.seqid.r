# Matches locus sequence names to the sequence names of a GFF, for
# gl.find.genes.for.loci and gl.find.loci.in.genes.
#
# DArT reports a reference sequence as "<accession>_<description>", for
# example NC_041728.1_chromosome_1, where the GFF names it by the accession
# alone (NC_041728.1). A name that occurs in the GFF is kept; a name that does
# not is replaced by the longest GFF name it starts with, followed by "_", so
# chr10_x matches chr10 but never chr1. Names with no such GFF name are kept
# and reported by the caller as absent from the GFF.
# @param chrom Character vector of locus sequence names (NA allowed).
# @param seqs Character vector of the GFF sequence names.
# @param verbose Verbosity; at 2 or more the renaming is reported.
# @return chrom with the matched names replaced.

utils.match.seqid <- function(chrom, seqs, verbose = 0) {
  seqs <- unique(stats::na.omit(as.character(seqs)))
  miss <- setdiff(unique(stats::na.omit(chrom)), seqs)
  if (length(miss) == 0L || length(seqs) == 0L) {
    return(chrom)
  }
  prefixes <- paste0(seqs, "_")
  to <- vapply(miss, function(u) {
    s <- seqs[startsWith(u, prefixes)]
    if (length(s)) s[which.max(nchar(s))] else NA_character_
  }, character(1))
  to <- to[!is.na(to)]
  if (length(to) == 0L) {
    return(chrom)
  }
  if (verbose >= 2) {
    cat(report(paste0("  ", length(to), " sequence name(s) matched to the GFF",
                      " sequence they start with (e.g. ", names(to)[1],
                      " as ", to[[1]], ")\n")))
  }
  hit <- chrom %in% names(to)
  chrom[hit] <- unname(to[chrom[hit]])
  chrom
}
