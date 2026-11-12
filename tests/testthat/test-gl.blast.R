# Characterization baseline for gl.blast (dartr-function-review, 2026-09-18).
# Snapshots what the function does today on testset.gl against a synthetic
# reference genome built from the object's own sequence tags: forty tags
# embedded exactly in contig1, thirty tags with two substitutions each in
# contig2, and a random contig3. The pre-review baseline pinned three
# defects (stale results after a failed run, ungated no-hit message, path
# string returned for a fasta query with no hit); those expectations were
# flipped to the approved behaviour recorded in
# function-review/reports/dartR.popgen/gl.blast.md (changes 1, 4, 5).
#
# The whole file is skipped when BLAST+ (blastn, makeblastdb) is not on PATH.

blast_available <- function() {
  nzchar(Sys.which("blastn")) && nzchar(Sys.which("makeblastdb"))
}

make_blast_fixture <- function(dir) {
  set.seed(1)
  gl <- dartR.data::testset.gl
  seqs <- as.character(gl@other$loc.metrics$TrimmedSequence)
  rnd <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE),
                           collapse = "")
  mut <- function(s, k) {
    v <- strsplit(s, "")[[1]]
    i <- sample(length(v), k)
    v[i] <- vapply(v[i], function(b) {
      sample(setdiff(c("A", "C", "G", "T"), b), 1)
    }, "")
    paste(v, collapse = "")
  }
  c1 <- paste0(rnd(200), paste0(seqs[1:40], rnd(200), collapse = ""))
  c2 <- paste0(rnd(200),
               paste0(vapply(seqs[41:70], mut, "", k = 2), rnd(200),
                      collapse = ""))
  c3 <- rnd(3000)
  genome <- file.path(dir, "genome.fasta")
  nohit <- file.path(dir, "genome_nohit.fasta")
  query <- file.path(dir, "query.fasta")
  writeLines(c(">contig1 chromosome 1", c1,
               ">contig2 plasmid", c2,
               ">contig3 unplaced", c3), genome)
  writeLines(c(">contig3 unplaced", c3), nohit)
  writeLines(c(rbind(paste0(">", locNames(gl)), seqs)), query)
  R.utils::gzip(genome, destname = paste0(genome, ".gz"),
                overwrite = TRUE, remove = FALSE)
  list(gl = gl, genome = genome, gz = paste0(genome, ".gz"),
       nohit = nohit, query = query, dir = dir)
}

blast_cols <- c("qseqid", "sacc", "stitle", "qseq", "sseq", "nident",
                "mismatch", "pident", "length", "evalue", "bitscore",
                "qstart", "qend", "sstart", "send", "gapopen", "gaps",
                "qlen", "slen", "PercentageOverlap")

test_that("genlight input: one hit per locus is merged into loc.metrics", {
  skip_if_not(blast_available(), "BLAST+ (blastn, makeblastdb) not on PATH")
  fx <- make_blast_fixture(withr::local_tempdir())
  expect_silent(res <- gl.blast(fx$gl, fx$genome, verbose = 0))
  lm <- res@other$loc.metrics
  lm0 <- fx$gl@other$loc.metrics

  expect_s4_class(res, "genlight")
  expect_equal(nLoc(res), 255)
  expect_equal(unique(ploidy(res)), 2)
  expect_equal(nrow(lm), 255)
  expect_equal(ncol(lm), ncol(lm0) + length(blast_cols))
  expect_setequal(setdiff(colnames(lm), colnames(lm0)), blast_cols)
  # locus order is preserved through the merge
  expect_equal(as.character(lm$AlleleID), as.character(lm0$AlleleID))
  expect_equal(lm$qseqid, 1:255)

  hit <- !is.na(lm$sacc)
  expect_equal(sum(hit), 59)
  expect_equal(as.vector(table(lm$sacc)), c(36, 23))
  expect_equal(names(table(lm$sacc)), c("contig1", "contig2"))
  # every kept hit satisfies the three documented thresholds
  expect_true(all(lm$pident[hit] >= 70))
  expect_true(all(lm$PercentageOverlap[hit] >= 0.8))
  expect_true(all(lm$bitscore[hit] >= 50))
  expect_equal(length(res@other$history),
               length(fx$gl@other$history) + 1)
  expect_output(gl.blast(fx$gl, fx$genome, verbose = 1),
                "59 of 255 sequences aligned after filtering")
})

test_that("a gzipped reference gives the same hits", {
  skip_if_not(blast_available(), "BLAST+ (blastn, makeblastdb) not on PATH")
  fx <- make_blast_fixture(withr::local_tempdir())
  a <- gl.blast(fx$gl, fx$genome, verbose = 0)
  b <- gl.blast(fx$gl, fx$gz, verbose = 0)
  expect_identical(b@other$loc.metrics$sacc, a@other$loc.metrics$sacc)
  expect_identical(b@other$loc.metrics$bitscore,
                   a@other$loc.metrics$bitscore)
})

test_that("paths with spaces are quoted for the shell", {
  skip_if_not(blast_available(), "BLAST+ (blastn, makeblastdb) not on PATH")
  dir <- file.path(withr::local_tempdir(), "dir with space")
  dir.create(dir)
  fx <- make_blast_fixture(dir)
  res <- gl.blast(fx$gl, fx$genome, verbose = 0)
  expect_equal(sum(!is.na(res@other$loc.metrics$sacc)), 59)
  res2 <- gl.blast(fx$query, fx$genome, verbose = 0)
  expect_equal(nrow(res2), 59)
})

test_that("fasta path input returns a data frame with one hit per query", {
  skip_if_not(blast_available(), "BLAST+ (blastn, makeblastdb) not on PATH")
  fx <- make_blast_fixture(withr::local_tempdir())
  res <- gl.blast(fx$query, fx$genome, verbose = 0)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 59)
  expect_setequal(colnames(res), blast_cols)
  expect_equal(anyDuplicated(res$qseqid), 0)
  expect_true(all(res$qseqid %in% locNames(fx$gl)))
})

test_that("no alignment: genlight comes back without BLAST columns, fasta path gives an empty data frame", {
  skip_if_not(blast_available(), "BLAST+ (blastn, makeblastdb) not on PATH")
  fx <- make_blast_fixture(withr::local_tempdir())
  expect_silent(res <- gl.blast(fx$gl, fx$nohit, verbose = 0))
  expect_equal(ncol(res@other$loc.metrics),
               ncol(fx$gl@other$loc.metrics))
  expect_equal(length(res@other$history),
               length(fx$gl@other$history) + 1)
  expect_output(gl.blast(fx$gl, fx$nohit, verbose = 1),
                "0 of 255 sequences aligned after filtering")
  expect_silent(res2 <- gl.blast(fx$query, fx$nohit, verbose = 0))
  expect_s3_class(res2, "data.frame")
  expect_equal(nrow(res2), 0)
  expect_setequal(colnames(res2), blast_cols)
})

test_that("nothing survives the R filter: no columns added, history appended, count reported", {
  skip_if_not(blast_available(), "BLAST+ (blastn, makeblastdb) not on PATH")
  fx <- make_blast_fixture(withr::local_tempdir())
  expect_output(res <- gl.blast(fx$gl, fx$genome, bitscore = 1e9,
                                verbose = 1),
                "0 of 255 sequences aligned after filtering")
  expect_equal(ncol(res@other$loc.metrics),
               ncol(fx$gl@other$loc.metrics))
  expect_equal(length(res@other$history),
               length(fx$gl@other$history) + 1)
})

test_that("second run on an annotated object suffixes the BLAST columns .x/.y", {
  skip_if_not(blast_available(), "BLAST+ (blastn, makeblastdb) not on PATH")
  fx <- make_blast_fixture(withr::local_tempdir())
  r1 <- gl.blast(fx$gl, fx$genome, verbose = 0)
  r2 <- gl.blast(r1, fx$genome, verbose = 0)
  lm <- r2@other$loc.metrics
  expect_equal(ncol(lm), ncol(fx$gl@other$loc.metrics) + 2 * 19 + 1)
  expect_true(all(c("sacc.x", "sacc.y") %in% colnames(lm)))
  expect_false("sacc" %in% colnames(lm))
})

test_that("invalid inputs and a failed tool stop instead of returning the previous run's hits", {
  skip_if_not(blast_available(), "BLAST+ (blastn, makeblastdb) not on PATH")
  fx <- make_blast_fixture(withr::local_tempdir())
  good <- gl.blast(fx$gl, fx$genome, verbose = 0)
  expect_equal(sum(!is.na(good@other$loc.metrics$sacc)), 59)
  expect_error(gl.blast(fx$gl, fx$genome, task = "blastx", verbose = 0),
               "should be one of")
  expect_error(gl.blast(file.path(fx$dir, "missing.fasta"), fx$genome,
                        verbose = 0),
               "query fasta file not found")
  expect_error(gl.blast(fx$gl, file.path(fx$dir, "missing_genome.fasta"),
                        verbose = 0),
               "reference genome file not found")
  expect_error(gl.blast(gl2gi(fx$gl, verbose = 0), fx$genome, verbose = 0),
               "must be a genlight object or the path")
  expect_error(gl.blast(fx$gl, fx$genome, Percentage_overlap = 2,
                        verbose = 0),
               "Percentage_overlap must be between 0 and")
  # a file makeblastdb cannot build a database from is reported as its failure
  empty <- file.path(fx$dir, "empty.fasta")
  file.create(empty)
  expect_error(gl.blast(fx$gl, empty, verbose = 0), "makeblastdb exited")
})
