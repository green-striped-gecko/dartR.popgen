# Regression tests for gl.find.loci.in.genes (function-review campaign).
#
# Synthetic GFF3 with genes whose matching text sits on the gene row, on a
# CDS, only on an mRNA, or on a pseudogene row. Expected loci are worked out
# by hand from the documented definition: a gene or pseudogene is selected
# when its own row or any descendant row matches the pattern.
#
# Before the review (dartR.popgen 1.2.2) text on mRNA rows and pseudogenes
# was ignored, verbose was not read, a missing GFF fell back to any object
# called 'gff' in the workspace, and sequence-name mismatches were silent.

gff_line <- function(...) paste(..., sep = "\t")

write_test_gff <- function(env = parent.frame()) {
  gff <- withr::local_tempfile(fileext = ".gff", .local_envir = env)
  writeLines(c(
    "##gff-version 3",
    # MHC1: pattern text on the gene row's children (CDS product)
    gff_line("chr1", "src", "gene", 100, 200, ".", "+", ".",
             "ID=gene-MHC1;Name=MHC1;gene=MHC1"),
    gff_line("chr1", "src", "mRNA", 100, 200, ".", "+", ".",
             "ID=rna-1;Parent=gene-MHC1;gene=MHC1;product=MHC class I antigen"),
    gff_line("chr1", "src", "CDS", 100, 200, ".", "+", "0",
             "ID=cds-1;Parent=rna-1;gene=MHC1;product=MHC class I antigen"),
    # LOC1: pattern text only on the mRNA row (as for NCBI lncRNA or when
    # the CDS product differs from the transcript product)
    gff_line("chr1", "src", "gene", 500, 800, ".", "-", ".",
             "ID=gene-LOC1;Name=LOC1;gene=LOC1"),
    gff_line("chr1", "src", "mRNA", 500, 800, ".", "-", ".",
             paste0("ID=rna-2;Parent=gene-LOC1;gene=LOC1;",
                    "product=major histocompatibility complex class II")),
    gff_line("chr1", "src", "gene", 1000, 1100, ".", "+", ".",
             "ID=gene-K;Name=K;gene=K;description=kinase"),
    # P: a pseudogene feature whose own row matches the pattern
    gff_line("chr2", "src", "pseudogene", 10, 90, ".", "+", ".",
             "ID=gene-P;Name=P;gene=P;description=MHC class I pseudogene"),
    gff_line("chr2", "src", "gene", 200, 300, ".", "+", ".",
             "ID=gene-TAP;Name=TAP1;gene=TAP1")
  ), gff)
  gff
}

# Loci 1-7: inside MHC1, inside LOC1, inside K, inside pseudogene P,
# inside TAP1, no position, on the MHC1 start coordinate.
make_test_gl <- function() {
  x <- dartR.data::testset.gl
  n <- nLoc(x)
  chrom <- rep(NA_character_, n)
  pos <- rep(NA_integer_, n)
  chrom[1:7] <- c("chr1", "chr1", "chr1", "chr2", "chr2", "chr1", "chr1")
  pos[1:7] <- c(150L, 600L, 1050L, 50L, 250L, NA_integer_, 100L)
  x@chromosome <- factor(chrom)
  x@position <- pos
  x
}

pat_mhc <- "(?i)major histocompatibility|\\bMHC\\b"

test_that("genes matched on any feature level, pseudogenes included", {
  gff <- write_test_gff()
  x <- make_test_gl()
  res <- gl.find.loci.in.genes(x, gff.file = gff, gene = pat_mhc, verbose = 0)
  expect_type(res, "character")
  # MHC1 via CDS/mRNA (loci 1 and 7, start coordinate is inside), LOC1 via
  # its mRNA only (locus 2), pseudogene P via its own row (locus 4)
  expect_setequal(res, locNames(x)[c(1, 2, 4, 7)])
  # genome order: chr1 100, 150, 600, then chr2 50
  expect_equal(res, locNames(x)[c(7, 1, 2, 4)])
})

test_that("pattern is matched against the whole attribute string", {
  gff <- write_test_gff()
  x <- make_test_gl()
  res <- gl.find.loci.in.genes(x, gff.file = gff, gene = "TAP", verbose = 0)
  expect_equal(res, locNames(x)[5])
  # 'ID' occurs in every attribute string, so every gene matches
  res_all <- gl.find.loci.in.genes(x, gff.file = gff, gene = "ID", verbose = 0)
  expect_setequal(res_all, locNames(x)[c(1, 2, 3, 4, 5, 7)])
})

test_that("genes linked only by the gene= key are found", {
  gff <- withr::local_tempfile(fileext = ".gff")
  writeLines(c(
    "##gff-version 3",
    gff_line("chr1", "src", "gene", 100, 200, ".", "+", ".",
             "ID=g1;Name=MHC1;gene=MHC1"),
    gff_line("chr1", "src", "CDS", 100, 200, ".", "+", "0",
             "ID=c1;gene=MHC1;product=MHC class I antigen")
  ), gff)
  x <- make_test_gl()
  expect_setequal(gl.find.loci.in.genes(x, gff, gene = "MHC", verbose = 0),
                  locNames(x)[c(1, 7)])
})

test_that("verbose = 0 is silent; verbose = 3 reports the summary", {
  gff <- write_test_gff()
  x <- make_test_gl()
  expect_silent(gl.find.loci.in.genes(x, gff.file = gff, gene = "TAP",
                                      verbose = 0))
  out <- capture.output(
    gl.find.loci.in.genes(x, gff.file = gff, gene = pat_mhc, verbose = 3)
  )
  expect_true(any(grepl("Genes matching the pattern: 3 \\(1 pseudogenes\\)",
                        out)))
  expect_true(any(grepl("Loci inside matching genes: 4", out)))
  expect_true(any(grepl("Completed:", out)))
})

test_that("no matching gene gives a generic message and an empty result", {
  gff <- write_test_gff()
  x <- make_test_gl()
  out <- capture.output(
    res <- gl.find.loci.in.genes(x, gff.file = gff, gene = "zzz", verbose = 1)
  )
  expect_identical(res, character(0))
  expect_true(any(grepl("no gene or pseudogene in the GFF matches", out)))
  expect_false(any(grepl("MHC", out)))
})

test_that("a missing GFF file errors, even with a 'gff' object in scope", {
  x <- make_test_gl()
  assign("gff", ape::read.gff(write_test_gff()), envir = globalenv())
  withr::defer(rm("gff", envir = globalenv()))
  expect_error(
    gl.find.loci.in.genes(x, gff.file = "no_such_file.gff", gene = "TAP",
                          verbose = 0),
    "Cannot find 'no_such_file.gff' or compressed 'no_such_file.gff.gz'"
  )
})

test_that("invalid arguments fail early with clear errors", {
  gff <- write_test_gff()
  x <- make_test_gl()
  expect_error(gl.find.loci.in.genes(x, gff, gene = c("MHC", "TAP"),
                                     verbose = 0), "single non-empty")
  expect_error(gl.find.loci.in.genes(x, gff, gene = NA_character_,
                                     verbose = 0), "single non-empty")
  expect_error(gl.find.loci.in.genes(x, gff, verbose = 0), "single non-empty")
  expect_error(gl.find.loci.in.genes(x, gene = "MHC", verbose = 0),
               "gff.file")
  y <- dartR.data::testset.gl
  y@chromosome <- NULL
  expect_error(gl.find.loci.in.genes(y, gff, gene = "MHC", verbose = 0),
               "chromosome")
})

test_that("SilicoDArT data are accepted", {
  gff <- write_test_gff()
  s <- dartR.data::testset.gs
  s@chromosome <- factor(c("chr1", rep(NA, nLoc(s) - 1)))
  s@position <- c(150L, rep(NA_integer_, nLoc(s) - 1))
  expect_equal(gl.find.loci.in.genes(s, gff, gene = "MHC", verbose = 0),
               locNames(s)[1])
})

test_that("unmapped loci and sequence-name mismatch are reported", {
  gff <- write_test_gff()
  x <- make_test_gl()
  out <- capture.output(
    gl.find.loci.in.genes(x, gff.file = gff, gene = "TAP", verbose = 1)
  )
  expect_true(any(grepl(paste0(nLoc(x) - 6, " of ", nLoc(x),
                               " loci have no chromosome or position"), out)))
  x@chromosome <- factor(sub("chr", "NC_0", as.character(x@chromosome)))
  out <- capture.output(
    r <- gl.find.loci.in.genes(x, gff.file = gff, gene = "TAP", verbose = 1)
  )
  expect_identical(r, character(0))
  expect_true(any(grepl("2 sequence name\\(s\\) of 6 loci do not occur", out)))
  expect_true(any(grepl("No locus sequence name matches the GFF", out)))
  # silent at verbose = 0
  expect_silent(gl.find.loci.in.genes(x, gff.file = gff, gene = "TAP",
                                      verbose = 0))
})

test_that("save2tmp writes the locus-gene table to tempdir()", {
  gff <- write_test_gff()
  x <- make_test_gl()
  before <- list.files(tempdir(), pattern = "^dartR_table_lociingenes_")
  gl.find.loci.in.genes(x, gff.file = gff, gene = pat_mhc, save2tmp = TRUE,
                        verbose = 0)
  new <- setdiff(list.files(tempdir(), pattern = "^dartR_table_lociingenes_"),
                 before)
  expect_length(new, 1)
  tab <- readRDS(file.path(tempdir(), new))
  withr::defer(unlink(file.path(tempdir(), new)))
  tab <- as.data.frame(tab)[order(tab$pos), ]
  expect_equal(tab$locus[tab$chrom == "chr1"], locNames(x)[c(7, 1, 2)])
  expect_equal(tab$gene_name, c("P", "MHC1", "MHC1", "LOC1"))
  expect_equal(tab$gene_type, c("pseudogene", "gene", "gene", "gene"))
})

test_that("gzip GFF works both as a direct path and as a companion file", {
  gff <- write_test_gff()
  x <- make_test_gl()
  gz <- withr::local_tempfile(fileext = ".gff.gz")
  con <- gzfile(gz, "w")
  writeLines(readLines(gff), con)
  close(con)
  expect_equal(
    gl.find.loci.in.genes(x, gff.file = gz, gene = "TAP", verbose = 0),
    locNames(x)[5]
  )
  expect_equal(
    gl.find.loci.in.genes(x, gff.file = sub("\\.gz$", "", gz), gene = "TAP",
                          verbose = 0),
    locNames(x)[5]
  )
})
