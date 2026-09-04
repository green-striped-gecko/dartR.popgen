# Regression tests for gl.find.genes.for.loci (function-review campaign).
#
# Synthetic GFF3 files with a few genes on two sequences and loci placed
# inside, left of, right of, and away from those genes. Expected values are
# computed by hand from the documented definition (distance to the nearest
# gene edge; 0 when inside; side "left" when locus < gene_start, "right"
# when locus > gene_end; upstream/downstream from the gene strand).
#
# Before the fix (dartR.popgen 1.2.2) the nearest-gene branch returned a
# single row with an empty locus for every locus outside a gene, because the
# data.table roll join does not produce the `i.locus`, `i.chrom` and `i.pos`
# columns the code read, and the NSE placeholders (`i.pos <- NA`) were used
# instead.

gff_line <- function(...) paste(..., sep = "\t")

write_test_gff <- function(env = parent.frame()) {
  gff <- withr::local_tempfile(fileext = ".gff", .local_envir = env)
  writeLines(c(
    "##gff-version 3",
    gff_line("chr1", "src", "gene", 100, 200, ".", "+", ".", "ID=gene-A;Name=A"),
    gff_line("chr1", "src", "mRNA", 100, 200, ".", "+", ".",
             "ID=rna-A;Parent=gene-A"),
    gff_line("chr1", "src", "CDS", 100, 200, ".", "+", "0",
             "ID=cds-A;Parent=rna-A;product=kinase A"),
    gff_line("chr1", "src", "gene", 500, 800, ".", "-", ".", "ID=gene-B;Name=B"),
    gff_line("chr1", "src", "mRNA", 500, 800, ".", "-", ".",
             "ID=rna-B;Parent=gene-B;product=phosphatase B"),
    gff_line("chr1", "src", "gene", 1500, 1600, ".", "+", ".",
             "ID=gene-C;Name=C;product=own product C"),
    gff_line("chr2", "src", "gene", 50, 60, ".", "+", ".", "ID=gene-D;Name=D")
  ), gff)
  gff
}

make_test_gl <- function() {
  x <- dartR.data::testset.gl
  n <- nLoc(x)
  chrom <- rep(NA_character_, n)
  pos <- rep(NA_integer_, n)
  chrom[1:6] <- c("chr1", "chr1", "chr1", "chr1", "chr3", "chr1")
  pos[1:6] <- c(150L, 300L, 1000L, 1450L, 10L, NA_integer_)
  x@chromosome <- factor(chrom)
  x@position <- pos
  x
}

test_that("loci inside, left of and right of genes get their own row", {
  gff <- write_test_gff()
  x <- make_test_gl()
  loci <- locNames(x)[1:6]

  res <- gl.find.genes.for.loci(x, gff.file = gff, loci = loci, verbose = 0)
  res <- as.data.frame(res)

  # locus 6 has no position and is dropped; the other five are reported once
  expect_equal(nrow(res), 5)
  expect_setequal(res$locus, loci[1:5])
  res <- res[match(loci[1:5], res$locus), ]

  expect_equal(res$chrom, c("chr1", "chr1", "chr1", "chr1", "chr3"))
  expect_equal(res$pos, c(150L, 300L, 1000L, 1450L, 10L))
  expect_equal(res$gene_name, c("A", "A", "B", "C", NA))
  expect_equal(res$distance_bp, c(0L, 100L, 200L, 50L, NA_integer_))
  expect_equal(res$nearest_side, c("inside", "right", "right", "left", NA))
  # strand-aware position: right of minus-strand B is upstream
  expect_equal(res$gene_strand, c("+", "+", "-", "+", NA))
  expect_equal(res$relative_position,
               c("inside", "downstream", "upstream", "upstream", NA))
})

test_that("gene_product comes from the gene row or its mRNA/CDS children", {
  gff <- write_test_gff()
  x <- make_test_gl()
  res <- as.data.frame(gl.find.genes.for.loci(x, gff.file = gff,
                                              loci = locNames(x)[1:4],
                                              verbose = 0))
  res <- res[match(locNames(x)[1:4], res$locus), ]
  expect_equal(res$gene_product,
               c("kinase A", "kinase A", "phosphatase B", "own product C"))
})

test_that("a locus inside two nested genes is assigned once, to the gene whose midpoint is closer", {
  gff <- withr::local_tempfile(fileext = ".gff")
  writeLines(c(
    "##gff-version 3",
    gff_line("chr1", "src", "gene", 100, 1000, ".", "+", ".", "ID=gene-outer;Name=outer"),
    gff_line("chr1", "src", "gene", 140, 160, ".", "+", ".", "ID=gene-inner;Name=inner")
  ), gff)
  x <- make_test_gl()
  res <- as.data.frame(gl.find.genes.for.loci(x, gff.file = gff,
                                              loci = locNames(x)[1],
                                              verbose = 0))
  expect_equal(nrow(res), 1)
  expect_equal(res$gene_name, "inner")
  expect_equal(res$distance_bp, 0L)
})

test_that("a locus at equal distance from two genes follows the documented tie-break", {
  gff <- withr::local_tempfile(fileext = ".gff")
  # locus at 300: A ends at 250 (50 bp away, midpoint 175, length 150),
  # B starts at 350 (50 bp away, midpoint 600, length 500): A wins on midpoint
  writeLines(c(
    "##gff-version 3",
    gff_line("chr1", "src", "gene", 100, 250, ".", "+", ".", "ID=gene-A;Name=A"),
    gff_line("chr1", "src", "gene", 350, 850, ".", "+", ".", "ID=gene-B;Name=B")
  ), gff)
  x <- make_test_gl()
  res <- as.data.frame(gl.find.genes.for.loci(x, gff.file = gff,
                                              loci = locNames(x)[2],
                                              verbose = 0))
  expect_equal(res$gene_name, "A")
  expect_equal(res$distance_bp, 50L)
  expect_equal(res$nearest_side, "right")
})

test_that("loci defaults to every locus with a position", {
  gff <- write_test_gff()
  x <- make_test_gl()
  res <- gl.find.genes.for.loci(x, gff.file = gff, verbose = 0)
  expect_equal(nrow(res), 5)
  expect_setequal(res$locus, locNames(x)[1:5])
})

test_that("unmapped loci and unknown sequence names are reported at verbose 1", {
  gff <- write_test_gff()
  x <- make_test_gl()
  expect_output(
    gl.find.genes.for.loci(x, gff.file = gff, loci = locNames(x)[5:6],
                           verbose = 1),
    "1 of 2 requested loci have no chromosome or position"
  )
  expect_output(
    gl.find.genes.for.loci(x, gff.file = gff, loci = locNames(x)[5:6],
                           verbose = 1),
    "do not occur in the GFF"
  )
  expect_output(
    gl.find.genes.for.loci(x, gff.file = gff, loci = locNames(x)[1],
                           verbose = 1),
    "Completed"
  )
})

test_that("a gzip-compressed GFF is read from the .gz path or the plain path", {
  gff <- write_test_gff()
  gz <- paste0(gff, ".gz")
  R.utils::gzip(gff, destname = gz, remove = TRUE)
  x <- make_test_gl()
  from_gz <- gl.find.genes.for.loci(x, gff.file = gz, loci = locNames(x)[1:2],
                                    verbose = 0)
  from_plain <- gl.find.genes.for.loci(x, gff.file = gff,
                                       loci = locNames(x)[1:2], verbose = 0)
  expect_equal(from_gz$gene_name, from_plain$gene_name)
  expect_equal(nrow(from_gz), 2)
})

test_that("unknown locus names warn and are skipped", {
  gff <- write_test_gff()
  x <- make_test_gl()
  expect_warning(
    res <- gl.find.genes.for.loci(x, gff.file = gff,
                                  loci = c(locNames(x)[1], "not_a_locus"),
                                  verbose = 0),
    "not found"
  )
  expect_equal(nrow(res), 1)
})
