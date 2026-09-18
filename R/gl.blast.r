#' @name gl.blast
#'
#' @title Aligns nucleotides sequences against those present in a target database
#' using blastn
#'
#' @description Basic Local Alignment Search Tool (BLAST; Altschul et al., 1990 &
#'  1997) is a sequence comparison algorithm optimized for speed used to search
#'  sequence databases for optimal local alignments to a query. This function
#'  creates fasta files, creates databases to run BLAST, runs blastn and filters
#'  these results to obtain the best hit per sequence.
#'
#'  This function can be used to run BLAST alignment of short-read (DArTseq
#'  data) and long-read sequences (Illumina, PacBio... etc). You can use
#'  reference genomes from NCBI, genomes from your private collection, contigs,
#'  scaffolds or any other genetic sequence that you would like to use as
#'  reference.
#'
#' @param x Either a genlight object containing a column named
#'  'TrimmedSequence' containing the sequence of the SNPs (the sequence tag)
#'  trimmed of adapters as provided by DArT; or a path to a fasta file with the
#'  query sequences [required].
#' @param ref_genome Path to a reference genome in fasta of fna format or in 
#' a compressed format ie "gz" extension [required].
#' @param task Four different tasks are supported: 1) "megablast", for very
#'  similar sequences (e.g, sequencing errors), 2) "dc-megablast", typically
#'  used for inter-species comparisons, 3) "blastn", the traditional program
#'  used for inter-species comparisons, 4) "blastn-short", optimized for
#'  sequences less than 30 nucleotides [default 'megablast'].
#' @param Percentage_identity Not a very sensitive or reliable measure of
#'  sequence similarity, however it is a reasonable proxy for evolutionary
#'  distance. The evolutionary distance associated with a 10 percent change in
#'  Percentage_identity is much greater at longer distances. Thus, a change from
#'  80 to 70 percent identity might reflect divergence 200 million years earlier
#'  in time, but the change from 30 percent to 20 percent might correspond to a
#'  billion year divergence time change [default 70].
#' @param Percentage_overlap Calculated as alignment length divided by the
#'  query length or subject length (whichever is shortest of the two lengths,
#'  i.e.  length / min(qlen,slen) ) [default 0.8].
#' @param bitscore A rule-of-thumb for inferring homology, a bit score of 50
#'  is almost always significant [default 50].
#' @param number_of_threads Number of threads (CPUs) to use in blastn search
#'  [default 2].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  brief progress messages; 3, progress and results summary; 5, full report
#'  [default 2, unless specified using gl.set.verbosity].
#'
#' @details \strong{Installing BLAST}
#'
#'  You can download the BLAST installs from:
#'  \url{https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/}
#'
#'  The executables blastn and makeblastdb must be findable on the system
#'  PATH. BLAST cannot open a database whose path contains spaces, so the R
#'  temporary directory (tempdir()) must not contain spaces; if it does, set
#'  the environment variable TMPDIR (TMP on Windows) to a path without spaces
#'  and restart R.
#'
#'  \strong{Running BLAST}
#'
#'  Four different tasks are supported: \itemize{ \item "megablast", for very
#'  similar sequences (e.g, sequencing errors) \item "dc-megablast", typically
#'  used for inter-species comparisons \item "blastn", the traditional program
#'  used for inter-species comparisons \item "blastn-short", optimized for
#'  sequences less than 30 nucleotides }
#'
#'  If  you  are  running  a  BLAST alignment of  similar  sequences,  for
#'  example  Turtle  Genome  Vs Turtle Sequences, the recommended parameters
#'  are: task = "megablast", Percentage_identity = 70, Percentage_overlap =  0.8
#'  and bitscore = 50.
#'
#'  If you are running a BLAST alignment of highly dissimilar sequences because
#'  you are probably looking for sex linked  hits in  a distantly  related
#'  species,  and  you  are aligning for example sequences of Chicken Genome Vs
#'  Bassiana, the recommended parameters are: task = "dc-megablast",
#'  Percentage_identity = 50, Percentage_overlap =  0.01 and bitscore = 30.
#'
#'  Be aware that running BLAST might take a long time (i.e. days) depending of
#'  the size of your query, the size of your database and the number of threads
#'  selected for your computer.
#'
#'  \strong{BLAST output}
#'
#'  The BLAST output is formatted as a table using output format 6, with columns
#'  defined in the following order: \itemize{ \item qseqid - Query Seq-id \item
#'  sacc - Subject accession \item stitle - Subject Title \item qseq - Aligned
#'  part of query sequence \item sseq - Aligned part of subject sequence \item
#'  nident - Number of identical matches \item mismatch - Number of mismatches
#'  \item pident - Percentage of identical matches \item length - Alignment
#'  length \item evalue - Expect value \item bitscore - Bit score \item qstart -
#'  Start of alignment in query \item qend - End of alignment in query \item
#'  sstart - Start of alignment in subject \item send - End of alignment in
#'  subject \item gapopen - Number of gap openings \item gaps - Total number of
#'  gaps \item qlen - Query sequence length \item slen - Subject sequence length
#'  \item PercentageOverlap - length / min(qlen,slen) }
#'
#'  Three tables (all aligned sequences, the aligned sequences that pass the
#'  filters, and one hit per sequence) are saved as RDS files in tempdir();
#'  their paths are printed when verbose >= 2. Each RDS holds a list of the
#'  call and the table.
#'
#'  \strong{BLAST filtering}
#'
#'  BLAST output is filtered by ordering the hits of each sequence first by the
#'  highest percentage identity, then the highest percentage overlap and then
#'  the highest bitscore. Only one hit per sequence is kept based on these
#'  selection criteria.
#'
#' @return If the input is a genlight object: returns a genlight object with one
#'  hit per sequence merged to the slot $other$loc.metrics; sequences without a
#'  hit have NA in the BLAST columns. If the input is a fasta file: returns a
#'  dataframe with one hit per sequence. When no hit passes the filters, the
#'  genlight object is returned without BLAST columns, or an empty dataframe
#'  with the BLAST columns is returned for a fasta file.
#'
#' @author Author(s): Berenice Talamantes Becerra & Luis Mijangos. Custodian:
#'  Luis Mijangos -- Post to \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \dontrun{
#' res <- gl.blast(x = testset.gl, ref_genome = "sequence.fasta")
#' # keep only loci with a confident hit
#' res <- gl.filter.locmetric(res, metric = "evalue", upper = 1e-5, lower = 0)
#' }
#'
#' @references
#' \itemize{
#' \item Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D.
#'  J. (1990). Basic local alignment search tool. Journal of molecular biology,
#'  215(3), 403-410.
#' \item Altschul, S. F., Madden, T. L., Schäffer, A. A., Zhang, J., Zhang,
#'  Z., Miller, W., & Lipman, D. J. (1997). Gapped BLAST and PSI-BLAST: a new
#'  generation of protein database search programs. Nucleic acids research,
#'  25(17), 3389-3402.
#' \item Pearson, W. R. (2013). An introduction to sequence similarity
#'  ("homology") searching. Current protocols in bioinformatics, 42(1), 3-1.
#'  }
#'
#' @seealso \code{\link[dartR.base]{gl.filter.locmetric}}
#'
#' @family reference genomes
#'
#' @export
#'

gl.blast <- function(x,
                     ref_genome,
                     task = "megablast",
                     Percentage_identity = 70,
                     Percentage_overlap = 0.8,
                     bitscore = 50,
                     number_of_threads = 2,
                     verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   verbose = verbose)
  
  # CHECK DATATYPE: a genlight object carrying sequence tags, or the path to a
  # fasta file of query sequences
  is_gl <- is(x, "genlight")
  if (is_gl) {
    if (is.null(x@other$loc.metrics$TrimmedSequence)) {
      stop(error(
        "\n\nFatal Error: TrimmedSequence column is required!.\n\n"
      ))
    }
  } else if (is.character(x) && length(x) == 1) {
    if (!file.exists(x)) {
      stop(error("\n\nFatal Error: query fasta file not found:", x, "\n\n"))
    }
  } else {
    stop(error(
      "\n\nFatal Error: x must be a genlight object or the path to a fasta 
      file.\n\n"
    ))
  }
  
  # FUNCTION SPECIFIC ERROR CHECKING
  if (!is.character(ref_genome) || length(ref_genome) != 1 ||
      !file.exists(ref_genome)) {
    stop(error(
      "\n\nFatal Error: reference genome file not found:", ref_genome, "\n\n"
    ))
  }
  task <- match.arg(task, c("megablast", "dc-megablast", "blastn",
                            "blastn-short"))
  if (!is.numeric(Percentage_identity) || Percentage_identity < 0 ||
      Percentage_identity > 100) {
    stop(error("\n\nFatal Error: Percentage_identity must be between 0 and 
               100.\n\n"))
  }
  if (!is.numeric(Percentage_overlap) || Percentage_overlap < 0 ||
      Percentage_overlap > 1) {
    stop(error("\n\nFatal Error: Percentage_overlap must be between 0 and 
               1.\n\n"))
  }
  if (!is.numeric(bitscore) || bitscore < 0) {
    stop(error("\n\nFatal Error: bitscore must be a non-negative number.\n\n"))
  }
  if (!is.numeric(number_of_threads) || number_of_threads < 1) {
    stop(error("\n\nFatal Error: number_of_threads must be at least 1.\n\n"))
  }
  # makeblastdb and blastn split a database path on whitespace, so the
  # database built in tempdir() cannot live under a path with spaces
  if (grepl("\\s", tempdir())) {
    stop(error(
      "\n\nFatal Error: the R temporary directory", tempdir(), "contains a 
      space and BLAST cannot open a database there. Set the environment 
      variable TMPDIR (TMP on Windows) to a path without spaces and restart 
      R.\n\n"
    ))
  }
  
  # Find the executables. Sys.which on unix; 'where' on windows.
  find_exec <- function(name) {
    if (grepl("unix", .Platform$OS.type, ignore.case = TRUE)) {
      path <- unname(Sys.which(name))
    } else {
      path <- tryCatch(
        system(sprintf("where %s", name), intern = TRUE)[1],
        warning = function(w) "",
        error = function(e) ""
      )
      if (is.na(path)) path <- ""
    }
    if (!nzchar(path)) {
      stop(error(
        "\n\nFatal Error: Executable for", name, "not found! Please make 
        sure that the software is correctly installed and on the PATH.\n\n"
      ))
    }
    path
  }
  path_makeblastdb <- find_exec("makeblastdb")
  path_blastn <- find_exec("blastn")
  
  # Run an external tool and stop on failure. Arguments are shell-quoted by
  # the caller, so paths with spaces are safe. The tool's stdout is shown at
  # verbose >= 3; its stderr always reaches the console so a failure is
  # explained before the stop.
  run_tool <- function(exec, args, what, stdin = "") {
    status <- system2(exec, args,
                      stdout = if (verbose >= 3) "" else FALSE,
                      stderr = "",
                      stdin = stdin)
    if (status != 0) {
      stop(error(
        "\n\nFatal Error:", what, "exited with status", status, ". See the 
        messages above.\n\n"
      ))
    }
  }
  
  # DO THE JOB
  
  # Intermediate files have fixed names in tempdir(). Remove anything left by
  # a previous run so that a failing step can never hand back stale results.
  fasta_input <- file.path(tempdir(), "fasta.input")
  db_prefix <- file.path(tempdir(), "db_blast")
  blast_out <- file.path(tempdir(), "output_blast.txt")
  unlink(c(fasta_input, blast_out,
           list.files(tempdir(), pattern = "^db_blast\\.",
                      full.names = TRUE)))
  
  # getting the query fasta files
  if (is_gl) {
    fasta.input <-
      c(rbind(
        paste0(">", seq_len(nLoc(x))),
        as.character(x@other$loc.metrics$TrimmedSequence)
      ))
    writeLines(fasta.input, fasta_input)
    n_query <- nLoc(x)
  } else {
    if (!file.copy(from = x, to = fasta_input, overwrite = TRUE)) {
      stop(error(
        "\n\nFatal Error: could not copy the query fasta file to 
        tempdir().\n\n"
      ))
    }
    n_query <- length(grep("^>", readLines(fasta_input)))
  }
  
  # if ref_genome is gzipped, decompress to a temp file
  if (grepl("\\.gz$", ref_genome, ignore.case = TRUE)) {
    # name of the uncompressed FASTA in tempdir()
    out_fa <- file.path(tempdir(),
                        basename(sub("\\.gz$", "", ref_genome)))
    R.utils::gunzip(ref_genome,
                    destname = out_fa,
                    overwrite = TRUE,
                    remove = FALSE)   # keep the original .gz
    ref_genome <- out_fa
  } 
  
  # creating BLAST databases. The genome is fed on stdin because makeblastdb
  # splits the value of -in on whitespace (it accepts a list of files), so a
  # genome path with spaces cannot be passed as an argument.
  if (verbose >= 2) {
    cat(report("  Building the BLAST database\n"))
  }
  run_tool(path_makeblastdb,
           c("-dbtype", "nucl",
             "-title", "ref_genome",
             "-out", shQuote(db_prefix)),
           what = "makeblastdb",
           stdin = ref_genome)
  
  blast_cols <- c("qseqid", "sacc", "stitle", "qseq", "sseq", "nident",
                  "mismatch", "pident", "length", "evalue", "bitscore",
                  "qstart", "qend", "sstart", "send", "gapopen", "gaps",
                  "qlen", "slen")
  
  # BLASTing
  if (verbose >= 2) {
    cat(report("  Starting BLASTing\n"))
  }
  run_tool(path_blastn,
           c("-task", task,
             "-db", shQuote(db_prefix),
             "-query", shQuote(fasta_input),
             "-out", shQuote(blast_out),
             "-perc_identity", Percentage_identity,
             "-num_threads", number_of_threads,
             "-outfmt", shQuote(paste("6", paste(blast_cols,
                                                 collapse = " ")))),
           what = "blastn")
  
  # reading file for filtering; an empty file means nothing aligned
  if (file.exists(blast_out) && file.info(blast_out)$size > 0) {
    blast_res_unfiltered <-
      read.table(
        file = blast_out,
        header = FALSE,
        sep = "\t",
        quote = "",
        dec = ".",
        fill = TRUE,
        comment.char = "",
        stringsAsFactors = FALSE
      )
  } else {
    blast_res_unfiltered <-
      as.data.frame(matrix(nrow = 0, ncol = length(blast_cols)))
  }
  colnames(blast_res_unfiltered) <- blast_cols
  
  if (verbose >= 2) {
    cat(report("  Starting filtering\n"))
  }
  
  # calculate percentage overlap ratio of the alignment length divided by the
  # query length or subject length (whichever is shortest of
  # the two lengths)
  blast_res_unfiltered$PercentageOverlap <-
    blast_res_unfiltered$length /
    pmin(blast_res_unfiltered$qlen, blast_res_unfiltered$slen)
  # filtering first by percentage overlap and bitscore
  blast_res_filtered <-
    blast_res_unfiltered[which(
      blast_res_unfiltered$PercentageOverlap >= Percentage_overlap &
        blast_res_unfiltered$bitscore >=
        bitscore
    ),]
  # splitting hits by sequence
  all_hits <-
    split(x = blast_res_filtered, f = blast_res_filtered$qseqid)
  # ordering by first considering the highest percentage identity, then the 
  # highest percentage overlap, then the highest bitscore. Only
  # one hit per sequence is kept based on these selection criteria.
  one_hit_temp <- lapply(all_hits, function(x) {
    x[order(x$pident,
            x$PercentageOverlap,
            x$bitscore,
            decreasing = TRUE),][1,]
  })
  
  if (length(one_hit_temp) > 0) {
    one_hit <- plyr::rbind.fill(one_hit_temp)
  } else {
    one_hit <- blast_res_filtered[0, ]
  }
  rownames(one_hit) <- NULL
  n_hit <- nrow(one_hit)
  
  # merging one hit per sequence with genlight object
  if (is_gl && n_hit > 0) {
    one_hit_temp <- x@other$loc.metrics
    one_hit_temp$qseqid <- seq_len(nLoc(x))
    x@other$loc.metrics <-
      merge(one_hit_temp,
            one_hit,
            by = "qseqid",
            all = TRUE)
  }
  
  # a result-affecting outcome: printed unless silent
  if (verbose >= 1) {
    if (n_hit > 0) {
      cat(report("  ", n_hit, "of", n_query,
                 "sequences aligned after filtering\n"))
    } else {
      cat(warn("  0 of", n_query, "sequences aligned after filtering;",
               "no BLAST columns were added\n"))
    }
  }
  
  match_call <-
    paste0(names(match.call()),
           "_",
           as.character(match.call()),
           collapse = "_")
  
  # creating file names
  temp_blast_unfiltered <- tempfile(pattern = "Blast_unfiltered_")
  temp_blast_filtered <- tempfile(pattern = "Blast_filtered_")
  temp_one_hit <- tempfile(pattern = "Blast_one_hit_")
  
  # saving to tempdir
  saveRDS(list(match_call, blast_res_unfiltered), 
          file = temp_blast_unfiltered)
  saveRDS(list(match_call, blast_res_filtered), file = temp_blast_filtered)
  saveRDS(list(match_call, one_hit), file = temp_one_hit)
  
  if (verbose >= 2) {
    cat(report("  BLAST tables saved as RDS files (list of call, table):\n"))
    cat(report("    all hits:     ", temp_blast_unfiltered, "\n"))
    cat(report("    filtered hits:", temp_blast_filtered, "\n"))
    cat(report("    one hit per sequence:", temp_one_hit, "\n"))
  }
  
  # ADD TO HISTORY
  if (is_gl) {
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()
  }
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n\n"))
  }
  
  if (is_gl) {
    return(x)
  } else {
    return(one_hit)
  }
  
}
