#' @name gl.read.structure
#'
#' @title Read output files produced by the program STRUCTURE
#'
#' @family population structure
#'
#' @description
#' Reads STRUCTURE output files from a folder and returns them as a
#' structure run object, the same object that \code{\link{gl.run.structure}}
#' returns, so it can be passed to \code{\link{gl.evanno}},
#' \code{\link{gl.plot.structure}} and \code{\link{gl.map.structure}}.
#' Optionally attaches population and individual names from the genlight
#' object that was analysed.
#'
#' @param folder.path Path to folder containing STRUCTURE output files
#' [required].
#' @param x The genlight object that was analysed, used to attach population
#' labels (and individual names) to the q-matrices [default NULL].
#' @param pattern Optional regular expression to select files in folder.path
#' (e.g. "_f$") [default NULL].
#' @param recursive Logical; search folder recursively [default FALSE].
#' @param rename_files Logical; if TRUE, renames the output files on disk to
#' <run name>_out. Stops, before renaming anything, if a file of that name
#' already exists [default FALSE].
#' @param prefix Optional prefix for the run names, e.g. "myrun" gives
#' "myrun.k2.r1" [default NULL, run names "k2.r1"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  brief progress messages; 3, progress and results summary; 5, full report
#'  [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' Only STRUCTURE output files are read: files containing both
#' "Estimated Ln Prob of Data" and "Estimated Allele Frequencies". Other files in the folder (data, params, log
#' files) are skipped, with a note at verbose >= 2. K is taken from each
#' file. Within each K, replicates are numbered in the order of the numbers
#' in their file names (so rep2 comes before rep10) and runs are named
#' k<K>.r<replicate>, as in \code{\link{gl.run.structure}}.
#'
#' Without x, the orig.pop column holds the population number used in the
#' STRUCTURE data file, and id the label from the file. With x, individuals
#' are matched to indNames(x) by name or, when the labels are the numbers
#' 1 to nInd(x) (as in files written by \code{\link{gl.run.structure}}), by
#' position; id then holds the individual names and orig.pop the population
#' names. Labels that match neither way stop the function with an error
#' naming them (STRUCTURE truncates labels longer than 11 characters).
#'
#' For runs with USEPOPINFO, individuals with a population prior are read
#' like other individuals: their q-matrix row holds the probability of
#' belonging to their given population and, for every other population, the
#' summed probability of ancestry from it over the GENSBACK generations,
#' which is also returned per generation in prior.anc.
#'
#' @return A list of class "structure.result", one element per run, named
#' k<K>.r<replicate>. Each element contains:
#' \itemize{
#'   \item summary: named numeric vector (k, est.ln.prob, mean.lnL, var.lnL)
#'   \item q.mat: data.frame (id, pct.miss, orig.pop, Group.1..Group.K)
#'   \item prior.anc: list of ancestry matrices for individuals with a
#'   population prior, or NULL
#'   \item files: path of the output file
#'   \item label: run name
#' }
#'
#' @author Author(s): Luis Mijangos & Bernd Gruber. Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \dontrun{
#' # read the output files kept by gl.run.structure
#' sr <- gl.run.structure(bc, k.range = 2:5, num.k.rep = 3,
#'                        exec = "./structure", delete.files = FALSE,
#'                        plot.dir = "structure_runs")
#' sr2 <- gl.read.structure(list.dirs("structure_runs",
#'                                    recursive = FALSE)[1], x = bc)
#' gl.evanno(sr2)
#' qmat <- gl.plot.structure(sr2, K = 3)
#' }
#'
#' @seealso \code{\link{gl.run.structure}}, \code{\link{gl.evanno}},
#' \code{\link{gl.plot.structure}}
#'
#' @export
gl.read.structure <- function(folder.path,
                              x = NULL,
                              pattern = NULL,
                              recursive = FALSE,
                              rename_files = FALSE,
                              prefix = NULL,
                              verbose = NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  if (!is.null(x)) {
    datatype <- utils.check.datatype(x, verbose = 0)
    if (!datatype %in% c("SNP", "SilicoDArT")) {
      stop(error(
        "The x parameter must be a genlight object containing SNP or",
        "SilicoDArT data.\n"
      ))
    }
  }
  
  if (!dir.exists(folder.path)) {
    stop(error(paste0("The folder '", folder.path, "' doesn't exist.\n")))
  }
  
  file_paths <- list.files(
    folder.path,
    full.names = TRUE,
    recursive = recursive
  )
  if (!is.null(pattern)) {
    file_paths <- file_paths[grepl(pattern, basename(file_paths))]
  }
  if (length(file_paths) == 0) {
    stop(error(paste0(
      "No files found in folder '", folder.path,
      "' with the requested criteria.\n"
    )))
  }
  
  # keep STRUCTURE output files only (data, params and log files are
  # skipped); the run log repeats the likelihood but has no allele
  # frequency section
  is_output <- vapply(file_paths, function(f) {
    lines <- tryCatch(readLines(f, warn = FALSE),
                      error = function(e) character(0))
    any(grepl("Estimated Ln Prob of Data", lines, fixed = TRUE)) &&
      any(grepl("Estimated Allele Frequencies", lines, fixed = TRUE))
  }, logical(1))
  if (any(!is_output) && verbose >= 2) {
    cat(report(
      "  Skipping", sum(!is_output), "file(s) that are not STRUCTURE output:",
      paste(basename(file_paths[!is_output]), collapse = ", "), "\n"
    ))
  }
  file_paths <- file_paths[is_output]
  if (length(file_paths) == 0) {
    stop(error(paste0(
      "No STRUCTURE output files found in '",
      folder.path, "'.\n"
    )))
  }
  
  # ---- helpers ----
  
  # Detect K: prefer MAXPOPS=, fallback to "populations assumed"
  detect_k <- function(file) {
    lines <- readLines(file, warn = FALSE, n = 300)
    
    # MAXPOPS= is typically present in STRUCTURE outputs
    mx <- grep("MAXPOPS\\s*=", lines, value = TRUE)
    if (length(mx) > 0) {
      k <- suppressWarnings(as.integer(sub(".*MAXPOPS\\s*=\\s*([0-9]+).*", "\\1", mx[1])))
      if (!is.na(k)) return(k)
    }
    
    # fallback
    pa <- grep("populations assumed", lines, ignore.case = TRUE, value = TRUE)
    if (length(pa) > 0) {
      k <- suppressWarnings(as.integer(gsub("[^0-9]", "", pa[1])))
      if (!is.na(k)) return(k)
    }
    
    NA_integer_
  }
  
  clean_q_line <- function(s) {
    s <- gsub("\\*+", "", s)
    s <- gsub("[()]", "", s)
    s <- sub("\\|\\s*$", "", s)
    s
  }
  
  .structureParseQmat2 <- function(q.mat.txt, pops = NULL, popdata = 1L) {
    q.mat.txt <- clean_q_line(q.mat.txt)
    
    # split by whitespace; keep stable columns
    cols_fixed <- c("row", "id", "pct.miss", "orig.pop")
    
    rows <- strsplit(q.mat.txt, "\\s+") |>
      purrr::map(function(q) {
        q <- q[q != "" & q != ":"]
        if (length(q) == 0) return(NULL)
        
        # If POPDATA==0 then orig.pop is absent; insert a dummy "1" after pct.miss.
        if (popdata == 0L) {
          if (length(q) < 4) return(NULL)
          q <- c(q[1:3], "1", q[4:length(q)])
        }
        
        df <- as.data.frame(t(q), stringsAsFactors = FALSE)
        n_groups <- max(0, ncol(df) - 4)
        colnames(df) <- c(cols_fixed, paste0("Group.", seq_len(n_groups)))
        df
      }) |>
      purrr::compact()
    
    if (length(rows) == 0) return(NULL)
    
    out <- dplyr::bind_rows(rows)
    
    # numeric coercion (suppress warnings from any stray text)
    num_cols <- c("row", "pct.miss", "orig.pop", grep("^Group\\.", names(out), value = TRUE))
    out[num_cols] <- lapply(out[num_cols], function(z) suppressWarnings(as.numeric(z)))
    
    if (!is.null(pops)) {
      out$orig.pop <- pops[out$orig.pop]
    }
    
    out
  }
  
  structureRead2 <- function(file, pops = NULL) {
    if (!file.exists(file)) stop(error("The file '", file, "' can't be found."))
    
    # read token stream for summary stats (fast enough, keeps your original logic)
    tokens <- scan(file, what = "character", quiet = TRUE)
    
    get_after <- function(pattern, offset, which_hit = 1L) {
      loc <- grep(pattern, tokens, fixed = FALSE)
      if (length(loc) < which_hit) return(NA_real_)
      suppressWarnings(as.numeric(tokens[loc[which_hit] + offset]))
    }
    
    est.ln.prob <- get_after("Estimated", 6, 1)
    mean.lnL    <- get_after("likelihood", 2, 1)
    var.lnL     <- get_after("likelihood", 2, 2)
    
    # Robust parameter extraction from tokens
    extract_param_int <- function(param) {
      loc <- grep(param, tokens, fixed = TRUE)
      if (length(loc) == 0) return(NA_integer_)
      v <- tokens[loc[1]]
      v <- sub(paste0(".*", param, "="), "", v)
      v <- sub(",", "", v, fixed = TRUE)
      suppressWarnings(as.integer(v))
    }
    
    maxpops  <- extract_param_int("MAXPOPS")
    gensback <- extract_param_int("GENSBACK")
    popdata  <- extract_param_int("POPDATA")
    if (is.na(popdata)) popdata <- 1L
    
    smry <- c(
      k = maxpops,
      est.ln.prob = est.ln.prob,
      mean.lnL = mean.lnL,
      var.lnL = var.lnL
    )
    
    # line-based parsing for Q-matrix
    lines <- scan(file, what = "character", sep = "\n", quiet = TRUE)
    
    first <- grep("\\(%Miss\\)", lines) + 1
    last  <- grep("Estimated Allele", lines) - 1
    if (length(first) == 0 || length(last) == 0 || first[1] > last[1]) {
      stop(error("Could not locate Q-matrix block in file: ", basename(file)))
    }
    
    tbl.txt <- lines[first[1]:last[1]]
    tbl.txt <- clean_q_line(tbl.txt)
    
    prior.lines <- grep("\\|", tbl.txt)
    
    no.prior <- NULL
    if (length(prior.lines) < length(tbl.txt)) {
      no.prior.q.txt <- if (length(prior.lines) == 0) tbl.txt else tbl.txt[-prior.lines]
      no.prior <- .structureParseQmat2(no.prior.q.txt, pops = pops, popdata = popdata)
    }
    
    # K=1 shortcut
    if (!is.na(maxpops) && maxpops == 1L) {
      if (!is.null(no.prior)) no.prior$row <- NULL
      return(list(summary = smry, q.mat = no.prior, prior.anc = NULL))
    }
    
    has.prior.df <- NULL
    prior.anc <- NULL
    
    if (length(prior.lines) > 0) {
      prior.parts <- strsplit(tbl.txt[prior.lines], "\\|", fixed = FALSE)
      
      prior.q.txt <- vapply(prior.parts, function(x) x[1], character(1))
      df <- .structureParseQmat2(prior.q.txt, pops = pops, popdata = popdata)
      
      # build prior ancestry matrices per individual (keyed by id)
      prior.anc <- purrr::imap(prior.parts, function(x, idx) {
        anc.mat <- matrix(NA_real_, nrow = maxpops, ncol = gensback + 1)
        rownames(anc.mat) <- paste0("Pop.", seq_len(nrow(anc.mat)))
        colnames(anc.mat) <- paste0("Gen.", 0:gensback)
        
        # x[-1] holds one "Pop j: g0 g1 ... gGENSBACK" block per other
        # population; drop the word "Pop" so the first token is j
        bits <- strsplit(x[-1], "\\s|:", perl = TRUE)
        bits <- lapply(bits, function(y) {
          y <- y[y != ""]
          if (length(y) > 0 && y[1] == "Pop") y <- y[-1]
          y
        })
        for (b in bits) {
          if (length(b) < 2) next
          pop_i <- suppressWarnings(as.integer(b[1]))
          if (is.na(pop_i) || pop_i < 1 || pop_i > maxpops) next
          vals <- suppressWarnings(as.numeric(b[-1]))
          if (length(vals) == ncol(anc.mat)) anc.mat[pop_i, ] <- vals
        }
        anc.mat
      })
      
      # name prior.anc by id if available
      if (!is.null(df) && "id" %in% names(df) && length(prior.anc) == nrow(df)) {
        names(prior.anc) <- df$id
      }
      
      # other populations: ancestry summed over generations; the given
      # population (no ancestry block, so NA) takes the value before '|'
      prob.mat <- t(vapply(seq_len(nrow(df)), function(i) {
        p <- rowSums(prior.anc[[i]])
        p[is.na(p)] <- df$Group.1[i]
        p
      }, numeric(maxpops)))
      
      colnames(prob.mat) <- paste0("Group.", seq_len(ncol(prob.mat)))
      df$Group.1 <- NULL
      has.prior.df <- cbind(df, prob.mat, stringsAsFactors = FALSE)
    }
    
    q.mat <- dplyr::bind_rows(no.prior, has.prior.df)
    if (is.null(q.mat) || nrow(q.mat) == 0) {
      stop(error("Parsed Q-matrix is empty for file: ", basename(file)))
    }
    
    # order by original row and drop helper
    if ("row" %in% names(q.mat)) {
      q.mat <- q.mat[order(q.mat$row), , drop = FALSE]
      q.mat$row <- NULL
    }
    rownames(q.mat) <- NULL
    
    # normalize Q rows to sum to 1 (Groups only)
    group_cols <- grep("^Group\\.", names(q.mat), value = TRUE)
    if (length(group_cols) > 0) {
      q <- as.matrix(q.mat[, group_cols, drop = FALSE])
      rs <- rowSums(q, na.rm = TRUE)
      rs[rs == 0] <- NA_real_
      q <- q / rs
      q.mat[, group_cols] <- q
    }
    
    list(summary = smry, q.mat = q.mat, prior.anc = prior.anc)
  }
  
  # ---- file inventory ----
  
  k_vals <- vapply(file_paths, detect_k, integer(1))
  file_info <- data.frame(
    f_name = file_paths,
    k = k_vals,
    stringsAsFactors = FALSE
  )
  
  if (all(is.na(file_info$k))) {
    stop(error("Could not detect K for any file in: ", folder.path))
  }
  
  # replicates in the order of the last number in the file name (rep2
  # before rep10), then by name
  last_num <- vapply(basename(file_info$f_name), function(b) {
    n <- regmatches(b, gregexpr("[0-9]+", b))[[1]]
    if (length(n) == 0) NA_real_ else as.numeric(n[length(n)])
  }, numeric(1))
  file_info <- file_info[order(file_info$k, last_num, file_info$f_name), ,
                         drop = FALSE]
  file_info <- dplyr::group_by(file_info, .data$k)
  file_info <- dplyr::mutate(file_info, rep = dplyr::row_number())
  file_info <- as.data.frame(file_info)
  
  file_info$label <- paste0("k", file_info$k, ".r", file_info$rep)
  if (!is.null(prefix)) {
    file_info$label <- paste0(prefix, ".", file_info$label)
  }
  
  # optionally rename files on disk
  if (isTRUE(rename_files)) {
    if (verbose >= 2) cat(report("  Renaming STRUCTURE files on disk.\n"))
    new_paths <- file.path(dirname(file_info$f_name), paste0(file_info$label, "_out"))
    clash <- new_paths[file.exists(new_paths) & new_paths != file_info$f_name]
    if (length(clash) > 0) {
      stop(error(
        "No files were renamed: these target files already exist:",
        paste(clash, collapse = ", "), "\n"
      ))
    }
    ok <- file.rename(file_info$f_name, new_paths)
    if (!all(ok)) warning(warn("Not all files were successfully renamed."))
    file_info$f_name <- ifelse(ok, new_paths, file_info$f_name)
  }
  
  # ---- parse all runs ----
  
  if (verbose >= 2) {
    cat(report("  Processing", nrow(file_info), "STRUCTURE output files.\n"))
  }
  
  run_results <- purrr::pmap(
    list(file_info$f_name, file_info$label),
    function(f, lab) {
      if (verbose >= 3) cat(report("  Processing:", basename(f), "\n"))
      out <- structureRead2(file = f)
      out$files <- f
      out$label <- lab
      out
    }
  )
  
  names(run_results) <- file_info$label
  
  # ---- attach pop data from genlight (optional) ----
  if (!is.null(x)) {
    ind_names <- indNames(x)
    ind_pops <- as.character(pop(x))
    
    run_results <- lapply(run_results, function(y) {
      if (is.null(y$q.mat) || !("id" %in% names(y$q.mat))) return(y)
      ids <- as.character(y$q.mat$id)
      if (all(ids %in% ind_names)) {
        idx <- match(ids, ind_names)
      } else if (all(grepl("^[0-9]+$", ids)) &&
                 setequal(as.integer(ids), seq_along(ind_names))) {
        # labels are positions in indNames(x), as written by
        # gl.run.structure: restore the names
        idx <- as.integer(ids)
        y$q.mat$id <- ind_names[idx]
        if (!is.null(y$prior.anc)) {
          names(y$prior.anc) <- ind_names[as.integer(names(y$prior.anc))]
        }
      } else {
        unmatched <- setdiff(ids, ind_names)
        stop(error(
          "Individuals in", basename(y$files),
          "match neither indNames(x) nor positions 1 to nInd(x):",
          paste(utils::head(unmatched, 10), collapse = ", "),
          if (length(unmatched) > 10) "..." else "",
          "(STRUCTURE truncates labels longer than 11 characters).\n"
        ))
      }
      y$q.mat$orig.pop <- ind_pops[idx]
      y
    })
  }
  
  class(run_results) <- c("structure.result", class(run_results))
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  run_results
}
