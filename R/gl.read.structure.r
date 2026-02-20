#' @name gl.read.structure
#'
#' @title Read output files produced by the program STRUCTURE
#'
#' @description
#' Reads and processes STRUCTURE output files, extracting run summaries and Q-matrices.
#' Optionally associates Q-matrices with population information from a genlight object.
#'
#' @param folder.path Path to folder containing STRUCTURE output files [required].
#' @param x Optional genlight object used to attach population labels to individuals [default NULL].
#' @param pattern Optional regex to filter files in folder.path (e.g. ".*_f$" or "out$") [default NULL].
#' @param recursive Logical; search folder recursively [default FALSE].
#' @param rename_files Logical; if TRUE, renames source files on disk based on K and replicate [default FALSE].
#' @param prefix Optional prefix used for renaming/labels. If NULL, uses longest common prefix of filenames [default NULL].
#' @param verbose Verbosity (as in dartR) [default 2 / gl.set.verbosity()].
#'
#' @return A list of class "structure.result". Each element contains:
#' \itemize{
#'   \item summary: named numeric vector (k, est.ln.prob, mean.lnL, var.lnL)
#'   \item q.mat: data.frame (id, pct.miss, orig.pop, Group.1..Group.K)
#'   \item prior.anc: optional list of ancestry matrices (if present)
#'   \item files: file path
#'   \item label: run label
#' }
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
  
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, build = "Jody", verbose = verbose)
  
  # Dependencies (keep explicit, but do not attach)
  req_packages <- c("purrr", "dplyr")
  for (pkg in req_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(error("Package ", pkg, " is needed for this function to work. Please install it."))
    }
  }
  
  if (!dir.exists(folder.path)) {
    stop(error("The folder '", folder.path, "' doesn't exist."))
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
    stop(error("No files found in folder '", folder.path, "' with the requested criteria."))
  }
  
  # ---- helpers ----
  
  longest_common_prefix <- function(strings) {
    strings <- as.character(strings)
    if (length(strings) == 0) return("")
    s0 <- strings[which.min(nchar(strings))]
    if (nchar(s0) == 0) return("")
    for (i in seq_len(nchar(s0))) {
      ch <- substr(s0, i, i)
      if (any(substr(strings, i, i) != ch)) {
        return(if (i == 1) "" else substr(s0, 1, i - 1))
      }
    }
    s0
  }
  
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
        
        # x[-1] contains ancestry strings after '|'
        bits <- strsplit(x[-1], "\\s|:", perl = TRUE)
        bits <- lapply(bits, function(y) y[y != ""])
        # expected: first token is pop, remaining are gens 0..gensback
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
      
      # reconstruct group probs (rowSums), fallback to first group when missing
      prob.mat <- t(vapply(seq_len(nrow(df)), function(i) {
        id_i <- df$id[i]
        anc <- prior.anc[[id_i]]
        p <- rowSums(anc, na.rm = TRUE)
        # if everything is NA, fallback using first group column (STRUCTURE behaviour varies)
        if (all(is.na(p)) || all(p == 0)) {
          # if Group.1 exists, replicate it; otherwise NA
          g1 <- if ("Group.1" %in% names(df)) df$Group.1[i] else NA_real_
          p <- rep(g1, maxpops)
        }
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
  
  file_info <- file_info[order(file_info$k, file_info$f_name), , drop = FALSE]
  file_info <- dplyr::group_by(file_info, .data$k)
  file_info <- dplyr::mutate(file_info, rep = dplyr::row_number())
  file_info <- as.data.frame(file_info)
  
  # label/prefix
  file_names <- basename(file_paths)
  if (is.null(prefix)) {
    prefix <- longest_common_prefix(file_names)
    if (prefix == "") prefix <- "structure"
  }
  
  file_info$label <- paste0(prefix, ".k", file_info$k, ".r", file_info$rep)
  
  # optionally rename files on disk
  if (isTRUE(rename_files)) {
    if (verbose >= 2) cat(report("Renaming STRUCTURE files on disk.\n"))
    new_paths <- file.path(dirname(file_info$f_name), paste0(file_info$label, "_out"))
    ok <- file.rename(file_info$f_name, new_paths)
    if (!all(ok)) warning(warn("Not all files were successfully renamed."))
    file_info$f_name <- ifelse(ok, new_paths, file_info$f_name)
  }
  
  # ---- parse all runs ----
  
  if (verbose >= 2) cat(report("Processing ", nrow(file_info), " STRUCTURE output files.\n"))
  
  run_results <- purrr::pmap(
    list(file_info$f_name, file_info$label),
    function(f, lab) {
      if (verbose >= 3) cat(report("Processing: ", basename(f), "\n"))
      out <- structureRead2(file = f)
      out$files <- f
      out$label <- lab
      out
    }
  )
  
  names(run_results) <- file_info$label
  
  # ---- attach pop data from genlight (optional) ----
  if (!is.null(x)) {
    datatype <- utils.check.datatype(x, verbose = 0)
    if (!datatype %in% c("SNP", "SilicoDArT")) {
      stop(error("The x parameter must be a genlight object containing SNP or SilicoDArT data."))
    }
    
    pop_tbl <- data.frame(id = indNames(x), pop = pop(x), stringsAsFactors = FALSE)
    
    run_results <- lapply(run_results, function(y) {
      if (is.null(y$q.mat) || !("id" %in% names(y$q.mat))) return(y)
      merged <- merge(y$q.mat, pop_tbl, by = "id", all.x = TRUE, sort = FALSE)
      # keep original ordering from y$q.mat
      merged <- merged[match(y$q.mat$id, merged$id), , drop = FALSE]
      merged$orig.pop <- merged$pop
      merged$pop <- NULL
      y$q.mat <- merged
      y
    })
  }
  
  class(run_results) <- c("structure.result", class(run_results))
  
  if (verbose >= 1) cat(report("\nCompleted: ", funname, "\n"))
  run_results
}
