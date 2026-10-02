#' @name utils.structure.read
#' @title Reads one STRUCTURE output file
#' @description
#' Parses a STRUCTURE output file into its run summary, q-matrix and, for
#' individuals with a population prior (USEPOPINFO), the ancestry matrices.
#' Shared by \code{\link{gl.read.structure}} and
#' \code{\link{utils.structure.run}} so that both read the file the same way.
#' @param file Path of the STRUCTURE output file [required].
#' @param pops Optional vector of population names, indexed by the population
#' numbers in the file, used to label orig.pop [default NULL].
#' @return A list with summary (named numeric vector: k, est.ln.prob,
#' mean.lnL, var.lnL), q.mat (data.frame: id, pct.miss, orig.pop,
#' Group.1..Group.K) and prior.anc (list of ancestry matrices named by id, or
#' NULL).
#' @author Author(s): Luis Mijangos & Bernd Gruber, adapted from strataG
#' (Eric Archer). Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @keywords internal
#' @noRd

utils.structure.read <- function(file, pops = NULL) {
  clean_q_line <- function(s) {
    s <- gsub("\\*+", "", s)
    s <- gsub("[()]", "", s)
    s <- sub("\\|\\s*$", "", s)
    s
  }
  
  .structureParseQmat <- function(q.mat.txt, pops = NULL, popdata = 1L) {
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
    no.prior <- .structureParseQmat(no.prior.q.txt, pops = pops, popdata = popdata)
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
    df <- .structureParseQmat(prior.q.txt, pops = pops, popdata = popdata)
    
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
