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
      out <- utils.structure.read(file = f)
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
