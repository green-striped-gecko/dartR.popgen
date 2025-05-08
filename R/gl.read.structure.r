#' @name gl.read.structure
#' 
#' @title Read output files produced by the program STRUCTURE
#' 
#' @description
#' This function reads and processes files produced by the population genetics
#' program STRUCTURE (Pritchard et al., 2000). It renames the files based on the
#' K value found in each file and consolidates the data into a standardized format.
#' 
#' @param folder.path Path to the folder containing STRUCTURE output files [required].
#' @param x Name of the genlight object to associate population information with
#'  [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log;
#'  3, progress and results summary; 5, full report [default 2 or as specified using 
#'  gl.set.verbosity].
#' 
#' @details
#' This function processes STRUCTURE output files by:
#' \itemize{
#'   \item Finding the K value (number of populations) in each file
#'   \item Renaming files based on their K value and replicate number
#'   \item Extracting Q-matrices and other relevant information
#'   \item Optionally associating the data with population information from a genlight object
#' }
#' 
#' STRUCTURE is a widely used software for inferring population structure from 
#' multilocus genotype data. This function helps to organize and consolidate the
#' output files for further analysis in R.
#' 
#' @return A list of class "structure.result" containing the processed results from
#'  each STRUCTURE output file, with names based on the file labels.
#' 
#' @author Luis Mijangos & Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' 
#' @examples
#' \dontrun{
#' # Assuming structure output files are in the "structure_results" folder
#' str_results <- gl.read.structure("structure_results")
#' 
#' # With population information from a genlight object
#' str_results <- gl.read.structure("structure_results", testset.gl)
#' }
#' 
#' @references
#' Pritchard, J.K., Stephens, M., and Donnelly, P. (2000). Inference of population 
#' structure using multilocus genotype data. Genetics 155, 945-959.
#' 
#' @seealso \code{\link{gl.plot.structure}}
#' 
#' @family structure functions
#' 
#' @importFrom purrr map
#' 
#' @export
#'

gl.read.structure <- function(folder.path,
                              x = NULL,
                              verbose = NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, build = "Jody", v = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  # Check if required packages are installed
  req_packages <- c("purrr", "dplyr")
  for (pkg in req_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(error("Package", pkg, " needed for this function to work. Please install it."))
    }
  }
  
  # Check if folder exists
  if (!dir.exists(folder.path)) {
    stop(error("Error: the folder '", folder.path, "' doesn't exist."))
  }
  
  # List all files in the directory
  file_paths <- list.files(folder.path, full.names = TRUE)
  file_names <- list.files(folder.path, full.names = FALSE)
  
  # Check if any files were found
  if (length(file_paths) == 0) {
    stop(error("Error: no files found in folder '", folder.path, "'"))
  }
  
  # DO THE JOB
  
  # Function to find the longest common prefix among strings
  longest_common_prefix <- function(x) {
    x <- as.character(x)  # Convert input to character vector
    if (length(x) == 0)
      return("")  # Return empty string if input is empty
    
    # Select the shortest string to avoid out-of-bounds errors
    s0 <- x[which.min(nchar(x))]
    
    # Iterate through each character position
    for (i in seq_len(nchar(s0))) {
      this_chr <- substr(s0, i, i)
      if (any(substr(x, i, i) != this_chr)) {
        # If mismatch found at position i, return prefix up to i-1
        return(if (i == 1) "" else substr(s0, 1, i - 1))
      }
    }
    
    # If no mismatch, the shortest string is the common prefix
    return(s0)
  }
  
  # Parse the Q-matrix text from STRUCTURE output
  .structureParseQmat2 <- function(q.mat.txt, pops, popdata) {
    # Clean up text by removing special characters
    q.mat.txt <- sub("[*]+", "", q.mat.txt)
    q.mat.txt <- sub("[(]", "", q.mat.txt)
    q.mat.txt <- sub("[)]", "", q.mat.txt)
    q.mat.txt <- sub("[|][ ]+$", "", q.mat.txt)
    
    # Define column names for the first four columns
    cols1to4 <- c("row", "id", "pct.miss", "orig.pop")
    
    # Split text by spaces and process each line
    result <- strsplit(q.mat.txt, " ") %>%
      purrr::map(function(q) {
        # Remove empty elements
        q <- q[!q %in% c("", " ", ":")]
        
        # Handle case when population data is not present
        if (popdata == 0) {
          q <- c(q[1:3], "1", q[4:length(q)])
        }
        
        # Convert to data frame
        q <- q %>%
          as.character() %>%
          rbind() %>%
          as.data.frame(stringsAsFactors = FALSE)
        
        # Assign column names (first 4 fixed columns + group columns)
        stats::setNames(q, c(cols1to4, paste("Group", 1:(ncol(q) - 4), sep = ".")))
      }) %>%
      dplyr::bind_rows() %>%
      # Convert numeric columns to actual numbers
      dplyr::mutate_at(dplyr::vars("row", "pct.miss", "orig.pop", dplyr::starts_with("Group.")),
                       as.numeric) 
    
    # Map population indices to names if provided
    if (!is.null(pops)) {
      result <- result %>% dplyr::mutate(orig.pop = pops[.data$orig.pop])
    }
    
    return(result)
  }
  
  # Read and parse STRUCTURE output file
  structureRead2 <- function(file, pops = NULL) {
    # Check if file exists
    if (!file.exists(file)) {
      stop(error("Error: the file '", file, "' can't be found."))
    }
    
    # Read the entire file as a character vector
    result <- scan(file, "character", quiet = TRUE)
    
    # Extract estimated log probability
    loc <- grep("Estimated", result, ignore.case = FALSE, value = FALSE)
    est.ln.prob <- as.numeric(result[loc[1] + 6])
    
    # Extract likelihood statistics
    loc <- grep("likelihood", result, ignore.case = FALSE, value = FALSE)
    mean.lnL <- as.numeric(result[loc[1] + 2])
    var.lnL <- as.numeric(result[loc[2] + 2])
    
    # Extract MAXPOPS parameter (number of populations, K)
    loc <- grep("MAXPOPS", result, value = FALSE)
    maxpops <- result[loc]
    maxpops <- sub("MAXPOPS=", "", maxpops)
    maxpops <- as.integer(sub(",", "", maxpops))
    
    # Extract GENSBACK parameter (number of generations)
    loc <- grep("GENSBACK", result, value = FALSE)
    gensback <- result[loc]
    gensback <- sub("GENSBACK=", "", gensback)
    gensback <- as.integer(sub(",", "", gensback))
    
    # Extract POPDATA parameter (population data flag)
    loc <- grep("POPDATA", result, value = FALSE)
    popdata <- result[loc]
    popdata <- sub("POPDATA=", "", popdata)
    popdata <- as.integer(sub(",", "", popdata))
    
    # Create summary statistics
    smry <- c(
      k = maxpops,
      est.ln.prob = est.ln.prob,
      mean.lnL = mean.lnL,
      var.lnL = var.lnL
    )
    
    # Read the file again, but by lines
    result <- scan(file, "character", sep = "\n", quiet = TRUE)
    
    # Extract the Q-matrix section
    first <- grep("(%Miss)", result, value = FALSE) + 1
    last <- grep("Estimated Allele", result, value = FALSE) - 1
    tbl.txt <- result[first:last]
    
    # Clean up text
    tbl.txt <- sub("[*]+", "", tbl.txt)
    tbl.txt <- sub("[(]", "", tbl.txt)
    tbl.txt <- sub("[)]", "", tbl.txt)
    tbl.txt <- sub("[|][ ]+$", "", tbl.txt)
    
    # Identify lines with prior information
    prior.lines <- grep("[|]", tbl.txt)
    
    # Process lines without prior information
    no.prior <- if (length(prior.lines) < length(tbl.txt)) {
      no.prior.q.txt <- if (length(prior.lines) == 0) {
        tbl.txt
      } else {
        tbl.txt[-prior.lines]
      }
      .structureParseQmat2(q.mat.txt = no.prior.q.txt, pops, popdata = popdata)
    } else {
      NULL
    }
    
    # For K=1, return simplified results
    if (maxpops == 1) {
      no.prior$row <- NULL
      return(list(
        summary = smry,
        q.mat = no.prior,
        prior.anc = NULL
      ))
    }
    
    # Process lines with prior information
    has.prior <- if (length(prior.lines) > 0) {
      prior.txt <- strsplit(tbl.txt[prior.lines], "[|]")
      prior.q.txt <- unlist(lapply(prior.txt, function(x) x[1]))
      df <- .structureParseQmat2(prior.q.txt, pops)
      
      # Extract prior ancestry information
      prior.anc <- purrr::map(prior.txt, function(x) {
        # Initialize ancestry matrix
        anc.mat <- matrix(NA, nrow = maxpops, ncol = gensback + 1)
        rownames(anc.mat) <- paste("Pop", 1:nrow(anc.mat), sep = ".")
        colnames(anc.mat) <- paste("Gen", 0:gensback, sep = ".")
        
        # Parse ancestry information
        x <- sapply(strsplit(x[-1], "\\s|[:]"), function(y) {
          y <- y[y != ""]
          y[-1]
        })
        
        # Fill in ancestry matrix
        for (i in 1:ncol(x)) {
          pop <- as.numeric(x[1, i])
          anc.mat[pop, ] <- as.numeric(x[-1, i])
        }
        anc.mat
      }) %>% stats::setNames(df$id)
      
      # Calculate probability matrix
      prob.mat <- t(sapply(1:nrow(df), function(i) {
        pop.probs <- rowSums(prior.anc[[i]])
        pop.probs[is.na(pop.probs)] <- df$Group.1[i]
        pop.probs
      }))
      colnames(prob.mat) <- paste("Group", 1:ncol(prob.mat), sep = ".")
      df$Group.1 <- NULL
      df <- cbind(df, prob.mat)
      list(df = df, prior.anc = prior.anc)
    } else {
      NULL
    }
    
    # Extract dataframe from has.prior results
    has.prior.df <- if (is.null(has.prior)) {
      NULL
    } else {
      has.prior$df
    }
    
    # Combine results and normalize Q-matrix values
    q.mat <- rbind(no.prior, has.prior.df)
    q.mat <- q.mat[order(q.mat$row), ]
    q.mat$row <- NULL
    rownames(q.mat) <- NULL
    
    # Normalize Q-matrix rows to sum to 1
    q.mat[, -(1:3)] <- t(apply(q.mat[, -(1:3)], 1, function(i) i / sum(i)))
    
    # Extract prior ancestry information
    prior.anc <- if (is.null(has.prior)) {
      NULL
    } else {
      has.prior$prior.anc
    }
    
    # Return results
    return(list(
      summary = smry,
      q.mat = q.mat,
      prior.anc = prior.anc
    ))
  }
  
  # Create dataframe to track files and their K values
  file_info <- data.frame(f_name = file_paths, k = NA, rep = NA)
  
  if (verbose >= 2) {
    cat(report("Processing", length(file_paths), "STRUCTURE output files\n"))
  }
  
  # Extract K value (number of populations) from each file
  for (i in 1:nrow(file_info)) {
    # Read first 100 lines to find K value
    file_lines <- readLines(file_info[i, "f_name"], 100)
    k_line <- grep("populations assumed", file_lines, value = FALSE)
    
    # Check if K value was found
    if (length(k_line) > 0) {
      file_info[i, "k"] <- as.numeric(gsub("[^0-9]", "", file_lines[k_line]))
      if (verbose >= 2) {
        cat(report("File", i, "has K =", file_info[i, "k"], "\n"))
      }
    } else {
      warning(warn("Could not find K value in file:", file_info[i, "f_name"]))
    }
  }
  
  # Sort by K value
  file_info <- file_info[order(file_info$k), ]
  
  # Assign replicate numbers within each K value
  file_info <- file_info %>%
    dplyr::group_by(k) %>%
    dplyr::mutate(rep = dplyr::row_number())
  
  # Find common prefix for renaming
  label <- longest_common_prefix(file_names)
  if (verbose >= 2) {
    cat(report("Common prefix for files:", label, "\n"))
  }
  
  # Create new filenames
  file_info$rename <- paste(folder.path, "/", label, ".k", file_info$k, 
                            ".r", file_info$rep, "_out", sep = "")
  
  file_info <- as.data.frame(file_info)
  rownames(file_info) <- file_info$rename
  
  # Rename files
  if (verbose >= 2) {
    cat(report("Renaming files based on K values and replicates\n"))
  }
  rename_status <- file.rename(file_info$f_name, file_info$rename)
  if (!all(rename_status)) {
    warning(warn("Not all files were successfully renamed."))
  }
  
  # Process each renamed file
  if (verbose >= 2) {
    cat(report("Processing STRUCTURE output files\n"))
  }
  
  result_files <- lapply(rownames(file_info), function(file_path) {
    # Read the STRUCTURE file
    tryCatch({
      if (verbose >= 3) {
        cat(report("Processing file:", basename(file_path), "\n"))
      }
      
      result <- structureRead2(file = file_path)
      
      # Add file information
      result <- c(result, list(
        files = file_path, 
        label = basename(file_path)
      ))
      
      # Save result to an R data file
      output_file <- paste(file_path, ".ws.rdata", sep = "")
      save(result, file = output_file)
      
      return(output_file)
    }, error = function(e) {
      warning(warn("Error processing file", file_path, ":", e$message))
      return(NULL)
    })
  })
  
  # Remove NULL results from failed processing
  result_files <- result_files[!sapply(result_files, is.null)]
  
  # Check if any results were successfully processed
  if (length(result_files) == 0) {
    stop(error("No files were successfully processed."))
  }
  
  # Load saved results
  if (verbose >= 2) {
    cat(report("Loading processed STRUCTURE results\n"))
  }
  
  run_results <- lapply(result_files, function(file_path) {
    result <- NULL
    load(file_path)
    return(result)
  })
  
  # Name the results by their labels
  names(run_results) <- sapply(run_results, function(x) x$label)
  
  # Associate with population data if genlight object provided
  if (!is.null(x)) {
    if (verbose >= 2) {
      cat(report("Associating results with population data from genlight object\n"))
    }
    
    # Check that x is a genlight object
    datatype <- utils.check.datatype(x, verbose = 0)
    if (!datatype %in% c("SNP", "SilicoDArT")) {
      stop(error("The x parameter must be a genlight object containing SNP or SilicoDArT data"))
    }
    
    pop_tbl <- data.frame(id = indNames(x),
                          pop = pop(x))
    
    run_results <- lapply(run_results, function(y) {
      ytmp <- y[["q.mat"]]
      ytmp2 <- merge(ytmp, pop_tbl, by = "id")
      ytmp$orig.pop <- ytmp2$pop
      y[["q.mat"]] <- ytmp
      return(y)
    })
  }
  
  # Set class for the results object
  class(run_results) <- c("structure.result", class(run_results))
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("\nCompleted:", funname, "\n"))
  }
  
  # RETURN
  return(run_results)
}