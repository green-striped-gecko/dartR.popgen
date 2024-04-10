#' Population size estimates using DNADot approach
#' 
#' @param name description
#' @param n.cores The number of cores to use. If "auto", it will 
#' use all but one available cores [default "auto"], if these are < number of 
#' datsets, otherwise it will be equal to the number of datesets.
#' @import readxl 
#' @import dplyr
#' @export
#' @author Carlo Pacioni, Adapted from Sherwin's MatLab code
gl.DNADot <- function(x=NULL, gen.file=NULL, header=FALSE, nonGenCols=NULL,
                      jj=0.7, minNtry, ppinc=0.05, validate=TRUE, pvalidate=0.5,
                      n.cores = "auto", outfile = "DNADot_out.csv",
                      outpath = tempdir(),
                      verbose=NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jackson",
                   verbose = verbose)
  
  # CHECK DATATYPE
  if (!is.null(x) & is(x, "genlight")) {
    dt <- utils.check.datatype(x, verbose = 0)
    if(dt != "SNP") stop(error("This function is SNP data only"))
    cat(report("  Analysis performed on the genlight object.\n"))
  }
  
  # set number of cores
  if(n.cores != 1) {
    if(n.cores == "auto"){
      n.cores <- parallel::detectCores() - 1
    }
  }
  
  # sort out input data
  #### Genotype files ####
  if(is.null(x)) {
    cat(report("  Analysis performed on the genotype file(s)\n"))
    extensions <- sapply(gen.file, grepl, pattern="xls$|xlsx$")
    if(all(extensions)) { # if excel
      InputData <- lapply(gen.file, read_excel, col_names = header)
      InputData <- lapply(InputData, data.table) # Ensure that data are a datable
      } else {
        if(any(extensions)) {
        stop(error("  It seems that some genotype files are in 
                                     excel, but others are not. Can process excel files only
                                     if all files are in excel"))
          } else { # any other data file
            InputData <- lapply(gen.file, data.table::fread, header = header)
          }
        names(InputData) <- tools::file_path_sans_ext(basename(gen.file))
      }
    
    if(!is.null(nonGenCols)) 
      InputData <- lapply(InputData, function(x) {x[, -nonGenCols]})
  }
  
  #### Genlight ####
  if (!is.null(x) & is(x, "genlight")) {
    pop_list <- seppop(x)
    formatDNADot.in <- function(x) {
      DNADotIn <- gl2related(x, save = FALSE, verbose = 0)
      DNADotIn <- DNADotIn[, -1]
    }
    
    InputData <- lapply(pop_list, formatDNADot.in)
  }
  
  #### Genind ####
  if (!is.null(x) & is(x, "genind")) {
    cat(report("  Analysis performed on the geninf object.\n"))
    pop_list <- seppop(x)
    InputData <- lapply(pop_list, genind2df, oneColPerAll = TRUE)
    InputData <- lapply(InputData, function(x) x[,-1])
  }
  
  # Adjust length of Ntry if there are multiple pops 
  if(length(minNtry)==1) {
    minNtry <- rep(minNtry, length(InputData))
  } else {
    # Check for error in Ntry's length
    if(length(minNtry)!=length(InputData)) 
      stop(error("The length of 'Ntry' is neither 1 or the number of the populations\n"))
    }
    
  # Check if Ntry values are adequate
  nSamples <- sapply(InputData, nrow)
  if(sum((minNtry - nSamples)>0)!=length(InputData)) 
    stop(error("The values of 'Ntry' have to be > than the sample size in the respective population\n"))
  
  # DO THE JOB #
  #### Randomise ####
  randomise <- function(x) {
    DNADotIn <- x[sample(seq_len(nrow(x)), 
                                size = nrow(x), 
                                replace = FALSE), ]
    return(DNADotIn)
  }
  InputData <- lapply(InputData, randomise)
  
  #### Validate ####
  if(validate) {
    InputDataV <- lapply(InputData, function(x) {
      x[rbinom(nrow(x), 1, pvalidate),]
    })
    names(InputDataV) <- paste(names(InputData), "V", sep = "_")
    InputData <- c(InputData, InputDataV)
    minNtry <- c(minNtry, minNtry)
  }
  
  # Update sample size
  nSamples <- sapply(InputData, nrow)
  
  process.input <- function(x) {
    Input2 <- t(x)
    RowCol <- dim(Input2)
    TwiceLoc <- RowCol[1]
    Inds <- RowCol[2]
    Input3 <- Input2
    Input4 <- Input2
    Input3 <- Input3[-seq(1, TwiceLoc, by = 2), ]
    Input4 <- Input4[-seq(2, TwiceLoc, by = 2), ]
    # preindit row=loci, col=inds*2 for alleles 1&2
    preindit <- cbind(Input3, Input4)
    return(preindit)
  }
  
  indit <- lapply(InputData, process.input)
  
  # Increments Ptry, Ntry (hypotheses for joint est. of Nc & p)
  minNtry <- 2 * minNtry # diploid data
  maxNtry <- 3 * minNtry
  Ninc <- round((maxNtry - minNtry) / 10)
  Ntry <- mapply(seq, minNtry, maxNtry, by = Ninc, SIMPLIFY = FALSE)
  ppinc <- ppinc
  Ptry <- seq(ppinc, 1 - ppinc, by = ppinc)
  
  #### Parallele execution ####
  if(n.cores > 1 & length(InputData) > 1) {
    if(length(InputData) < n.cores) n.cores <- length(InputData)
    cl <- parallel::makeCluster(n.cores)
    on.exit(expr=parallel::stopCluster(cl))
    catch <- parallel::clusterEvalQ(cl, library("dartR.popgen"))
    
    parallel::clusterExport(cl, 
                            varlist=c("indit", "Ptry", "minNtry", "maxNtry", 
                                      "Ntry", "jj"), 
                            envir=environment()) 
    
    res <- parallel::clusterMap(cl = cl, fun = utils.DNADot, indit,  
                                minNtry=minNtry, maxNtry, Ntry=Ntry, 
                               MoreArgs = list(Ptry=Ptry, jj = jj)
                               )
    res <- do.call(cbind, res)
  } else {
    #### Serial execution ####
    res <- mapply(utils.DNADot, indit, minNtry, maxNtry, Ntry, 
                  MoreArgs = list(Ptry=Ptry, jj=jj))
  }
  res <- cbind(Param=c("minNtry", "maxNtry", "jj", "Nest", "SD", "SE", "n"), 
               rbind(data.frame(round(res, 1)), nSamples))

  write.csv(res, file(file.path(outpath, outfile)), row.names = FALSE)
  
  return(res)
  }
  