#' Population size estimates using DNADot approach
#' 
#' @param name description
#' @param n.cores The number of cores to use. If "auto", it will 
#' use all but one available cores [default "auto"], if these are < number of 
#' datsets, otherwise it will be equal to the number of datesets.
#' @import readxl 
#' @import dplyr
#' @author Carlo Pacioni, Adapted from Sherwin's MatLab code
gl.DNADot <- function(x=NULL, gen.file=NULL, header=FALSE, nonGenCols=NULL,
                      jj=0.7, minNtry, ppinc=0.05, validate=TRUE, pvalidate=0.5,
                      n.cores = "auto",
                      verbose) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jackson",
                   verbose = verbose)
  
  # CHECK DATATYPE
  if (!is.null(x)) {
    dt <- utils.check.datatype(x, verbose = 0)
    if(st != "SNP") stop(error("This function is SNP data only"))
    cat(report("  Analysis performed on the genlight object.\n"))
  }
  
  # avoid global binding error
  Bin <-
    r <-
    L.r <-
    U.r <- L.r.null <- U.r.null <- Freq <- Var1 <-  NULL 
  
  # set number of cores
  if(n.cores != 1) {
    if(n.cores == "auto"){
      n.cores <- parallel::detectCores() - 1
    }
  
  # sort out input data
  # If genotype fiels are provided
  if(is.null(x)) {
    extensions <- sapply(gen.file, grepl, pattern="xls$|xlsx$")
    if(all(extensions)) {
      InputData <- lapply(gen.file, read_excel, col_names = header)
      InputData <- lapply(InputData, data.table) # Ensure that data are a datable
      } else {
        if(any(extensions)) {
        stop(error("  It seems that some genotype files are in 
                                     excel, but others are not. Can process excel files only
                                     if all files are in excel"))
          } else {
            InputData <- lapply(gen.file, data.table::fread, header = header)
          }
        names(InputData) <- tools::file_path_sans_ext(basename(gen.file))
      }
    
    if(!is.null(nonGenCols)) 
      InputData <- lapply(InputData, function(x) {x[, -nonGenCols]})
  }
  
  # if a genlight object is provided #
  if (!is.null(x) & is(x, "genlight")) {
    pop_list <- seppop(x)
    formatDNADot.in <- function(x) {
      DNADotIn <- gl2related(glNoSecLD, save = FALSE, verbose = 0)
      DNADotIn <- DNADotIn[, -1]
    }
    
    InputData <- lapply(pop_list, formatDNADot.in)
  }
  
  # Adjust length of Ntry if there are multiple pops 
  if(length(minNtry)==1) {
    Ntry <- rep(minNtry, length(InputData))
  } else {
    # Check for error in Ntry's length
    if(length(minNtry)!=length(InputData)) 
      stop(error("The length of 'Ntry' and input data is not the same"))
    }
    
  # Check if Ntry values are adequate
  nSamples <- sapply(InputData, nrow)
  if(sum((minNtry - nSamples)>0)!=length(InputData)) 
    stop(error("The values of 'Ntry' have to be > than the sample size in the respective population"))
  
  # DO THE JOB #
  # Randomise
  randomise <- function(x) {
    DNADotIn <- x[sample(seq_len(nrow(x)), 
                                size = nrow(x), 
                                replace = FALSE), ]
    return(DNADotIn)
  }
  InputData <- lapply(InputData, randomise)
  
  if(validate) {
    InputDataV <- lapply(InputData, function(x) {
      x[rbinom(nrow(x), 1, pvalidate),]
    })
    names(InputDataV) <- paste0(names(InputData), "V")
    InputData <- c(InputData, InputDataV)
  }
  
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
  Ntry <- seq(minNtry, maxNtry, by = Ninc)
  ppinc <- ppinc
  Ptry <- seq(ppinc, 1 - ppinc, by = ppinc)
  
  if(n.cores > 1 & length(InputData) > 1) {
    if(length(InputData) < n.cores) n.cores <- length(InputData)
    cl <- parallel::makeCluster(n.cores)
    on.exit(expr=parallel::stopCluster(cl))
    catch <- parallel::clusterEvalQ(cl, library("dartR.popgen"))
    
    parallel::clusterExport(cl, 
                            varlist=c("indit", "Ptry", "Ntry", "jj"), 
                            envir=environment()) 
    
    res <- parallel::parLapply(cl = cl, X = indit, fun = utils.DNADot, 
                               Ptry=Ptry, Ntry=Ntry, jj = jj)
  } else {
    res <- utils.DNADot(indit, Ptry, Ntry, jj)
  }
  return(res)
}