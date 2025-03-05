#' Select Loci Panel Based on Various Methods
#'
#' @description
#' This function selects a panel of loci from a genomic dataset (`genlight` object)
#' based on various selection methods.
#'
#' @param x A `dartR\genlight` object containing the genomic data.
#' @param method A character string specifying the selection method. Options include:
#'   \itemize{
#'     \item `"dapc"`: Select loci contributing most to discrimination between populations using DAPC (Discriminant Analysis of Principal Components).
#'     \item `"pahigh"`: Select loci with private alleles having high frequency.
#'     \item `"random"`: Randomly select loci.
#'     \item `"monopop"`: Select monomorphic loci within populations.
#'     \item `"stratified"`: Stratified sampling of loci based on allele frequencies.
#'     \item `"hafall"`: Select loci with the highest allele frequencies across all populations.
#'     \item `"hafpop"`: Select loci with the highest allele frequencies within each population.
#'   }
#' @param nl An integer specifying the number of loci to select.
#' @param plot.out Logical. If `TRUE`, generates plots summarizing selected loci.
#' @param plot.file A character string specifying the file name for saving plots. If `NULL`, plots are not saved.
#' @param plot.dir A character string specifying the directory to save plots. Defaults to the working directory.
#' @param verbose Integer level of verbosity for reporting progress and information.
#'
#' @details
#' The function applies various methods to select loci based on the input `genlight` object.
#' Each method has specific criteria for selecting loci:
#' \itemize{
#'   \item `dapc`: Performs DAPC and identifies loci with the highest contributions to discrimination between population pairs.
#'   \item `pahigh`: Identifies loci with private alleles that have high frequency differences between populations.
#'   \item `random`: Selects loci randomly.
#'   \item `monopop`: Selects loci that are monomorphic within populations.
#'   \item `stratified`: Uses stratified sampling to select loci based on allele frequencies.
#'   \item `hafall`: Selects loci with the highest allele frequencies across the dataset.
#'   \item `hafpop`: Selects loci with the highest allele frequencies within individual populations.
#' }
#'
#' @return A `genlight` object containing the selected loci.
#'
#' @examples
#' # Example usage:
#'
#' # Select 20 loci randomly
#' selected <- gl.select.panel(possums.gl, method = "random", nl = 50)
#'
#' # Select loci based on DAPC
#' selected <- gl.select.panel(possums.gl, method = "dapc", nl = 5)
#'
#' @export

gl.select.panel <- function(x, 
                            method="random", 
                            nl=10,
                            plot.out = TRUE,
                            plot.file = NULL,
                            plot.dir = NULL,
                            verbose = NULL) {
  
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # FLAG SCRIPT START
  
  funname <- match.call()[[1]]
  utils.flag.start(
    func = funname,
    build = "Jody",
    verbose = verbose
  )
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  #1. dapc
  #2. private allele high frequency
  #3. random
  #4. monomorph by pops
  #5. stratified 
  #6. highest allele frequencies overall
  #7. highest allele frequencies per pop
  
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  ### DO THE JOB
  
  x <- x[order(pop(x)),]
  res <- list()
  
  if (method=="dapc"){ 
    com <- t(combn(nPop(x), 2))
    pops <- seppop(x)
    nl <- ceiling(nl/nrow(com)) 
        
    for (i in 1:nrow(com)){
      dummy <- pops[c(com[i,1],com[i,2])] 
      dummy <-do.call(rbind,dummy)
      dd <- dapc(dummy, n.pca=20, n.da=5)
      ll <- sort(dd$var.contr[,1], decreasing=TRUE)
      index <-names(ll)[1:nl]
      res[[i]] <- index
      
    }
    
    #selpos <- paste0("X",unlist(res))
    selloc <- unique(unlist(res))
    #get rid of .locall
    selloc <- strsplit(selloc, "\\.")
    selloc <- unlist(lapply(selloc, function(x) x[1]))
    
    
  }
  
  if (method=="pahigh"){
    prxx <- gl.report.pa(x, loc.names = TRUE, verbose = 0)
    qq <-prxx$names_loci
    panxx <- lapply(qq, function(x) list(pa1=x$pop1_pop2, pa2=x$pop2_pop1))
    
    com <- t(combn(nPop(x), 2))
    pops <- seppop(x)
    nl <- ceiling(nl/(nrow(com)*2))
        
        
    res <- list()
    for (i in 1:nrow(com)){
      pas <- panxx[[i]]$pa1
      p1 <- pops[[com[i,1]]]
      p2 <- pops[[com[i,2]]]
      p1p <- gl.keep.loc(p1, loc.list = pas, verbose = 0)
      p2p <- gl.keep.loc(p2, loc.list = pas, verbose = 0) 
      res1<- names(sort(rowSums(cbind(gl.alf(p1p)* !gl.alf(p2p))), decreasing=TRUE)[1:nl])
      
      pas <- panxx[[i]]$pa2
      p1 <- pops[[com[i,2]]]
      p2 <- pops[[com[i,1]]]
      p1p <- gl.keep.loc(p1, loc.list = pas, verbose = 0)
      p2p <- gl.keep.loc(p2, loc.list = pas, verbose = 0) 
      res2<- names(sort(rowSums(cbind(gl.alf(p1p)* !gl.alf(p2p))), decreasing=TRUE)[1:nl]) 
      res[[i]] <- c(res1,res2)
      
    }
    
    selloc <- unique(unlist(res))
    
    
    
  }
  if (method=="random"){
    #random selection
    selloc <- locNames(x)[sample(nLoc(x), nl, replace = FALSE)]
  }
  if (method=="monopop"){ 
    res <- list()
    pops <- seppop(x)
    nl <- nl/length(pops)
    
    for (i in 1:nPop(x)){
      dummy <- pops[[i]]
      index <- abs(colMeans(as.matrix(dummy), na.rm=TRUE)-1)==1
      dl<-  locNames(dummy)[index]
      mons <- sample(dl, nl, replace=FALSE)
      res[[i]] <- sample(mons, nl)
    }
    
    
    selloc <- unique(unlist(res))
    
    }
  if (method=="stratified"){
    res <- list()
    pops <- seppop(x)
    nl <- nl/length(pops)
    for (i in 1:nPop(x)){
      
      dummy <- pops[[i]]
      
      df <- data.frame(id=locNames(dummy), freq=0.5-(abs(gl.alf(dummy)[,1]-0.5)))
      df <- df[order(df$freq),]
      
      self <- seq(0,0.5, length=nl+1)
      
      cf <- cut(df$freq, breaks=self, include.lowest=TRUE)
      
      scf <- split(df$id, cf)
      dres <- list()
      for (ii in 1:length(scf)){
        dres[[ii]] <- sample(scf[[ii]], 1)
      }
      res[[i]] <- unlist(dres)
    }
    
    
    selloc <- unique(unlist(res))
  }
  if (method=="hafall"){
    
   index <- order(0.5-abs(0.5-gl.alf(x)[,1]), decreasing = TRUE)
    selloc <- locNames(x)[index[1:nl]]
    
    
  }
  if (method=="hafpop"){ 
    
    res <- list()
    pops <- seppop(x)
    nl <- nl/length(pops)
    for (i in 1:nPop(x)){
      
      dummy <- pops[[i]]
      index <- order(0.5-abs(0.5-gl.alf(pops[[i]])[,1]), decreasing = TRUE)
      res[[i]] <- locNames(pops[[i]])[index[1:nl]]
    }
    
    
    selloc <- unique(unlist(res))
    
    }
  
  
  
  #filter object to keep only selected loci
  xx <- gl.keep.loc(x, selloc, verbose = verbose)
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(xx)
  
  
  
}

