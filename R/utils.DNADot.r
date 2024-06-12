#' Run DNADot analysis
#' 
#' This function is used internally to run the DNADot analysis
#'
#'
#'
#' @param indit Input data for 
#' @param maxNtry Maximum value of N to be tested
#' @param Ntry Values of N to be trialed
#' @param Ptry Allele frequencies values to be trialed
#' @inheritParams gl.DNADot
#' @inheritParams utils.HyperGeometricP
#' @importFrom stats sd
#' @return A vector of probabilities
#' @author Carlo Pacioni, Adapted from Sherwin's MatLab code
utils.DNADot <- function(indit, minNtry, maxNtry, Ntry, Ptry, jj) {
  
  locind <- dim(indit)
  L1 <- locind[1] # Number of loci
  n <- locind[2] / 2 # Number of samples
  
  # prep for loop
  L <- 0
  alcodelen <- vector("numeric", length = L1)                   # new
  alcodelenvector <- apply(indit, 1, FUN = function(x) length(unique(x)))
  indal <- matrix(NA, ncol = ncol(indit), nrow = sum(alcodelenvector))
  for (l1 in 1:L1) {
    alcode <- sort(unique(indit[l1, ]))
    alcodelen[l1] <- length(alcode)
    for (al in 1:alcodelen[l1]) {
      LTEMP <- L + al
      indal[LTEMP, ] <- indit[l1, ] == alcode[al]
    }
    L <- LTEMP
  }
  
  indalTEMP1 <- indal
  indalq <- (indalTEMP1 - 1) * (-1)
  indalfold <- rbind(indal, indalq)
  p <- vector("numeric", length = L)
  for (l in 1:L) {
    sumTEMP <- sum(indalfold[l, ])
    nTEMP <- 2 * n
    p[l] <- sumTEMP / nTEMP
  }
  
  j <- round(2 * n * jj) # number of alleles
  HypExp <- array(0, dim = c(length(Ntry), length(Ptry), j + 1))
  propdetNit <- matrix(0, nrow = L, ncol = j + 1)
  BIAS1 <- array(0, dim = c(length(Ntry), length(Ptry), L))
  
  for (Nindex in 1:length(Ntry)) {
    for (pindex in 1:length(Ptry)) { # loop Hyp Exp in sample of 'j',
      parg <- Ptry[pindex]
      Narg <- Ntry[Nindex]
      narg <- j
      HypExpTEMP1 <- utils.HyperGeometricP(parg, Narg, narg)
      HypExp[Nindex, pindex, 1:(j + 1)] <- HypExpTEMP1
    }
  }
  
  
  counter1 <- 1
  edges <- numeric(2*n-j+2)
  for (e in seq(-0.5, j+0.5, by=1)) {
    edges[counter1] <- e
    counter1 <- counter1 + 1
  }
  
  indalTEMP <-  matrix(NA, nrow = L, ncol=j)
  detNit <- matrix(NA, nrow = L, ncol = (2*n-j+1))
  for (l in 1:L) {
    for (jndx in 1:(2*n-j+1)) {
      indalTEMP[l, 1:j] <- indalfold[l, jndx:(jndx+j-1)] # takes a long time to execute
      detNit[l, jndx] <- sum(indalTEMP[l, 1:j]) # takes a long time to execute
    }
  }
  
  for (l in 1:L) {
    binnedFrq <- hist(detNit[l,], breaks=edges,plot=FALSE)$counts
    #binnedFrq <- table(cut(detNit[l,], breaks=edges)) # Thought this was faster but it seems not
    PropNit <- binnedFrq / sum(binnedFrq)
    propdetNit[l, 1:(j+1)] <- PropNit
  }
  
  for (l in 1:L) {
    for (Nindex in 1:length(Ntry)) {
      for (pindex in 1:length(Ptry)) {
        HypExpTEMP2 <- rep(0, j+1)
        aa <- HypExp[Nindex, pindex, 1:(j+1)]
        HypExpTEMP2 <- as.numeric(aa)
        PropTEMP <- propdetNit[l,]
        BiasTEMP1 <- abs(PropTEMP - HypExpTEMP2)  
        BiasTEMP2 <- sum(BiasTEMP1)
        BIAS1[Nindex, pindex, l] <- BiasTEMP2
      }
    }
  }
  


  Nest <- vector("numeric", length = L)
  Pest <- vector("numeric", length = L)
  for (l in 1:L) {
    locBIAS <- BIAS1[, , l]
    Min <- min(locBIAS)
    index <- which(locBIAS == Min, arr.ind = TRUE) #
    Nest[l] <- mean(Ntry[index[1]])
    Pest[l] <- mean(Ptry[index[2]])
    print(c(l,Pest[l],Nest[l]))
  }
  

  AveNestDiploid <- mean(Nest)
  SeNestDiploid <- sd(Nest) / sqrt(L)
  AveNest <- AveNestDiploid / 2
  SeNest <- SeNestDiploid / 2
  
  SDdetNit <- apply(detNit, 1, sd)
  AvedetNit <- apply(detNit, 1, mean)
  CVdetNit <- (SDdetNit / AvedetNit)
  
  LocNbr <- 1:L
  CVNest <- cbind(CVdetNit[1:L], Nest[1:L], LocNbr)
  CVNestSort <- t(cbind(CVNest[order(CVNest[,1]), ]))
  NestLoCV <- CVNestSort[2, 1:round(L/10)]
  CVLoLoci <- CVNestSort[3, 1:round(L/10)]


  AveNestLoCVdiploid <- mean(NestLoCV)
  AveNestLoCV <- AveNestLoCVdiploid / 2


  SDNestLoCVdiploid <- sd(NestLoCV)
  SDNestLoCV <- SDNestLoCVdiploid / 2
  SENestLoCV <- SDNestLoCV / sqrt(L)

  
  Output <- c(minNtry/2, maxNtry/2, jj, AveNestLoCV, SDNestLoCV, SENestLoCV)
  # write.csv(Output, file="CensusOutputTest.csv")
  return(Output)
  
  
}