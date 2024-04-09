#' Hypergeometric  probability for DNADot
#' 
#' This function is used internally to calculate the hypergeometric expectation 
#' needed for DNADot
#'
#'
#'
#' @param parg allele freq being trialled (Ptry)
#' @param Narg Pop size being trialled (Ntry)
#' @param narg Total number of alleles: j <- round(2 * n * jj) where n is the 
#'    sample size and JJ is the proportion of the subsample being trialled
#' @return A vector of probabilities
#' @author Adapted from Sherwin's MatLab code
HypExactExp230607 <- function(parg, Narg, narg) {
  q <- 1 - parg
  qN <- round(q * Narg)
  pN <- round(parg * Narg)
  
  HypExp <- rep(0, narg + 1)
  
  for (x1 in 1:(narg + 1)) {
    x <- x1 - 1
    notx <- narg - x
    notx1 <- notx + 1
    
    qNpNdetect <- rep(1, narg)
    qNundetect <- rep(0, notx)
    pNdetect <- rep(0, ifelse(x == 0, 1, x))
    nN <- rep(0, narg)
    nNDenom <- rep(1, narg)
    
    for (det1 in 1:x) {
      det <- det1 - 1
      if(det1 == 0) det1 <- 1 # to prevent error, but this is actually ignore in the results when x=0
      # pNdetect[det1] <- (pN - det) / (x - det)
      if (is.infinite((pN - det) / (x - det))) pNdetect[det1] <- 1 else
        pNdetect[det1] <- (pN - det) / (x - det)
    }
    
    for (undet1 in 1:notx) {
      undet <- undet1 - 1
      # qNundetect[undet1] <- (qN - undet) / (notx - undet)
      if (is.infinite((qN - undet) / (notx - undet))) qNundetect[undet1] <- 1 else
        qNundetect[undet1] <- (qN - undet) / (notx - undet)
    }
    
    if (notx == narg) {
      qNpNdetect[1:notx] <- rev(qNundetect)
    } else if (x == narg) {
      qNpNdetect[1:x] <- rev(pNdetect)
    } else {
      qNpNdetect[1:notx] <- rev(qNundetect)
      qNpNdetect[(notx + 1):narg] <- rev(pNdetect)
    }
    
    for (nn in 1:narg) {
      nN[nn] <- (Narg - nn + 1) / (narg - nn + 1)
    }
    
    nNDenom[1:narg] <- rev(nN[1:narg])
    
    HypElements <- qNpNdetect / nNDenom
    HypExp[x1] <- prod(HypElements)
  }
  
  return(HypExp)
}