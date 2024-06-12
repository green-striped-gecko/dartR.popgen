#################################################################################################
### 1. MAIN FUNCTION. Employs gradient descent to find optimal values for
###                   allele proportion (p) and census population (N)
#################################################################################################
Find.Population <- function(allele_frequencies, Niter = 100, learn_rate = 0.2, zero_level = 5e-4, plot=FALSE){
    # Finds the minimum in therms of the number of each allele
    x1 <- 0
    x2 <- 0
    for(i in 1:length(allele_frequencies)){
      x1 <- x1 + (i-1)*allele_frequencies[i]
      x2 <- x2 + allele_frequencies[i]*(i-1)**2
    }

    s <- length(allele_frequencies)
    p0 <- x1/s
    a <- (x2-x1)/x1
    N0 <- (a-s+1)/(a+p0-x1)
    N1 <- s + (N0-s)/2.
           
    pN0 <- p0*N0
    qN0 <- (1-p0)*N0

    pN_start <- p0*N1
    qN_start <- (1-p0)*N1

    # Likelihood function
    SE <- function(pN, qN){
            s <- length(allele_frequencies) - 1
            hyp <- new.HyperGeom(s, pN, qN)
            # Can be changed to desired metric
            SE <- abs(allele_frequencies - hyp)**2.
            return( log(sum(SE)) )
        }

    # Calls Gradient Descent
    v <- Gradient.Descent(  SE,
                            x0 = pN_start, 
                            y0 = qN_start, 
                            Niter = Niter, 
                            learn_rate = learn_rate,
                            zero_level = zero_level,
                            plot=plot
                        )
    # Obtains the allele proportion and population estimates
    Min <- which.min(c(SE(pN0,qN0), SE(v[1],v[2])))
    if(Min==2){
      p <- v[1]/(v[1]+v[2])
      N <- v[1] + v[2]
    }else{
      p <- p0
      N <- N0
    }
    return( c(p,N) )
}

###########################################################################################
######  2. Find Minimum point of a 2D-function using Gradient Descent
###########################################################################################
Gradient.Descent <- function(ff, x0, y0, Niter = 100, hx = .1, hy = .1, learn_rate = .2, zero_level = 5.e-4, plot=TRUE){
    # Employs matrix operations for easier readability 
    vk <- matrix(c(x0,y0),2,1)
    loss_function <- c()
    Nest <- c()
    count <- 0
    for(i in 1:Niter){
        df <- grad(ff,vk,hx,hy)
        #---------------------------------
        # Gradient descent
        vk <- vk - learn_rate*df/sqrt(norm(df))
        #--------------------------------
        loss_function[i] <- ff(vk[1],vk[2])
        Nest[i] <- vk[1] + vk[2]
        count <- count +1
        if(i>1){
          if(!is.null(loss_function[i-1]) && abs((loss_function[i] - loss_function[i-1])/loss_function[1]) < zero_level) break
        }
        
    }
    if(plot){
        dfl <- data.frame('Niter'=1:count,"Loss"=1 - abs((loss_function - loss_function[1])/loss_function[1]), 'DNest'=1-abs((Nest - tail(Nest,1))/tail(Nest,1))) 
        p <- ggplot(dfl,aes(Niter)) + 
             geom_line(aes(y = Loss,color="Loss")) + 
             geom_line(aes(y = DNest,color="DNest"))
        print(p)
    }
    return(c(vk[1],vk[2]))
}

###########################################################################################
# 2.1 Supplementary. Computes the gradient of a 2d function
grad <- function(ff,v,hx,hy){
    x <- v[1]
    y <- v[2]
    f <- ff(x,y)
    # First derivative X
    fi <- ff(x+hx,y)
    f_i <- ff(x-hx,y)
    fx <- (fi - f_i)/(2*hx)
    # First derivative Y
    fk <- ff(x,y+hy)
    f_k <- ff(x,y-hy)
    fy <- (fk - f_k)/(2*hy)
    df <- matrix(c(fx,fy),2,1)
    return(df)
}


#########################################################################################
### 3. Suplementary. Modified Hypergeometric function with number of allele-types as input
new.HyperGeom <- function(narg, pN, qN) {
  Narg <- pN + qN
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
  
  return(abs(HypExp))
}


###########################################################################################
### 4. This method its unstable for this particular problem (negative Hessian).
###    It could be used in other projects (converges very fast to the solution)
###########################################################################################
####### Find Maximal point of a 2D-function using Newton-Rapson Method (Main)
Newton.Rapson <- function(ff, x0, y0, Niter = 10, hx = .1, hy = .1, zero_level = 1.e-12){
    # Employs matrix operations for easier readability (overkill in 2D)
    vk <- matrix(c(x0,y0),2,1)
    for(i in 1:Niter){
        derivatives_f <- derivate(ff,vk,hx,hy)
        df <- matrix(head(derivatives_f,2),2,1)
        H  <- matrix(tail(derivatives_f,4),2,2)
        # Check if Hessian is invertible
        if(abs(det(H)) < zero_level) break
        # Check for non-zero gradient
        #print(derivatives_f)
        if(norm(df)<zero_level) break
        #--------------------------
        # Newton-Rapson Method
        vk <- vk - solve(H) %*% df
        #print(vk)
    }
    return(c(vk[1],vk[2]))
}

###########################################################################################
# 4.1   Calculates finite differences for NR method (suplementary)
derivate <- function(ff,v,hx,hy){
    x <- v[1]
    y <- v[2]
    f <- ff(x,y)
    # First derivative X
    fi <- ff(x+hx,y)
    f_i <- ff(x-hx,y)
    fx <- (fi - f_i)/(2*hx)
    # First derivative Y
    fk <- ff(x,y+hy)
    f_k <- ff(x,y-hy)
    fy <- (fk - f_k)/(2*hy)
    # Second derivative
    fxx <- (fi - 2.*f + f_i)/(hx**2.) # X
    fyy <- (fk - 2.*f + f_k)/(hy**2.) # Y
    # Crossed second derivative
    fik <- ff(x+hx,y+hy)
    f_ik <- ff(x-hx,y-hy)
    fxy <- (fik - fi - fk + 2.*f - f_i - f_k + f_ik)/(2.*hx*hy)
    # Store values in list to easily define Hessian and gradient
    return(c(fx,fy,fxx,fxy,fxy,fyy))
}
###############################################################################