#' @name gl.select.panel
#'
#' @title Select a panel of loci based on various methods
#'
#' @family panel selection
#'
#' @description
#' This function selects a panel of loci from a SNP genlight object based on
#' various selection methods, for example to design a targeted genotyping
#' panel. Use \code{\link{gl.check.panel}} to compare the panel with the
#' full data.
#'
#' @param x A genlight object with SNP data [required].
#' @param method The selection method, one of "random", "dapc", "pahigh",
#'   "monopop", "stratified", "hafall", "hafpop", "pic" or "picdart" (see
#'   Details) [default "random"].
#' @param nl The number of loci to select, between 1 and nLoc(x)
#'   [default 10].
#' @param exact Logical. If TRUE, the panel has exactly \code{nl} loci:
#'   surplus loci are dropped at random, and missing ones are added at random
#'   from the remaining loci. If FALSE, the panel holds what the method
#'   selected, which can be more or fewer than \code{nl}; the function stops
#'   if the method selected no loci [default TRUE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'   brief progress messages; 3, progress and results summary; 5, full report
#'   [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' Each method has specific criteria for selecting loci:
#' \itemize{
#'   \item "random": selects loci at random.
#'   \item "dapc": for each pair of populations, runs a DAPC and takes the
#'   ceiling(nl / number of pairs) loci contributing most to the first
#'   discriminant function.
#'   \item "pahigh": for each pair of populations and each direction, takes
#'   the ceiling(nl / (2 * number of pairs)) private-allele loci whose private
#'   alleles have the highest frequency.
#'   \item "monopop": for each population, takes up to ceiling(nl / number of
#'   populations) loci at random among those monomorphic in that population.
#'   \item "stratified": for each population, bins loci by minor allele
#'   frequency into ceiling(nl / number of populations) classes between 0 and
#'   0.5 and takes one locus at random from each class.
#'   \item "hafall": takes the loci with the highest minor allele frequency
#'   (closest to 0.5) across all individuals.
#'   \item "hafpop": for each population, takes the ceiling(nl / number of
#'   populations) loci with the highest minor allele frequency in that
#'   population.
#'   \item "pic": takes the loci with the highest polymorphic information
#'   content, 1 - (p^2 + q^2) - 2 p^2 q^2.
#'   \item "picdart": takes the loci with the highest AvgPIC, recalculated
#'   with gl.recalc.metrics.
#' }
#' The per-population and per-pair methods can select more than \code{nl}
#' loci (rounding up, then combining); with \code{exact = TRUE} the surplus
#' is dropped at random, so the ranking is not used to choose which loci go.
#' "random", "monopop", "stratified" and the \code{exact} adjustment use
#' R's random number generator; call \code{set.seed()} first for a
#' reproducible panel.
#'
#' All methods except "random", "hafall", "pic" and "picdart" need
#' populations; "dapc" and "pahigh" need at least two. "dapc" leaves out,
#' for each pair, the loci with no calls in that pair.
#'
#' @return A genlight object with the selected loci; individuals are in the
#'   same order as in \code{x}.
#'
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#'   \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # Select 50 loci randomly
#' selected <- gl.select.panel(possums.gl, method = "random", nl = 50)
#'
#' # Select 5 loci based on DAPC
#' selected <- gl.select.panel(possums.gl, method = "dapc", nl = 5)
#'
#' @export
#' @importFrom utils combn

gl.select.panel<-
  function(x,
           method="random",
           nl=10,
           exact=TRUE,
           verbose = NULL) {


    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    # FLAG SCRIPT START

    funname <- match.call()[[1]]
    utils.flag.start(func = funname, verbose = verbose)

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

    # FUNCTION SPECIFIC ERROR CHECKING

    methods <- c("random", "dapc", "pahigh", "monopop", "stratified",
                 "hafall", "hafpop", "pic", "picdart")
    if (!is.character(method) || length(method) != 1 ||
        !method %in% methods) {
      stop(error("'method' must be one of:", paste(methods, collapse = ", "),
                 "\n"))
    }
    if (!is.numeric(nl) || length(nl) != 1 || is.na(nl) || nl < 1 ||
        nl != round(nl) || nl > nLoc(x)) {
      stop(error(paste0("'nl' must be a whole number between 1 and nLoc(x) (",
                        nLoc(x), ").\n")))
    }
    pop_methods <- c("dapc", "pahigh", "monopop", "stratified", "hafpop")
    if (method %in% pop_methods) {
      if (is.null(pop(x)) || any(is.na(pop(x)))) {
        stop(error(paste0("Method '", method, "' needs every individual ",
                          "assigned to a population.\n")))
      }
      if (method %in% c("dapc", "pahigh") && nPop(x) < 2) {
        stop(error(paste0("Method '", method, "' needs at least two ",
                          "populations.\n")))
      }
    }

    # Allele frequencies (reference, alternate) per locus from the genotype
    # matrix. gl.alf() runs the dartR.base datatype check, which stops on
    # population subsets whose genotypes are all 0/1 when the object carries
    # no SNP metadata.
    alf <- function(g) {
      q <- colMeans(as.matrix(g), na.rm = TRUE) / 2
      q[is.nan(q)] <- NA
      cbind(alf1 = 1 - q, alf2 = q)
    }

    ### DO THE JOB

    res <- list()

    if (method=="dapc"){
      com <- t(combn(nPop(x), 2))
      pops <- seppop(x)
      nl2 <- ceiling(nl/nrow(com))

      for (i in 1:nrow(com)){
        dummy <- pops[c(com[i,1],com[i,2])]
        dummy <-do.call(rbind,dummy)
        # glPca (inside dapc) stops on loci with no calls in the pair
        called <- colSums(!is.na(as.matrix(dummy))) > 0
        dummy <- dummy[, called]
        dd <- dapc(dummy, n.pca=20, n.da=5)
        # rows of var.contr follow the locus order of dummy; their names do
        # not always equal locNames (e.g. "1" for locus "X1")
        ll <- order(dd$var.contr[,1], decreasing=TRUE)
        res[[i]] <- locNames(dummy)[ll[seq_len(min(nl2, length(ll)))]]

      }

      selloc <- unique(unlist(res))

    }

    if (method=="pahigh"){
      com <- t(combn(nPop(x), 2))
      pops <- seppop(x)
      nl2 <- ceiling(nl/(nrow(com)*2))

      # Private-allele loci of population a against b, with the rule of
      # gl.report.pa(): an allele present in a and absent from b. Computed
      # here because gl.report.pa() needs networkD3, tibble and tidyr (for
      # its plot) on every call.
      pa_loci <- function(a, b) {
        qa <- alf(pops[[a]])[, 2]
        qb <- alf(pops[[b]])[, 2]
        locNames(x)[which((qb == 0 & qa != 0) | (qb == 1 & qa != 1))]
      }
      panxx <- lapply(seq_len(nrow(com)), function(i) {
        list(pa1 = pa_loci(com[i, 1], com[i, 2]),
             pa2 = pa_loci(com[i, 2], com[i, 1]))
      })

      # highest-frequency private alleles of population a against b
      top_pa <- function(pas, a, b) {
        if (length(pas) == 0) return(character(0))
        p1p <- gl.keep.loc(pops[[a]], loc.list = pas, verbose = 0)
        p2p <- gl.keep.loc(pops[[b]], loc.list = pas, verbose = 0)
        score <- rowSums(alf(p1p) * !alf(p2p))
        names(score) <- locNames(p1p)
        utils::head(names(sort(score, decreasing = TRUE)), nl2)
      }

      res <- list()
      for (i in 1:nrow(com)){
        res1 <- top_pa(panxx[[i]]$pa1, com[i,1], com[i,2])
        res2 <- top_pa(panxx[[i]]$pa2, com[i,2], com[i,1])
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
      nl2 <- ceiling(nl/length(pops))

      for (i in 1:nPop(x)){
        dummy <- pops[[i]]
        cm <- colMeans(as.matrix(dummy), na.rm = TRUE)
        index <- !is.na(cm) & (cm == 0 | cm == 2)
        dl<-  locNames(dummy)[index]
        # a population can hold fewer monomorphic loci than its share
        mons <- dl[sample.int(length(dl), min(nl2, length(dl)))]
        res[[i]] <- mons
      }


      selloc <- unique(unlist(res))

    }
    if (method=="stratified"){
      res <- list()
      pops <- seppop(x)
      nl2 <- ceiling(nl/length(pops))
      for (i in 1:nPop(x)){

        dummy <- pops[[i]]

        df <- data.frame(id=locNames(dummy), freq=0.5-(abs(alf(dummy)[,1]-0.5)))
        df <- df[order(df$freq),]

        self <- seq(0,0.5, length=nl2+1)

        cf <- cut(df$freq, breaks=self, include.lowest=TRUE)

        scf <- split(df$id, cf)
        dres <- list()
        for (ii in 1:length(scf)){

          if (length(scf[[ii]]) < 1) next
          dres[[ii]] <- scf[[ii]][sample.int(length(scf[[ii]]), 1)]
        }
        res[[i]] <- unlist(dres)
      }


      selloc <- unique(unlist(res))
    }
    if (method=="hafall"){

      index <- order(0.5-abs(0.5-alf(x)[,1]), decreasing = TRUE)
      selloc <- locNames(x)[index[1:nl]]


    }
    if (method=="hafpop"){

      res <- list()
      pops <- seppop(x)
      nl2 <- ceiling(nl/length(pops))
      for (i in 1:nPop(x)){

        index <- order(0.5-abs(0.5-alf(pops[[i]])[,1]), decreasing = TRUE)
        res[[i]] <- locNames(pops[[i]])[index[seq_len(min(nl2, length(index)))]]
      }


      selloc <- unique(unlist(res))

    }

    if (method=="pic"){

      res <- list()
      af <- alf(x)
      p <- af[,1]
      q <- af[,2]
      pic <- order(1-(p^2+q^2)-2*p^2*q^2, decreasing = TRUE)
      index <- pic[1:nl]
      selloc<- locNames(x)[index]
    }


    if (method=="picdart"){

      res <- list()
      x <- gl.recalc.metrics(x,verbose=0)
      pic <- order(x@other$loc.metrics$AvgPIC, decreasing = TRUE)
      index <- pic[1:nl]
      selloc<- locNames(x)[index]
    }

    selloc <- selloc[!is.na(selloc)]

    # gl.keep.loc() returns the whole object when given no loci
    if (length(selloc) == 0 && !exact) {
      stop(error(paste0("Method '", method, "' selected no loci (none ",
                        "monomorphic in a population, or no private ",
                        "alleles). Use exact = TRUE to fill the panel with ",
                        "random loci, or another method.\n")))
    }

    if (exact) {  #add/remove random loci in case exact is wanted
      if (length(selloc) > nl)
      {
        selloc <- selloc[sample.int(length(selloc), nl)]
      }
      if (length(selloc)< nl)
      {
        if (verbose >= 1) {
          cat(warn(paste0("  Warning: method '", method, "' selected ",
                          length(selloc), " loci; adding ",
                          nl - length(selloc), " random loci to reach nl = ",
                          nl, ".\n")))
        }
        #add random loci
        lcs <- locNames(x)
        rest <- lcs[!lcs %in% selloc]
        selloc <- c(selloc, rest[sample.int(length(rest), nl - length(selloc))])
      }

    }



    #filter object to keep only selected loci
    xx <- gl.keep.loc(x, selloc, verbose = 0)

    if (verbose >= 3) {
      cat(report("  Method:", method, "; loci selected:", nLoc(xx), "of",
                 nLoc(x), "\n"))
    }

    # FLAG SCRIPT END

    if (verbose >= 1) {
      cat(report("Completed:", funname, "\n"))
    }

    return(xx)



  }
