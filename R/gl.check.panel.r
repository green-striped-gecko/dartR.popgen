#' Check a snp panel for a specified parameter
#'
#' @description
#' This function checks a panel how good it is to recreate the specified parameter of conservation concern (Ne, Fst, Ho etc.)
#'
#' @param x A 'dartR or genlight' object containing the SNP panel genomic data.
#' @param xorig A 'dartR or genlight' object containing the original genomic data for comparison.
#' @param parameter A character string specifying the parameter to check. Options include: Fst, He, Ho, Ne, Nall, Fis.
#' @param neest.path  Path to neestimator (see gl.LDNe) 

#' @param plot.out Logical. If `TRUE`, generates plots summarizing selected loci.
#' @param plot.file A character string specifying the file name for saving plots. If `NULL`, plots are not saved.
#' @param plot.dir A character string specifying the directory to save plots. Defaults to the working directory.
#' @param verbose Integer level of verbosity for reporting progress and information.
#'
#' @details
#' The function applies various methods to select loci based on the input 'dartR or genlight' object.
#' Each method has specific criteria for selecting loci:

#'
#' @return A plot and the result of the linear regression
#'
#' @examples
#' # Example usage:
#'
#' # Select 20 loci randomly
#' selected <- gl.select.panel(possums.gl, method = "random", nl = 50)
#' gl.check.panel(selected, possums.gl, parameter="Fst")
#'
#' @export
#' @importFrom ggpmisc stat_poly_eq use_label

gl.check.panel <- function(x, xorig, parameter="Fst", neest.path = NULL, 
                           plot.out = TRUE, plot.file = NULL, plot.dir = NULL, 
                           verbose = NULL) {

  
#order populatios after pop to make sure they have the same order


    
x<- x[order(pop(x)),]
xorig <- xorig[order(pop(xorig)),]

if (!sum(pop(x) == pop(xorig))==nInd(x)) {
  stop("The populations in x and xorig do not match. Please ensure that the panel is a subset of the original data.")
}

gg<- NULL

fis_orig <- fis_panel <- fst_orig <- fst_panel <- 
  hetso_orig <- hetso_panel <- het_orig <- het_panel <- 
  nall_orig <- nall_panel <- nes_orig <- nes_panel <- NULL

#check fsts

if (parameter == "Fst") {
  

fst_x <- as.numeric(gl.fst.pop(xorig, verbose = 0, nboots=1))

fst_test <- as.numeric(gl.fst.pop(x, verbose = 0, nboots = 1))

fsts <- data.frame(fst_orig=fst_x, fst_panel=fst_test)
fsts <- fsts[complete.cases(fsts),]

gg<- ggplot(fsts, aes(x=fst_orig, y=fst_panel)) + geom_point() + geom_smooth(aes(x=fst_orig, y=fst_panel), method="lm") +stat_poly_eq(use_label(c("eq", "R2")))
res <- fsts
}
#check expected heterozygosity

if (parameter=="He") {
het_x <- gl.report.heterozygosity(xorig, verbose=0)
het_test <- gl.report.heterozygosity(x, verbose=0)
hets <- data.frame(het_orig=het_x$He, het_panel=het_test$He)
gg<- ggplot(hets, aes(x=het_orig, y=het_panel)) + geom_point() + geom_smooth(aes(x=het_orig, y=het_panel), method="lm") +stat_poly_eq(use_label(c("eq", "R2")))
res <- hets
}

#check number of alleles
if (parameter=="Na") {
nall_x <- gl.report.allelerich(xorig, verbose=0,plot.display = F)
nall_test <-  gl.report.allelerich(x, verbose=0,plot.display = F)
# check pcoas
nalls <- data.frame(nall_orig=nall_x$`Allelic Richness per population`$`mean_corrected_richness`, nall_panel= nall_test$`Allelic Richness per population`$`mean_corrected_richness`)

gg<- ggplot(nalls, aes(x=nall_orig, y=nall_panel)) + geom_point() + geom_smooth(aes(x=nall_orig, y=nall_panel), method="lm") +stat_poly_eq(use_label(c("eq", "R2")))
res <- nalls
}


#check FIS
if (parameter=="Fis") {


fis_x <- gl.report.heterozygosity(xorig, verbose=0)
fis_test <- gl.report.heterozygosity(x, verbose=0)

fiss <- data.frame(fis_orig=fis_x$FIS, fis_panel= fis_test$FIS)

gg<- ggplot(fiss, aes(x=fis_orig, y=fis_panel)) + geom_point() + geom_smooth(aes(x=fis_orig, y=fis_panel), method="lm") + stat_poly_eq(use_label(c("eq", "R2")))
res <- fiss
}


#check observed heterozygosity
if (parameter=="Ho") {
hetso_x <- gl.report.heterozygosity(xorig, verbose=0)
hetso_test <- gl.report.heterozygosity(x, verbose=0)

hetso <- data.frame(hetso_orig=hetso_x$Ho, hetso_panel =hetso_test$Ho)

gg<- ggplot(hetso, aes(x=hetso_orig, y=hetso_panel)) + geom_point() + geom_smooth(aes(x=hetso_orig, y=hetso_panel), method="lm") +stat_poly_eq(use_label(c("eq", "R2")))
res <- hetso
}

#check Ne
if (parameter=="Ne") {
ne_x <- gl.LDNe(xorig, neest.path = neest.path, critical = 0.05, verbose=0,mating = "random")

ne_test <- gl.LDNe(x, neest.path = neest.path, critical = c(0.05), verbose=0,mating = "random", singleton.rm = F) 

nes_x <- as.numeric(unlist(lapply(ne_x, function(x) x$`Frequency 1`[6])))

nes_test <- as.numeric(unlist(lapply(ne_test, function(x) x$`Frequency 1`[6])))

nes<- data.frame(nes_orig=nes_x,  nes_panel = nes_test)

nes[sapply(nes, function(x) x=="Inf")] <- NA  # replace Inf with NA]


gg<- ggplot(nes, aes(x=nes_orig, y=nes_panel)) + geom_point() + geom_smooth(aes(x=nes_orig, y=nes_panel), method="lm")+stat_poly_eq(use_label(c("eq", "R2")))


res <- nes
}
print(gg)
return(res)
}



