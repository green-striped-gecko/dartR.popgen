% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.select.panel.R
\name{gl.select.panel}
\alias{gl.select.panel}
\title{Select Loci Panel Based on Various Methods}
\usage{
gl.select.panel(
  x,
  method = "random",
  nl = 10,
  plot.out = TRUE,
  plot.file = NULL,
  plot.dir = NULL,
  verbose = NULL
)
}
\arguments{
\item{x}{A 'dartR or genlight' object containing the genomic data.}

\item{method}{A character string specifying the selection method. Options include:
\itemize{
  \item `"dapc"`: Select loci contributing most to discrimination between populations using DAPC (Discriminant Analysis of Principal Components).
  \item `"pahigh"`: Select loci with private alleles having high frequency.
  \item `"random"`: Randomly select loci.
  \item `"monopop"`: Select monomorphic loci within populations.
  \item `"stratified"`: Stratified sampling of loci based on allele frequencies.
  \item `"hafall"`: Select loci with the highest allele frequencies across all populations.
  \item `"hafpop"`: Select loci with the highest allele frequencies within each population.
}}

\item{nl}{An integer specifying the number of loci to select.}

\item{plot.out}{Logical. If `TRUE`, generates plots summarizing selected loci.}

\item{plot.file}{A character string specifying the file name for saving plots. If `NULL`, plots are not saved.}

\item{plot.dir}{A character string specifying the directory to save plots. Defaults to the working directory.}

\item{verbose}{Integer level of verbosity for reporting progress and information.}
}
\value{
A 'dartR or genlight' object containing the selected loci.
}
\description{
This function selects a panel of loci from a genomic dataset ('dartR or genlight' object)
based on various selection methods.
}
\details{
The function applies various methods to select loci based on the input 'dartR or genlight' object.
Each method has specific criteria for selecting loci:
\itemize{
  \item `dapc`: Performs DAPC and identifies loci with the highest contributions to discrimination between population pairs.
  \item `pahigh`: Identifies loci with private alleles that have high frequency differences between populations.
  \item `random`: Selects loci randomly.
  \item `monopop`: Selects loci that are monomorphic within populations.
  \item `stratified`: Uses stratified sampling to select loci based on allele frequencies.
  \item `hafall`: Selects loci with the highest allele frequencies across the dataset.
  \item `hafpop`: Selects loci with the highest allele frequencies within individual populations.
}
}
\examples{
# Example usage:

# Select 20 loci randomly
selected <- gl.select.panel(possums.gl, method = "random", nl = 50)

# Select loci based on DAPC
selected <- gl.select.panel(possums.gl, method = "dapc", nl = 5)

}
