% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.get.allele.freq.r
\name{utils.get.allele.freq}
\alias{utils.get.allele.freq}
\title{utils.get.allele.freq}
\usage{
utils.get.allele.freq(x, verbose = 2)
}
\arguments{
\item{x}{Name of the genlight object containing the SNP or Tag P/A
(SilicoDArT) data [required].}

\item{verbose}{Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
progress log; 3, progress and results summary; 5, full report
[default 2 or as specified using gl.set.verbosity].}
}
\value{
A matrix with allele (SNP data) or presence/absence frequencies
(Tag P/A data) broken down by population and locus
}
\description{
Generates percentage allele frequencies by locus and population -- This is copy from package dartR function gl.percent.freq
}
\details{
This is a support script, to take SNP data or SilicoDArT presence/absence
data grouped into populations in a genlight object \{adegenet\} and generate
a table of allele frequencies for each population and locus
}
\examples{
m <-  utils.get.allele.freq(testset.gl)
}
\author{
Custodian: Arthur Georges (Post to
\url{https://groups.google.com/d/forum/dartr})
}
