% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.structure.run.r
\name{utils.structure.run}
\alias{utils.structure.run}
\title{Utility function to run Structure}
\usage{
utils.structure.run(
  g,
  k.range,
  num.k.rep,
  label,
  delete.files = TRUE,
  exec,
  burnin,
  numreps,
  noadmix,
  freqscorr,
  randomize,
  seed,
  pop.prior,
  locpriorinit,
  maxlocprior,
  gensback,
  migrprior,
  pfrompopflagonly,
  popflag,
  inferalpha,
  alpha,
  unifprioralpha,
  alphamax,
  alphapriora,
  alphapriorb
)
}
\arguments{
\item{g}{a gtypes object [see \code{strataG}].}

\item{k.range}{vector of values to for \code{maxpop} in multiple runs. If set
to \code{NULL}, a single STRUCTURE run is conducted with \code{maxpops}
groups. If specified, do not also specify \code{maxpops}.}

\item{num.k.rep}{number of replicates for each value in \code{k.range}.}

\item{label}{label to use for input and output files}

\item{delete.files}{logical. Delete all files when STRUCTURE is finished?}

\item{exec}{name of executable for STRUCTURE. Defaults to "structure".}

\item{burnin}{Number of burnin reps [default 10000].}

\item{numreps}{Number of MCMC replicates [default 1000].}

\item{noadmix}{Logical. No admixture? [default TRUE].}

\item{freqscorr}{Logical. Correlated frequencies? [default FALSE].}

\item{randomize}{Randomize [default TRUE].}

\item{seed}{Set random seed [default 0].}

\item{pop.prior}{A character specifying which population prior model to use:
"locprior" or "usepopinfo" [default NULL].}

\item{locpriorinit}{Parameterizes locprior parameter r - how informative the
populations are. Only used when pop.prior = "locprior" [default 1].}

\item{maxlocprior}{Specifies range of locprior parameter r. Only used when 
pop.prior = "locprior" [default 20].}

\item{gensback}{Integer defining the number of generations back to test for 
immigrant ancestry. Only used when pop.prior = "usepopinfo" [default 2].}

\item{migrprior}{Numeric between 0 and 1 listing migration prior. Only used 
when pop.prior = "usepopinfo" [default 0.05].}

\item{pfrompopflagonly}{Logical. update allele frequencies from individuals 
specified by popflag. Only used when pop.prior = "usepopinfo" [default TRUE].}

\item{popflag}{A vector of integers (0, 1) or logicals identifiying whether 
or not to use strata information. Only used when pop.prior = "usepopinfo"
 [default NULL].}

\item{inferalpha}{Logical. Infer the value of the model parameter # from the 
data; otherwise is fixed at the value alpha which is chosen by the user. 
This option is ignored under the NOADMIX model. Small alpha implies that 
most individuals are essentially from one population or another, while 
alpha > 1 implies that most individuals are admixed [default FALSE].}

\item{alpha}{Dirichlet parameter for degree of admixture. This is the 
initial value if inferalpha = TRUE [default 1].}

\item{unifprioralpha}{Logical. Assume a uniform prior for alpha which runs 
between 0 and alphamax. This model seems to work fine; the alternative model
 (when unfprioralpha = 0) is to take alpha as having a Gamma prior, with 
 mean alphapriora × alphapriorb, and variance alphapriora × alphapriorb^2 
 [default TRUE].}

\item{alphamax}{Maximum for uniform prior on alpha when 
unifprioralpha = TRUE [default 20].}

\item{alphapriora}{Parameters of Gamma prior on alpha when 
unifprioralpha = FALSE [default 0.05].}

\item{alphapriorb}{Parameters of Gamma prior on alpha when 
unifprioralpha = FALSE [default 0.001].}
}
\value{
\describe{ \item{\code{structureRun}}{a list where each element is a
list with results from \code{structureRead} and a vector of the filenames
used} \item{\code{structureWrite}}{a vector of the filenames used by
STRUCTURE} \item{\code{structureRead}}{a list containing: \describe{
\item{\code{summary}}{new locus name, which is a combination of loci in
group} \item{\code{q.mat}}{data.frame of assignment probabilities for each
id} \item{\code{prior.anc}}{list of prior ancestry estimates for each
individual where population priors were used} \item{\code{files}}{vector of
input and output files used by STRUCTURE} \item{\code{label}}{label for the
run} } } }
}
\description{
These functions were copied from package strataG, which is no longer on CRAN
 (maintained by Eric Archer)
}
\author{
Bernd Gruber (bugs? Post to
 \url{https://groups.google.com/d/forum/dartr}); original implementation of
 Eric Archer \url{https://github.com/EricArcher/strataG}
}
