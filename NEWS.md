# dartR.popgen (development version)

## Bug fixes

* `gl.collapse` builds population groups as connected components. The old
  single pass could leave a population in two groups when similar
  populations formed a chain, and `gl.merge.pop()` then stopped ("not
  present in the dataset"; 24 of 363 random test matrices). When all
  populations fall into one group, the one-population result is returned
  instead of an error, and when none merge the input `fd` is returned
  unchanged. `tloc` must now match the value used in `gl.fixed.diff()`:
  previously the returned matrices were recomputed with `gl.collapse`'s
  `tloc` (default 0) while the grouping used the input's, so the help
  example mixed 0.05 and 0. Arguments are checked; NA distances do not join
  populations. **Calls whose `tloc` differs from the one used to build `fd`
  now stop with an error; runtime about doubles (one extra
  `gl.fixed.diff()` run for the check).**
* `gl.ld.distance` no longer fails when `ld.resolution` exceeds the
  largest distance ("subscript out of bounds" inside `fields::stats.bin`);
  bins are now computed in base R, with the same means, so `fields` is no
  longer needed. No empty bin is added when the bins end exactly on the
  largest distance. The returned table gains `n.pairs`, the number of
  pairs behind each bin mean. The table prints only at `verbose >= 3`;
  `pop.colors` accepts a palette function; the threshold line is labelled
  in the legend and documented as an R.squared threshold. **The returned
  table gains a column.**
* `gl.outflank` analyses each SNP once, from the 0/1/2 genotype matrix, and
  now reproduces the OutFLANK package exactly. Previously the genlight was
  converted to genind and both allele columns of every SNP were analysed,
  so `numberHighFstOutliers` and `numberLowFstOutliers` were doubled (26
  reported for 13 flagged loci in a simulation) and `dfInferred` and the
  q-values departed from OutFLANK; locus names containing a dot failed.
  Loci missing in every individual stay in the output as NA, so `index` and
  `results` have one row per input locus in input order (previously they
  were dropped and later loci shifted). New `verbose` argument; `...`
  removed (it was ignored); clear errors for fewer than two populations,
  SilicoDArT input or a missing qvalue package; the plot uses `Hmin`; the
  help page states that `index` is TRUE for loci that are not outliers.
  About 6 times faster. **Outlier counts, `dfInferred`, q-values and
  `meanAlleleFreq` (now the reference-allele frequency) change; calls that
  pass extra arguments now error.**
* `gl.find.loci.in.genes` selects a gene when its own row or any descendant
  feature (mRNA, lnc_RNA, CDS, exon, followed through `Parent`) matches
  `gene`, and treats pseudogenes as genes. Previously only gene rows and
  CDS rows were searched, so genes whose matching text sits on a transcript
  row were missed (platypus NCBI annotation: 33 of 38 loci found for
  "receptor", 6 of 50 for "uncharacterized"). A missing `gff.file` now
  errors; previously the function silently used any object called `gff` in
  the workspace. `verbose` is honoured (`verbose = 0` is silent), loci
  without a position and sequence names absent from the GFF are reported
  at `verbose >= 1`, SilicoDArT data are accepted, and `save2tmp` saves
  the locus-gene table. **The loci returned change for genes matched
  through transcript rows and for pseudogenes.**
* `gl.TajimasD`: `sim_pval` now tests D against a neutral null of unlinked
  SNPs, simulated in R (derived counts drawn with P(k) proportional to 1/k
  from 2N sequences); `rep` alone runs it, and ms/sample_stats are no longer
  needed (`ms.path` is accepted and ignored). Previously ms simulated all
  sites on one non-recombining locus from N (not 2N) sequences, which for
  unlinked DArT SNPs gives far too wide a null. The documentation now says
  that `Pval.normal` and `Pval.beta` assume one non-recombining locus and
  are conservative for unlinked SNPs, and that ascertainment of polymorphic
  loci shifts D upwards. D is computed from exact allele counts with
  per-site sample sizes in Watterson's term. **`sim_pval` changes; D changes
  in the 5th decimal on complete data and by a few percent with missing
  calls; SilicoDArT input now errors; `simulation.out` holds the simulated
  D values.** `utils.get.allele.freq` works on genlight objects without
  dartR flags and honours `verbose`; `verbose = 0` is silent.
* `gl.sfs` builds the spectrum only from loci scored in every individual and
  reports how many loci with missing calls were excluded (verbose >= 1).
  Previously such loci were counted against the full sample size, so they
  fell into lower classes and a locus fixed in all called individuals
  appeared polymorphic. A folded multi-population spectrum is now folded on
  the minor allele of all populations combined (fastsimcoal2 joint MAF
  spectrum); previously each population was folded on its own. On a
  multi-population spectrum, `minbinsize` sets to zero the cells whose total
  count is below it and keeps the full array; previously it removed the
  first classes of every population, dropping polymorphisms private to one
  population. **Spectra change for data with missing calls (and so do
  `gl.run.stairway2` and `gl.run.epos` results built on them), for folded
  multi-population spectra, and for multi-population spectra with
  `minbinsize > 0`. SilicoDArT input now errors.** New `plot.theme`
  argument; `plot.file` works with `plot.out = FALSE`; bars sit at their
  class numbers.
* `gl.run.epos` sends EPOS the SFS with its true class numbers. Previously
  `folded = FALSE` was run as a folded SFS of twice the sample size (no
  `-U`), `minbinsize = 2` relabelled doubletons as singletons, and
  `minbinsize = 0` (zero class instead of `L`) always failed. `upper` and
  `lower` are now passed to epos2plot. **Estimates change for
  `folded = FALSE`, `minbinsize >= 2` and non-default `upper`/`lower`;
  default calls are unchanged. SilicoDArT input now errors.** New `seed`
  argument (passed to epos and bootSfs); `u` can be left NULL as documented;
  `L` is required unless `minbinsize = 0`; an unnamed `sfs` vector is
  accepted; missing binaries and EPOS failures stop with a clear message;
  `verbose = 0` is silent.
* `gl.run.stairway2`: when `L` is not given and the SFS is computed from
  `x`, the default `L` counts only the loci without missing calls (the loci
  `gl.sfs` uses) x 69, instead of all loci x 69. **Results change for data
  with missing calls when `L` is left at its default.**
* `gl.run.stairway2` runs Stairway Plot 2 in a new subfolder of `tempdir()`
  and `cleanup = TRUE` removes only that subfolder. Previously it deleted
  the whole session `tempdir()`, including binaries downloaded there by
  `gl.download.binary()` and files saved by other functions, and a
  `plot.file` save then failed after the run. `run = FALSE` now returns
  instead of stopping with "object 'res' not found". **The returned list
  gains `run.dir`, the run folder (NULL once removed); the history column
  `" low75"` is renamed `low75`; SilicoDArT input now errors at every
  verbosity.** The plot uses `plot.theme`; Stairway Plot 2 output prints
  only at `verbose >= 3`; the user's `future` plan is restored after
  `parallel > 1`; a missing binary folder, missing Java or a failed run
  stop with a clear message. Estimates are unchanged.
* `gl.run.popcluster`: the likelihood table (`best_run`) is numeric (NA
  where PopCluster writes "-"). It was kept as text, so the LogL(K) plot
  had an alphabetically sorted axis and was drawn upside down (the best K
  lowest), and the DLK plots were scrambled. **`best_run` columns become
  numeric and the plots change.** `output.path` defaults to `tempdir()`
  (input files were written to the working directory), `cleanup` works,
  `plot_theme` is applied, the panels are combined with patchwork instead
  of gridExtra, PopCluster's output is shown only at `verbose >= 3`, and
  arguments and the executable are checked before running.
* `gl.map.popcluster` now draws with `gl.map.snmf` (and `gl.map.structure`),
  so each population's bars are drawn at its own centre, `movepops` is
  added by name and extra populations in `x` are ignored. **Maps and
  default colours change as for `gl.map.snmf`.**
* `gl.plot.popcluster`: `verbose` follows `gl.set.verbosity()`, `plot.K`
  must be one K of the run, a palette function is accepted, and `aes()`
  uses `.data` (removing the R CMD check NOTE about global variables).

* `gl.map.snmf` now draws with `gl.map.structure`, so it gets the fixes
  made there: each population's bars are drawn at its own centre (with
  population levels not in alphabetical order they were drawn at another
  population's centre), `movepops` is added to `lon`/`lat` by name, and
  populations of `x` that are not in `qmat` are ignored instead of causing
  an error. **Maps change in those cases, and default colours change from
  `rainbow()` to the dartR palette.** New `plot.out` and `verbose`.
* `gl.run.snmf`: `cleanup` works (its code came after `return()`, so the
  LEA run files were never removed); with `cleanup = TRUE` (default)
  `best_run` holds run names such as `"K2/run1"` instead of paths.
  **`best_run` changes with the default.** `plot.out` is respected (the
  cross-entropy plot is still returned), `plot.file` is saved in
  `plot.dir` (default `tempdir()`) instead of the working directory, LEA's
  output is shown only at `verbose >= 3`, and `minK`, `maxK`, `rep` are
  checked.
* `gl.plot.snmf`: the dendrogram (`den = TRUE`) clusters the distance
  matrix itself instead of Euclidean distances between its rows. **The
  individual order changes with `den = TRUE`.** `verbose` follows
  `gl.set.verbosity()` and is silent at 0, `plot.K` must be one K of the
  run (clear error), a palette function is accepted in `color.clusters`,
  and `aes_()` is replaced.
* `gl.run.faststructure`: results are stored by K, so any `k.range` works.
  Previously `k.range` had to start at 2 without gaps: `3:4` or `c(2, 4)`
  stopped with "subscript out of bounds" after all runs had finished, and
  `1:2` dropped K = 1 and returned an empty element. **Calls with other
  `k.range` values now return results.**
* `gl.run.faststructure`: each call writes to a new subfolder of `output`,
  whose default is now `tempdir()` instead of the working directory, and
  reads only that subfolder; previously files of earlier runs in the same
  folder were read back. With `seed`, replicate r uses `seed + r - 1`;
  previously every replicate was identical. **Seeded results change after
  the first replicate.** The executables are checked before running, a
  failed run stops with the K and replicate, progress and program output
  follow `verbose`, and the likelihood plot uses `theme_dartR()` with new
  `plot.out`, `plot.theme`, `plot.dir`, `plot.file` arguments. `gsubfn`
  is no longer needed.
* `gl.plot.faststructure` now draws with `gl.plot.structure`, so it gets
  the fixes made there: each mode is the average of its replicates (it was
  one replicate), K = 1 after another K is not duplicated, the dendrogram
  clusters the distance matrix itself. **Returned q-matrices change when a
  K has several modes, and are data.tables.** `k.range = NULL` plots every
  K; new `dis.mat`, `plot.out`, `plot.dir`, `plot.file`, `verbose`.
  `gl.plot.structure` gains `label.size` (default 12, the previous fixed
  size).

* `gl.read.structure`: individuals with a population prior (USEPOPINFO)
  are read correctly. Previously the ancestry blocks after `|` were not
  parsed, so these individuals got equal membership in every cluster (e.g.
  0.5 / 0.5 instead of 0.042 / 0.958) and an all-NA `prior.anc`. **q
  values and `prior.anc` change for every individual with a prior.**
* `gl.read.structure`: with `x`, individuals are matched by name or, for
  files written by `gl.run.structure`, by position, and their names are
  restored. Previously such files got `orig.pop = NA` for every
  individual without a message; ids that match neither way now stop with
  an error.
* `gl.read.structure`: only STRUCTURE output files are read, so a folder
  kept by `gl.run.structure(delete.files = FALSE)` no longer fails on its
  log and params files. Runs are named `k<K>.r<replicate>` (as in
  `gl.run.structure`), with replicates in file-name number order (rep2
  before rep10); `prefix` is added only when set. **Run names change.**
  `rename_files = TRUE` stops instead of overwriting existing files. Help
  page completed.
* `gl.read.structure` and `gl.run.structure` now read STRUCTURE output
  files with one shared internal reader, `utils.structure.read`, instead of
  two separate copies of the parser (the `usepopinfo` bug above was a
  difference between them). `gl.run.structure` output is unchanged.

* `gl.evanno` / `utils.structure.evanno`: delta K is `NA`, with a warning,
  where all replicates of a K report the same LnP(K). Previously it was
  `Inf` (sd = 0), as on a real STRUCTURE run of three testset.gl
  populations. **`delta.k` changes from `Inf` to `NA` in that case.**
* `gl.evanno` / `utils.structure.evanno`: LnP'(K), |LnP''(K)| and delta K
  are computed only from K - 1 and K + 1. Previously a `k.range` with gaps
  was differenced as if consecutive (delta K at K = 3 from K = 1 and
  K = 5). **These values become `NA` where they were computed across a
  gap.**
* `gl.evanno`: new arguments `plot.theme` (default `theme_dartR()`),
  `plot.dir`, `plot.file` and `verbose`; clear errors for a wrong input or
  fewer than three K; a warning when delta K needs more replicates. The
  combined figure is built with patchwork and returned as
  `plots$combined` (`gridExtra` is no longer used); `aes_string()`
  replaced. Help pages rewritten, including the delta K method.
* `gl.run.structure`: the Evanno panel passes `verbose` to
  `utils.structure.evanno`, so its warnings (delta K undefined because
  replicates agree, K values not consecutive, one replicate per K) are
  shown at `verbose >= 1` when the plot is drawn or saved.
* `gl.map.structure`: each population's bars are drawn at its own centre.
  Previously the centres were matched to the q-matrix in the wrong
  direction, so when the population levels of `x` were not in alphabetical
  order, bars were drawn at another population's centre under that
  population's label. **Maps change for any `x` whose population levels are
  not alphabetical.**
* `gl.map.structure`: `movepops` is added to the `lon` and `lat` columns by
  name (previously by position, which moved latitude on genlights storing
  `lat` first), and its rows can be matched by population name.
  **`movepops` shifts change for genlights whose `latlon` stores `lat`
  first.**
* `gl.map.structure`: populations of `x` that are not in `qmat` are ignored
  (previously an error after many leaflet warnings); a `qmat` population
  without coordinates, a `latlon` without `lon`/`lat` columns and malformed
  `qmat` or `K` give informative errors. One population and K = 1 work.
  `K` accepts a mode label such as `"2.2"`. New arguments `plot.colors`
  (default: the `gl.plot.structure` palette instead of `rainbow()`, **so
  default map colours change**), `plot.out` and `verbose`; `leaflet` is
  guarded; the help example is corrected.
* `gl.map.structure`: individuals whose `orig.pop` is NA are dropped with
  a warning at `verbose >= 1`. Previously they added empty rows to the
  returned tables and many leaflet warnings, and a q-matrix with no
  populations at all stopped with "attempt to select less than one
  element"; it now stops with a message pointing to `x` in
  `gl.read.structure`.
* `gl.plot.structure`: when Clumpak finds more than one mode at a K, each
  mode is now the average of its replicates, as documented. Previously every
  mode showed only its first replicate. **Returned q-matrices and bar
  heights change whenever a K has more than one mode.**
* `gl.plot.structure`: requesting K = 1 after another K (e.g. `K = c(2, 1)`)
  no longer duplicates the earlier K panels, which were relabelled as extra
  modes (`2.1`, `2.2`). **The returned list and the plot lose the duplicate
  panels.**
* `gl.plot.structure`: with `den = TRUE` the returned tables keep the
  population names (only the plot hides them), so they can be passed to
  `gl.map.structure()`. The dendrogram now clusters `dis.mat` itself;
  previously it clustered Euclidean distances between the rows of the
  distance matrix. **Individual order in `den = TRUE` plots changes.**
* `gl.plot.structure`: arguments are checked before any work: `den = TRUE`
  needs `x` or `dis.mat` (and `dis.mat` alone now works), `met_clumpp`,
  missing K values and too few colours give informative errors;
  `color_clusters` accepts a palette function as documented; `proxy` and
  `reshape2` are guarded; `aes_()` is replaced by `aes()`, removing the
  ggplot2 deprecation warning; help page corrected (`border_ind` default,
  `plot.dir`, return value, `k_name`).

* `gl.run.structure`: the STRUCTURE runs are now always returned. Previously
  the Evanno step ran unconditionally and stopped with "must have at least
  two values of k" whenever `k.range` had fewer than three values, so the
  finished runs were lost even with `plot.out = FALSE`. **Calls with one or
  two K values now succeed.**
* `gl.run.structure`: individuals are passed to STRUCTURE by index and the
  names in `indNames(x)` are restored in `q.mat` and `prior.anc`. Previously
  names longer than 11 characters were truncated by STRUCTURE and collapsed
  to one id without an error, and names with spaces produced an NA popflag
  and a STRUCTURE failure. **`q.mat` rows are now in the order of
  `indNames(x)` instead of alphabetical by id, and runs are named
  `k<K>.r<replicate>` instead of carrying a timestamp label.** An unnamed
  `popflag` is matched to individuals in `indNames(x)` order.
* `gl.run.structure`: STRUCTURE runs in a temporary folder that is removed
  whether or not the run succeeds. Previously a time-stamped folder was
  created in the working directory, left behind on failure, and any folder
  of the same name was deleted first. **With `delete.files = FALSE` the
  files are now kept in a time-stamped folder under `plot.dir`.** The
  command line is quoted, so an executable path with a space works;
  `plot.dir = NULL` by default so `gl.set.wd()` is honoured; STRUCTURE's own
  output is shown only at `verbose >= 3`, with one progress line per run at
  `verbose >= 2`; the Evanno layout shows mean LnP(K) and LnP'(K) instead of
  the first panel twice; `tidyr` is guarded and the superseded
  `gather`/`spread` calls replaced.
* `gl.ld.haplotype`: haplotypes are now identified from the pairwise LD of
  adjacent SNPs in the LD matrix. Previously they were read off a rotated
  copy of the matrix with a misaligned start index, so one adjacent pair in
  LD could be reported as a block spanning the whole chromosome, and start
  and end positions were rounded to three significant digits. **The
  haplotype table changes for every call with `haplo_id = TRUE`.**
* `gl.ld.haplotype`: `plot.save` is honoured. Previously a PDF was written
  to `plot.dir` on every call when no haplotypes were drawn, and never when
  they were. **Default calls no longer write files.**
* `gl.ld.haplotype`: `ind.limit` is a true minimum; a population with
  exactly `ind.limit` individuals is analysed instead of skipped.
* `gl.ld.haplotype`: chromosomes with fewer than four SNPs after filtering
  are skipped with a warning instead of aborting the run; unknown
  `chrom_name` / `pop_name` values, missing `@chromosome` / `@position`
  slots and SilicoDArT input now stop with a clear message; the
  heterozygosity and SNP-position tracks, title and axis of each plot use
  only that chromosome's SNPs; heatmap cells with LD of 0 or negative LD are
  drawn instead of left as holes; PLINK intermediates are written and read
  from the same temporary path (the function failed after `gl.set.wd()`);
  missing packages stop with an error instead of returning `-1`; output is
  silent at `verbose = 0`.

* `gl.find.genes.for.loci`: loci that do not overlap a gene were returned as
  a single row with an empty locus name, no gene and no distance, because
  the nearest-gene join did not carry the locus identifiers. Each such locus
  now gets its own row with the nearest gene, the distance to its closest
  edge, and the side. **Results change for every locus outside a gene.**
* `gl.nhybrids`: the results table `aa-PofZ.csv` no longer mislabels the
  posterior-probability columns. Previously the NewHybrids `IndivName`
  column was kept, so `P0` held the string "NoName" and every probability
  sat one column to the right of its label. **Values read from the labelled
  columns of `aa-PofZ.csv` change for all users.**
* `gl.nhybrids`: `nhyb.directory = NULL` (write the input file and stop)
  no longer errors on macOS/Linux.
* `gl.nhybrids`: a NewHybrids run that produces no output now stops with a
  clear error instead of failing cryptically or silently reusing a previous
  run's results; paths of 100+ characters, which crash the NewHybrids
  binary, are rejected with an explanation; `aa-Pi.hist` is now delivered
  to `outpath`; the caller's working directory is restored on failure; and
  `loc.metrics` in the returned object now track the selected loci
  positionally.

* `gl.blast`: a run whose external step failed (unknown `task`, missing
  query fasta or reference genome, makeblastdb or blastn error) no longer
  returns the previous run's hits as its result. Inputs are validated up
  front, temporary files from earlier runs are removed, and a non-zero exit
  of makeblastdb or blastn stops with the tool's message. A reference genome
  path containing spaces now works (the genome is passed to makeblastdb on
  stdin), and a BLAST install under a path with spaces, such as Program
  Files on Windows, is no longer refused; an R temporary directory with
  spaces is refused with an explanation, because BLAST cannot open a
  database there. **A fasta-file query with no surviving hit now returns an
  empty data frame with the BLAST columns instead of the input path.**

* `gl.LDNe`: `mating = "monogamy"`, the documented value, errored with
  "missing value where TRUE/FALSE needed"; only the undocumented
  abbreviation `"mono"` ran. Both are now accepted.
* `gl.LDNe`: when the first population had fewer than three individuals it
  was given the *next* population's jackknife confidence limits, and that
  population received NA. NeEstimator prints no jackknife line for such a
  population, and the placeholder was inserted one position too late.
  **Jackknife CI values change for datasets whose first population has
  fewer than three individuals.**
* `gl.LDNe`: an unrecognised `Waples.correction` no longer applies the
  genome-length formula silently (it produced negative Ne); an invalid
  `pairing` no longer fails on an internal object; non-SNP data, a missing
  or non-scalar `Waples.correction.value`, and an unsupported operating
  system now stop with an explanation. **Callers passing the string
  `"NULL"` as `Waples.correction` must pass `NULL` instead.**
* `gl.LDNe`: `plot.file` without `plot.out` no longer errors with "object
  'p3' not found" after NeEstimator has run; `plot_colors_pop` accepts a
  palette function or a vector longer than the number of populations, and
  reports too few colours as such; `naive = TRUE` keeps the population
  names on the returned list; the output file in `outpath` is refreshed by
  a second run instead of silently keeping the first run's results; a
  relative `outpath` such as `"."` now means the caller's working
  directory, as documented.
* `gl.LDNe`: each call runs NeEstimator in its own directory under
  `tempdir()` and removes it afterwards, so concurrent calls (forked
  parallel runs share the parent's `tempdir()`) no longer overwrite each
  other's input and output files, and a failed run stops instead of reading
  a file left by an earlier one. NeEstimator's own byproduct files no
  longer persist in `tempdir()` after the call.

## Improvements

* `gl.find.genes.for.loci`: two new output columns, `gene_strand` and
  `relative_position` ("inside", "upstream", "downstream"); `gene_product`
  is taken from the mRNA, CDS or exon children when the gene row has none
  (the NCBI and Ensembl GFF3 layout); `loci` is optional and defaults to
  every locus with a position; a locus at exactly equal distance from two
  genes now follows the documented tie-break (closest midpoint, shorter
  gene, gene_id); loci without a position and sequence names absent from
  the GFF are reported at `verbose >= 1`; a locus on a sequence without
  genes is returned with NA gene columns; the roxygen example runs.
* `gl.blast`: console output honours `verbose` (makeblastdb output only at
  `verbose >= 3`, the "k of N sequences aligned" count at `verbose >= 1`,
  history entry and "Completed" on every exit); the note pointing to
  `gl.list.reports()`/`gl.print.reports()`, which do not exist, is replaced
  by the paths of the three saved tables; documentation corrected.
* `gl.LDNe`: console output honours `verbose` (`gl2genepop` silent below
  3, NeEstimator's log only at `verbose >= 3`, the result tables printed at
  `verbose >= 2`); `@return` now describes the named list of data frames
  that the function returns; the `plot.dir`, `plot.file` and
  `plot_colors_pop` defaults, the `mating` values and the second example
  are documented as they behave.
