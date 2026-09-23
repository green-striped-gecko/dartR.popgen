# dartR.popgen (development version)

## Bug fixes

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
