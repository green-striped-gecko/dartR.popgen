# dartR.popgen (development version)

## Bug fixes

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
