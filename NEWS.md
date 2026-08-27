# dartR.popgen (development version)

## Bug fixes

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
