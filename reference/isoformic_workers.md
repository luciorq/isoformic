# Start or stop the enrichment workers

`isoformic_workers()` starts `n` `mirai` daemons ahead for the
enrichment functions
([`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md)
and the other functions of that topic), stopping those of the
`"isoformic"` compute profile that were running. Every enrichment call
then uses them as they are (replacing those that died), whatever its
`workers`, as long as `workers` is above 1. The daemons get the
session's library paths and one BLAS and OpenMP thread each, and
starting them leaves the session's random number stream where it was.
`isoformic_workers(0)` stops the daemons of the `"isoformic"` profile,
whoever started them; they are also stopped when the session ends or
isoformic is unloaded.

The section *Workers* of
[`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md)
describes how the calls use the daemons, their memory and the
reproducibility of their p-values.

## Usage

``` r
isoformic_workers(n)
```

## Arguments

- n:

  Number of daemons to start: a whole number of at least 1, `"auto"` for
  [`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html),
  or 0 to stop the daemons.

## Value

The number of daemons now running, invisibly.

## See also

[`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md);
the `isoformic.workers` option in
[isoformic-options](https://luciorq.github.io/isoformic/reference/isoformic-options.md).

## Examples

``` r
if (FALSE) { # \dontrun{
isoformic_workers(4)
# ... calls of run_enrichment() with workers = 4 ...
isoformic_workers(0)
} # }
```
