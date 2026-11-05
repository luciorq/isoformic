# Run Terminus on Salmon output from R

`run_terminus()` builds the consensus groups of Terminus for a set of
samples without leaving R: it runs `terminus group` on every Salmon
output directory, then `terminus collapse` on all of them, and returns
the groups read by
[`read_terminus_groups()`](https://luciorq.github.io/isoformic/reference/read_terminus_groups.md),
ready for `add_uncertainty(groups = )`. Terminus is installed from
Bioconda into a conda environment managed by `condathis` the first time
it is needed.

## Usage

``` r
run_terminus(
  paths,
  output_dir,
  consensus_thresh = 0.5,
  min_spread = 0.1,
  threads = 4L,
  workers = 1L,
  keep_logs = FALSE,
  env_name = "terminus-env",
  verbose = "output"
)
```

## Arguments

- paths:

  Salmon output directories (or paths to their `quant.sf`).

- output_dir:

  Folder for the Terminus output; one subfolder per sample holds
  `groups.txt` and, after the consensus, `clusters.txt` and the
  collapsed quantification.

- consensus_thresh:

  Share of samples in which two transcripts must be grouped for the
  consensus to group them (`--consensus-thresh`).

- min_spread:

  Minimum spread of a transcript for its edges to be collapse candidates
  (`--min-spread`).

- threads:

  Threads for writing the collapsed quantifications.

- workers:

  Number of samples grouped at the same time (in parallel processes on
  Unix-like systems). Each one needs several gigabytes of memory for a
  human transcriptome.

- keep_logs:

  Logical, keep the logs of `terminus group`.

- env_name:

  Name of the conda environment with Terminus.

- verbose:

  Passed to
  [`condathis::run()`](https://luciorq.github.io/condathis/reference/run.html).

## Value

The consensus groups, as returned by
[`read_terminus_groups()`](https://luciorq.github.io/isoformic/reference/read_terminus_groups.md),
with the output folder in the `output_dir` attribute.

## Details

Every directory needs inferential replicates (`--numGibbsSamples` or
`--numBootstraps`) and equivalence classes with their weights
(`--dumpEq --dumpEqWeights`);
[`salmon_quant()`](https://luciorq.github.io/isoformic/reference/salmon_quant.md)
writes both. Without the weights Terminus stops with an unhelpful error,
so the directories are checked first. Terminus names the folder of each
sample after the basename of its Salmon directory, which must therefore
be unique.

`terminus group` writes a log of several gigabytes per human sample
(`delta.log`) that `collapse` does not use; it is deleted as soon as
each sample is grouped, and the other logs after the consensus, unless
`keep_logs = TRUE`.

## See also

[`read_terminus_groups()`](https://luciorq.github.io/isoformic/reference/read_terminus_groups.md),
[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md),
[`salmon_quant()`](https://luciorq.github.io/isoformic/reference/salmon_quant.md)

## Examples

``` r
if (FALSE) { # \dontrun{
groups <- run_terminus(quant_dirs, output_dir = "terminus")
iso <- add_uncertainty(iso, paths = quant_dirs, groups = groups)
tx_uncertainty(iso)
} # }
```
