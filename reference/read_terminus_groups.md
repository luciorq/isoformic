# Read the transcript groups that Terminus builds

Terminus groups the transcripts whose abundances the reads cannot
estimate separately and reports one consensus grouping for a set of
samples. Run it on Salmon output made with `--numGibbsSamples` (or
`--numBootstraps`), `--dumpEq` and `--dumpEqWeights`:

## Usage

``` r
read_terminus_groups(path, ignore_after_bar = TRUE)
```

## Arguments

- path:

  The `clusters.txt` file, or a folder that holds it (the output folder
  of one sample, or the folder given to `-o` when it holds one sample
  folder).

- ignore_after_bar:

  Logical. Keep only the part of transcript names before the first `|`.

## Value

A tibble with `group_id`, `transcript_id` and `group_size`.

## Details

    terminus group -d quant/sample1 -o terminus/out
    terminus collapse -d quant/sample1 quant/sample2 -o terminus/out

`collapse` writes `clusters.txt` in the output folder of every sample:
one line per group, the name of the group followed by its transcripts,
separated by commas. `read_terminus_groups()` reads that file into a
long table that
[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
accepts in its `groups` argument, so that the groups of Terminus replace
the groups built from the correlation of the replicates.

## See also

[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md),
[`isoform_ambiguity()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)

## Examples

``` r
clusters <- tempfile(fileext = ".txt")
writeLines(
  c(
    "NewTr1,ENST00000831217.1,ENST00000831165.1",
    "NewTr2,ENST00000831423.1,ENST00000831157.1,ENST00000831246.1"
  ),
  clusters
)
read_terminus_groups(clusters)
#> # A tibble: 5 × 3
#>   group_id transcript_id     group_size
#>   <chr>    <chr>                  <int>
#> 1 NewTr1   ENST00000831217.1          2
#> 2 NewTr1   ENST00000831165.1          2
#> 3 NewTr2   ENST00000831423.1          3
#> 4 NewTr2   ENST00000831157.1          3
#> 5 NewTr2   ENST00000831246.1          3
```
