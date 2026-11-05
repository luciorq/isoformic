# Read and cut the transcript trees of TreeTerminus

TreeTerminus builds, from Salmon replicates, a forest of transcript
trees: the leaves are transcripts and every internal node is the sum of
the transcripts below it, built so that the inferential uncertainty
falls as the tree is ascended. Unlike the flat groups of Terminus, a
tree lets the resolution change from one isoform to the next: an isoform
the reads determine stays a leaf, one they do not is replaced by the
smallest node that they do determine.

- `read_treeterminus()` reads the trees (Newick, one per line):
  `cluster_nwk.txt` for the trees of all samples
  (`treeterminus group --mean_inf true`), `group_nwk.txt` for one
  sample, or the output of `treeterminus consensus`. TreeTerminus labels
  the leaves with the 0-based position of the transcript in `quant.sf`;
  `salmon_dir` (any of the Salmon directories, which share the order) or
  `transcripts` turns them into transcript IDs.

- `cut_transcript_trees()` computes the inferential relative variance of
  every node from the replicates stored in an
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
  (the draws of the transcripts below the node are summed), replaces
  every transcript whose own variance is at least `cutoff` by its lowest
  ancestor below `cutoff`, and returns the groups, which do not overlap
  (when two chosen nodes are nested, the larger one is kept). The result
  goes to `add_uncertainty(groups = )`, so that the confidence classes
  and `collapse_groups` work at the resolution the data support.

Transcripts whose whole path to the root stays above the cutoff are
unresolved: by default they are left as they are
(`unresolved = "transcript"`), where the confidence classes still flag
them as uncertain; `unresolved = "root"` merges every tree that has one
into a single group.

## Usage

``` r
read_treeterminus(
  path,
  salmon_dir = NULL,
  transcripts = NULL,
  ignore_after_bar = TRUE
)

cut_transcript_trees(
  trees,
  x,
  cutoff = NULL,
  unresolved = c("transcript", "root")
)
```

## Arguments

- path:

  A Newick file of TreeTerminus, or a folder holding `cluster_nwk.txt`
  (or a single `group_nwk.txt`).

- salmon_dir:

  A Salmon output directory quantified against the same index, whose
  `quant.sf` gives the order of the transcripts.

- transcripts:

  Alternatively, the transcript names in `quant.sf` order.

- ignore_after_bar:

  Logical. Keep only the part of transcript names before the first `|`.

- trees:

  The result of `read_treeterminus()`.

- x:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
  with inferential replicates (see
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)).

- cutoff:

  Inferential relative variance below which a node is well determined
  (the `inf_rv_cutoff` of
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  by default).

- unresolved:

  `"transcript"` or `"root"`, see above.

## Value

`read_treeterminus()`: a tibble with one row per node and leaf:
`tree_id`, `node_id` (internal nodes; `NA` for leaves), `transcript_id`
(leaves; `NA` for nodes) and `parent_id` (`NA` for the roots).
`cut_transcript_trees()`: a tibble with `transcript_id`, `group_id`,
`group_size`, `group_inf_rv` and `tree_id` for the transcripts that are
grouped, with the cutoff in the `cutoff` attribute.

## References

Singh NP, Love MI, Patro R (2023). TreeTerminus: creating transcript
trees using inferential replicate counts. *iScience* 26(6):106961.

## See also

[`read_terminus_groups()`](https://luciorq.github.io/isoformic/reference/read_terminus_groups.md),
[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)

## Examples

``` r
nwk <- tempfile(fileext = ".txt")
writeLines(c("((0,1),2);", "(3,4);"), nwk)
read_treeterminus(nwk, transcripts = c("tA1", "tA2", "tA3", "tB1", "tB2"))
#> # A tibble: 8 × 4
#>   tree_id node_id transcript_id parent_id
#>     <int>   <int> <chr>             <int>
#> 1       1       1 NA                   NA
#> 2       1       2 NA                    1
#> 3       1      NA tA1                   2
#> 4       1      NA tA2                   2
#> 5       1      NA tA3                   1
#> 6       2       3 NA                   NA
#> 7       2      NA tB1                   3
#> 8       2      NA tB2                   3
```
