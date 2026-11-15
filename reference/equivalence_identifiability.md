# Which transcripts can the equivalence classes tell apart?

Salmon summarises the reads as equivalence classes: sets of transcripts
with the number of fragments compatible with exactly that set. The
abundance of a transcript is determined by these counts only when the
classes separate it from the other transcripts. Two isoforms that always
appear together in the same classes, for example, can exchange fragments
without any change to the class counts.

`equivalence_identifiability()` reads the classes
(`salmon quant --dumpEq`) and reports, for every transcript, whether it
is identifiable from them: whether its abundance is fixed once the total
of every class is fixed. This is a property of the classes observed in a
sample, not of the posterior draws: a transcript can be non-identifiable
and still show a small inferential variance (see Details).

## Usage

``` r
equivalence_identifiability(
  paths,
  transcripts = NULL,
  ignore_after_bar = TRUE,
  max_dense = 400L,
  max_iterations = 50000L
)
```

## Arguments

- paths:

  Character vector of Salmon output directories (containing `quant.sf`)
  or paths to `quant.sf` files. Names, when present, are the sample IDs;
  otherwise the directory name is used.

- transcripts:

  Optional transcript IDs to keep. Reading every transcript of a human
  annotation for many samples takes gigabytes of memory, so
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  only reads the expressed ones.

- ignore_after_bar:

  Logical. Keep only the part of transcript names before the first `|`
  (GENCODE FASTA headers used as Salmon references).

- max_dense:

  Largest block (number of transcripts) solved with a dense
  decomposition. Larger blocks use the iterative solver.

- max_iterations:

  Maximum number of iterations of the iterative solver per block. The
  transcripts of a block that does not converge get `NA`.

## Value

A tibble with one row per transcript that has fragments in at least one
sample: `transcript_id`, `n_samples` (samples in which it has fragments
and was resolved), `n_identifiable` (in how many of them it is
identifiable), `identifiable_fraction` and `identifiable` (`TRUE` when
identifiable in all of them, `NA` when it was resolved in none).

## Details

With `M` the class by transcript membership matrix of one sample, a
transcript is identifiable when its unit vector lies in the row space of
`M`, that is, when no vector of the null space of `M` changes its
abundance. The function avoids the decomposition of the full matrix:

1.  A class with a single unresolved transcript resolves it (a class of
    one transcript does so directly); repeated until nothing changes.

2.  Transcripts that are members of exactly the same classes cannot be
    told apart: they are marked non-identifiable, replaced by their sum,
    and step 1 is repeated.

3.  What is left splits into independent blocks. Blocks of up to
    `max_dense` transcripts are solved with a singular value
    decomposition. Larger blocks need the `Matrix` package: two random
    vectors are projected on the null space with conjugate gradients,
    and a transcript with a non-zero coordinate is non-identifiable.

Limits:

- Only class membership is used, not the conditional probabilities of
  `--dumpEqWeights`. Isoforms that share all their fragments have
  proportional weights, so they stay non-identifiable with the weights;
  other cases may be resolved by them.

- The annotation is taken as complete, and the classes are those
  observed in the sample: a transcript can be identifiable in a deep
  sample and not in a shallow one. Transcripts without fragments in a
  sample are not assessed in it.

- This is not the range of optima of Zheng et al. (2022), which needs
  their sufficient statistics and a linear program per transcript.

- Human samples have one block of tens of thousands of transcripts
  linked by multi-mapping fragments; expect minutes per sample.

In the six GSE129540 samples of the quantification stability article,
the estimates of non-identifiable transcripts differed more between two
halves of the reads of a sample than those of identifiable transcripts
of the same expression: a median absolute log2 ratio of 2.4 against 0.7
at 20 to 100 reads. The posterior draws call three quarters of both
kinds uncertain, so the two measures are close to independent.

## References

Zheng H, Ma C, Kingsford C (2022). Deriving ranges of optimal estimated
transcript expression due to nonidentifiability. Journal of
Computational Biology 29:121-139.
[doi:10.1089/cmb.2021.0444](https://doi.org/10.1089/cmb.2021.0444)

## See also

[`equivalence_support()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
for the share of unique fragments,
[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
for the uncertainty measured by the posterior draws.

## Examples

``` r
# Two transcripts that always share their classes, and one with its own
dir <- file.path(tempdir(), "identifiability_example")
dir.create(file.path(dir, "aux_info"), recursive = TRUE, showWarnings = FALSE)
writeLines(
  c("3", "2", "tx1", "tx2", "tx3", "2\t0\t1\t10", "1\t2\t5"),
  file.path(dir, "aux_info", "eq_classes.txt")
)
equivalence_identifiability(c(sample1 = dir))
#> # A tibble: 3 × 5
#>   transcript_id n_samples n_identifiable identifiable_fraction identifiable
#>   <chr>             <int>          <int>                 <dbl> <lgl>       
#> 1 tx1                   1              0                     0 FALSE       
#> 2 tx2                   1              0                     0 FALSE       
#> 3 tx3                   1              1                     1 TRUE        
```
