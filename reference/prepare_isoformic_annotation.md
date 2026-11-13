# Convert a GFF3/GTF Annotation to Parquet Files

Parses a GENCODE/Ensembl style GFF3 or GTF annotation file with DuckDB
in a single pass and writes one Parquet table per feature level (genes,
transcripts, exons) plus metadata and sequence-region tables into an
*annotation directory*. This directory is what
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
objects point to through their `annot_path` property.

## Usage

``` r
prepare_isoformic_annotation(
  input_path,
  output_path = NULL,
  file_type = c("auto", "gff", "gtf"),
  annot_name = NULL,
  overwrite = FALSE,
  haplotype = NULL
)
```

## Arguments

- input_path:

  Path to a GFF3 or GTF file, optionally gzip-compressed.

- output_path:

  Directory under which the annotation directory is created. `NULL`
  (default) uses
  [`get_isoformic_cache()`](https://luciorq.github.io/isoformic/reference/get_isoformic_cache.md).

- file_type:

  `"auto"` (default) detects the format from the file extension;
  otherwise `"gff"` or `"gtf"`.

- annot_name:

  Name of the annotation directory. Defaults to the input file name,
  lower-cased, with dots replaced by underscores.

- overwrite:

  Logical. Re-create the annotation directory even if a complete one
  already exists.

- haplotype:

  Optional label (for example `"maternal"` or `"hap1"`) stored in the
  `haplotype` column of every table, for annotations of one haplotype of
  a phased assembly. See
  [`merge_annotations()`](https://luciorq.github.io/isoformic/reference/merge_annotations.md).

## Value

The path to the annotation directory (character string).

## Details

Reading a full human GENCODE release (about 3 million features,
gzip-compressed) takes a few seconds; the result is cached so repeated
calls with the same file and output location return immediately.

The parser is vocabulary-agnostic so that annotations from different
providers and assemblies work the same way:

- GENCODE GFF3/GTF (GRCh38, GRCm39): `gene`/`transcript`/`exon` features
  with `gene_id`, `transcript_id`, `gene_type`, `transcript_type`.

- Ensembl and Ensembl Rapid Release GFF3 (e.g. T2T-CHM13
  `GCA_009914755.4`): `ncRNA_gene`, `mRNA`, `lnc_RNA`, ... features,
  `gene:`/`transcript:` prefixed IDs, `biotype`, `rank`.

- NCBI RefSeq GFF3 (e.g. `GCF_009914755.1_T2T-CHM13v2.0`):
  `gene-`/`rna-` prefixed IDs, `gene_biotype`, `mRNA`/`lnc_RNA`/`miRNA`
  features and pseudogenes carrying exons directly.

- T2T consortium Liftoff and CAT/Liftoff GFF3 (RefSeq or GENCODE lifted
  to CHM13) and HPRC per-haplotype CAT annotations:
  `transcript_biotype`, `gbkey`, and the `source_gene` /
  `source_transcript` attributes, kept as `source_gene_id` /
  `source_transcript_id` so that lifted annotations can be mapped back
  to the GRCh38 GENCODE or RefSeq identifiers.

- MANE GFF3/GTF and CHESS GFF3 (same vocabulary as GENCODE/RefSeq).

- FANTOM CAT (FANTOM5 lv1-lv4) and FANTOM6 CAT GTF: `__na` placeholders
  become `NA`, `geneClass` is used as the gene (and, when the transcript
  has no type, the transcript) type, and `geneCategory`, `DHS_type`,
  `coding_status`, `cumulative_support` and `TIEScore` are kept. CAGE
  cluster identifiers come from the FANTOM ID mapping table, see
  [`add_fantom_cat_tables()`](https://luciorq.github.io/isoformic/reference/add_fantom_cat_tables.md).

- StringTie, Bambu, FLAIR, IsoQuant and SQANTI3 GTF files from long-read
  or assembled transcriptomes: `reference_id`/`cmp_ref`, `ref_gene_id`
  and `ref_gene_name` populate the source identifiers and gene names;
  novel transcripts without a biotype get `transcript_type = NA` until
  [`add_sqanti_classification()`](https://luciorq.github.io/isoformic/reference/add_sqanti_classification.md)
  assigns SQANTI3 structural categories.

Genes are features of a gene-like type or parents of transcripts;
transcripts are features of a transcript-like type or parents of exons.
Transcript types are normalised (`mRNA` to `protein_coding`, `lnc_RNA`
to `lncRNA`) and fall back to the feature type (`miRNA`, `snRNA`, ...)
when no biotype attribute exists.

- genes: `gene_id`, `gene_name`, `gene_type`, `source_gene_id`,
  `gene_class`, `gene_category`, `dhs_type`, `coding_status`,
  `cumulative_support` (FANTOM CAT attributes, `NA` elsewhere),
  `haplotype`, `seqid`, `start_pos`, `end_pos`, `strand`.

- transcripts: `transcript_id`, `transcript_name`, `transcript_type`,
  `gene_id`, `source_gene_id`, `source_transcript_id`, `tags`
  (comma-separated `tag` attributes such as `MANE_Select`, `basic`,
  `Ensembl_canonical`), `is_canonical` (`TRUE` for MANE Select, Ensembl
  canonical or RefSeq Select transcripts), `coding_status` (FANTOM CAT),
  `transcript_support_level` (GENCODE TSL), `tie_score` (FANTOM CAT
  TIEScore), `cds_start`, `cds_end`, `cds_length` (aggregated from the
  `CDS` features; `NA` for non-coding transcripts), `haplotype`,
  `seqid`, `start_pos`, `end_pos`, `strand`, `n_exons`, `tx_length`.

- exons: `exon_id`, `transcript_id`, `exon_number`, `haplotype`,
  `seqid`, `start_pos`, `end_pos`, `strand`.

- metadata: `key`, `value` (source file, format, header pragmas,
  creation time).

- seq_regions: `seqid`, `start_pos`, `end_pos`.

## Examples

``` r
gff <- system.file("extdata", "example.gff3", package = "isoformic")
annot_dir <- prepare_isoformic_annotation(gff, output_path = tempdir())
#> ℹ Using cached annotation directory /tmp/RtmpmIi7db/example.
list.files(annot_dir)
#> [1] "annot_data_exons.parquet"       "annot_data_genes.parquet"      
#> [3] "annot_data_transcripts.parquet" "annot_metadata.parquet"        
#> [5] "annot_seq_regions.parquet"     
head(read_parquet(file.path(annot_dir, "annot_data_transcripts.parquet")))
#> # A tibble: 6 × 21
#>   transcript_id     transcript_name transcript_type gene_id       source_gene_id
#>   <chr>             <chr>           <chr>           <chr>         <chr>         
#> 1 ENST00000282397.9 FLT1-201        protein_coding  ENSG00000102… NA            
#> 2 ENST00000615611.4 FLT1-206        protein_coding  ENSG00000102… NA            
#> 3 ENST00000617835.4 FLT1-208        protein_coding  ENSG00000102… NA            
#> 4 ENST00000543394.2 FLT1-205        protein_coding  ENSG00000102… NA            
#> 5 ENST00000540678.2 FLT1-203        protein_coding  ENSG00000102… NA            
#> 6 ENST00000541932.5 FLT1-204        protein_coding  ENSG00000102… NA            
#> # ℹ 16 more variables: source_transcript_id <chr>, tags <chr>,
#> #   coding_status <chr>, transcript_support_level <chr>, tie_score <dbl>,
#> #   cds_start <int>, cds_end <int>, cds_length <int>, is_canonical <lgl>,
#> #   haplotype <chr>, seqid <chr>, start_pos <int>, end_pos <int>, strand <chr>,
#> #   n_exons <int>, tx_length <int>
```
