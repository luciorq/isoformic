# Annotation directories ------------------------------------------------------
#
# An "annotation directory" is the on-disk, Parquet-backed representation of a
# GFF3/GTF annotation used by `IsoformicExperiment`. It holds:
#
#   annot_data_genes.parquet        one row per gene
#   annot_data_transcripts.parquet  one row per transcript (+ n_exons, tx_length)
#   annot_data_exons.parquet        one row per exon
#   annot_metadata.parquet          key/value table describing the source
#   annot_seq_regions.parquet       one row per sequence (chromosome) length
#   annot_data_sequences.parquet    optional, one row per transcript sequence
#                                   (see extract_transcript_sequences())
#   annot_data_orfs.parquet         optional, the ORFs computed so far (written
#                                   by transcript_orfs(), see R/orfs.R)
#
# All readers are DuckDB based and read only the file requested.

annot_file_names <- function() {
  c(
    genes = "annot_data_genes.parquet",
    transcripts = "annot_data_transcripts.parquet",
    exons = "annot_data_exons.parquet",
    metadata = "annot_metadata.parquet",
    seq_regions = "annot_seq_regions.parquet",
    sequences = "annot_data_sequences.parquet",
    splice_sites = "annot_data_splice_sites.parquet",
    orfs = "annot_data_orfs.parquet"
  )
}

annot_file_path <- function(annot_path, type) {
  fs::path(annot_path, annot_file_names()[[type]])
}

#' Is `path` a complete annotation directory?
#' @keywords internal
#' @noRd
is_annotation_dir <- function(path) {
  if (!rlang::is_string(path) || !isTRUE(fs::dir_exists(path))) {
    return(FALSE)
  }
  required <- annot_file_names()[c("genes", "transcripts", "exons")]
  all(fs::file_exists(fs::path(path, required)))
}

#' Resolve an annotation source into an annotation directory path
#'
#' Accepts either an existing annotation directory or a GFF3/GTF file.
#' Files are converted (once, cached) with [prepare_isoformic_annotation()].
#'
#' @keywords internal
#' @noRd
resolve_annot_path <- function(annot_path, output_path = NULL) {
  check_string(annot_path)
  if (is_annotation_dir(annot_path)) {
    return(as.character(fs::path_abs(annot_path)))
  }
  if (isTRUE(fs::file_exists(annot_path)) && !isTRUE(fs::is_dir(annot_path))) {
    return(prepare_isoformic_annotation(
      input_path = annot_path,
      output_path = output_path
    ))
  }
  cli::cli_abort(
    message = c(
      x = "{.arg annot_path} must be an annotation directory or a GFF3/GTF file.",
      i = "Received: {.path {annot_path}}.",
      i = "Create an annotation directory with {.fn prepare_isoformic_annotation}."
    ),
    class = "isoformic_annot_path_dont_exist"
  )
}

#' Fingerprint of an annotation file for cache invalidation
#' @keywords internal
#' @noRd
annot_signature <- function(annot_path, type) {
  file <- annot_file_path(annot_path, type)
  if (!isTRUE(fs::file_exists(file))) {
    return(paste0("missing:", file))
  }
  info <- fs::file_info(file)
  paste(
    as.character(file),
    as.character(info$size),
    format(as.numeric(info$modification_time), digits = 17),
    sep = "|"
  )
}

#' Read one table of an annotation directory
#' @keywords internal
#' @noRd
read_annotation_table <- function(annot_path, type, columns = NULL) {
  type <- rlang::arg_match(
    type,
    c(
      "genes",
      "transcripts",
      "exons",
      "metadata",
      "seq_regions",
      "sequences",
      "splice_sites"
    )
  )
  file <- annot_file_path(annot_path, type)
  if (!isTRUE(fs::file_exists(file))) {
    if (type %in% c("metadata", "seq_regions")) {
      return(empty_annotation_table(type))
    }
    if (identical(type, "sequences")) {
      cli::cli_abort(
        message = c(
          x = "The annotation in {.path {annot_path}} has no transcript sequences.",
          i = "Add them with {.fn extract_transcript_sequences} from a genome or transcript FASTA."
        ),
        class = "isoformic_no_sequences"
      )
    }
    if (identical(type, "splice_sites")) {
      cli::cli_abort(
        message = c(
          x = "The annotation in {.path {annot_path}} has no splice-site sequences.",
          i = "Add them with {.fn extract_splice_sites} from a genome FASTA."
        ),
        class = "isoformic_no_splice_sites"
      )
    }
    cli::cli_abort(
      message = c(
        x = "Annotation file {.path {file}} does not exist.",
        i = "Re-create the annotation with {.fn prepare_isoformic_annotation}."
      ),
      class = "isoformic_annot_path_dont_exist"
    )
  }
  read_parquet(file, columns = columns)
}

empty_annotation_table <- function(type) {
  switch(
    type,
    genes = tibble::tibble(
      gene_id = character(),
      gene_name = character(),
      gene_type = character(),
      source_gene_id = character(),
      gene_class = character(),
      gene_category = character(),
      dhs_type = character(),
      coding_status = character(),
      cumulative_support = character(),
      haplotype = character(),
      seqid = character(),
      start_pos = integer(),
      end_pos = integer(),
      strand = character()
    ),
    transcripts = tibble::tibble(
      transcript_id = character(),
      transcript_name = character(),
      transcript_type = character(),
      gene_id = character(),
      source_gene_id = character(),
      source_transcript_id = character(),
      tags = character(),
      coding_status = character(),
      transcript_support_level = character(),
      tie_score = numeric(),
      cds_start = integer(),
      cds_end = integer(),
      cds_length = integer(),
      is_canonical = logical(),
      haplotype = character(),
      seqid = character(),
      start_pos = integer(),
      end_pos = integer(),
      strand = character(),
      n_exons = integer(),
      tx_length = integer()
    ),
    exons = tibble::tibble(
      exon_id = character(),
      transcript_id = character(),
      exon_number = integer(),
      haplotype = character(),
      seqid = character(),
      start_pos = integer(),
      end_pos = integer(),
      strand = character()
    ),
    metadata = tibble::tibble(key = character(), value = character()),
    seq_regions = tibble::tibble(
      seqid = character(),
      start_pos = integer(),
      end_pos = integer()
    ),
    sequences = tibble::tibble(
      transcript_id = character(),
      sequence = character(),
      seq_length = integer()
    )
  )
}

#' Convert a key/value metadata table to a named list (and back)
#' @keywords internal
#' @noRd
kv_to_list <- function(kv) {
  if (is.null(kv) || !nrow(kv)) {
    return(list())
  }
  out <- as.list(as.character(kv$value))
  names(out) <- as.character(kv$key)
  out
}

list_to_kv <- function(x) {
  if (!length(x)) {
    return(tibble::tibble(key = character(), value = character()))
  }
  scalar <- vapply(
    x,
    function(v) is.atomic(v) && length(v) == 1L && !is.null(v),
    logical(1L)
  )
  if (!all(scalar)) {
    cli::cli_warn(c(
      "!" = "Dropping non-scalar metadata entries: {.val {names(x)[!scalar]}}."
    ))
  }
  x <- x[scalar]
  tibble::tibble(
    key = names(x),
    value = vapply(x, function(v) as.character(v), character(1L))
  )
}

# Preparation -----------------------------------------------------------------

#' Convert a GFF3/GTF Annotation to Parquet Files
#'
#' Parses a GENCODE/Ensembl style GFF3 or GTF annotation file with DuckDB in
#' a single pass and writes one Parquet table per feature level (genes,
#' transcripts, exons) plus metadata and sequence-region tables into an
#' *annotation directory*. This directory is what [IsoformicExperiment]
#' objects point to through their `annot_path` property.
#'
#' Reading a full human GENCODE release (3 to 4 million features,
#' gzip-compressed) takes under a minute; the result is cached so repeated
#' calls with the same file and output location return immediately.
#'
#' @param input_path Path to a GFF3 or GTF file, optionally gzip-compressed.
#' @param output_path Directory under which the annotation directory is
#'   created. `NULL` (default) uses [get_isoformic_cache()].
#' @param file_type `"auto"` (default) detects the format from the file
#'   extension; otherwise `"gff"` or `"gtf"`.
#' @param annot_name Name of the annotation directory. Defaults to the
#'   input file name, lower-cased, with dots replaced by underscores.
#' @param overwrite Logical. Re-create the annotation directory even if a
#'   complete one already exists.
#' @param haplotype Optional label (for example `"maternal"` or `"hap1"`)
#'   stored in the `haplotype` column of every table, for annotations of one
#'   haplotype of a phased assembly. See [merge_annotations()].
#'
#' @returns The path to the annotation directory (character string).
#'
#' @details
#' The parser is vocabulary-agnostic so that annotations from different
#' providers and assemblies work the same way:
#'
#' * GENCODE GFF3/GTF (GRCh38, GRCm39): `gene`/`transcript`/`exon` features
#'   with `gene_id`, `transcript_id`, `gene_type`, `transcript_type`.
#' * Ensembl and Ensembl Rapid Release GFF3 (e.g. T2T-CHM13
#'   `GCA_009914755.4`): `ncRNA_gene`, `mRNA`, `lnc_RNA`, ... features,
#'   `gene:`/`transcript:` prefixed IDs, `biotype`, `rank`.
#' * NCBI RefSeq GFF3 (e.g. `GCF_009914755.1_T2T-CHM13v2.0`): `gene-`/`rna-`
#'   prefixed IDs, `gene_biotype`, `mRNA`/`lnc_RNA`/`miRNA` features and
#'   pseudogenes carrying exons directly.
#' * T2T consortium Liftoff and CAT/Liftoff GFF3 (RefSeq or GENCODE lifted
#'   to CHM13) and HPRC per-haplotype CAT annotations: `transcript_biotype`,
#'   `gbkey`, and the `source_gene` / `source_transcript` attributes, kept as
#'   `source_gene_id` / `source_transcript_id` so that lifted annotations can
#'   be mapped back to the GRCh38 GENCODE or RefSeq identifiers.
#' * MANE GFF3/GTF and CHESS GFF3 (same vocabulary as GENCODE/RefSeq).
#' * FANTOM CAT (FANTOM5 lv1-lv4) and FANTOM6 CAT GTF: `__na` placeholders
#'   become `NA`, `geneClass` is used as the gene (and, when the transcript
#'   has no type, the transcript) type, and `geneCategory`, `DHS_type`,
#'   `coding_status`, `cumulative_support` and `TIEScore` are kept. CAGE
#'   cluster identifiers come from the FANTOM ID mapping table, see
#'   [add_fantom_cat_tables()].
#' * StringTie, Bambu, FLAIR, IsoQuant and SQANTI3 GTF files from long-read
#'   or assembled transcriptomes: `reference_id`/`cmp_ref`, `ref_gene_id` and
#'   `ref_gene_name` populate the source identifiers and gene names; novel
#'   transcripts without a biotype get `transcript_type = NA` until
#'   [add_sqanti_classification()] assigns SQANTI3 structural categories.
#'
#' Genes are features of a gene-like type or parents of transcripts;
#' transcripts are features of a transcript-like type or parents of exons.
#' Transcript types are normalised (`mRNA` to `protein_coding`, `lnc_RNA`
#' to `lncRNA`) and fall back to the feature type (`miRNA`, `snRNA`, ...)
#' when no biotype attribute exists.
#'
#' * genes: `gene_id`, `gene_name`, `gene_type`, `source_gene_id`,
#'   `gene_class`, `gene_category`, `dhs_type`, `coding_status`,
#'   `cumulative_support` (FANTOM CAT attributes, `NA` elsewhere),
#'   `haplotype`, `seqid`, `start_pos`, `end_pos`, `strand`.
#' * transcripts: `transcript_id`, `transcript_name`, `transcript_type`,
#'   `gene_id`, `source_gene_id`, `source_transcript_id`, `tags`
#'   (comma-separated `tag` attributes such as `MANE_Select`, `basic`,
#'   `Ensembl_canonical`), `is_canonical` (`TRUE` for MANE Select, Ensembl
#'   canonical or RefSeq Select transcripts), `coding_status` (FANTOM CAT),
#'   `transcript_support_level` (GENCODE TSL), `tie_score` (FANTOM CAT
#'   TIEScore), `cds_start`, `cds_end`, `cds_length` (aggregated from the
#'   `CDS` features; `NA` for non-coding transcripts), `haplotype`,
#'   `seqid`, `start_pos`, `end_pos`, `strand`, `n_exons`, `tx_length`.
#' * exons: `exon_id`, `transcript_id`, `exon_number`, `haplotype`, `seqid`,
#'   `start_pos`, `end_pos`, `strand`.
#' * metadata: `key`, `value` (source file, format, header pragmas,
#'   creation time).
#' * seq_regions: `seqid`, `start_pos`, `end_pos`.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "isoformic")
#' annot_dir <- prepare_isoformic_annotation(gff, output_path = tempdir())
#' list.files(annot_dir)
#' head(read_parquet(file.path(annot_dir, "annot_data_transcripts.parquet")))
#'
#' @export
prepare_isoformic_annotation <- function(
  input_path,
  output_path = NULL,
  file_type = c("auto", "gff", "gtf"),
  annot_name = NULL,
  overwrite = FALSE,
  haplotype = NULL
) {
  assert_file_exists(input_path)
  check_string(haplotype, allow_null = TRUE)
  file_type <- rlang::arg_match(file_type)
  if (identical(file_type, "auto")) {
    file_type <- detect_annotation_type(input_path)
  }
  if (is.null(output_path) || identical(output_path, "")) {
    output_path <- get_isoformic_cache()
  }
  if (is.null(annot_name)) {
    annot_name <- annotation_name_from_path(input_path)
  }
  check_string(annot_name)
  annot_dir <- as.character(fs::path_abs(fs::path(output_path, annot_name)))

  if (isFALSE(overwrite) && is_annotation_dir(annot_dir)) {
    cli::cli_inform(c(
      i = "Using cached annotation directory {.path {annot_dir}}."
    ))
    return(annot_dir)
  }
  fs::dir_create(annot_dir)

  input_path <- as.character(fs::path_abs(input_path))
  genes_file <- as.character(annot_file_path(annot_dir, "genes"))
  transcripts_file <- as.character(annot_file_path(annot_dir, "transcripts"))
  exons_file <- as.character(annot_file_path(annot_dir, "exons"))
  seq_regions_file <- as.character(annot_file_path(annot_dir, "seq_regions"))

  cli::cli_inform(c(
    i = "Parsing {.path {fs::path_file(input_path)}} ({toupper(file_type)}) into {.path {annot_dir}}."
  ))

  duckdb_run(
    sql_string = build_annotation_sql(file_type),
    input_path = input_path,
    genes_file = genes_file,
    transcripts_file = transcripts_file,
    exons_file = exons_file,
    seq_regions_file = seq_regions_file,
    haplotype = haplotype %||% NA_character_
  )

  header <- read_annotation_header(input_path)
  metadata <- c(
    list(
      annot_name = annot_name,
      source_path = input_path,
      source_file = fs::path_file(input_path),
      file_type = file_type,
      created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      isoformic_version = as.character(utils::packageVersion("isoformic"))
    ),
    header$metadata
  )
  write_parquet(list_to_kv(metadata), annot_file_path(annot_dir, "metadata"))
  # The SQL wrote one row per sequence with the largest feature end. Declared
  # lengths (`##sequence-region`) replace them where the file has them.
  declared <- read_sequence_regions(input_path)
  if (nrow(declared) > 0L) {
    observed <- read_parquet(seq_regions_file)
    observed <- observed[!observed$seqid %in% declared$seqid, ]
    write_parquet(dplyr::bind_rows(declared, observed), seq_regions_file)
  }
  return(annot_dir)
}

#' Build the DuckDB SQL that parses a GFF3/GTF file into Parquet tables
#'
#' The parser does not depend on one vocabulary. Genes are features of a
#' gene-like type (`gene`, `pseudogene`, `ncRNA_gene`, ...) or features that
#' are the `Parent` of a transcript; transcripts are features of a
#' transcript-like type (`transcript`, `mRNA`, `lnc_RNA`, `miRNA`, ...) or
#' features that are the `Parent` of an exon. Identifiers come from
#' `gene_id`/`transcript_id` attributes when present and from `ID`/`Parent`
#' otherwise, with the `gene:`, `transcript:`, `gene-`, `rna-` prefixes used
#' by Ensembl and NCBI removed. This covers GENCODE, Ensembl (including Rapid
#' Release), NCBI RefSeq and the T2T-CHM13 Liftoff/CAT annotations.
#'
#' @keywords internal
#' @noRd
build_annotation_sql <- function(file_type) {
  attr_expr <- function(keys) {
    exprs <- vapply(
      keys,
      function(key) {
        pattern <- if (identical(file_type, "gff")) {
          sprintf("(^|;)%s=([^;]+)", key)
        } else {
          sprintf("(^|; *)%s \"?([^\";]+)\"?", key)
        }
        # '__na' (FANTOM CAT), 'NA', 'N/A' and '.' are placeholders for missing
        sprintf(
          "nullif(nullif(nullif(nullif(nullif(regexp_extract(attributes, '%s', 2), ''), '__na'), 'NA'), 'N/A'), '.')",
          pattern
        )
      },
      character(1L)
    )
    paste0("coalesce(", paste(exprs, collapse = ", "), ")")
  }
  strip_prefix <- function(expr) {
    sprintf(
      "regexp_replace(%s, '^(gene|transcript|rna|exon|cds|mrna)[:-]', '', 'i')",
      expr
    )
  }
  gene_types <- c("gene", "pseudogene", "ncRNA_gene", "protein_coding_gene")
  transcript_types <- c(
    "transcript",
    "mRNA",
    "lnc_RNA",
    "ncRNA",
    "pseudogenic_transcript",
    "primary_transcript",
    "unconfirmed_transcript",
    "miRNA",
    "snRNA",
    "snoRNA",
    "scaRNA",
    "scRNA",
    "rRNA",
    "tRNA",
    "guide_RNA",
    "antisense_RNA",
    "Y_RNA",
    "vault_RNA",
    "telomerase_RNA",
    "RNase_MRP_RNA",
    "RNase_P_RNA",
    "SRP_RNA",
    "V_gene_segment",
    "D_gene_segment",
    "J_gene_segment",
    "C_gene_segment"
  )
  sql_list <- function(x) paste0("('", paste(x, collapse = "', '"), "')")
  # Feature types of transcripts map to a transcript type when the file
  # carries no biotype attribute.
  type_to_biotype <- paste0(
    "CASE type ",
    "WHEN 'mRNA' THEN 'protein_coding' ",
    "WHEN 'lnc_RNA' THEN 'lncRNA' ",
    "WHEN 'pseudogenic_transcript' THEN 'pseudogene' ",
    "WHEN 'transcript' THEN NULL ",
    "WHEN 'exon' THEN NULL ",
    "ELSE type END"
  )
  biotype_normalise <- function(expr) {
    sprintf(
      "CASE %s WHEN 'mRNA' THEN 'protein_coding' WHEN 'lnc_RNA' THEN 'lncRNA' ELSE %s END",
      expr,
      expr
    )
  }
  # GTF files repeat `tag "..."` several times per feature; GFF3 files use
  # one comma-separated `tag=` attribute.
  tags_expr <- if (identical(file_type, "gff")) {
    attr_expr(c("tag"))
  } else {
    "nullif(array_to_string(regexp_extract_all(attributes, '(^|; *)tag \"([^\"]+)\"', 2), ','), '')"
  }
  tx_type_expr <- biotype_normalise(paste0(
    "coalesce(",
    attr_expr(c("transcript_type", "transcript_biotype", "biotype")),
    ", ",
    type_to_biotype,
    ", ",
    attr_expr(c("gbkey")),
    ")"
  ))
  gene_type_expr <- biotype_normalise(
    attr_expr(c("gene_type", "gene_biotype", "biotype", "geneClass"))
  )

  paste0(
    "CREATE TEMP TABLE parsed AS
SELECT
  type, seqid, start_pos, end_pos, strand,
  ",
    strip_prefix(attr_expr(c("ID"))),
    " AS id,
  ",
    strip_prefix(attr_expr(c("Parent"))),
    " AS parent,
  ",
    attr_expr(c("gene_id")),
    " AS gene_id_attr,
  ",
    attr_expr(c("transcript_id")),
    " AS transcript_id_attr,
  ",
    attr_expr(c("gene_name", "Name", "gene", "ref_gene_name")),
    " AS gene_name,
  ",
    gene_type_expr,
    " AS gene_type,
  ",
    attr_expr(c("transcript_name", "source_transcript_name")),
    " AS transcript_name,
  ",
    tx_type_expr,
    " AS transcript_type,
  ",
    attr_expr(c("source_gene", "ref_gene_id")),
    " AS source_gene_id,
  ",
    attr_expr(c("source_transcript", "reference_id", "cmp_ref")),
    " AS source_transcript_id,
  ",
    tags_expr,
    " AS tags,
  ",
    attr_expr(c("coding_status")),
    " AS coding_status,
  ",
    attr_expr(c("transcript_support_level")),
    " AS transcript_support_level,
  TRY_CAST(",
    attr_expr(c("TIEScore")),
    " AS DOUBLE) AS tie_score,
  ",
    attr_expr(c("geneClass")),
    " AS gene_class,
  ",
    attr_expr(c("geneCategory")),
    " AS gene_category,
  ",
    attr_expr(c("DHS_type")),
    " AS dhs_type,
  ",
    attr_expr(c("cumulative_support")),
    " AS cumulative_support,
  ",
    attr_expr(c("exon_id")),
    " AS exon_id,
  TRY_CAST(",
    attr_expr(c("exon_number", "rank")),
    " AS INTEGER) AS exon_number
FROM read_csv(
  {input_path},
  delim = '\t', header = false, comment = '#',
  auto_detect = false, quote = '', escape = '',
  columns = {{
    'seqid': 'VARCHAR', 'source': 'VARCHAR', 'type': 'VARCHAR',
    'start_pos': 'INTEGER', 'end_pos': 'INTEGER', 'score': 'VARCHAR',
    'strand': 'VARCHAR', 'phase': 'VARCHAR', 'attributes': 'VARCHAR'
  }}
)
WHERE type NOT IN ('five_prime_UTR', 'three_prime_UTR', 'start_codon',
                   'stop_codon', 'intron', 'region', 'chromosome', 'biological_region',
                   'enhancer', 'silencer', 'match', 'cDNA_match');

CREATE TEMP TABLE cds_rows AS
  SELECT
    coalesce(transcript_id_attr, parent) AS transcript_id,
    start_pos, end_pos
  FROM parsed
  WHERE type = 'CDS';

CREATE TEMP TABLE cds_stats AS
  SELECT
    transcript_id,
    min(start_pos)::INTEGER AS cds_start,
    max(end_pos)::INTEGER AS cds_end,
    sum(end_pos - start_pos + 1)::INTEGER AS cds_length
  FROM cds_rows
  WHERE transcript_id IS NOT NULL
  GROUP BY transcript_id;

CREATE TEMP TABLE exon_rows AS
  SELECT
    exon_id,
    coalesce(transcript_id_attr, parent) AS transcript_id,
    exon_number, seqid, start_pos, end_pos, strand
  FROM parsed
  WHERE type = 'exon';

CREATE TEMP TABLE tx_rows AS
  SELECT
    coalesce(transcript_id_attr, id) AS transcript_id,
    transcript_name, transcript_type,
    coalesce(gene_id_attr, parent, id) AS gene_id,
    source_gene_id, source_transcript_id, tags,
    coding_status, transcript_support_level, tie_score,
    seqid, start_pos, end_pos, strand
  FROM parsed
  WHERE type NOT IN ('exon', 'CDS')
    AND (
      type IN ",
    sql_list(transcript_types),
    "
      OR id IN (SELECT DISTINCT transcript_id FROM exon_rows WHERE transcript_id IS NOT NULL)
    );

CREATE TEMP TABLE gene_rows AS
  SELECT
    coalesce(gene_id_attr, id) AS gene_id,
    gene_name, gene_type, source_gene_id,
    gene_class, gene_category, dhs_type, coding_status, cumulative_support,
    seqid, start_pos, end_pos, strand
  FROM parsed
  WHERE type NOT IN ('exon', 'CDS')
    AND (
      type IN ",
    sql_list(gene_types),
    "
      OR id IN (SELECT DISTINCT gene_id FROM tx_rows WHERE gene_id IS NOT NULL)
    );

-- Files without gene features (StringTie, FLAIR, some GTFs): derive one gene
-- row per gene_id from its transcripts.
INSERT INTO gene_rows
  SELECT
    t.gene_id,
    any_value(p.gene_name) AS gene_name,
    any_value(p.gene_type) AS gene_type,
    any_value(t.source_gene_id) AS source_gene_id,
    any_value(p.gene_class) AS gene_class,
    any_value(p.gene_category) AS gene_category,
    any_value(p.dhs_type) AS dhs_type,
    any_value(p.coding_status) AS coding_status,
    any_value(p.cumulative_support) AS cumulative_support,
    any_value(t.seqid) AS seqid,
    min(t.start_pos) AS start_pos,
    max(t.end_pos) AS end_pos,
    any_value(t.strand) AS strand
  FROM tx_rows t
  LEFT JOIN parsed p
    ON p.type NOT IN ('exon', 'CDS') AND coalesce(p.transcript_id_attr, p.id) = t.transcript_id
  WHERE t.gene_id IS NOT NULL
    AND t.gene_id NOT IN (SELECT gene_id FROM gene_rows WHERE gene_id IS NOT NULL)
  GROUP BY t.gene_id;

COPY (
  SELECT
    gene_id,
    coalesce(gene_name, gene_id) AS gene_name,
    gene_type, source_gene_id,
    gene_class, gene_category, dhs_type, coding_status, cumulative_support,
    {haplotype}::VARCHAR AS haplotype,
    seqid, start_pos, end_pos, strand
  FROM gene_rows
  WHERE gene_id IS NOT NULL
  QUALIFY row_number() OVER (PARTITION BY gene_id ORDER BY start_pos) = 1
  ORDER BY seqid, start_pos, gene_id
) TO {genes_file} (FORMAT parquet, COMPRESSION zstd);

CREATE TEMP TABLE exon_stats AS
  SELECT
    transcript_id,
    count(*)::INTEGER AS n_exons,
    sum(end_pos - start_pos + 1)::INTEGER AS tx_length
  FROM exon_rows
  GROUP BY transcript_id;

COPY (
  SELECT
    t.transcript_id,
    coalesce(t.transcript_name, t.transcript_id) AS transcript_name,
    coalesce(t.transcript_type, g.gene_type) AS transcript_type,
    t.gene_id, t.source_gene_id, t.source_transcript_id,
    t.tags,
    t.coding_status, t.transcript_support_level, t.tie_score,
    c.cds_start, c.cds_end, c.cds_length,
    coalesce(regexp_matches(t.tags, 'MANE[_ ]Select|Ensembl_canonical|RefSeq Select'), false) AS is_canonical,
    {haplotype}::VARCHAR AS haplotype,
    t.seqid, t.start_pos, t.end_pos, t.strand,
    coalesce(e.n_exons, 0)::INTEGER AS n_exons,
    coalesce(e.tx_length, t.end_pos - t.start_pos + 1)::INTEGER AS tx_length
  FROM tx_rows t
  LEFT JOIN exon_stats e USING (transcript_id)
  LEFT JOIN cds_stats c USING (transcript_id)
  LEFT JOIN (
    SELECT gene_id, any_value(gene_type) AS gene_type FROM gene_rows GROUP BY gene_id
  ) g ON g.gene_id = t.gene_id
  WHERE t.transcript_id IS NOT NULL
  QUALIFY row_number() OVER (PARTITION BY t.transcript_id ORDER BY t.start_pos) = 1
  ORDER BY t.seqid, t.start_pos, t.transcript_id
) TO {transcripts_file} (FORMAT parquet, COMPRESSION zstd);

COPY (
  SELECT exon_id, transcript_id, exon_number, {haplotype}::VARCHAR AS haplotype,
         seqid, start_pos, end_pos, strand
  FROM exon_rows
  WHERE transcript_id IS NOT NULL
  ORDER BY transcript_id, exon_number, start_pos
) TO {exons_file} (FORMAT parquet, COMPRESSION zstd);

COPY (
  SELECT seqid, 1::INTEGER AS start_pos, max(end_pos)::INTEGER AS end_pos
  FROM parsed
  GROUP BY seqid
  ORDER BY seqid
) TO {seq_regions_file} (FORMAT parquet, COMPRESSION zstd);
"
  )
}

#' Detect annotation file type from its extension
#' @keywords internal
#' @noRd
detect_annotation_type <- function(path) {
  ext <- tolower(fs::path_ext(sub("\\.(gz|bgz|zst|bz2)$", "", path)))
  if (ext %in% c("gff", "gff3")) {
    return("gff")
  }
  if (ext %in% c("gtf")) {
    return("gtf")
  }
  cli::cli_abort(
    message = c(
      x = "Could not detect the annotation format of {.path {path}}.",
      i = "Set {.arg file_type} to {.val gff} or {.val gtf}."
    ),
    class = "isoformic_unsupported_file_type"
  )
}

#' Derive a file-system friendly annotation name from a file path
#' @keywords internal
#' @noRd
annotation_name_from_path <- function(path) {
  path |>
    fs::path_ext_remove() |>
    fs::path_file() |>
    tolower() |>
    gsub(pattern = ".", replacement = "_", fixed = TRUE)
}

#' Read every `##sequence-region` pragma of a GFF3 file
#'
#' GENCODE writes one pragma before the features of each sequence, not all of
#' them in the header, so the whole file is scanned (DuckDB, plain or gzip).
#' A pragma line has no tab and lands whole in the first column.
#'
#' @returns A tibble with `seqid`, `start_pos` and `end_pos`, in file order,
#'   one row per sequence (the first pragma of a sequence wins).
#' @keywords internal
#' @noRd
read_sequence_regions <- function(path) {
  path <- as.character(fs::path_abs(path))
  con <- duckdb_connect()
  cols <- paste0("'c", seq_len(9L), "': 'VARCHAR'", collapse = ", ")
  sql <- glue::glue_sql(
    "SELECT c1 AS line
     FROM read_csv({path}, delim = '\t', header = false, auto_detect = false,
       quote = '', escape = '', columns = {{{DBI::SQL(cols)}}},
       null_padding = true, strict_mode = false, max_line_size = 10000000)
     WHERE starts_with(c1, '##sequence-region')",
    .con = con
  )
  lines <- DBI::dbGetQuery(con, sql)$line
  res <- empty_annotation_table("seq_regions")
  if (length(lines)) {
    parts <- strsplit(
      trimws(sub("^##sequence-region\\s+", "", lines)),
      "\\s+"
    )
    parts <- parts[lengths(parts) >= 3L]
    res <- tibble::tibble(
      seqid = vapply(parts, function(p) p[1], character(1L)),
      start_pos = as.integer(vapply(parts, function(p) p[2], character(1L))),
      end_pos = as.integer(vapply(parts, function(p) p[3], character(1L)))
    )
    res <- res[!duplicated(res$seqid), ]
  }
  res
}

#' Read the comment header of a GFF3/GTF file
#'
#' Reads leading `#` lines (works for gzip-compressed files) and returns
#' the `##key value` pragmas and `#key: value` comments as a named list plus
#' the `##sequence-region` entries as a tibble.
#'
#' @keywords internal
#' @noRd
read_annotation_header <- function(path, chunk_size = 512L, max_lines = 1e6) {
  con <- gzfile(path, open = "rt")
  withr::defer(close(con))
  header <- character()
  total <- 0L
  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (!length(lines)) {
      break
    }
    is_comment <- startsWith(lines, "#")
    first_data <- match(FALSE, is_comment)
    if (!is.na(first_data)) {
      header <- c(header, lines[seq_len(first_data - 1L)])
      break
    }
    header <- c(header, lines)
    total <- total + length(lines)
    if (total >= max_lines) {
      break
    }
  }
  pragmas <- header[startsWith(header, "##")]
  pragmas <- sub("^##", "", pragmas)
  comments <- header[startsWith(header, "#") & !startsWith(header, "##")]
  comments <- sub("^#!?", "", comments)

  metadata <- list()
  for (line in comments) {
    key <- sub("[:=\\s].*$", "", line, perl = TRUE)
    value <- trimws(sub("^[^:=\\s]+[:=\\s]\\s*", "", line, perl = TRUE))
    if (nzchar(key) && is.null(metadata[[key]])) {
      metadata[[key]] <- value
    }
  }
  for (line in pragmas[!startsWith(pragmas, "sequence-region")]) {
    key <- sub("\\s.*$", "", line)
    value <- trimws(sub("^[^\\s]+\\s*", "", line, perl = TRUE))
    if (nzchar(key) && is.null(metadata[[key]])) {
      metadata[[key]] <- value
    }
  }
  seq_lines <- pragmas[startsWith(pragmas, "sequence-region")]
  seq_regions <- empty_annotation_table("seq_regions")
  if (length(seq_lines)) {
    parts <- strsplit(
      trimws(sub("^sequence-region\\s+", "", seq_lines)),
      "\\s+"
    )
    seq_regions <- tibble::tibble(
      seqid = vapply(parts, function(p) p[1], character(1L)),
      start_pos = as.integer(vapply(parts, function(p) p[2], character(1L))),
      end_pos = as.integer(vapply(parts, function(p) p[3], character(1L)))
    )
  }
  list(metadata = metadata, seq_regions = seq_regions)
}

# Convenience readers ---------------------------------------------------------

#' Read a Prepared Annotation into Memory
#'
#' Reads the Parquet tables of an annotation directory (or converts a
#' GFF3/GTF file first, see [prepare_isoformic_annotation()]) and returns
#' them as a list of tibbles. This is the table-based counterpart of the
#' `annot_data_*()` accessors of [IsoformicExperiment].
#'
#' @param file_path Path to an annotation directory or to a GFF3/GTF file.
#' @param file_type Passed to [prepare_isoformic_annotation()] when
#'   `file_path` is a file.
#' @param output_path Passed to [prepare_isoformic_annotation()].
#'
#' @returns A named list with tibbles `genes`, `transcripts`, `exons`,
#'   `seq_regions` and a named list `metadata`.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "isoformic")
#' annot <- prepare_annotation(gff, output_path = tempdir())
#' names(annot)
#' head(annot$genes)
#'
#' @export
prepare_annotation <- function(
  file_path,
  file_type = c("auto", "gff", "gtf"),
  output_path = NULL
) {
  file_type <- rlang::arg_match(file_type)
  if (is_annotation_dir(file_path)) {
    annot_dir <- file_path
  } else {
    annot_dir <- prepare_isoformic_annotation(
      input_path = file_path,
      output_path = output_path,
      file_type = file_type
    )
  }
  list(
    genes = read_annotation_table(annot_dir, "genes"),
    transcripts = read_annotation_table(annot_dir, "transcripts"),
    exons = read_annotation_table(annot_dir, "exons"),
    seq_regions = read_annotation_table(annot_dir, "seq_regions"),
    metadata = kv_to_list(read_annotation_table(annot_dir, "metadata"))
  )
}

#' Exon Table for Selected Genes From an Annotation
#'
#' Table-based helper returning one row per exon for the transcripts of the
#' requested genes, joined with transcript and gene attributes. This is the
#' input format of the data-frame method of [plot_genomic_context()].
#' With an [IsoformicExperiment] use [exon_table()] instead.
#'
#' @param gene_name Character vector of gene names (or gene IDs).
#' @param file_path Path to an annotation directory or a GFF3/GTF file.
#' @param file_type Passed to [prepare_isoformic_annotation()].
#' @param output_path Passed to [prepare_isoformic_annotation()].
#'
#' @returns A tibble with columns `gene_id`, `gene_name`, `transcript_id`,
#'   `transcript_name`, `transcript_type`, `exon_id`, `exon_number`, `seqid`,
#'   `start_pos`, `end_pos`, `strand`.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "isoformic")
#' exons <- prepare_exon_annotation("RBPJ", gff, output_path = tempdir())
#' head(exons)
#'
#' @export
prepare_exon_annotation <- function(
  gene_name,
  file_path,
  file_type = c("auto", "gff", "gtf"),
  output_path = NULL
) {
  annot <- prepare_annotation(
    file_path = file_path,
    file_type = file_type,
    output_path = output_path
  )
  build_exon_table(
    genes = annot$genes,
    transcripts = annot$transcripts,
    exons = annot$exons,
    gene_name = gene_name
  )
}

#' Join exons with transcript and gene attributes for selected genes
#' @keywords internal
#' @noRd
build_exon_table <- function(
  genes,
  transcripts,
  exons,
  gene_name = NULL,
  transcript_id = NULL
) {
  if (!is.null(gene_name)) {
    gene_ids <- genes$gene_id[
      genes$gene_name %in% gene_name | genes$gene_id %in% gene_name
    ]
    missing <- gene_name[
      !(gene_name %in% genes$gene_name | gene_name %in% genes$gene_id)
    ]
    if (length(missing)) {
      cli::cli_warn(c(
        "!" = "Gene{?s} not found in the annotation: {.val {missing}}."
      ))
    }
    transcripts <- transcripts[transcripts$gene_id %in% gene_ids, ]
  }
  if (!is.null(transcript_id)) {
    transcripts <- transcripts[transcripts$transcript_id %in% transcript_id, ]
  }
  tx_cols <- intersect(
    c(
      "transcript_id",
      "transcript_name",
      "transcript_type",
      "gene_id",
      "haplotype"
    ),
    colnames(transcripts)
  )
  exons <- exons[exons$transcript_id %in% transcripts$transcript_id, ]
  exons$haplotype <- NULL
  exons |>
    dplyr::left_join(
      transcripts[, tx_cols],
      by = "transcript_id"
    ) |>
    dplyr::left_join(
      genes[, c("gene_id", "gene_name")],
      by = "gene_id"
    ) |>
    dplyr::select(
      "gene_id",
      "gene_name",
      "transcript_id",
      "transcript_name",
      "transcript_type",
      dplyr::any_of("haplotype"),
      "exon_id",
      "exon_number",
      "seqid",
      "start_pos",
      "end_pos",
      "strand"
    ) |>
    dplyr::arrange(.data$gene_id, .data$transcript_id, .data$start_pos)
}
