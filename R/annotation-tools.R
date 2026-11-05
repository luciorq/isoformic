# Tools that modify or combine annotation directories ------------------------

#' Merge Per-Haplotype Annotations Into One Annotation Directory
#'
#' @description
#' Phased (trio-binned or pangenome) assemblies come with one annotation per
#' haplotype, for example the HPRC CAT files `HG00438.1.chm13.gff3.gz` and
#' `HG00438.2.chm13.gff3.gz`. `merge_annotations()` converts each file (or
#' reuses an existing annotation directory), tags every table with the
#' haplotype label, makes identifiers unique by appending `_<haplotype>`, and
#' writes one combined annotation directory. Transcript, gene and exon IDs
#' therefore become `<id>_<haplotype>`, the convention expected for the row
#' names of a haplotype-resolved [IsoformicExperiment] (rename the sequences
#' of the per-haplotype transcript FASTA files the same way before building
#' the Salmon index).
#'
#' The `source_gene_id` / `source_transcript_id` columns (GRCh38 identifiers
#' kept by Liftoff/CAT) are left untouched, which is what
#' [haplotype_summary()] uses to compare haplotypes of the same gene.
#'
#' @param annot_paths Named character vector: names are haplotype labels,
#'   values are GFF3/GTF files or annotation directories.
#' @param output_path Directory under which the merged annotation directory
#'   is created. `NULL` uses [get_isoformic_cache()].
#' @param annot_name Name of the merged annotation directory.
#' @param id_suffix Logical. Append `_<haplotype>` to gene, transcript and
#'   exon identifiers (default `TRUE`). Set to `FALSE` when the identifiers
#'   are already unique across haplotypes.
#' @param overwrite Logical. Replace an existing merged directory.
#'
#' @returns The path to the merged annotation directory.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "isoformic")
#' merged <- merge_annotations(
#'   c(hap1 = gff, hap2 = gff),
#'   output_path = tempdir(),
#'   annot_name = "example_diploid",
#'   overwrite = TRUE
#' )
#' tx <- read_parquet(file.path(merged, "annot_data_transcripts.parquet"))
#' table(tx$haplotype)
#' head(tx$transcript_id)
#'
#' @export
merge_annotations <- function(
  annot_paths,
  output_path = NULL,
  annot_name = "merged_haplotypes",
  id_suffix = TRUE,
  overwrite = FALSE
) {
  if (
    !is.character(annot_paths) ||
      length(annot_paths) < 2L ||
      is.null(names(annot_paths)) ||
      !all(nzchar(names(annot_paths))) ||
      anyDuplicated(names(annot_paths))
  ) {
    cli::cli_abort(
      message = c(
        x = "{.arg annot_paths} must be a named character vector with at least two unique haplotype labels."
      ),
      class = "isoformic_invalid_haplotypes"
    )
  }
  check_string(annot_name)
  if (is.null(output_path)) {
    output_path <- get_isoformic_cache()
  }
  merged_dir <- as.character(fs::path_abs(fs::path(output_path, annot_name)))
  if (isFALSE(overwrite) && is_annotation_dir(merged_dir)) {
    cli::cli_inform(c(
      i = "Using cached annotation directory {.path {merged_dir}}."
    ))
    return(merged_dir)
  }
  parts <- lapply(names(annot_paths), function(hap) {
    path <- annot_paths[[hap]]
    dir <- if (is_annotation_dir(path)) {
      path
    } else {
      prepare_isoformic_annotation(
        path,
        output_path = output_path,
        haplotype = hap
      )
    }
    annot <- prepare_annotation(dir)
    suffix <- if (isTRUE(id_suffix)) paste0("_", hap) else ""
    annot$genes$haplotype <- hap
    annot$transcripts$haplotype <- hap
    annot$exons$haplotype <- hap
    annot$genes$gene_id <- paste0(annot$genes$gene_id, suffix)
    annot$transcripts$gene_id <- paste0(annot$transcripts$gene_id, suffix)
    annot$transcripts$transcript_id <- paste0(
      annot$transcripts$transcript_id,
      suffix
    )
    annot$exons$transcript_id <- paste0(annot$exons$transcript_id, suffix)
    annot$exons$exon_id <- ifelse(
      is.na(annot$exons$exon_id),
      NA_character_,
      paste0(annot$exons$exon_id, suffix)
    )
    annot$seq_regions$haplotype <- hap
    if (has_sequences(dir)) {
      seqs <- read_annotation_table(dir, "sequences")
      seqs$transcript_id <- paste0(seqs$transcript_id, suffix)
      annot$sequences <- seqs
    }
    annot$metadata <- stats::setNames(
      as.list(unlist(annot$metadata)),
      paste0(hap, ".", names(annot$metadata))
    )
    annot
  })
  fs::dir_create(merged_dir)
  bind <- function(name) dplyr::bind_rows(lapply(parts, function(p) p[[name]]))
  write_parquet(bind("genes"), annot_file_path(merged_dir, "genes"))
  write_parquet(bind("transcripts"), annot_file_path(merged_dir, "transcripts"))
  write_parquet(bind("exons"), annot_file_path(merged_dir, "exons"))
  write_parquet(bind("seq_regions"), annot_file_path(merged_dir, "seq_regions"))
  if (all(vapply(parts, function(p) !is.null(p$sequences), logical(1L)))) {
    write_parquet(bind("sequences"), annot_file_path(merged_dir, "sequences"))
  }
  metadata <- c(
    list(
      annot_name = annot_name,
      haplotypes = paste(names(annot_paths), collapse = ";"),
      id_suffix = as.character(isTRUE(id_suffix)),
      created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      isoformic_version = as.character(utils::packageVersion("isoformic"))
    ),
    unlist(lapply(parts, function(p) p$metadata), recursive = FALSE)
  )
  write_parquet(list_to_kv(metadata), annot_file_path(merged_dir, "metadata"))
  merged_dir
}

#' Haplotype labels of an annotation (empty when not haplotype-resolved)
#' @keywords internal
#' @noRd
annotation_haplotypes <- function(transcripts) {
  if (!"haplotype" %in% colnames(transcripts)) {
    return(character())
  }
  unique(stats::na.omit(transcripts$haplotype))
}

#' Add SQANTI3 Structural Categories to an Annotation
#'
#' @description
#' Long-read transcript catalogs (IsoSeq, FLAIR, Bambu, StringTie) carry no
#' biotype for novel transcripts. SQANTI3 classifies every isoform against a
#' reference (`*_classification.txt`: `isoform`, `structural_category`,
#' `associated_gene`, `associated_transcript`, `subcategory`, ...).
#' `add_sqanti_classification()` writes those categories into the
#' `transcript_type` column of an annotation directory (and the associated
#' reference transcript / gene into `source_transcript_id`, `source_gene_id`
#' and `gene_name` when they are missing), so the SQANTI categories behave
#' like transcript types everywhere in the package (plots, enrichment by
#' category, [tx_diversity()]). The categories are part of
#' [tx_type_palette()].
#'
#' @param x An annotation directory path or an [IsoformicExperiment] (its
#'   annotation directory is modified in place and caches refresh
#'   automatically).
#' @param classification Path to the SQANTI3 classification file (tab
#'   separated, optionally gzipped) or a data frame with at least the
#'   `isoform` and `structural_category` columns.
#' @param overwrite_types Logical. Replace existing transcript types
#'   (default `TRUE`); with `FALSE` only missing types are filled.
#'
#' @returns The annotation directory path, invisibly.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "isoformic")
#' annot <- prepare_isoformic_annotation(gff, output_path = tempdir(), annot_name = "sqanti_demo")
#' tx <- read_parquet(file.path(annot, "annot_data_transcripts.parquet"))
#' classification <- data.frame(
#'   isoform = tx$transcript_id[1:3],
#'   structural_category = c("full-splice_match", "novel_in_catalog", "novel_not_in_catalog"),
#'   associated_gene = "RBPJ",
#'   associated_transcript = c(tx$transcript_id[1], "novel", "novel")
#' )
#' add_sqanti_classification(annot, classification)
#' tx <- read_parquet(file.path(annot, "annot_data_transcripts.parquet"))
#' head(tx[, c("transcript_id", "transcript_type")], 4)
#'
#' @export
add_sqanti_classification <- function(
  x,
  classification,
  overwrite_types = TRUE
) {
  annot_dir <- if (S7::S7_inherits(x, IsoformicExperiment)) x@annot_path else x
  if (!is_annotation_dir(annot_dir)) {
    cli::cli_abort(
      message = c(
        x = "{.arg x} must be an annotation directory or an {.cls IsoformicExperiment}."
      ),
      class = "isoformic_annot_path_dont_exist"
    )
  }
  if (!is.data.frame(classification)) {
    classification <- read_tabular(classification)
  }
  required <- c("isoform", "structural_category")
  missing_cols <- setdiff(required, colnames(classification))
  if (length(missing_cols)) {
    cli::cli_abort(
      message = c(
        x = "SQANTI classification is missing column{?s}: {.field {missing_cols}}."
      ),
      class = "isoformic_invalid_sqanti"
    )
  }
  classification <- tibble::as_tibble(classification)
  classification$isoform <- as.character(classification$isoform)
  classification <- dplyr::distinct(
    classification,
    .data$isoform,
    .keep_all = TRUE
  )

  tx <- read_annotation_table(annot_dir, "transcripts")
  idx <- match(tx$transcript_id, classification$isoform)
  hit <- !is.na(idx)
  n_hit <- sum(hit)
  if (!n_hit) {
    cli::cli_abort(
      message = c(
        x = "None of the {nrow(classification)} classified isoforms match the annotation's transcript IDs."
      ),
      class = "isoformic_invalid_sqanti"
    )
  }
  category <- as.character(classification$structural_category[idx[hit]])
  if (isTRUE(overwrite_types)) {
    tx$transcript_type[hit] <- category
  } else {
    fill <- hit & is.na(tx$transcript_type)
    tx$transcript_type[
      fill
    ] <- as.character(classification$structural_category[idx[fill]])
  }
  if ("associated_transcript" %in% colnames(classification)) {
    assoc <- as.character(classification$associated_transcript[idx[hit]])
    assoc[assoc %in% c("novel", "NA", "")] <- NA_character_
    tx$source_transcript_id[hit] <- ifelse(
      is.na(assoc),
      tx$source_transcript_id[hit],
      assoc
    )
  }
  genes <- read_annotation_table(annot_dir, "genes")
  if ("associated_gene" %in% colnames(classification)) {
    assoc_gene <- as.character(classification$associated_gene[idx[hit]])
    assoc_gene[
      assoc_gene %in% c("novel", "NA", "") | grepl("^novelGene", assoc_gene)
    ] <- NA_character_
    gene_ids <- tx$gene_id[hit]
    lookup <- stats::setNames(assoc_gene, gene_ids)
    lookup <- lookup[!is.na(lookup) & !duplicated(names(lookup))]
    no_name <- is.na(genes$gene_name) | genes$gene_name == genes$gene_id
    replace <- no_name & genes$gene_id %in% names(lookup)
    genes$gene_name[replace] <- unname(lookup[genes$gene_id[replace]])
    write_parquet(genes, annot_file_path(annot_dir, "genes"))
  }
  write_parquet(tx, annot_file_path(annot_dir, "transcripts"))
  cli::cli_inform(c(
    v = "Assigned SQANTI categories to {n_hit} of {nrow(tx)} transcripts in {.path {annot_dir}}."
  ))
  invisible(annot_dir)
}
