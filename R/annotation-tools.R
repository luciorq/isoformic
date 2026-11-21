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
#' The tables are copied by DuckDB from the Parquet files of each haplotype
#' to the merged files, so they are not read into R: merging two annotations
#' of a whole human release with their transcript sequences takes seconds
#' and less than 2 GB of memory. The sequences are copied by one DuckDB
#' thread, which keeps their memory to a few hundred megabytes.
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
#' gff <- system.file("extdata", "example.gff3.gz", package = "isoformic")
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
  dirs <- vapply(
    names(annot_paths),
    function(hap) {
      path <- annot_paths[[hap]]
      if (is_annotation_dir(path)) {
        return(as.character(fs::path_abs(path)))
      }
      prepare_isoformic_annotation(
        path,
        output_path = output_path,
        haplotype = hap
      )
    },
    character(1L)
  )
  fs::dir_create(merged_dir)
  suffixes <- if (isTRUE(id_suffix)) paste0("_", names(dirs)) else NULL
  # Each table is copied by DuckDB from the Parquet files of the haplotypes,
  # so that no table (the sequences of a whole release take gigabytes) is
  # read into R
  con <- duckdb_connect()
  for (type in c("genes", "transcripts", "exons", "seq_regions", "sequences")) {
    files <- as.character(annot_file_path(dirs, type))
    present <- fs::file_exists(files)
    # Sequences are kept only when every haplotype has them; a haplotype
    # without the (optional) sequence region table is left out of it
    if (identical(type, "sequences") && !all(present)) {
      next
    }
    if (!any(present)) {
      empty <- empty_annotation_table(type)
      empty$haplotype <- character()
      write_parquet(empty, annot_file_path(merged_dir, type))
      next
    }
    merge_annotation_table(
      con,
      files = files[present],
      haplotypes = names(dirs)[present],
      suffixes = suffixes[present],
      type = type,
      output_file = as.character(annot_file_path(merged_dir, type))
    )
  }
  part_metadata <- lapply(names(dirs), function(hap) {
    meta <- kv_to_list(read_annotation_table(dirs[[hap]], "metadata"))
    stats::setNames(as.list(unlist(meta)), paste0(hap, ".", names(meta)))
  })
  metadata <- c(
    list(
      annot_name = annot_name,
      haplotypes = paste(names(annot_paths), collapse = ";"),
      id_suffix = as.character(isTRUE(id_suffix)),
      created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      isoformic_version = as.character(utils::packageVersion("isoformic"))
    ),
    unlist(part_metadata, recursive = FALSE)
  )
  write_parquet(list_to_kv(metadata), annot_file_path(merged_dir, "metadata"))
  merged_dir
}

#' Copy one table of several haplotype annotations into one Parquet file
#'
#' One `COPY` of a `UNION ALL BY NAME` of the haplotype files, in the order
#' of `files`: each row keeps its position (DuckDB preserves the insertion
#' order), the `haplotype` column is set (added at the end when a file has
#' none, as for the sequence regions) and the identifiers get the suffix of
#' their haplotype. Columns missing from some files are filled with `NULL`,
#' like [dplyr::bind_rows()].
#' @keywords internal
#' @noRd
merge_annotation_table <- function(
  con,
  files,
  haplotypes,
  suffixes,
  type,
  output_file
) {
  id_columns <- switch(
    type,
    genes = "gene_id",
    transcripts = c("transcript_id", "gene_id"),
    exons = c("exon_id", "transcript_id"),
    sequences = "transcript_id",
    character()
  )
  with_haplotype <- !identical(type, "sequences")
  selects <- vapply(
    seq_along(files),
    function(i) {
      file <- DBI::dbQuoteString(con, files[[i]])
      columns <- colnames(DBI::dbGetQuery(
        con,
        paste0("SELECT * FROM read_parquet(", file, ") LIMIT 0")
      ))
      if (with_haplotype && !"haplotype" %in% columns) {
        columns <- c(columns, "haplotype")
      }
      exprs <- vapply(
        columns,
        function(column) {
          name <- as.character(DBI::dbQuoteIdentifier(con, column))
          if (with_haplotype && identical(column, "haplotype")) {
            return(paste(DBI::dbQuoteString(con, haplotypes[[i]]), "AS", name))
          }
          if (column %in% id_columns && length(suffixes)) {
            # A missing identifier stays missing
            return(paste0(
              "CAST(",
              name,
              " AS VARCHAR) || ",
              DBI::dbQuoteString(con, suffixes[[i]]),
              " AS ",
              name
            ))
          }
          name
        },
        character(1L)
      )
      paste0(
        "SELECT ",
        paste(exprs, collapse = ", "),
        " FROM read_parquet(",
        file,
        ")"
      )
    },
    character(1L)
  )
  # An order-preserving COPY buffers whole row groups; with the default of
  # 122,880 rows a row group of sequences takes hundreds of megabytes.
  # Several threads also buffer the batches that are read out of order: for
  # the sequences of two GENCODE 50 annotations (3 GB of strings) that is
  # 6 GB with four threads and 0.5 GB with one, which is about as fast as
  # the I/O allows
  row_group <- ""
  if (identical(type, "sequences")) {
    row_group <- ", ROW_GROUP_SIZE 20000"
    threads <- DBI::dbGetQuery(con, "SELECT current_setting('threads') AS n")$n
    DBI::dbExecute(con, "SET threads = 1")
    on.exit(
      DBI::dbExecute(con, paste0("SET threads = ", as.integer(threads))),
      add = TRUE
    )
  }
  DBI::dbExecute(
    con,
    paste0(
      "COPY (",
      paste(selects, collapse = " UNION ALL BY NAME "),
      ") TO ",
      DBI::dbQuoteString(con, output_file),
      " (FORMAT parquet, COMPRESSION zstd",
      row_group,
      ")"
    )
  )
  invisible(output_file)
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
#' `transcript_type` column of an annotation directory (and, when the
#' columns are present, the associated reference transcript into
#' `source_transcript_id` and the associated gene into the `gene_name` of
#' genes without a name), so the SQANTI categories behave
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
#' gff <- system.file("extdata", "example.gff3.gz", package = "isoformic")
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
