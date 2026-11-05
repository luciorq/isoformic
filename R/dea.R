# Differential expression results ---------------------------------------------

#' Differential Expression Results in an IsoformicExperiment
#'
#' @description
#' An [IsoformicExperiment] stores differential expression results for
#' transcripts (`det`) and genes (`deg`) in its `dea` property, together
#' with the significance cutoffs used to flag features. Results can come
#' from any tool: on assignment the tables are standardised to the column
#' names `log2FC`, `pvalue` and `qvalue` (see Details), checked against the
#' annotation, and stored as tibbles.
#'
#' * `de_tx(x)` / `de_gene(x)`: results with an added `det_sig` / `deg_sig`
#'   column (`"yes"`/`"no"`) computed from the cutoffs, sorted by absolute
#'   `log2FC`. Cutoffs can be overridden per call.
#' * `de_tx(x) <- value` / `de_gene(x) <- value`: set results.
#' * `de_cutoffs(x)` / `de_cutoffs(x) <- value`: the stored cutoffs, a list
#'   with `log2fc` (default 1), `pvalue` (default 0.05) and `use_fdr`
#'   (default `TRUE`, i.e. compare `qvalue` rather than `pvalue`).
#' * `de_table(x)`: a long table combining gene and transcript results,
#'   one row per feature, with `feature_id`, `feature_name`,
#'   `feature_type` (`"gene"` or the transcript type), `gene_id`,
#'   `gene_name`, `log2FC`, `pvalue`, `qvalue` and `is_de`. This is the input
#'   of [plot_log2fc()].
#'
#' @details
#' Recognised column synonyms: `log2FoldChange`/`logFC` for `log2FC`;
#' `padj`/`FDR`/`adj.P.Val`/`svalue` for `qvalue`; `PValue`/`P.Value` for
#' `pvalue`. Row names are used as IDs when no ID column is present, and
#' GENCODE-style `|`-suffixed transcript IDs are truncated. When `qvalue`
#' is absent it is computed with Benjamini-Hochberg adjustment of `pvalue`.
#'
#' @param self An [IsoformicExperiment].
#' @param value A data frame of results (`de_tx<-`, `de_gene<-`) or a list
#'   of cutoffs (`de_cutoffs<-`).
#' @param log2fc_cutoff,pvalue_cutoff,use_fdr Optional per-call overrides
#'   of the stored cutoffs.
#' @param sig_only Logical, return only significant features.
#' @param ... Passed on to methods.
#'
#' @returns Tibbles, or the modified object for replacement functions.
#'
#' @examples
#' iso <- example_isoformic()
#' de_cutoffs(iso)
#' head(de_tx(iso))
#' head(de_gene(iso, sig_only = TRUE))
#' head(de_table(iso))
#'
#' de_cutoffs(iso) <- list(log2fc = 0.5)
#' sum(de_tx(iso)$det_sig == "yes")
#'
#' @name dea
#' @rdname dea
NULL

#' @rdname dea
#' @export
de_tx <- S7::new_generic(
  "de_tx",
  "self",
  function(
    self,
    log2fc_cutoff = NULL,
    pvalue_cutoff = NULL,
    use_fdr = NULL,
    sig_only = FALSE,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(de_tx, IsoformicExperiment) <- function(
  self,
  log2fc_cutoff = NULL,
  pvalue_cutoff = NULL,
  use_fdr = NULL,
  sig_only = FALSE,
  ...
) {
  flag_de_results(
    self,
    level = "det",
    log2fc_cutoff = log2fc_cutoff,
    pvalue_cutoff = pvalue_cutoff,
    use_fdr = use_fdr,
    sig_only = sig_only
  )
}

#' @rdname dea
#' @export
de_gene <- S7::new_generic(
  "de_gene",
  "self",
  function(
    self,
    log2fc_cutoff = NULL,
    pvalue_cutoff = NULL,
    use_fdr = NULL,
    sig_only = FALSE,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(de_gene, IsoformicExperiment) <- function(
  self,
  log2fc_cutoff = NULL,
  pvalue_cutoff = NULL,
  use_fdr = NULL,
  sig_only = FALSE,
  ...
) {
  flag_de_results(
    self,
    level = "deg",
    log2fc_cutoff = log2fc_cutoff,
    pvalue_cutoff = pvalue_cutoff,
    use_fdr = use_fdr,
    sig_only = sig_only
  )
}

#' @rdname dea
#' @export
`de_tx<-` <- S7::new_generic("de_tx<-", "self", function(self, value) {
  S7::S7_dispatch()
})

S7::method(`de_tx<-`, IsoformicExperiment) <- function(self, value) {
  dea <- self@dea
  dea[["det"]] <- standardize_de_table(value, level = "transcript")
  self@dea <- dea
  self
}

#' @rdname dea
#' @export
`de_gene<-` <- S7::new_generic("de_gene<-", "self", function(self, value) {
  S7::S7_dispatch()
})

S7::method(`de_gene<-`, IsoformicExperiment) <- function(self, value) {
  dea <- self@dea
  dea[["deg"]] <- standardize_de_table(value, level = "gene")
  self@dea <- dea
  self
}

#' @rdname dea
#' @export
de_cutoffs <- S7::new_generic("de_cutoffs", "self")

S7::method(de_cutoffs, IsoformicExperiment) <- function(self) {
  resolve_cutoffs(self@dea[["cutoffs"]])
}

#' @rdname dea
#' @export
`de_cutoffs<-` <- S7::new_generic(
  "de_cutoffs<-",
  "self",
  function(self, value) {
    S7::S7_dispatch()
  }
)

S7::method(`de_cutoffs<-`, IsoformicExperiment) <- function(self, value) {
  dea <- self@dea
  dea[["cutoffs"]] <- resolve_cutoffs(value)
  self@dea <- dea
  self
}

#' @rdname dea
#' @export
de_table <- S7::new_generic("de_table", "self")

S7::method(de_table, IsoformicExperiment) <- function(self) {
  combine_deg_det_longer(self)
}

default_cutoffs <- function() {
  list(log2fc = 1, pvalue = 0.05, use_fdr = TRUE)
}

#' Merge user cutoffs with defaults and validate them
#' @keywords internal
#' @noRd
resolve_cutoffs <- function(cutoffs = NULL) {
  res <- default_cutoffs()
  if (is.null(cutoffs)) {
    return(res)
  }
  if (!is.list(cutoffs)) {
    cli::cli_abort(
      message = c(x = "Cutoffs must be a named list."),
      class = "isoformic_invalid_cutoffs"
    )
  }
  # Backwards compatible aliases
  aliases <- c(
    log2FC = "log2fc",
    lfc = "log2fc",
    fdr = "use_fdr",
    padj = "pvalue"
  )
  nms <- names(cutoffs)
  nms[nms %in% names(aliases)] <- aliases[nms[nms %in% names(aliases)]]
  names(cutoffs) <- nms
  unknown <- setdiff(nms, names(res))
  if (length(unknown)) {
    cli::cli_abort(
      message = c(
        x = "Unknown cutoff{?s}: {.val {unknown}}.",
        i = "Valid names: {.val {names(res)}}."
      ),
      class = "isoformic_invalid_cutoffs"
    )
  }
  res[nms] <- cutoffs
  if (
    !rlang::is_scalar_double(res$log2fc) &&
      !rlang::is_scalar_integer(res$log2fc)
  ) {
    cli::cli_abort(
      c(x = "{.field log2fc} cutoff must be a single number."),
      class = "isoformic_invalid_cutoffs"
    )
  }
  if (
    !is.numeric(res$pvalue) ||
      length(res$pvalue) != 1L ||
      res$pvalue < 0 ||
      res$pvalue > 1
  ) {
    cli::cli_abort(
      c(x = "{.field pvalue} cutoff must be a single number in [0, 1]."),
      class = "isoformic_invalid_cutoffs"
    )
  }
  if (!rlang::is_scalar_logical(res$use_fdr)) {
    cli::cli_abort(
      c(x = "{.field use_fdr} must be TRUE or FALSE."),
      class = "isoformic_invalid_cutoffs"
    )
  }
  res
}

#' Standardise the `dea` list passed to the constructor
#' @keywords internal
#' @noRd
standardize_dea <- function(dea, annot_path = NULL) {
  if (is.null(dea) || !length(dea)) {
    return(list())
  }
  if (!is.list(dea)) {
    cli::cli_abort(
      message = c(
        x = "{.arg dea} must be a list with elements {.field det}, {.field deg} and/or {.field cutoffs}."
      ),
      class = "isoformic_invalid_dea"
    )
  }
  out <- list()
  if (!is.null(dea[["det"]])) {
    out[["det"]] <- standardize_de_table(dea[["det"]], level = "transcript")
  }
  if (!is.null(dea[["deg"]])) {
    out[["deg"]] <- standardize_de_table(dea[["deg"]], level = "gene")
  }
  if (!is.null(dea[["cutoffs"]])) {
    out[["cutoffs"]] <- resolve_cutoffs(dea[["cutoffs"]])
  }
  for (nm in uncertainty_table_names()) {
    if (is.data.frame(dea[[nm]])) {
      out[[nm]] <- dea[[nm]]
    }
  }
  out
}

#' Names of the inferential uncertainty tables kept in `dea`
#' @keywords internal
#' @noRd
uncertainty_table_names <- function() {
  c("uncertainty", "ambiguity", "confidence")
}

#' Standardise a differential expression table
#'
#' Harmonises column names from common tools, ensures an ID column and a
#' `qvalue` column, and returns a tibble.
#'
#' @keywords internal
#' @noRd
standardize_de_table <- function(x, level = c("transcript", "gene")) {
  level <- rlang::arg_match(level)
  id_col <- paste0(level, "_id")
  if (!is.data.frame(x)) {
    x <- rlang::try_fetch(
      as.data.frame(x, stringsAsFactors = FALSE, optional = TRUE),
      error = function(cnd) {
        cli::cli_abort(
          message = c(
            x = "Differential expression results must be a data frame."
          ),
          class = "isoformic_invalid_dea",
          parent = cnd
        )
      }
    )
  }
  nms <- colnames(x)
  synonyms <- list(
    log2FC = c("log2FoldChange", "logFC", "log2fc", "log2_fc"),
    pvalue = c("PValue", "P.Value", "pval", "p_value"),
    qvalue = c("padj", "FDR", "adj.P.Val", "svalue", "fdr", "q_value")
  )
  for (target in names(synonyms)) {
    if (!target %in% nms) {
      hit <- intersect(synonyms[[target]], nms)
      if (length(hit)) {
        nms[nms == hit[1L]] <- target
      }
    }
  }
  if (!id_col %in% nms) {
    alt <- intersect(
      c(
        "feature_id",
        "id",
        "ID",
        if (level == "transcript") c("tx_id", "transcript") else "gene"
      ),
      nms
    )
    if (length(alt)) {
      nms[nms == alt[1L]] <- id_col
    }
  }
  colnames(x) <- nms
  if (!id_col %in% colnames(x)) {
    rn <- rownames(x)
    if (is.null(rn) || identical(rn, as.character(seq_len(nrow(x))))) {
      cli::cli_abort(
        message = c(
          x = "Could not find a {.field {id_col}} column (or row names) in the {level}-level results."
        ),
        class = "isoformic_invalid_dea"
      )
    }
    x[[id_col]] <- rn
  }
  x <- tibble::as_tibble(x, .name_repair = "minimal")
  x[[id_col]] <- as.character(x[[id_col]])
  if (identical(level, "transcript")) {
    x[[id_col]] <- sub("\\|.*$", "", x[[id_col]])
  }
  check_de_columns(x, id_col, if (level == "transcript") "det" else "deg")
  if (!"qvalue" %in% colnames(x)) {
    x$qvalue <- stats::p.adjust(x$pvalue, method = "BH")
  }
  dplyr::relocate(x, dplyr::all_of(id_col))
}

#' Add a significance flag column to stored DE results
#' @keywords internal
#' @noRd
flag_de_results <- function(
  self,
  level = c("det", "deg"),
  log2fc_cutoff = NULL,
  pvalue_cutoff = NULL,
  use_fdr = NULL,
  sig_only = FALSE
) {
  level <- rlang::arg_match(level)
  res <- self@dea[[level]]
  if (is.null(res)) {
    cli::cli_abort(
      message = c(
        x = "No {if (level == 'det') 'transcript' else 'gene'}-level differential expression results in the object.",
        i = "Set them with {.code {if (level == 'det') 'de_tx' else 'de_gene'}(x) <- results}."
      ),
      class = "isoformic_missing_dea_results"
    )
  }
  cutoffs <- de_cutoffs(self)
  cutoffs$log2fc <- log2fc_cutoff %||% cutoffs$log2fc
  cutoffs$pvalue <- pvalue_cutoff %||% cutoffs$pvalue
  cutoffs$use_fdr <- use_fdr %||% cutoffs$use_fdr
  cutoffs <- resolve_cutoffs(cutoffs)
  p_col <- if (isTRUE(cutoffs$use_fdr) && "qvalue" %in% colnames(res)) {
    "qvalue"
  } else {
    "pvalue"
  }
  sig <- !is.na(res$log2FC) &
    !is.na(res[[p_col]]) &
    abs(res$log2FC) >= cutoffs$log2fc &
    res[[p_col]] <= cutoffs$pvalue
  res[[paste0(level, "_sig")]] <- ifelse(sig, "yes", "no")
  res <- res[order(-abs(res$log2FC), na.last = TRUE), ]
  if (isTRUE(sig_only)) {
    res <- res[res[[paste0(level, "_sig")]] == "yes", ]
  }
  res
}

#' Combine gene and transcript DE results in one long table
#' @keywords internal
#' @noRd
combine_deg_det_longer <- function(self) {
  dea <- self@dea
  if (is.null(dea[["det"]]) && is.null(dea[["deg"]])) {
    cli::cli_abort(
      message = c(
        x = "No differential expression results in the object.",
        i = "Set them with {.fn de_tx<-} and {.fn de_gene<-}."
      ),
      class = "isoformic_missing_dea_results"
    )
  }
  genes <- self@annot_data_genes[, c("gene_id", "gene_name")]
  parts <- list()
  if (!is.null(dea[["deg"]])) {
    parts$deg <- join_annotation(de_gene(self), genes, "gene_id") |>
      dplyr::mutate(
        feature_id = .data$gene_id,
        feature_name = .data$gene_name,
        feature_type = "gene",
        is_de = .data$deg_sig
      ) |>
      select_deg_det_cols()
  }
  if (!is.null(dea[["det"]])) {
    tx <- self@annot_data_transcripts[, c(
      "transcript_id",
      "transcript_name",
      "transcript_type",
      "gene_id"
    )]
    parts$det <- join_annotation(de_tx(self), tx, "transcript_id") |>
      join_annotation(genes, "gene_id") |>
      dplyr::mutate(
        feature_id = .data$transcript_id,
        feature_name = .data$transcript_name,
        feature_type = .data$transcript_type,
        is_de = .data$det_sig
      ) |>
      select_deg_det_cols()
  }
  dplyr::bind_rows(parts)
}

# Left join of annotation columns onto a results table; the annotation wins
# over columns of the same name that the table already carries, which would
# otherwise be suffixed with .x and .y.
join_annotation <- function(data, annotation, by) {
  overlap <- setdiff(intersect(colnames(data), colnames(annotation)), by)
  data <- data[, setdiff(colnames(data), overlap), drop = FALSE]
  dplyr::left_join(data, annotation, by = by)
}

select_deg_det_cols <- function(data) {
  dplyr::select(
    data,
    dplyr::any_of(c(
      "feature_id",
      "feature_name",
      "feature_type",
      "gene_id",
      "gene_name",
      "log2FC",
      "pvalue",
      "qvalue",
      "is_de"
    ))
  )
}

# Table-based helpers ---------------------------------------------------------

#' Merge Gene and Transcript Level Differential Expression Tables
#'
#' Table-based helper combining gene-level and transcript-level results into
#' the long format used by [plot_log2fc()]. With an [IsoformicExperiment],
#' [de_table()] produces the same table from the stored results.
#'
#' @param deg_table Gene-level results with `gene_id`, `gene_name`,
#'   `log2FC` and `pvalue` columns.
#' @param det_table Transcript-level results with `transcript_id`,
#'   `transcript_name`, `transcript_type`, `log2FC` and `pvalue` columns.
#' @param logfc_cut Absolute log2 fold-change cutoff for significance.
#' @param pval_cut P-value cutoff for significance. Like the cutoffs of
#'   [de_tx()], a feature is significant when `abs(log2FC) >= logfc_cut` and
#'   the p-value is `<= pval_cut`.
#' @param use_fdr Logical, test the `qvalue` column instead of `pvalue` (when
#'   present), as the stored cutoffs of an [IsoformicExperiment] do.
#' @param transcript_types Transcript types to keep. `NULL` keeps all.
#'
#' @returns A tibble with `feature_id`, `feature_name`, `feature_type`,
#'   `gene_name`, `log2FC`, `pvalue` (and `qvalue` when present) and `is_de`.
#'
#' @examples
#' deg_table <- data.frame(
#'   gene_id = c("gene1", "gene2"),
#'   gene_name = c("GeneA", "GeneB"),
#'   log2FC = c(1.5, -2.0),
#'   pvalue = c(0.01, 0.04)
#' )
#' det_table <- data.frame(
#'   transcript_id = c("tx1", "tx2", "tx3"),
#'   transcript_name = c("Transcript1", "Transcript2", "Transcript3"),
#'   transcript_type = c("protein_coding", "lncRNA", "processed_transcript"),
#'   gene_name = c("GeneA", "GeneA", "GeneB"),
#'   log2FC = c(1.2, -1.8, 0.5),
#'   pvalue = c(0.02, 0.03, 0.2)
#' )
#' join_DEG_DET(deg_table, det_table, logfc_cut = 1, pval_cut = 0.05)
#'
#' @export
# nolint start: object_name_linter.
join_DEG_DET <- function(
  deg_table,
  det_table,
  logfc_cut = 1,
  pval_cut = 0.05,
  use_fdr = FALSE,
  transcript_types = default_transcript_types()
) {
  deg <- tibble::as_tibble(deg_table) |>
    dplyr::mutate(
      feature_id = as.character(.data$gene_id),
      feature_name = .data$gene_name,
      feature_type = "gene"
    )
  det <- tibble::as_tibble(det_table)
  if (!is.null(transcript_types)) {
    det <- dplyr::filter(det, .data$transcript_type %in% transcript_types)
  }
  det <- det |>
    dplyr::mutate(
      feature_id = as.character(.data$transcript_id),
      feature_name = .data$transcript_name,
      feature_type = .data$transcript_type
    )
  keep <- c(
    "feature_id",
    "feature_name",
    "feature_type",
    "gene_id",
    "gene_name",
    "log2FC",
    "pvalue",
    "qvalue"
  )
  res <- dplyr::bind_rows(
    dplyr::select(deg, dplyr::any_of(keep)),
    dplyr::select(det, dplyr::any_of(keep))
  )
  p_col <- if (isTRUE(use_fdr) && "qvalue" %in% colnames(res)) {
    "qvalue"
  } else {
    "pvalue"
  }
  res |>
    dplyr::mutate(
      is_de = dplyr::if_else(
        !is.na(.data$log2FC) &
          abs(.data$log2FC) >= .env$logfc_cut &
          !is.na(.data[[p_col]]) &
          .data[[p_col]] <= .env$pval_cut,
        "yes",
        "no"
      )
    )
}
# nolint end

default_transcript_types <- function() {
  c(
    "retained_intron",
    "protein_coding_CDS_not_defined",
    "processed_transcript",
    "nonsense_mediated_decay",
    "lncRNA",
    "protein_coding",
    "pseudogene",
    "non_stop_decay",
    "processed_pseudogene",
    "transcribed_unprocessed_pseudogene",
    "transcribed_unitary_pseudogene",
    "unprocessed_pseudogene",
    "unitary_pseudogene"
  )
}

#' Flag Transcripts Whose Gene Is Differentially Expressed
#'
#' Table-based helper adding a `DEG_sig` column (`"YES"`/`"NO"`) to a
#' transcript table. Row order is preserved. With an [IsoformicExperiment],
#' [isoform_switches()] relates transcript and gene significance directly.
#'
#' @param deg_sig_vector Character vector of transcript names (or IDs)
#'   belonging to significantly differentially expressed genes.
#' @param det_table Transcript-level table with a `transcript_name` (or
#'   `transcript_id`) column.
#'
#' @returns `det_table` with an extra `DEG_sig` column.
#'
#' @examples
#' det_table <- data.frame(
#'   transcript_name = c("tx1", "tx2", "tx3"),
#'   log2FC = c(2.5, -1.2, 0.8)
#' )
#' is_deg_sig(c("tx1", "tx3"), det_table)
#'
#' @export
is_deg_sig <- function(deg_sig_vector, det_table) {
  key <- if ("transcript_name" %in% colnames(det_table)) {
    det_table$transcript_name
  } else {
    det_table$transcript_id
  }
  det_table$DEG_sig <- ifelse(key %in% deg_sig_vector, "YES", "NO")
  det_table
}
