#' Prepare Data for Gene and Transcript Expression Profile Plot
#'
#' This function processes gene and transcript-level expression data, along with differential
#' expression results, to prepare a tidy data frame suitable for plotting expression profiles
#' across different sample groups.
#'
#' @param txi_gene A `tibble` or `tximport` output containing gene-level expression abundances.
#' If `NULL`, gene-level abundances will be summarized from `txi_transcript`. Default is `NULL`.
#' @param txi_transcript A `tibble` or `tximport` output containing transcript-level expression abundances.
#' @param sample_metadata A `data.frame` or `tibble` containing sample metadata. The first column
#' should contain sample names matching the column names in `txi_gene` and `txi_transcript`.
#' @param tx_to_gene A `data.frame` or `tibble` containing transcript-to-gene mapping information.
#' Must include columns specified by `gene_col` and `tx_col`.
#' @param de_result_gene A `data.frame` or `tibble` containing differential expression results at the gene level.
#' Must include `gene_name`, `log2FC`, and `qvalue` columns.
#' @param de_result_transcript A `data.frame` or `tibble` containing differential expression results at the transcript level.
#' Must include `transcript_name`, `log2FC`, and `qvalue` columns.
#' @param var A string specifying the column name in `sample_metadata` that indicates the grouping variable (e.g., treatment, condition).
#' @param var_levels A character vector specifying the levels of `var` to include in the contrasts.
#' @param gene_col A string specifying the column name in `tx_to_gene` that contains gene names. Default is `"gene_name"`.
#' @param tx_col A string specifying the column name in `tx_to_gene` that contains transcript names. Default is `"transcript_name"`.
#' @param pvalue_cutoff A numeric value specifying the p-value cutoff for determining significant differential expression. Default is `0.05`.
#' @param lfc_cutoff A numeric value specifying the log2 fold-change cutoff for determining significant differential expression. Default is `1`.
#' @param use_fdr A logical value indicating whether to use the false discovery rate (`qvalue`) instead of p-value for significance cutoff. Default is `TRUE`.
#'
#' @return A `tibble` containing processed expression data and differential expression flags, ready for plotting.
#'
#' @details The function combines gene and transcript expression data with differential expression results to generate a tidy data frame. It filters significant genes and transcripts based on specified cutoffs and prepares the data for plotting expression profiles across specified sample groups.
#'
#' @examples
#' \dontrun{
#' # Assuming txi_gene, txi_transcript, sample_metadata, tx_to_gene, de_result_gene,
#' # and de_result_transcript are pre-loaded data frames:
#'
#' # Prepare data for plotting
#' expr_df <- prepare_profile_data(
#'   txi_gene = txi_gene,
#'   txi_transcript = txi_transcript,
#'   sample_metadata = sample_metadata,
#'   tx_to_gene = tx_to_gene,
#'   de_result_gene = de_result_gene,
#'   de_result_transcript = de_result_transcript,
#'   var = "condition",
#'   var_levels = c("control", "treatment"),
#'   gene_col = "gene_name",
#'   tx_col = "transcript_name",
#'   pvalue_cutoff = 0.05,
#'   lfc_cutoff = 1,
#'   use_fdr = TRUE
#' )
#'
#' # View the prepared data
#' head(expr_df)
#'
#' # Plotting example (assuming ggplot2 is installed)
#' library(ggplot2)
#' ggplot(expr_df, aes(x = condition, y = mean_TPM, fill = DE)) +
#'   geom_bar(stat = "identity", position = position_dodge()) +
#'   facet_wrap(~ parent_gene + transcript_type)
#' }
#'
#' @export
prepare_profile_data <- function(
    txi_gene = NULL, # txi abundance genes
    txi_transcript, # txi abundance transcripts
    sample_metadata, # metadata
    tx_to_gene, # tx2gene table (id_dictionary)
    de_result_gene, # gene level DE analysis results
    de_result_transcript, # transcript level DE analysis results
    var, # column name in sample metadata to differ groups
    var_levels, # vector of the contrasts
    gene_col = "gene_name",
    tx_col = "transcript_name",
    pvalue_cutoff = 0.05,
    lfc_cutoff = 1,
    use_fdr = TRUE) {
  # dependencies
  .data <- rlang::.data
  .env <- rlang::.env
  `:=` <- rlang::`:=`
  txi_transcript <- convert_to_isoformic_tibble(txi_transcript)
  # Extract tx2gene from gene annotation table
  # renamed gene annotation to gene_metadata
  # TODO: should gene_metadata be renamed to transcript_metadata?
  # gene_metadata <- tibble::as_tibble(tx_to_gene)
  gene_metadata <- tx_to_gene
  # tx2gene_df <- gene_metadata |>
  #   dplyr::filter(.data[[gene_col]] %in% genes_to_plot) |>
  #   dplyr::select({{ tx_col }}, {{ gene_col }}) |>
  #   dplyr::rename(TXNAME = {{ tx_col }}, GENEID = {{ gene_col }})

  if (isFALSE("transcript_name" %in% colnames(txi_transcript))) {
    txi_transcript <- txi_transcript |>
      dplyr::left_join(
        y = dplyr::select(tx_to_gene, "transcript_id", "transcript_name"),
        by = "transcript_id"
      )
  }

  # Format sample metadata
  # filter only samples included in var_levels
  # First is used as sample names
  sample_col <- colnames(dplyr::select(sample_metadata, 1))
  condition_df <- sample_metadata |>
    dplyr::select({{ sample_col }}, {{ var }}) |>
    dplyr::filter(.data[[var]] %in% var_levels) |>
    tibble::as_tibble()

  # if txi_gene not provided summarize from txi_transcriopt
  if (is.null(txi_gene)) {
    txi_gene <- summarize_to_gene(
      txi_transcript = txi_transcript,
      tx_to_gene = tx_to_gene
    )
  }
  # Format expression data
  # TODO should be replaced with a function to check if it is `tximport` output
  # only used with tximport direct output
  if (!tibble::is_tibble(txi_gene)) {
    if (names(txi_gene) %in% "abundance") {
      txi_gene <- txi_gene$abundance
    }
    txi_gene <- txi_gene |>
      base::as.data.frame() |>
      tibble::rownames_to_column("genename") |>
      tibble::as_tibble()
  }
  # TODO: should be replaced with a function to check if it is `tximport` output
  # only used with tximport direct output
  if (!tibble::is_tibble(txi_transcript)) {
    if (names(txi_transcript) %in% "abundance") {
      txi_transcript <- txi_transcript$abundance
    }
    txi_transcript <- txi_transcript |>
      base::as.data.frame() |>
      tibble::rownames_to_column("genename") |>
      tibble::as_tibble()
  }



  # Format gene names and transcript names in expression tables
  geneid_col <- colnames(gene_metadata)[2]
  # txid_col <- colnames(gene_metadata)[1]
  first_col <- colnames(txi_gene)[1]
  genename_conversion_df <- tibble::as_tibble(
    dplyr::select(
      gene_metadata, {{ geneid_col }}, {{ gene_col }}
    )
  )

  # txname_conversion_df <- tibble::as_tibble(
  #  dplyr::select(
  #    gene_metadata, {{ txid_col }}, {{ tx_col }}
  #  )
  # )

  txi_gene <- txi_gene |>
    dplyr::left_join(genename_conversion_df, by = c(gene_id = geneid_col)) |>
    dplyr::select(-c({{ first_col }})) |>
    dplyr::rename(genename = {{ gene_col }}) |>
    dplyr::relocate("genename") |>
    dplyr::filter(!is.na(.data$genename)) |>
    dplyr::distinct()

  first_col <- colnames(txi_transcript)[1]
  # names(txi_transcript)
  # txi_transcript <- txname_conversion_df |>
  #   dplyr::left_join(txi_transcript, by = c("transcript_id")) |>
  txi_transcript <- txi_transcript |>
    dplyr::select(-c({{ first_col }})) |>
    dplyr::rename(genename = "transcript_name") |>
    dplyr::relocate("genename") |>
    dplyr::filter(!is.na(.data$genename)) |>
    dplyr::distinct()

  # Prepare data to plot
  gene_expr_df <- txi_gene |>
    tidyr::pivot_longer(
      -c("genename"),
      names_to = "sample",
      values_to = "TPM"
    ) |>
    dplyr::left_join(condition_df, by = c("sample" = sample_col)) |>
    dplyr::group_by(.data$genename, .data[[var]]) |>
    dplyr::summarise(
      mean_TPM = base::mean(.data$TPM),
      SD = stats::sd(.data$TPM)
    ) |>
    dplyr::ungroup()

  # Add differential expression information
  de_genes_df <- de_result_gene |>
    tibble::as_tibble() |>
    dplyr::filter(.data$qvalue <= .env$pvalue_cutoff) |>
    dplyr::filter(abs(.data$log2FC) >= .env$lfc_cutoff) |>
    dplyr::select("gene_name") |>
    dplyr::mutate(`DE` = TRUE)

  # filter DE genes
  gene_expr_df <- gene_expr_df |>
    dplyr::mutate(transcript_type = "gene") |>
    dplyr::left_join(de_genes_df, by = c(genename = "gene_name")) |>
    dplyr::distinct()

  gene_expr_df <- gene_expr_df |>
    dplyr::mutate(parent_gene = "genename")

  # Add differential expression information for tx
  transcript_expr_df <- txi_transcript |>
    tidyr::pivot_longer(
      -c("genename"),
      names_to = "sample",
      values_to = "TPM"
    ) |>
    dplyr::left_join(condition_df, by = c("sample" = sample_col)) |>
    dplyr::group_by(.data$genename, .data[[var]]) |>
    dplyr::summarise(
      mean_TPM = base::mean(.data$TPM),
      SD = stats::sd(.data$TPM)
    ) |>
    dplyr::ungroup()

  de_tx_df <- de_result_transcript |>
    tibble::as_tibble() |>
    dplyr::select("transcript_name", "qvalue", "log2FC") |>
    dplyr::filter(.data$qvalue <= .env$pvalue_cutoff) |>
    dplyr::filter(abs(.data$log2FC) >= .env$lfc_cutoff) |>
    dplyr::select(1) |>
    dplyr::mutate(DE = TRUE)

  ## Add transcript type information
  tx_type_df <- tibble::as_tibble(
    dplyr::select(
      gene_metadata,
      {{ gene_col }},
      {{ tx_col }},
      "transcript_type"
    )
  )

  transcript_expr_df <- transcript_expr_df |>
    dplyr::left_join(tx_type_df, by = c(genename = tx_col))
  transcript_expr_df <- transcript_expr_df |>
    dplyr::rename(parent_gene = {{ gene_col }})
  # filter DE transcripts
  transcript_expr_df <- transcript_expr_df |>
    dplyr::left_join(de_tx_df, by = c("genename" = "transcript_name")) |>
    dplyr::distinct()

  # unite gene and transcript expression summaries by condition
  expr_df <- dplyr::bind_rows(gene_expr_df, transcript_expr_df)

  # relevel for x axis Plot order
  expr_df <- expr_df |>
    dplyr::mutate({{ var }} := forcats::as_factor(.data[[var]])) |>
    dplyr::mutate({{ var }} := forcats::fct_relevel(.data[[var]], var_levels))

  expr_df <- expr_df |>
    dplyr::mutate(DE = dplyr::if_else(is.na(.data$DE), "No", "Yes")) |>
    dplyr::mutate(DE = factor(.data$DE, levels = c("Yes", "No")))
  return(expr_df)
}
