#' Prepare Data for Profile Plot
#' @return a `tibble` with data to plot.
#' @export
prepare_profile_data <- function(
  txi_gene = NULL, # txi abundance genes
  txi_transcript, # txi abundance transcripts
  sample_metadata, # metadata
  tx_to_gene, # tx2gene table (id_dictionary)
  de_result_gene, # gene level DE analysis results
  de_result_transcript, #transcript level DE analysis results
  var, # column name in sample metadata to differ groups
  var_levels, # vector of the contrasts
  gene_col = "gene_name",
  tx_col = "transcript_name",
  pvalue_cutoff = 0.05,
  lfc_cutoff = 1,
  use_fdr = TRUE
) {
  # dependencies
  .data <- rlang::.data
  `%>%` <- dplyr::`%>%`

  # Extract tx2gene from gene annotation table
  # renamed gene annotation to gene_metadata
  # TODO: should gene_metadata be renamed to transcript_metadata?
  # gene_metadata <- tibble::as_tibble(tx_to_gene)
  gene_metadata <- tx_to_gene
  # tx2gene_df <- gene_metadata %>%
  #   dplyr::filter(.data[[gene_col]] %in% genes_to_plot) %>%
  #   dplyr::select({{ tx_col }}, {{ gene_col }}) %>%
  #   dplyr::rename(TXNAME = {{ tx_col }}, GENEID = {{ gene_col }})

  # Format sample metadata
  # filter only samples included in var_levels
  # First is used as sample names
  sample_col <- colnames(dplyr::select(sample_metadata, 1))
  condition_df <- sample_metadata %>%
    dplyr::select({{ sample_col }}, {{ var }}) %>%
    dplyr::filter(.data[[var]] %in% var_levels) %>%
    tibble::as_tibble()

  # if txi_gene not provided summarize from txi_transcriopt
  if (is.null(txi_gene)) {
    txi_gene = summarize_to_gene(
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
    txi_gene <- txi_gene %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column("genename") %>%
      tibble::as_tibble()
  }
  # TODO: should be replaced with a function to check if it is `tximport` output
  # only used with tximport direct output
  if (!tibble::is_tibble(txi_transcript)) {
    if (names(txi_transcript) %in% "abundance") {
      txi_transcript <- txi_transcript$abundance
    }
    txi_transcript <- txi_transcript %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column("genename") %>%
      tibble::as_tibble()
  }



  # Format gene names and transcript names in expression tables
  geneid_col <- colnames(gene_metadata)[2]
  txid_col <- colnames(gene_metadata)[1]
  first_col <- colnames(txi_gene)[1]
  genename_conversion_df <- tibble::as_tibble(
    dplyr::select(
      gene_metadata, {{ geneid_col }}, {{ gene_col }}
    )
  )

  txname_conversion_df <- tibble::as_tibble(
    dplyr::select(
      gene_metadata, {{ txid_col }}, {{ tx_col }})
  )

  txi_gene <- txi_gene %>%
    dplyr::left_join(genename_conversion_df, by = c(gene_id = geneid_col)) %>%
    dplyr::select(-c({{ first_col }})) %>%
    dplyr::rename(genename = {{ gene_col }}) %>%
    dplyr::relocate(genename) %>%
    dplyr::filter(!is.na(genename)) %>%
    dplyr::distinct()

  first_col <- colnames(txi_transcript)[1]
  #names(txi_transcript)
  # txi_transcript <- txname_conversion_df %>%
  #   dplyr::left_join(txi_transcript, by = c("transcript_id")) %>%
  txi_transcript <- txi_transcript |>
    dplyr::select(-c({{ first_col }})) |>
    dplyr::rename(genename = transcript_name) |>
    dplyr::relocate(genename) |>
    dplyr::filter(!is.na(genename)) |>
    dplyr::distinct()

  # Prepare data to plot
  gene_expr_df <- txi_gene %>%
    tidyr::pivot_longer(-genename, names_to = "sample", values_to = "TPM") %>%
    dplyr::left_join(condition_df, by = c("sample" = sample_col)) %>%
    dplyr::group_by(genename, .data[[var]]) %>%
    dplyr::summarise(mean_TPM = mean(TPM), SD = sd(TPM)) %>%
    dplyr::ungroup()

  # Add differential expression information
  de_genes_df <- de_result_gene %>%
    tibble::as_tibble() %>%
    dplyr::filter(qvalue <= pvalue_cutoff) %>%
    dplyr::filter(abs(log2FC) >= lfc_cutoff) %>%
    dplyr::select(gene_name) %>%
    dplyr::mutate(DE = TRUE)



  # filter DE genes
  gene_expr_df <- gene_expr_df %>%
    dplyr::mutate(transcript_type = "gene") %>%
    dplyr::left_join(de_genes_df, by = c(genename = "gene_name")) %>%
    dplyr::distinct()

  gene_expr_df <- gene_expr_df %>%
    dplyr::mutate(parent_gene = genename)

  # Add differential expression information for tx
  transcript_expr_df <- txi_transcript %>%
    tidyr::pivot_longer(-genename, names_to = "sample", values_to = "TPM") %>%
    dplyr::left_join(condition_df, by = c("sample" = sample_col)) %>%
    dplyr::group_by(genename, .data[[var]]) %>%
    dplyr::summarise(mean_TPM = mean(TPM), SD = sd(TPM)) %>%
    dplyr::ungroup()

  de_tx_df <- de_result_transcript %>%
    tibble::as_tibble() %>%
    dplyr::select(transcript_name, qvalue, log2FC) %>%
    dplyr::filter(qvalue <= pvalue_cutoff) %>%
    dplyr::filter(abs(log2FC) >= lfc_cutoff) %>%
    dplyr::select(1) %>%
    dplyr::mutate(DE = TRUE)

  ## Add transcript type information
  tx_type_df <- tibble::as_tibble(dplyr::select(gene_metadata, {{ gene_col }}, {{ tx_col }}, transcript_type))

  transcript_expr_df <- transcript_expr_df %>%
    dplyr::left_join(tx_type_df, by = c(genename = tx_col))
  transcript_expr_df <- transcript_expr_df |>
    dplyr::rename(parent_gene = {{ gene_col }})
  # filter DE transcripts
  transcript_expr_df <- transcript_expr_df %>%
    dplyr::left_join(de_tx_df, by = c("genename" = "transcript_name")) %>%
    dplyr::distinct()

  # unite gene and transcript expression summaries by condition
  expr_df <- dplyr::bind_rows(gene_expr_df, transcript_expr_df)

  # relevel for x axis Plot order
  expr_df <- expr_df %>%
    dplyr::mutate({{ var }} := forcats::as_factor(.data[[var]])) %>%
    dplyr::mutate({{ var }} := forcats::fct_relevel(.data[[var]], var_levels))

  expr_df <- expr_df %>%
    dplyr::mutate(DE = dplyr::if_else(is.na(DE), "No", "Yes")) %>%
    dplyr::mutate(DE = factor(DE, levels = c("Yes", "No")))
  return(expr_df)
}


# Utils
summarize_to_gene <- function(txi_transcript, tx_to_gene) {
  id_df <- tx_to_gene |>
    dplyr::select(transcript_id, gene_id) |>
    dplyr::distinct()

  txi_gene <- id_df |>
    dplyr::left_join(txi_transcript, by = "transcript_id") |>
    dplyr::select(-c(transcript_id, transcript_name)) |>
    tidyr::pivot_longer(
      -gene_id,
      names_to = "samples",
      values_to = "expr"
    ) |>
    dplyr::group_by(gene_id, samples) |>
    dplyr::summarise(
      mean_expr = base::mean(expr, na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(mean_expr))
  txi_gene <- txi_gene |>
    tidyr::pivot_wider(
      names_from = "samples",
      values_from = "mean_expr"
    )
  return(txi_gene)
}
