# This is the big script with all the functions
# TODO: Temporary file, move functions to their own files

#' Create Transcript to Gene Relationship Table
#'
#' This function can extract a transcript to gene relationship table
#'   from GENCODE annotation files, such as the transcriptome FASTA and
#'   GFF3 or GTF annotation files.
#'
#' @param file_path Path to file containing the reference.
#'
#' @inheritParams download_reference
#' @export
make_tx_to_gene <- function(file_path, file_type = "fasta") {
  if (!isTRUE(fs::file_exists(file_path))) {
    cli::cli_abort(
      c(`x` = "{.path {file_path}} do not exist.")
    )
  }
  if (isTRUE(file_type == "fasta")) {
    fasta_lines <- readr::read_lines(file_path)
    vector_detect <- stringr::str_detect(fasta_lines, "^>.")
    header_table <- fasta_lines[vector_detect]
    header_table <- header_table |>
      tibble::as_tibble() |>
      tidyr::separate(
        value,
        into = paste0("col_", 1:9),
        sep = "\\|"
      ) |>
      dplyr::select(-c("col_9"))
    header_table$col_1 <- stringr::str_remove(header_table$col_1, ">")
  } else if (isTRUE(file_type == "gff")) {
    # TODO: Actually implement parsing GFF3 and GTF files
    cli::cli_abort(
      c(
        `x` = "{.var GTF} and {.var GFF3} files are not supported yet.",
        `!` = "Use the GENCODE {.var FASTA} file."
      )
    )
  }
  names(header_table) <- c(
    "transcript_id", "gene_id", "havanna_gene_id",
    "havanna_transcript_id", "transcript_name",
    "gene_name", "entrez_id", "transcript_type"
  )
  return(header_table)
}

#' Add Differential Gene Expression Results to Transcript Table
#' @export
is_deg_sig <- function(DegsigVector, DET_table) {
  DETs_DEGs <- DET_table |>
    dplyr::filter(transcript_name %in% DegsigVector)
  DETs_DEGs$DEG_sig <- "YES"
  DETs_notDEGs <- DET_table |>
    dplyr::filter(!transcript_name %in% DegsigVector)
  DETs_notDEGs$DEG_sig <- "NO"
  DET_table_final <- dplyr::bind_rows(DETs_DEGs, DETs_notDEGs)
  return(DET_table_final)
}

#' Merge Gene and Transcript Level Differential Expression Tables
#' @export
join_DEG_DET <- function(DEG_tab, DET_final_tab, logfc_cut, pval_cut) {
  DEG_tab_mod <- DEG_tab |>
    dplyr::rename(id = gene_id)
  DEG_tab_mod <- DEG_tab_mod |>
    dplyr::rename(name = gene_name)
  DEG_tab_mod <- DEG_tab_mod |>
    dplyr::mutate(transcript_type = "gene")
  DEG_tab_mod <- DEG_tab_mod |>
    dplyr::mutate(gene_name = DEG_tab_mod$name)
  drop_columns <- c("DEG_sig")
  DET_final_tab <- DET_final_tab |>
    dplyr::filter(
      transcript_type %in% c(
        "retained_intron",
        "processed_transcript",
        "nonsense_mediated_decay",
        "lncRNA", "protein_coding",
        "pseudogene", "non_stop_decay",
        "processed_pseudogene",
        "transcribed_unprocessed_pseudogene",
        "transcribed_unitary_pseudogene",
        "unprocessed_pseudogene",
        "unitary_pseudogene"
      )
    )
  DET_final_tab_mod <- DET_final_tab |>
    dplyr::select(-dplyr::one_of(drop_columns))
  DET_final_tab_mod <- DET_final_tab_mod |>
    dplyr::rename(id = transcript_id)
  DET_final_tab_mod <- DET_final_tab_mod |>
    dplyr::rename(name = transcript_name)
  DEG_tab_mod <- DEG_tab_mod[
    colnames(DEG_tab_mod)[colnames(DEG_tab_mod)
    %in% colnames(DET_final_tab_mod)]
  ]
  DEGs_DETs_table <- dplyr::bind_rows(DEG_tab_mod, DET_final_tab_mod)
  DEGs_DETs_table$significance <- c()
  DEGs_DETs_table$abs_log2FC <- base::abs(DEGs_DETs_table$log2FC)
  DEGs_DETs_table$significance <- "not_sig"
  DEGs_DETs_table$significance[DEGs_DETs_table$abs_log2FC > logfc_cut &
    DEGs_DETs_table$pvalue < pval_cut] <- "sig"
  return(DEGs_DETs_table)
}



#' Run Functional Analysis
#' @export
run_enrichment <- function(det_df,
                           genesets_list,
                           pval_cutoff = 0.05,
                           lfc_cutoff = 1) {
  tx_types_list <- list(
    "protein_coding",
    c(
      "retained_intron", "processed_transcript",
      "nonsense_mediated_decay"
    ),
    "retained_intron",
    "processed_transcript",
    "nonsense_mediated_decay"
  )
  # TODO: @luciorq Better define the non-coding category
  # + maybe something that translates:
  # + "non-coding isoform from protein coding gene + lncRNAs"
  # @iza editei aqui e coloquei unproductive no lugar de non-coding
  # tambem removi os lncRNAs ali em cima e aqui
  tx_type_names <- c(
    "protein_coding",
    "unproductive",
    "retained_intron",
    "processed_transcript",
    "nonsense_mediated_decay"
  )
  fgsea_results_df <- base::seq_along(tx_type_names) |>
    purrr::map_dfr(
      .f = \(x) {
        type_vec <- tx_types_list[[x]]
        type_name <- tx_type_names[x]
        res_fgsea <- det_df |>
          dplyr::filter(transcript_type %in% type_vec) |>
          dplyr::filter(!is.na(log2FC)) |>
          dplyr::arrange(log2FC)
        # TODO: @luciorq add option for log2FC * log10pvalue
        ranks <- res_fgsea$log2FC
        names(ranks) <- res_fgsea$gene_name
        fgsea_results <- fgsea::fgseaMultilevel(
          pathways = genesets_list,
          stats = ranks,
          eps = 0
        ) |>
          tibble::as_tibble() |>
          dplyr::filter(pval < 0.05) |>
          dplyr::mutate(experiment = type_name)
        return(fgsea_results)
      }
    )
  return(fgsea_results_df)
}

#' Plot Log2 Fold-Change Results
#' @export
plot_log2FC <- function(DEG_DET_table, selected_gene, custom_colors = NULL) {
  if (all(DEG_DET_table$significance[DEG_DET_table$gene_name %in% selected_gene] == "sig")) {
    ggplot2::ggplot(
      data = DEG_DET_table[DEG_DET_table$gene_name %in% selected_gene, ],
      mapping = ggplot2::aes(
        x = name,
        y = log2FC,
        fill = transcript_type
      )
    ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = tx_type_color_names) +
      ggplot2::theme_bw()
  } else {
    ggplot2::ggplot(
      data = DEG_DET_table[DEG_DET_table$gene_name %in% selected_gene, ],
      mapping = ggplot2::aes(
        x = name,
        y = log2FC,
        fill = transcript_type,
        alpha = significance
      )
    ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = tx_type_color_names) +
      ggplot2::theme_bw()
  }
  # ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5))
}
