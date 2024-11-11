#' Create Transcript-to-Gene Relationship Table
#'
#' Extracts a transcript-to-gene mapping table from GENCODE annotation files, such as the transcriptome FASTA file.
#' Currently, only FASTA files are supported.
#'
#' @param file_path A character string specifying the path to the reference file (e.g., GENCODE FASTA file).
#' @param file_type A character string specifying the type of the reference file. Currently, only `"fasta"` is supported.
#'   Default is `"fasta"`.
#'
#' @return A `tibble` containing the transcript-to-gene mapping information, including transcript IDs, gene IDs,
#'   transcript names, gene names, and transcript types.
#'
#' @details
#' The function reads the headers of the FASTA file and extracts relevant information to create a mapping table.
#' For GTF or GFF3 files, support is not yet implemented.
#'
#' @examples
#' # Assuming you have downloaded the GENCODE transcriptome FASTA file:
#' fasta_file <- download_reference(
#'   version = "43",
#'   organism = "human",
#'   file_type = "fasta",
#'   output_path = "data-raw"
#' )
#'
#' # Create the transcript-to-gene mapping table
#' tx_to_gene <- make_tx_to_gene(file_path = fasta_file, file_type = "fasta")
#'
#' # View the first few rows
#' head(tx_to_gene)
#'
#' @export
make_tx_to_gene <- function(file_path, file_type = c("fasta", "gff")) {
  file_type <- rlang::arg_match(file_type)
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

#' Annotate Transcripts with Differential Gene Expression Significance
#'
#' Adds a column to a transcript-level differential expression table indicating whether each transcript
#' originates from a gene that is significantly differentially expressed.
#'
#' @param DegsigVector A character vector containing the names of transcripts from significantly differentially expressed genes.
#' @param DET_table A `data.frame` or `tibble` containing transcript-level differential expression results,
#'   including a `transcript_name` column.
#'
#' @return A `tibble` with an additional column `DEG_sig` indicating whether the transcript is from a significantly
#'   differentially expressed gene (`"YES"` or `"NO"`).
#'
#' @examples
#' # Sample data
#' significant_transcripts <- c("transcript1", "transcript3")
#' DET_table <- data.frame(
#'   transcript_name = c("transcript1", "transcript2", "transcript3", "transcript4"),
#'   log2FC = c(2.5, -1.2, 0.8, -0.5),
#'   pvalue = c(0.01, 0.2, 0.03, 0.6)
#' )
#'
#' # Annotate transcripts with DEG significance
#' DET_table_annotated <- is_deg_sig(DegsigVector = significant_transcripts, DET_table = DET_table)
#'
#' # View the result
#' print(DET_table_annotated)
#'
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
#'
#' Combines gene-level and transcript-level differential expression results into a single table,
#' annotates the combined data with significance labels based on specified cutoffs, and filters
#' transcripts based on their types.
#'
#' @param DEG_tab A `data.frame` or `tibble` containing gene-level differential expression results,
#'   including `gene_id`, `gene_name`, `log2FC`, and `pvalue` columns.
#' @param DET_final_tab A `data.frame` or `tibble` containing transcript-level differential expression results,
#'   including `transcript_id`, `transcript_name`, `transcript_type`, `log2FC`, and `pvalue` columns.
#' @param logfc_cut A numeric value specifying the absolute log2 fold-change cutoff for significance.
#' @param pval_cut A numeric value specifying the p-value cutoff for significance.
#'
#' @return A `tibble` combining gene and transcript differential expression results, with additional columns:
#'   - `id`: gene or transcript ID.
#'   - `name`: gene or transcript name.
#'   - `transcript_type`: type of transcript or `"gene"` for gene-level entries.
#'   - `abs_log2FC`: absolute value of log2 fold-change.
#'   - `significance`: `"sig"` if significant based on cutoffs, `"not_sig"` otherwise.
#'
#' @examples
#' # Sample gene-level data
#' DEG_tab <- data.frame(
#'   gene_id = c("gene1", "gene2"),
#'   gene_name = c("GeneA", "GeneB"),
#'   log2FC = c(1.5, -2.0),
#'   pvalue = c(0.01, 0.04)
#' )
#'
#' # Sample transcript-level data
#' DET_final_tab <- data.frame(
#'   transcript_id = c("tx1", "tx2", "tx3"),
#'   transcript_name = c("Transcript1", "Transcript2", "Transcript3"),
#'   transcript_type = c("protein_coding", "lncRNA", "processed_transcript"),
#'   log2FC = c(1.2, -1.8, 0.5),
#'   pvalue = c(0.02, 0.03, 0.2)
#' )
#'
#' # Merge and annotate differential expression results
#' DEGs_DETs_table <- join_DEG_DET(
#'   DEG_tab = DEG_tab,
#'   DET_final_tab = DET_final_tab,
#'   logfc_cut = 1,
#'   pval_cut = 0.05
#' )
#'
#' # View the result
#' print(DEGs_DETs_table)
#'
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




#' Run Gene Set Enrichment Analysis for Different Transcript Types
#'
#' Performs gene set enrichment analysis (GSEA) on differential expression results for various transcript types,
#' using the `fgsea` package. The function iterates over specified transcript types, filters the data accordingly,
#' and runs GSEA for each type.
#'
#' @param det_df A `data.frame` or `tibble` containing transcript-level differential expression results,
#'   including `transcript_type`, `log2FC`, and `gene_name` columns.
#' @param genesets_list A list of gene sets to be used in the enrichment analysis.
#' @param pval_cutoff A numeric value specifying the p-value cutoff for the enrichment results. Default is `0.05`.
#' @param lfc_cutoff A numeric value specifying the log2 fold-change cutoff for filtering transcripts. Default is `1`.
#'
#' @return A `tibble` containing the enrichment analysis results for each transcript type, including pathway names,
#'   p-values, adjusted p-values, and the transcript type (experiment).
#'
#' @details
#' The function defines a list of transcript types and their corresponding labels.
#' It then filters the input differential expression data for each transcript type, ranks the genes by log2 fold-change,
#' and performs GSEA using the `fgsea` package.
#'
#' @examples
#' # Sample differential expression data
#' det_df <- data.frame(
#'   gene_name = c("GeneA", "GeneB", "GeneC", "GeneD"),
#'   transcript_type = c(
#'     "protein_coding", "retained_intron",
#'     "processed_transcript", "nonsense_mediated_decay"
#'   ),
#'   log2FC = c(1.5, -2.0, 0.8, -1.2)
#' )
#'
#' # Sample gene sets
#' genesets_list <- list(
#'   Pathway1 = c("GeneA", "GeneC"),
#'   Pathway2 = c("GeneB", "GeneD")
#' )
#'
#' # Run enrichment analysis
#' fgsea_results_df <- run_enrichment(
#'   det_df = det_df,
#'   genesets_list = genesets_list,
#'   pval_cutoff = 0.05,
#'   lfc_cutoff = 1
#' )
#'
#' # View the results
#' print(fgsea_results_df)
#'
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

  # TODO: @iza editei aqui e coloquei unproductive no lugar de non-coding
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


#' Plot Log2 Fold-Change Results for Selected Genes
#'
#' Creates a bar plot of log2 fold-change values for transcripts of a selected gene,
#' differentiating transcript types and significance levels.
#'
#' @param DEG_DET_table A `data.frame` or `tibble` containing combined gene and transcript differential expression results,
#'   including `name`, `log2FC`, `transcript_type`, `significance`, and `gene_name` columns.
#' @param selected_gene A character string specifying the gene name to plot.
#' @param custom_colors An optional named vector of colors for different transcript types.
#'
#' @return A `ggplot2` object representing the bar plot.
#'
#' @details
#' The function filters the input table for the selected gene and creates a bar plot of log2 fold-change values.
#' If all transcripts are significant, it plots without adjusting alpha transparency; otherwise, it adjusts alpha
#' based on significance. The function uses predefined colors for transcript types, which can be overridden
#' by providing `custom_colors`.
#'
#' @examples
#' # Sample data
#' DEGs_DETs_table <- data.frame(
#'   name = c("Transcript1", "Transcript2", "GeneA"),
#'   log2FC = c(1.5, -2.0, 0.8),
#'   transcript_type = c("protein_coding", "lncRNA", "gene"),
#'   significance = c("sig", "not_sig", "sig"),
#'   gene_name = c("GeneA", "GeneA", "GeneA")
#' )
#'
#' # Plot log2 fold-change for the selected gene
#' plot_obj <- plot_log2FC(
#'   DEG_DET_table = DEGs_DETs_table,
#'   selected_gene = "GeneA"
#' )
#'
#' # Display the plot
#' print(plot_obj)
#'
#' @export
plot_log2FC <- function(DEG_DET_table, selected_gene, custom_colors = NULL) {
  if (all(DEG_DET_table$significance[DEG_DET_table$gene_name %in% selected_gene] == "sig")) {
    plot_obj <- ggplot2::ggplot(
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
    plot_obj <- ggplot2::ggplot(
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
  return(plot_obj)
  # ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5))
}
