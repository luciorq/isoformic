#' Plot Transcripts Genomic Context
#'
#' Generate a genomic context plot for a specified gene, displaying its transcripts
#' along with their types and annotations.
#' The function utilizes the `plotgardener` package to create a detailed
#' visualization of the genomic context, including an ideogram, chromosome
#' highlight, and transcript structures.
#' It requires a `ContextData` object that contains the necessary genomic
#' information and annotations.
#' The plot can be customized with various parameters such as offsets,
#' label limits, and ideogram references.
#'
#' @param gene_name Character string specifying the name of the gene to plot.
#' @param context_data A `ContextData` object containing genomic context information.
#' @param limit_label Logical indicating whether to limit the length of
#' transcript labels to avoid overlap (default is `TRUE`).
#' @param show_guides Logical indicating whether to show guide lines on the plot
#' (default is `FALSE`).
#' @param y_offset Numeric value to adjust the vertical position of the plot
#' (default is `0`).
#' @param height_offset Numeric value to adjust the height of the plot
#' (default is `0`).
#' @param downstream_offset Numeric value to extend the downstream region
#' beyond the gene's end position (default is `0`).
#' @param upstream_offset Numeric value to extend the upstream region
#' beyond the gene's start position (default is `0`).
#' @param ideogram_reference Character string specifying the reference genome
#' for the ideogram. Options include "hg38", "hg19", "mm11",
#' "mm10", or "none" (default is "hg38").
#' @return A `plotgardener` object representing the genomic context plot.
#' @export
plot_genomic_context <- function(
    gene_name,
    context_data,
    limit_label = TRUE,
    show_guides = FALSE,
    y_offset = 0,
    height_offset = 0,
    downstream_offset = 0,
    upstream_offset = 0,
    ideogram_reference = c("hg38", "hg19", "mm11", "mm10", "none")) {
  .data <- rlang::.data

  rlang::check_required(gene_name)
  rlang::check_required(context_data)
  ideogram_reference <- rlang::arg_match(ideogram_reference)

  if (rlang::is_character(gene_name)) {
    context_data@gene_name <- gene_name
  } else {
    cli::cli_abort("`gene_name` must be a character vector.")
  }

  custom_assembly <- plotgardener::assembly(
    # Genome = "hg38_GENCODE34",
    Genome = paste0(
      context_data@assembly_name,
      context_data@annotation_name,
      sep = "_"
    ),
    TxDb = context_data@txdb,
    OrgDb = context_data@orgdb_package,
    BSgenome = context_data@bsgenome_package,
    gene.id.column = "TXSYMBOL",
    display.column = "TXSYMBOL"
  )

  tx_biotype_df <- context_data@txdb$conn |>
    dplyr::tbl("transcript") |>
    dplyr::collect() |>
    dplyr::mutate(
      gene_symbol = stringr::str_replace(.data$tx_symbol, "-\\d+$", "")
    ) |>
    dplyr::filter(.data$gene_symbol %in% context_data@gene_name) |>
    dplyr::select("tx_symbol", "tx_biotype", "tx_ensemblid")

  tx_type_color_list <- context_data@tx_type_palette
  tx_type_color_df <- tibble::tibble(
    tx_biotype = names(tx_type_color_list),
    color = unname(tx_type_color_list)
  )

  tx_color_df <- tx_biotype_df |>
    dplyr::left_join(tx_type_color_df, by = "tx_biotype") |>
    dplyr::rename(
      transcript = "tx_symbol",
      tx_name = "tx_ensemblid"
    )

  transcript_highlights <- data.frame(
    transcript = tx_color_df$transcript,
    color = tx_color_df$color
  )
  exons_gr_list <- GenomicFeatures::exonsBy(
    context_data@txdb,
    by = "tx",
    use.names = TRUE
  )
  selected_exons_gr_list <- exons_gr_list[unique(tx_color_df$transcript)]
  start_value <- min(unlist(lapply(selected_exons_gr_list, \(x) {
    IRanges::start(IRanges::ranges(x))
  })))
  end_value <- max(unlist(lapply(selected_exons_gr_list, \(x) {
    IRanges::end(IRanges::ranges(x))
  })))

  start_value <- as.integer(
    round(start_value - (0.1 * ((end_value - start_value) + upstream_offset)))
  )
  end_value <- as.integer(
    round(end_value + (0.1 * ((end_value - start_value) + downstream_offset)))
  )
  chr_string <- as.character(selected_exons_gr_list[[1]]@seqnames@values)

  # ===========================================================================
  # Plot start here
  # ===========================================================================
  # rlang::is_bare_numeric(y_offset)

  plotgardener::pageCreate(
    width = 7.5,
    height = 5.0 + height_offset,
    default.units = "inches",
    showGuides = show_guides
  )

  plotgardener::plotText(
    label = stringr::str_replace(chr_string, "chr", "Chromosome "),
    fontcolor = "dark grey",
    x = 6.5 + 0.5,
    y = 0.9,
    just = "right"
  )

  if (end_value - start_value < 100000) {
    end_value_highlight <- start_value + 150000
  } else {
    end_value_highlight <- end_value
  }
  region_param <- plotgardener::pgParams(
    chrom = chr_string,
    chromstart = start_value,
    chromend = end_value_highlight
  )

  # TODO: @luciorq - Check which other assemblies would be supported
  # + by plotgardener
  if (ideogram_reference == "none") {
    chr_length <- context_data@txdb$conn |>
      dplyr::tbl("chrominfo") |>
      dplyr::filter(.data$chrom == chr_string) |>
      dplyr::pull(length)

    if (isTRUE(is.na(chr_length))) {
      chr_length <- context_data@txdb$conn |>
        dplyr::tbl("transcript") |>
        dplyr::filter(.data$tx_chrom == chr_string) |>
        dplyr::slice_max(order_by = .data$tx_end) |>
        dplyr::pull("tx_end")
    }
    chrom_plot <- plotgardener::plotGenomeLabel(
      chrom = chr_string,
      chromstart = 0L,
      chromend = chr_length,
      assembly = custom_assembly,
      fontcolor = "dark grey",
      x = 0.5,
      y = 0.5,
      length = 6.5,
      default.units = "inches",
      scale = "Mb"
    )
  } else {
    chrom_plot <- plotgardener::plotIdeogram(
      chrom = chr_string,
      assembly = ideogram_reference,
      orientation = "h",
      x = 0.5,
      y = 0.5, # grid::unit(0.25, "npc"),
      width = 6.5,
      height = 0.3,
      just = c("left", "top"),
      default.units = "inches"
    )
  }
  plotgardener::annoHighlight(
    plot = chrom_plot,
    params = region_param,
    fill = "red",
    y = 0.3,
    height = 0.7,
    just = c("left", "top"),
    default.units = "inches"
  )

  plotgardener::annoZoomLines(
    plot = chrom_plot,
    params = region_param,
    y0 = 1,
    x1 = c(0.5, 6.5),
    y1 = 1.6,
    default.units = "inches",
    just = c("left", "top")
  )

  transcripts_res <- plotgardener::plotTranscripts(
    chrom = chr_string,
    chromstart = start_value,
    chromend = end_value,
    assembly = custom_assembly,
    labels = "transcript",
    just = c("center", "center"),
    x = 7.5 / 2,
    y = ((4.5 + y_offset) / 2),
    width = 6.5,
    height = 3.6 + (height_offset * 2),
    draw = TRUE,
    limitLabel = limit_label,
    fill = "grey",
    colorbyStrand = FALSE,
    transcriptHighlights = transcript_highlights,
    strandSplit = TRUE,
    fontsize = 10,
    boxHeight = grid::unit(3.5, "mm"),
    stroke = 0.05
  )

  plotgardener::plotLegend(
    legend = unique(tx_color_df$tx_biotype),
    fill = unique(tx_color_df$color),
    border = FALSE,
    x = 0.6,
    y = 4.0 + height_offset,
    width = 1.3 * length(unique(tx_color_df$tx_biotype)),
    height = 0.5,
    just = c("left", "top"),
    orientation = "h"
  )

  plotgardener::plotGenomeLabel(
    chrom = chr_string,
    chromstart = start_value,
    chromend = end_value,
    assembly = custom_assembly,
    x = 0.5,
    y = 4.5 + height_offset,
    length = 6.5,
    default.units = "inches"
  )

  return(invisible(transcripts_res))
}
