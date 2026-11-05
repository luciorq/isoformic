#' Transcript Type Color Palette
#'
#' Default colour palette for transcript types: the GENCODE/Ensembl biotypes
#' used by the workflow, the `gene` level, generic labels for annotations
#' without biotypes (`unknown`, `novel`) and the SQANTI3 structural
#' categories assigned by [add_sqanti_classification()] to long-read
#' transcript catalogs.
#'
#' @returns A named character vector of hexadecimal colours.
#'
#' @examples
#' head(tx_type_palette())
#'
#' @export
tx_type_palette <- function() {
  extra <- c(
    unknown = "#bdbdbd",
    novel = "#8c6bb1",
    misc_RNA = "#d9d9d9",
    miRNA = "#d9d9d9",
    snRNA = "#d9d9d9",
    snoRNA = "#d9d9d9",
    `full-splice_match` = "#a6d854",
    `incomplete-splice_match` = "#66c2a5",
    novel_in_catalog = "#fdb462",
    novel_not_in_catalog = "#bc80bd",
    genic = "#d9d9d9",
    antisense = "#80b1d3",
    fusion = "#fccde5",
    intergenic = "#d9d9d9",
    genic_intron = "#d9d9d9",
    # FANTOM CAT gene classes
    coding_mRNA = "#a6d854",
    lncRNA_intergenic = "#e78ac3",
    lncRNA_antisense = "#8da0cb",
    lncRNA_divergent = "#fc8d62",
    lncRNA_sense_intronic = "#66c2a5",
    sense_overlap_RNA = "#ffd92f",
    short_ncRNA = "#e5c494",
    uncertain_coding = "#b3b3b3"
  )
  fixed_tx_biotypes <- c(
    "gene",
    "protein_coding",
    "retained_intron",
    "protein_coding_CDS_not_defined",
    "nonsense_mediated_decay",
    "lncRNA",
    "processed_pseudogene",
    "transcribed_unprocessed_pseudogene",
    "unprocessed_pseudogene",
    "non_stop_decay",
    "transcribed_unitary_pseudogene",
    "pseudogene",
    "unitary_pseudogene",
    "processed_transcript"
  )
  tx_type_color_names <- c(
    "#fb8072",
    "#a6d854",
    "#8da0cb",
    "#fc8d62",
    "#66c2a5",
    "#e78ac3",
    "#ffd92f",
    "#e5c494",
    "#d9d9d9",
    "#d9d9d9",
    "#d9d9d9",
    "#ffffb3",
    "#d9d9d9",
    "#d9d9d9"
  )
  names(tx_type_color_names) <- fixed_tx_biotypes
  return(c(tx_type_color_names, extra))
}
