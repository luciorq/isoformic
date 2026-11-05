# Accessor generics for IsoformicExperiment -----------------------------------

#' Accessors for IsoformicExperiment Objects
#'
#' @description
#' S7 generics giving `SummarizedExperiment`-style access to the three
#' levels of an [IsoformicExperiment]. All of them return in-memory tibbles
#' or matrices; the annotation tables are read from Parquet on first use
#' and cached (see the class documentation).
#'
#' **Samples**
#' * `col_data()`, `col_data<-`: sample metadata (`sample_id` first).
#' * `col_names()`: sample IDs of the assays.
#'
#' **Transcript level**
#' * `assays()`, `assay(x, name)`, `assay_names()`, `assay<-`: transcript
#'   matrices. `assay()` without a name returns the first assay.
#' * `row_names()`: transcript IDs of the assays.
#' * `row_data()`: transcript annotation (with gene columns) for the assayed
#'   transcripts, in assay row order.
#' * `tx_annot()`: the same join for *all* annotated transcripts.
#' * `tx_to_gene()`: two-column transcript/gene mapping.
#'
#' **Gene level**
#' * `gene_assays()`, `gene_assay(x, name)`, `gene_assays<-`: gene matrices,
#'   derived by summing isoforms unless set explicitly.
#' * `gene_data()`: gene annotation for the genes represented in the assays.
#' * `gene_names()`: gene IDs of the gene assays.
#'
#' **Annotation**
#' * `annot_data_transcripts()`, `annot_data_genes()`, `annot_data_exons()`:
#'   the raw annotation tables.
#' * `annot_data()`: transcript-centric join of transcripts and genes
#'   (exons optional).
#' * `annot_metadata()`: named list describing the annotation source.
#' * `exon_table(x, gene_name)`: exons of the transcripts of selected genes.
#'
#' **Housekeeping**
#' * `clear_cache()`: drop cached derived data for the object (and, with
#'   `annotation = TRUE`, the shared annotation tables).
#'
#' @param self An [IsoformicExperiment].
#' @param value Replacement value.
#' @param name Name of the assay. `NULL` selects the first assay.
#' @param gene_name Character vector of gene names or gene IDs.
#' @param with_exons Logical, also join the exon table (one row per exon).
#' @param annotation Logical, also clear the shared annotation cache.
#' @param ... Passed on to methods.
#'
#' @returns See the individual descriptions above.
#'
#' @examples
#' iso <- example_isoformic()
#' col_data(iso)
#' assay_names(iso)
#' assay(iso, "abundance")[1:3, 1:3]
#' gene_assay(iso, "abundance")[1:3, 1:3]
#' head(row_data(iso))
#' head(gene_data(iso))
#' head(exon_table(iso, "RBPJ"))
#'
#' @name IsoformicExperiment-accessors
#' @rdname IsoformicExperiment-accessors
NULL

# Samples ---------------------------------------------------------------------

#' @rdname IsoformicExperiment-accessors
#' @export
col_data <- S7::new_generic("col_data", "self")

S7::method(col_data, IsoformicExperiment) <- function(self) {
  self@col_data
}

#' @rdname IsoformicExperiment-accessors
#' @export
`col_data<-` <- S7::new_generic("col_data<-", "self", function(self, value) {
  S7::S7_dispatch()
})

S7::method(`col_data<-`, IsoformicExperiment) <- function(self, value) {
  self@col_data <- value
  self
}

#' @rdname IsoformicExperiment-accessors
#' @export
col_names <- S7::new_generic("col_names", "self")

S7::method(col_names, IsoformicExperiment) <- function(self) {
  dimnames(self)[[2L]]
}

# Transcript level ------------------------------------------------------------

#' @rdname IsoformicExperiment-accessors
#' @export
assays <- S7::new_generic("assays", "self")

S7::method(assays, IsoformicExperiment) <- function(self) {
  self@assays
}

#' @rdname IsoformicExperiment-accessors
#' @export
assay_names <- S7::new_generic("assay_names", "self")

S7::method(assay_names, IsoformicExperiment) <- function(self) {
  names(self@assays)
}

#' @rdname IsoformicExperiment-accessors
#' @export
assay <- S7::new_generic("assay", "self", function(self, name = NULL, ...) {
  S7::S7_dispatch()
})

S7::method(assay, IsoformicExperiment) <- function(self, name = NULL, ...) {
  pick_assay(self@assays, name, level = "transcript")
}

#' @rdname IsoformicExperiment-accessors
#' @export
`assay<-` <- S7::new_generic("assay<-", "self", function(self, name, value) {
  S7::S7_dispatch()
})

S7::method(`assay<-`, IsoformicExperiment) <- function(self, name, value) {
  check_string(name)
  assays <- self@assays
  assays[[name]] <- value
  self@assays <- assays
  self
}

#' @rdname IsoformicExperiment-accessors
#' @export
row_names <- S7::new_generic("row_names", "self")

S7::method(row_names, IsoformicExperiment) <- function(self) {
  dimnames(self)[[1L]]
}

#' @rdname IsoformicExperiment-accessors
#' @export
row_data <- S7::new_generic("row_data", "self")

S7::method(row_data, IsoformicExperiment) <- function(self) {
  ids <- row_names(self)
  if (is.null(ids)) {
    return(tx_gene_annotation(self)[0L, ])
  }
  signature <- paste(rlang::hash(ids), annot_full_signature(self))
  object_cache_fetch(self, "row_data", signature, function() {
    annot <- tx_gene_annotation(self)
    annot[match(ids, annot$transcript_id), ]
  })
}

#' @rdname IsoformicExperiment-accessors
#' @export
tx_annot <- S7::new_generic("tx_annot", "self")

S7::method(tx_annot, IsoformicExperiment) <- function(self) {
  tx_gene_annotation(self)
}

#' @rdname IsoformicExperiment-accessors
#' @export
tx_to_gene <- S7::new_generic("tx_to_gene", "self")

S7::method(tx_to_gene, IsoformicExperiment) <- function(self) {
  self@annot_data_transcripts[, c("transcript_id", "gene_id")]
}

#' Transcript annotation joined with gene attributes (cached)
#' @keywords internal
#' @noRd
tx_gene_annotation <- function(self) {
  object_cache_fetch(
    self,
    "tx_gene_annotation",
    annot_full_signature(self),
    function() {
      tx <- self@annot_data_transcripts
      genes <- self@annot_data_genes
      tx |>
        dplyr::select(
          -dplyr::any_of(c("seqid", "start_pos", "end_pos", "strand"))
        ) |>
        dplyr::left_join(
          genes[, intersect(
            c("gene_id", "gene_name", "gene_type"),
            colnames(genes)
          )],
          by = "gene_id"
        ) |>
        dplyr::relocate("gene_name", "gene_type", .after = "gene_id")
    }
  )
}

# Gene level ------------------------------------------------------------------

#' @rdname IsoformicExperiment-accessors
#' @export
gene_assays <- S7::new_generic("gene_assays", "self")

S7::method(gene_assays, IsoformicExperiment) <- function(self) {
  self@gene_assays
}

#' @rdname IsoformicExperiment-accessors
#' @export
`gene_assays<-` <- S7::new_generic(
  "gene_assays<-",
  "self",
  function(self, value) {
    S7::S7_dispatch()
  }
)

S7::method(`gene_assays<-`, IsoformicExperiment) <- function(self, value) {
  self@gene_assays <- value
  self
}

#' @rdname IsoformicExperiment-accessors
#' @export
gene_assay <- S7::new_generic(
  "gene_assay",
  "self",
  function(self, name = NULL, ...) {
    S7::S7_dispatch()
  }
)

S7::method(gene_assay, IsoformicExperiment) <- function(
  self,
  name = NULL,
  ...
) {
  pick_assay(self@gene_assays, name, level = "gene")
}

#' @rdname IsoformicExperiment-accessors
#' @export
gene_names <- S7::new_generic("gene_names", "self")

S7::method(gene_names, IsoformicExperiment) <- function(self) {
  ga <- self@gene_assays
  if (!length(ga)) {
    return(NULL)
  }
  rownames(ga[[1L]])
}

#' @rdname IsoformicExperiment-accessors
#' @export
gene_data <- S7::new_generic("gene_data", "self")

S7::method(gene_data, IsoformicExperiment) <- function(self) {
  genes <- self@annot_data_genes
  ids <- gene_names(self)
  if (is.null(ids)) {
    return(genes[0L, ])
  }
  genes[match(ids, genes$gene_id), ]
}

pick_assay <- function(assays, name, level) {
  if (!length(assays)) {
    cli::cli_abort(
      message = c(x = "The object has no {level}-level assays."),
      class = "isoformic_no_assays"
    )
  }
  if (is.null(name)) {
    return(assays[[1L]])
  }
  if (is.numeric(name)) {
    name <- names(assays)[name]
  }
  if (!rlang::is_string(name) || !name %in% names(assays)) {
    cli::cli_abort(
      message = c(
        x = "Assay {.val {name}} not found.",
        i = "Available {level} assays: {.val {names(assays)}}."
      ),
      class = "isoformic_assay_not_found"
    )
  }
  assays[[name]]
}

# Annotation ------------------------------------------------------------------

#' @rdname IsoformicExperiment-accessors
#' @export
annot_data_transcripts <- S7::new_generic("annot_data_transcripts", "self")

S7::method(annot_data_transcripts, IsoformicExperiment) <- function(self) {
  self@annot_data_transcripts
}

#' @rdname IsoformicExperiment-accessors
#' @export
annot_data_genes <- S7::new_generic("annot_data_genes", "self")

S7::method(annot_data_genes, IsoformicExperiment) <- function(self) {
  self@annot_data_genes
}

#' @rdname IsoformicExperiment-accessors
#' @export
annot_data_exons <- S7::new_generic("annot_data_exons", "self")

S7::method(annot_data_exons, IsoformicExperiment) <- function(self) {
  self@annot_data_exons
}

#' @rdname IsoformicExperiment-accessors
#' @export
annot_data <- S7::new_generic(
  "annot_data",
  "self",
  function(self, with_exons = FALSE, ...) {
    S7::S7_dispatch()
  }
)

S7::method(annot_data, IsoformicExperiment) <- function(
  self,
  with_exons = FALSE,
  ...
) {
  genes <- self@annot_data_genes |>
    dplyr::select(-dplyr::any_of(c("seqid", "strand"))) |>
    dplyr::rename(gene_start = "start_pos", gene_end = "end_pos")
  res <- dplyr::left_join(self@annot_data_transcripts, genes, by = "gene_id")
  if (isTRUE(with_exons)) {
    exons <- self@annot_data_exons |>
      dplyr::select(-dplyr::any_of(c("seqid", "strand"))) |>
      dplyr::rename(exon_start = "start_pos", exon_end = "end_pos")
    res <- dplyr::left_join(res, exons, by = "transcript_id")
  }
  res
}

#' @rdname IsoformicExperiment-accessors
#' @export
annot_metadata <- S7::new_generic("annot_metadata", "self")

S7::method(annot_metadata, IsoformicExperiment) <- function(self) {
  self@annot_metadata
}

#' @rdname IsoformicExperiment-accessors
#' @export
exon_table <- S7::new_generic(
  "exon_table",
  "self",
  function(self, gene_name, ...) {
    S7::S7_dispatch()
  }
)

S7::method(exon_table, IsoformicExperiment) <- function(self, gene_name, ...) {
  rlang::check_required(gene_name)
  build_exon_table(
    genes = self@annot_data_genes,
    transcripts = self@annot_data_transcripts,
    exons = self@annot_data_exons,
    gene_name = gene_name
  )
}

# Housekeeping ----------------------------------------------------------------

#' @rdname IsoformicExperiment-accessors
#' @export
clear_cache <- S7::new_generic(
  "clear_cache",
  "self",
  function(self, annotation = FALSE, ...) {
    S7::S7_dispatch()
  }
)

S7::method(clear_cache, IsoformicExperiment) <- function(
  self,
  annotation = FALSE,
  ...
) {
  cache <- attr(self, ".cache", exact = TRUE)
  rm(list = ls(cache, all.names = TRUE), envir = cache)
  if (isTRUE(annotation)) {
    rm(list = ls(annotation_cache, all.names = TRUE), envir = annotation_cache)
  }
  invisible(self)
}
