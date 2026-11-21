# Enrichment ------------------------------------------------------------------

#' Transcript-Type Aware Gene Set Enrichment
#'
#' @description
#' Runs gene set enrichment analysis (`fgsea::fgseaMultilevel()`) separately
#' for transcripts of each type, ranking transcripts by the statistic chosen
#' with `rank_by` (log2 fold-change by default, see *Ranking*). Gene sets
#' defined on gene symbols are expanded to all transcripts of their genes,
#' so each transcript category (protein coding, retained intron,
#' nonsense-mediated decay, ...) is tested on its own and in the combined
#' `"unproductive"` category. A category with fewer than two ranked
#' transcripts is skipped.
#'
#' The [IsoformicExperiment] method uses the stored transcript-level
#' results and the annotation; the data frame method takes a transcript
#' table with `transcript_name` and `transcript_type` columns and the
#' columns `rank_by` needs (`log2FC`; `log2FC` and `pvalue`; or
#' `lfc_median` and `lfc_sd`), plus a `tx_to_gene` table.
#'
#' [run_productivity_enrichment()] and [run_apa_enrichment()] run the same
#' analysis on genes ranked by their change in productive fraction or in
#' 3' end usage, and take gene set collections in the same way.
#'
#' `read_gmt()` reads gene sets in GMT format into a named list, and
#' `plot_enrichment()` draws a dot plot of enrichment results.
#'
#' @section Ranking:
#' `rank_by` chooses the statistic that ranks the transcripts:
#'
#' * `"log2FC"` (default): the log2 fold-change of the stored results, a
#'   point estimate that ignores how well it is determined.
#' * `"signed_p"`: `-log10(pvalue)` with the sign of the fold-change, so that
#'   the test of the differential expression weighs the ranking. A p-value
#'   of zero gets the largest finite value of the table.
#' * `"lfc_z"`: the median fold-change over the inferential replicates
#'   divided by its standard deviation (`lfc_median / lfc_sd` of
#'   [de_confidence()]), so that a fold-change the quantification cannot
#'   determine ranks low. It needs the fold-change table that
#'   `add_uncertainty(x, group = )` stores, or these two columns in a data
#'   frame (the table of [de_confidence()] works). The ranking has the
#'   contrast of that table, the last level of `add_uncertainty()` against
#'   the first; a warning (`isoformic_contrast_mismatch`) says so when most
#'   of its fold-changes have the sign opposite to the stored results.
#'   Transcripts that were not assessed, that were collapsed into a group
#'   or, for an object, that are not in its differential expression results
#'   are not ranked. A fold-change that is the same in every replicate
#'   gets the largest finite value of the table.
#'
#' @section Collections of gene sets:
#' `genesets_list` is either one collection (a named list of gene symbol
#' vectors) or several (a named list of such lists, for example
#' `list(hallmark = read_gmt(h), reactome = read_gmt(r))`). Each collection is
#' tested on its own, so that the adjusted p-values are those of the
#' collection, and the result has a `collection` column, by which
#' `plot_enrichment()` splits its panels. Every gene set and every collection
#' needs a name.
#'
#' @section Workers:
#' Each collection and transcript category is one `fgsea` call. With
#' `workers` above 1 (or the `isoformic.workers` option), these calls run in
#' parallel on that many `mirai` daemons (at most one per call), which need
#' the `mirai` package. `workers = "auto"` uses
#' `parallelly::availableCores()`, which follows the CPUs a job scheduler
#' (Slurm, the CPU affinity mask, a cgroup quota) allocated; without
#' `parallelly` it runs serially.
#'
#' The daemons are fresh R processes that receive only the ranking of their
#' call and the gene sets cut down to the IDs of that ranking (which gives
#' the same result), never the session's objects, so their memory does not
#' grow with the session. The session serializes all calls at once when it
#' sends them, which costs it about the size of each ranking and of the gene
#' sets cut down to it, once per call.
#'
#' The first call with `workers` above 1 starts the daemons (one per
#' `fgsea` call, up to `workers`), with the session's library paths and one
#' BLAS and OpenMP thread each, and keeps them for the session: later calls
#' reuse them, so a loop of small calls (one contrast and a few collections
#' each) does not pay their start-up and their loading of `fgsea` every
#' time. A call that needs more daemons launches the missing ones; a call
#' with a smaller `workers` than the daemons running starts them again for
#' that call (one per `fgsea` call, up to its `workers`), so that no call
#' runs more `fgsea` calls at once than its `workers`. Daemons that died
#' between calls (killed for their memory, for example) are replaced by the
#' next call; when all of them die during a call, it fails after 30 seconds
#' without a connected daemon instead of waiting forever, and the next call
#' starts them again. An idle daemon holds the memory of its R process and of
#' `fgsea` (about 0.3 GB once it has run a call). Stop the daemons with
#' [`isoformic_workers(0)`][isoformic_workers()]; they are also stopped when
#' the session ends or isoformic is unloaded, and an interrupted call cancels
#' its `fgsea` calls on them. [`isoformic_workers(n)`][isoformic_workers()]
#' starts `n` daemons ahead
#' (stopping those running), with the same settings,
#' which every call then uses as they are (replacing those that died), and
#' so do daemons started with `mirai::daemons(n, .compute = "isoformic")`
#' (`workers` then only decides whether to use them; they keep their own
#' thread settings, but every call puts the session's library paths first
#' while it runs; a call fails when none of them is connected, because
#' isoformic does not restart them).
#'
#' Warnings of `fgsea` (ties in the ranking, imprecise p-values) raised on a
#' worker are raised again in the session when all calls have returned, once
#' per distinct message and call, naming the collection and the category.
#'
#' `fgsea` estimates p-values by sampling. With one worker the calls use the
#' session's random number stream, in order, as in earlier versions, so
#' `set.seed()` before the call reproduces them. With more than one worker
#' every call gets its own seed, drawn from the session's stream before the
#' calls start: the result is the same for any number of workers above one,
#' and whether or not the daemons were already running (and reproducible
#' with `set.seed()`), but its sampled p-values differ slightly from those
#' of a serial run. Starting daemons with `mirai::daemons()` draws from the
#' session's stream, so daemons started that way after `set.seed()` change
#' the p-values of the calls that follow; start them before `set.seed()`,
#' or with [isoformic_workers()], which leaves the stream where it was. A
#' single call always runs in the session.
#'
#' `BiocParallel::MulticoreParam()` (and `nproc` on Unix-like systems) forks
#' the session: the forked workers share its memory, but tools that add up
#' the resident memory of processes (`top`, Slurm accounting without
#' `UsePss`) count the whole session once per worker.
#'
#' @param self An [IsoformicExperiment] or a transcript-level DE table.
#' @param genesets_list Named list of character vectors (gene symbols), as
#'   returned by [read_gmt()], or a named list of such lists (see
#'   *Collections of gene sets*).
#' @param rank_by Statistic that ranks the transcripts: `"log2FC"`
#'   (default), `"signed_p"` or `"lfc_z"` (see *Ranking*).
#' @param tx_to_gene Data frame with `transcript_name` and `gene_name`
#'   columns (and optionally `transcript_type`, used to choose the default
#'   `tx_types`); required by the data frame method, ignored otherwise.
#' @param pval_cutoff Keep pathways with `pval` at or below this value;
#'   `NULL` keeps all.
#' @param tx_types Named list defining the transcript categories to test.
#'   Defaults to [default_enrichment_types()], which picks
#'   `processed_transcript` or `protein_coding_CDS_not_defined` depending
#'   on the annotation version.
#' @param min_size,max_size Minimum and maximum gene set size: sets with
#'   fewer or more genes are not tested. `run_enrichment()` counts, in each
#'   transcript category, the genes of the set that have at least one ranked
#'   transcript of that category (the `size` column of the result counts
#'   transcripts); the gene-level enrichments count the ranked genes of the
#'   set.
#' @param workers Number of worker processes that run the `fgsea` calls
#'   (one per collection and transcript category) in parallel: a whole
#'   number, or `"auto"` for the CPUs available to the session (see
#'   *Workers*). Defaults to the `isoformic.workers` option, or 1 (serial).
#' @param ... Passed on to `fgsea::fgseaMultilevel()`. Each `fgsea` call runs
#'   with `BiocParallel::SerialParam()` unless `nproc` or `BPPARAM` is given.
#'   Prefer `workers`: `nproc` parallelises inside one call, which was slower
#'   than a serial run in our benchmarks, and on Unix-like systems it forks
#'   the session (see *Workers*).
#' @param path Path to a GMT file.
#' @param enrichment_df Result of `run_enrichment()`.
#' @param n Number of top pathways (by adjusted p-value) to show.
#' @param experiments Categories to show; `NULL` shows all.
#' @param nes_cutoff Minimum absolute NES to show.
#' @param strip_prefix Regular expression removed from pathway names before
#'   plotting (default drops MSigDB collection prefixes such as `REACTOME_`).
#' @param label_width Pathway labels are wrapped to this many characters.
#' @param annotation_types Transcript types present in the annotation,
#'   used to choose between `processed_transcript` and
#'   `protein_coding_CDS_not_defined`.
#'
#' @returns `run_enrichment()` returns a tibble with the `fgsea` columns
#'   plus `experiment` (the transcript category), `rank_by` and, for several
#'   collections, `collection`. `read_gmt()` returns a named list of
#'   character vectors (the unique gene symbols of each set).
#'   `plot_enrichment()` returns a `ggplot` object, with one panel per
#'   collection when there is a `collection` column.
#'   `default_enrichment_types()` returns a named list of transcript types
#'   per category.
#'
#' @seealso [isoformic_workers()] to start or stop the daemons ahead.
#'
#' @examples
#' iso <- example_isoformic()
#' gmt <- system.file("extdata", "reactome_subset.gmt", package = "isoformic")
#' genesets <- read_gmt(gmt)
#' length(genesets)
#'
#' if (requireNamespace("fgsea", quietly = TRUE)) {
#'   res <- run_enrichment(iso, genesets, pval_cutoff = NULL, min_size = 2)
#'   head(res)
#'   plot_enrichment(res, n = 10)
#' }
#'
#' @name enrichment
#' @rdname enrichment
NULL

#' @rdname enrichment
#' @export
run_enrichment <- S7::new_generic(
  "run_enrichment",
  "self",
  function(
    self,
    genesets_list,
    tx_to_gene = NULL,
    pval_cutoff = 0.05,
    tx_types = NULL,
    min_size = 15,
    max_size = 500,
    rank_by = c("log2FC", "signed_p", "lfc_z"),
    workers = getOption("isoformic.workers", 1L),
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(run_enrichment, IsoformicExperiment) <- function(
  self,
  genesets_list,
  tx_to_gene = NULL,
  pval_cutoff = 0.05,
  tx_types = NULL,
  min_size = 15,
  max_size = 500,
  rank_by = c("log2FC", "signed_p", "lfc_z"),
  workers = getOption("isoformic.workers", 1L),
  ...
) {
  rlang::check_required(genesets_list)
  rank_by <- rlang::arg_match(rank_by)
  det <- de_tx(self)
  annot <- tx_annot(self)
  # The annotation wins over columns of the same name in the stored results
  det <- join_annotation(
    det,
    annot[, c(
      "transcript_id",
      "transcript_name",
      "transcript_type",
      "gene_name"
    )],
    by = "transcript_id"
  )
  if (identical(rank_by, "lfc_z")) {
    conf <- self@dea[["confidence"]]
    if (is.null(conf)) {
      cli::cli_abort(
        message = c(
          x = "{.code rank_by = \"lfc_z\"} needs the fold-change table of {.code add_uncertainty(group = )}.",
          i = "Store it with {.code add_uncertainty(x, group = )}."
        ),
        class = "isoformic_no_inferential_replicates"
      )
    }
    conf <- conf[conf$feature_type %in% "transcript", ]
    idx <- match(det$transcript_id, conf$feature_id)
    det$lfc_median <- conf$lfc_median[idx]
    det$lfc_sd <- conf$lfc_sd[idx]
    # The fold-change table has the contrast of `add_uncertainty(group = , levels = )`,
    # which may not be the one of the stored results
    same_sign <- sign(det$lfc_median) == sign(det$log2FC)
    same_sign <- same_sign[
      !is.na(same_sign) & det$lfc_median != 0 & det$log2FC != 0
    ]
    if (length(same_sign) >= 2L && mean(same_sign) < 0.5) {
      cli::cli_warn(
        message = c(
          "!" = paste0(
            "Most fold-changes of the fold-change table of ",
            "{.code add_uncertainty(group = )} have the sign ",
            "opposite to the stored differential expression results."
          ),
          i = paste0(
            "{.code rank_by = \"lfc_z\"} has the contrast of ",
            "{.code add_uncertainty(x, group = , levels = )}: check its levels."
          )
        ),
        class = "isoformic_contrast_mismatch"
      )
    }
  }
  tx_types <- tx_types %||%
    default_enrichment_types(unique(annot$transcript_type))
  run_fgsea_by_type(
    det = det,
    rank_id = "transcript_id",
    genesets_list = genesets_list,
    tx_to_gene = annot[, c("transcript_id", "gene_name")],
    map_id = "transcript_id",
    pval_cutoff = pval_cutoff,
    tx_types = tx_types,
    min_size = min_size,
    max_size = max_size,
    rank_by = rank_by,
    workers = workers,
    ...
  )
}

# nolint start: object_name_linter.
S7::method(run_enrichment, S7::class_data.frame) <- function(
  self,
  genesets_list,
  tx_to_gene = NULL,
  pval_cutoff = 0.05,
  tx_types = NULL,
  min_size = 15,
  max_size = 500,
  rank_by = c("log2FC", "signed_p", "lfc_z"),
  workers = getOption("isoformic.workers", 1L),
  ...
) {
  rlang::check_required(genesets_list)
  rank_by <- rlang::arg_match(rank_by)
  if (is.null(tx_to_gene)) {
    cli::cli_abort(
      message = c(
        x = "{.arg tx_to_gene} is required when {.arg self} is a data frame."
      ),
      class = "isoformic_invalid_tx_to_gene"
    )
  }
  det <- tibble::as_tibble(self)
  required <- c(
    "transcript_name",
    "transcript_type",
    if (!identical(rank_by, "lfc_z")) "log2FC"
  )
  missing_cols <- setdiff(required, colnames(det))
  if (length(missing_cols)) {
    cli::cli_abort(
      message = c(
        x = "{.arg self} is missing column{?s}: {.field {missing_cols}}."
      ),
      class = "isoformic_invalid_de_table"
    )
  }
  types_in_annot <- if ("transcript_type" %in% colnames(tx_to_gene)) {
    unique(tx_to_gene$transcript_type)
  } else {
    unique(det$transcript_type)
  }
  tx_types <- tx_types %||% default_enrichment_types(types_in_annot)
  run_fgsea_by_type(
    det = det,
    rank_id = "transcript_name",
    genesets_list = genesets_list,
    tx_to_gene = tibble::as_tibble(tx_to_gene)[, c(
      "transcript_name",
      "gene_name"
    )],
    map_id = "transcript_name",
    pval_cutoff = pval_cutoff,
    tx_types = tx_types,
    min_size = min_size,
    max_size = max_size,
    rank_by = rank_by,
    workers = workers,
    ...
  )
}
# nolint end

#' @rdname enrichment
#' @export
default_enrichment_types <- function(annotation_types = NULL) {
  processed_or_cds <- "protein_coding_CDS_not_defined"
  if (
    !is.null(annotation_types) &&
      "processed_transcript" %in% annotation_types &&
      !"protein_coding_CDS_not_defined" %in% annotation_types
  ) {
    processed_or_cds <- "processed_transcript"
  }
  list(
    protein_coding = "protein_coding",
    unproductive = c(
      "retained_intron",
      processed_or_cds,
      "nonsense_mediated_decay"
    ),
    retained_intron = "retained_intron",
    processed_or_cds = processed_or_cds,
    nonsense_mediated_decay = "nonsense_mediated_decay"
  ) |>
    stats::setNames(c(
      "protein_coding",
      "unproductive",
      "retained_intron",
      processed_or_cds,
      "nonsense_mediated_decay"
    ))
}

#' Abort unless `genesets_list` is a named list of gene symbol vectors
#'
#' Or a named list of such lists, see `geneset_collections()`.
#' @keywords internal
#' @noRd
check_genesets <- function(genesets_list) {
  geneset_collections(genesets_list)
  invisible(genesets_list)
}

#' The collections of `genesets_list`
#'
#' One collection is a named list of gene symbol vectors; several are a named
#' list of such lists. Returns a list of collections, without names for a
#' single collection given as such.
#' @keywords internal
#' @noRd
geneset_collections <- function(genesets_list) {
  named <- function(x) {
    !is.null(names(x)) && !anyNA(names(x)) && all(nzchar(names(x)))
  }
  # A gene set: a vector, an empty set, or a list of single values (as JSON
  # and YAML readers give)
  is_set <- function(set) {
    is.null(set) ||
      (is.atomic(set) && !is.list(set)) ||
      (is.list(set) &&
        is.null(names(set)) &&
        all(vapply(
          set,
          function(v) is.atomic(v) && length(v) == 1L,
          logical(1L)
        )))
  }
  is_collection <- function(x) {
    is.list(x) && named(x) && all(vapply(x, is_set, logical(1L)))
  }
  as_collection <- function(x) {
    lapply(x, function(set) as.character(unlist(set, use.names = FALSE)))
  }
  if (is_collection(genesets_list)) {
    return(list(as_collection(genesets_list)))
  }
  if (
    is.list(genesets_list) &&
      length(genesets_list) > 0L &&
      named(genesets_list) &&
      !anyDuplicated(names(genesets_list)) &&
      all(vapply(genesets_list, is_collection, logical(1L)))
  ) {
    return(lapply(genesets_list, as_collection))
  }
  cli::cli_abort(
    message = c(
      x = "{.arg genesets_list} must be a named list of gene symbol vectors.",
      i = "Every gene set needs a name. Several collections are a named list of such lists."
    ),
    class = "isoformic_invalid_genesets"
  )
}

#' Stack the results of every collection
#'
#' `results` has one table per collection of `geneset_collections()`. The
#' `collection` column is added only when several collections were given.
#' @keywords internal
#' @noRd
stack_collections <- function(collections, results) {
  if (is.null(names(collections))) {
    return(results[[1L]])
  }
  res <- lapply(seq_along(collections), function(i) {
    out <- results[[i]]
    out$collection <- rep(names(collections)[[i]], nrow(out))
    out
  })
  dplyr::bind_rows(res)
}

#' Run the fgsea tasks of an enrichment, serially or on workers
#'
#' Each task is a list of `pathways` and `stats`; the other arguments of
#' `fgsea::fgseaMultilevel()` are common to all tasks.
#' @keywords internal
#' @noRd
run_fgsea_tasks <- function(tasks, labels, min_size, max_size, workers, ...) {
  res <- run_tasks(
    tasks,
    fun = fgsea_multilevel,
    args = list(minSize = min_size, maxSize = max_size, eps = 0, ...),
    workers = workers,
    labels = labels
  )
  lapply(res, function(r) tibble::as_tibble(as.data.frame(r)))
}

#' Gene sets cut down to the IDs of a ranking
#'
#' `fgsea` matches every set against the names of the ranking and drops the
#' other members before it tests, so the result is the same; the sets keep
#' their names and order, and a set with no member left is kept empty.
#' @keywords internal
#' @noRd
sets_in_ranking <- function(sets, ids) {
  if (!length(sets)) {
    return(sets)
  }
  members <- unlist(sets, use.names = FALSE)
  set_index <- rep.int(seq_along(sets), lengths(sets))
  keep <- members %in% ids
  out <- split(
    members[keep],
    factor(set_index[keep], levels = seq_along(sets))
  )
  names(out) <- names(sets)
  out
}

#' Which transcript sets have a number of genes within the size limits
#'
#' Counts, for each term, the distinct genes with at least one transcript in
#' `ids` (the ranking of one category).
#'
#' @param sets_long Data frame with `term`, `gene_name` and `map_id` columns,
#'   one row per transcript of a set.
#' @returns A logical vector along `terms`.
#' @keywords internal
#' @noRd
sets_with_gene_count <- function(
  sets_long,
  map_id,
  ids,
  terms,
  min_size,
  max_size
) {
  ranked <- sets_long[[map_id]] %in% ids
  term <- match(sets_long$term[ranked], terms)
  gene <- sets_long$gene_name[ranked]
  first <- !duplicated(cbind(term, match(gene, gene)))
  n_genes <- tabulate(term[first], nbins = length(terms))
  n_genes >= min_size & n_genes <= max_size
}

#' The label of a task, for error messages
#' @keywords internal
#' @noRd
collection_label <- function(collections, i) {
  if (is.null(names(collections))) {
    "the gene set collection"
  } else {
    paste0("collection \"", names(collections)[[i]], "\"")
  }
}

#' A `fgsea` result without rows
#' @keywords internal
#' @noRd
empty_fgsea <- function() {
  tibble::tibble(
    pathway = character(),
    pval = numeric(),
    padj = numeric(),
    log2err = numeric(),
    ES = numeric(),
    NES = numeric(),
    size = integer(),
    leadingEdge = list(),
    experiment = character()
  )
}

#' Estimate over its standard deviation, for a ranking
#'
#' An estimate that is the same in every replicate (a standard deviation of
#' zero) is the most certain one: it gets the largest finite value of the
#' vector. No change in any replicate is zero.
#' @keywords internal
#' @noRd
certainty_z <- function(estimate, sd) {
  z <- estimate / sd
  z[!is.finite(z)] <- NA_real_
  certain <- sd %in% 0 & !estimate %in% c(0, NA)
  finite_z <- abs(z[!is.na(z)])
  z[certain] <- sign(estimate[certain]) * max(finite_z, 1)
  z[sd %in% 0 & estimate %in% 0] <- 0
  z
}

#' The statistic that ranks the rows of a differential expression table
#' @keywords internal
#' @noRd
enrichment_rank <- function(det, rank_by) {
  needed <- switch(
    rank_by,
    log2FC = "log2FC",
    signed_p = c("log2FC", "pvalue"),
    lfc_z = c("lfc_median", "lfc_sd")
  )
  missing_cols <- setdiff(needed, colnames(det))
  if (length(missing_cols)) {
    cli::cli_abort(
      message = c(
        x = "{.code rank_by = \"{rank_by}\"} needs the column{?s} {.field {missing_cols}}.",
        i = if (identical(rank_by, "lfc_z")) {
          "They come from {.fn de_confidence}."
        }
      ),
      class = if (identical(rank_by, "lfc_z")) {
        "isoformic_no_inferential_replicates"
      } else {
        "isoformic_invalid_de_table"
      }
    )
  }
  if (identical(rank_by, "log2FC")) {
    return(as.numeric(det$log2FC))
  }
  if (identical(rank_by, "lfc_z")) {
    return(certainty_z(as.numeric(det$lfc_median), as.numeric(det$lfc_sd)))
  }
  score <- -log10(as.numeric(det$pvalue))
  # A p-value of zero is the smallest of the table, not an infinite score
  finite <- score[is.finite(score)]
  score[is.infinite(score)] <- if (length(finite)) max(finite) else 1
  sign(as.numeric(det$log2FC)) * score
}

#' Run fgsea once per transcript category
#' @keywords internal
#' @noRd
run_fgsea_by_type <- function(
  det,
  rank_id,
  genesets_list,
  tx_to_gene,
  map_id,
  pval_cutoff,
  tx_types,
  min_size,
  max_size,
  rank_by = "log2FC",
  workers = 1L,
  ...
) {
  rlang::check_installed("fgsea")
  check_genesets(genesets_list)
  det$rank_stat <- enrichment_rank(det, rank_by)
  if (identical(rank_by, "lfc_z") && !any(is.finite(det$rank_stat))) {
    cli::cli_abort(
      message = c(
        x = "{.code rank_by = \"lfc_z\"} needs at least two inferential replicates.",
        i = "No transcript has a standard deviation of its fold-change."
      ),
      class = "isoformic_no_inferential_replicates"
    )
  }
  det <- det[!is.na(det$rank_stat) & !is.na(det[[rank_id]]), ]
  gene_map <- dplyr::distinct(tx_to_gene[, c(map_id, "gene_name")])
  collections <- geneset_collections(genesets_list)
  # One ranking per transcript category with at least two transcripts
  ranks_by_type <- lapply(names(tx_types), function(type_name) {
    sub <- det[det$transcript_type %in% tx_types[[type_name]], ]
    if (nrow(sub) < 2L) {
      return(NULL)
    }
    ranks <- stats::setNames(sub$rank_stat, sub[[rank_id]])
    ranks <- ranks[!duplicated(names(ranks))]
    sort(ranks)
  })
  names(ranks_by_type) <- names(tx_types)
  types <- names(tx_types)[!vapply(ranks_by_type, is.null, logical(1L))]
  # The requested count, not capped at the number of tasks: kept daemons
  # are sized on it (see run_tasks())
  workers <- resolve_workers(workers, Inf)
  on_workers <- min(workers, length(collections) * length(types)) > 1L
  # Tasks in the order of the serial loop: collections, then categories
  tasks <- list()
  task_collection <- integer()
  task_type <- character()
  labels <- character()
  for (i in seq_along(collections)) {
    sets <- collections[[i]]
    # Expand gene-symbol sets to transcript IDs once (vectorised join)
    sets_long <- tibble::tibble(
      term = rep(names(sets), lengths(sets)),
      gene_name = as.character(unlist(sets, use.names = FALSE))
    ) |>
      dplyr::inner_join(
        gene_map,
        by = "gene_name",
        relationship = "many-to-many"
      )
    tx_sets <- split(sets_long[[map_id]], sets_long$term)
    for (type_name in types) {
      # The size limits count the genes of a set that have a ranked
      # transcript of the category, not its transcripts
      in_size <- sets_with_gene_count(
        sets_long,
        map_id = map_id,
        ids = names(ranks_by_type[[type_name]]),
        terms = names(tx_sets),
        min_size = min_size,
        max_size = max_size
      )
      # On workers, only the transcripts of the ranking: fgsea drops the
      # others, and every task is serialized at dispatch, so it carries the
      # category's share only.
      pathways <- if (on_workers) {
        sets_in_ranking(tx_sets[in_size], names(ranks_by_type[[type_name]]))
      } else {
        tx_sets[in_size]
      }
      tasks[[length(tasks) + 1L]] <- list(
        pathways = pathways,
        stats = ranks_by_type[[type_name]]
      )
      task_collection <- c(task_collection, i)
      task_type <- c(task_type, type_name)
      labels <- c(
        labels,
        paste0(
          collection_label(collections, i),
          ", transcript category \"",
          type_name,
          "\""
        )
      )
    }
  }
  # The sets are already filtered on their gene counts
  fgsea_res <- run_fgsea_tasks(
    tasks,
    labels = labels,
    min_size = 1,
    max_size = Inf,
    workers = workers,
    ...
  )
  results <- lapply(seq_along(collections), function(i) {
    out <- lapply(which(task_collection == i), function(k) {
      res <- fgsea_res[[k]]
      if (!is.null(pval_cutoff)) {
        res <- res[!is.na(res$pval) & res$pval <= pval_cutoff, ]
      }
      res$experiment <- rep(task_type[[k]], nrow(res))
      res
    })
    out <- dplyr::bind_rows(out)
    if (ncol(out)) out else empty_fgsea()
  })
  res <- stack_collections(collections, results)
  res$rank_by <- rep(rank_by, nrow(res))
  dplyr::relocate(res, dplyr::any_of("collection"), .after = "rank_by")
}

#' Gene set enrichment on a gene table ranked by one of its columns
#'
#' The largest absolute value is kept for a gene symbol that is repeated.
#' @keywords internal
#' @noRd
gene_rank_enrichment <- function(
  table,
  rank_col,
  genesets_list,
  experiment,
  rank_by,
  min_size,
  max_size,
  workers = 1L,
  ...
) {
  table <- table[!is.na(table[[rank_col]]) & !is.na(table$gene_name), ]
  table <- table[order(-abs(table[[rank_col]])), ]
  table <- table[!duplicated(table$gene_name), ]
  ranks <- sort(stats::setNames(table[[rank_col]], table$gene_name))
  collections <- geneset_collections(genesets_list)
  workers <- resolve_workers(workers, Inf)
  on_workers <- min(workers, length(collections)) > 1L
  tasks <- lapply(collections, function(sets) {
    if (on_workers) {
      sets <- sets_in_ranking(sets, names(ranks))
    }
    list(pathways = sets, stats = ranks)
  })
  names(tasks) <- NULL
  results <- run_fgsea_tasks(
    tasks,
    labels = vapply(
      seq_along(collections),
      function(i) collection_label(collections, i),
      character(1L)
    ),
    min_size = min_size,
    max_size = max_size,
    workers = workers,
    ...
  )
  res <- stack_collections(collections, results)
  res$experiment <- rep(experiment, nrow(res))
  res$rank_by <- rep(rank_by, nrow(res))
  dplyr::relocate(res, dplyr::any_of("collection"), .after = "rank_by")
}

#' Run `fgsea::fgseaMultilevel()` serially unless asked otherwise
#'
#' Without `BPPARAM` or `nproc`, fgsea uses `BiocParallel::bpparam()`, which
#' starts socket workers (`SnowParam`) on Windows; they can hang R CMD check.
#' @keywords internal
#' @noRd
fgsea_multilevel <- function(...) {
  args <- list(...)
  if (is.null(args[["BPPARAM"]]) && is.null(args[["nproc"]])) {
    args[["BPPARAM"]] <- BiocParallel::SerialParam()
  }
  do.call(fgsea::fgseaMultilevel, args)
}

#' @rdname enrichment
#' @export
read_gmt <- function(path) {
  assert_file_exists(path)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  sets <- lapply(parts, function(p) unique(p[-(1:2)][nzchar(p[-(1:2)])]))
  names(sets) <- vapply(parts, function(p) p[[1L]], character(1L))
  sets
}

#' @rdname enrichment
#' @export
plot_enrichment <- function(
  enrichment_df,
  n = 30,
  experiments = NULL,
  nes_cutoff = 0,
  strip_prefix = "^[A-Z]+_",
  label_width = 50
) {
  data <- tibble::as_tibble(enrichment_df)
  required <- c("pathway", "NES", "padj", "size", "experiment")
  missing_cols <- setdiff(required, colnames(data))
  if (length(missing_cols)) {
    cli::cli_abort(
      message = c(
        x = "{.arg enrichment_df} is missing column{?s}: {.field {missing_cols}}."
      ),
      class = "isoformic_invalid_enrichment"
    )
  }
  if (!is.null(experiments)) {
    data <- data[data$experiment %in% experiments, ]
  }
  data <- data[!is.na(data$NES) & abs(data$NES) >= nes_cutoff, ]
  if (!nrow(data)) {
    cli::cli_abort(
      c(x = "No pathways to plot."),
      class = "isoformic_empty_plot_data"
    )
  }
  top <- data |>
    dplyr::group_by(.data$pathway) |>
    dplyr::summarise(
      min_padj = min(.data$padj, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$min_padj) |>
    dplyr::slice_head(n = n)
  data <- data[data$pathway %in% top$pathway, ]
  order_df <- data |>
    dplyr::group_by(.data$pathway) |>
    dplyr::summarise(mean_nes = mean(.data$NES), .groups = "drop") |>
    dplyr::arrange(.data$mean_nes)
  data$pathway <- factor(data$pathway, levels = order_df$pathway)
  data$experiment <- factor(
    data$experiment,
    levels = unique(enrichment_df$experiment)
  )
  limit <- max(abs(data$NES))
  pretty_label <- function(x) {
    x <- as.character(x)
    if (!is.null(strip_prefix)) {
      x <- sub(strip_prefix, "", x)
    }
    x <- gsub("_", " ", x, fixed = TRUE)
    vapply(
      strwrap(x, width = label_width, simplify = FALSE),
      paste,
      character(1L),
      collapse = "\n"
    )
  }
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$experiment,
      y = .data$pathway,
      color = .data$NES,
      size = .data$size
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient2(
      low = "#2166ac",
      mid = "grey85",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-limit, limit)
    ) +
    ggplot2::scale_size_area(max_size = 7) +
    ggplot2::scale_y_discrete(labels = pretty_label) +
    ggplot2::scale_x_discrete(labels = function(x) {
      gsub("_", " ", x, fixed = TRUE)
    }) +
    ggplot2::labs(x = NULL, y = NULL, color = "NES", size = "Set size") +
    isoformic_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8)
    ) +
    # The same pathway name in two collections is two results
    if ("collection" %in% colnames(data)) {
      ggplot2::facet_wrap(ggplot2::vars(.data$collection), scales = "free_y")
    }
}
