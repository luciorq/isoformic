#' Read and cut the transcript trees of TreeTerminus
#'
#' @description
#' TreeTerminus builds, from Salmon replicates, a forest of transcript trees:
#' the leaves are transcripts and every internal node is the sum of the
#' transcripts below it, built so that the inferential uncertainty falls as
#' the tree is ascended. Unlike the flat groups of Terminus, a tree lets the
#' resolution change from one isoform to the next: an isoform the reads
#' determine stays a leaf, one they do not is replaced by the smallest node
#' that they do determine.
#'
#' * `read_treeterminus()` reads the trees (Newick, one per line):
#'   `cluster_nwk.txt` for the trees of all samples (`treeterminus group
#'   --mean_inf true`), `group_nwk.txt` for one sample, or the output of
#'   `treeterminus consensus`. TreeTerminus labels the leaves with the 0-based
#'   position of the transcript in `quant.sf`; `salmon_dir` (any of the Salmon
#'   directories, which share the order) or `transcripts` turns them into
#'   transcript IDs.
#' * `cut_transcript_trees()` computes the inferential relative variance of
#'   every node from the replicates stored in an [IsoformicExperiment] (the
#'   draws of the transcripts below the node are summed), replaces every
#'   transcript whose own variance is at least `cutoff` by its lowest ancestor
#'   below `cutoff`, and returns the groups, which do not overlap (when two
#'   chosen nodes are nested, the larger one is kept). The result goes to
#'   `add_uncertainty(groups = )`, so that the confidence classes and
#'   `collapse_groups` work at the resolution the data support.
#'
#' Transcripts whose whole path to the root stays above the cutoff are
#' unresolved: by default they are left as they are (`unresolved =
#' "transcript"`), where the confidence classes still flag them as uncertain;
#' `unresolved = "root"` merges every tree that has one into a single group.
#'
#' @param path A Newick file of TreeTerminus, or a folder holding
#'   `cluster_nwk.txt` (or a single `group_nwk.txt`).
#' @param salmon_dir A Salmon output directory quantified against the same
#'   index, whose `quant.sf` gives the order of the transcripts.
#' @param transcripts Alternatively, the transcript names in `quant.sf` order.
#' @param ignore_after_bar Logical. Keep only the part of transcript names
#'   before the first `|`.
#' @param trees The result of `read_treeterminus()`.
#' @param x An [IsoformicExperiment] with inferential replicates (see
#'   [add_uncertainty()]).
#' @param cutoff Inferential relative variance below which a node is well
#'   determined (the `inf_rv_cutoff` of [add_uncertainty()] by default).
#' @param unresolved `"transcript"` or `"root"`, see above.
#'
#' @returns `read_treeterminus()`: a tibble with one row per node and leaf:
#'   `tree_id`, `node_id` (internal nodes; `NA` for leaves), `transcript_id`
#'   (leaves; `NA` for nodes) and `parent_id` (`NA` for the roots).
#'   `cut_transcript_trees()`: a tibble with `transcript_id`, `group_id`,
#'   `group_size`, `group_inf_rv` and `tree_id` for the transcripts that are
#'   grouped, with the cutoff in the `cutoff` attribute.
#'
#' @references
#' Singh NP, Love MI, Patro R (2023). TreeTerminus: creating transcript trees
#' using inferential replicate counts. *iScience* 26(6):106961.
#'
#' @seealso [read_terminus_groups()], [add_uncertainty()]
#'
#' @examples
#' nwk <- tempfile(fileext = ".txt")
#' writeLines(c("((0,1),2);", "(3,4);"), nwk)
#' read_treeterminus(nwk, transcripts = c("tA1", "tA2", "tA3", "tB1", "tB2"))
#'
#' @name treeterminus
#' @aliases read_treeterminus cut_transcript_trees
NULL

#' @rdname treeterminus
#' @export
read_treeterminus <- function(
  path,
  salmon_dir = NULL,
  transcripts = NULL,
  ignore_after_bar = TRUE
) {
  rlang::check_required(path)
  file <- treeterminus_file(as.character(path))
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) {
    cli::cli_abort(
      message = c(x = "{.file {file}} has no trees."),
      class = "isoformic_invalid_treeterminus"
    )
  }
  if (is.null(transcripts) && !is.null(salmon_dir)) {
    quant <- fs::path(salmon_dir, "quant.sf")
    assert_file_exists(quant, arg = "salmon_dir")
    transcripts <- read_tabular(quant, columns = "Name")$Name
  }
  if (!is.null(transcripts) && isTRUE(ignore_after_bar)) {
    transcripts <- sub("\\|.*$", "", transcripts)
  }
  parsed <- parse_newick_forest(lines)
  labels <- parsed$label[!is.na(parsed$label)]
  numeric_labels <- all(grepl("^[0-9]+$", labels))
  if (numeric_labels) {
    if (is.null(transcripts)) {
      cli::cli_abort(
        message = c(
          x = "The leaves of {.file {fs::path_file(file)}} are transcript positions.",
          i = "Give {.arg salmon_dir} or {.arg transcripts} to name them."
        ),
        class = "isoformic_invalid_treeterminus"
      )
    }
    position <- as.integer(parsed$label) + 1L
    if (any(position > length(transcripts), na.rm = TRUE)) {
      cli::cli_abort(
        message = c(
          x = paste0(
            "Leaf positions go up to {max(position, na.rm = TRUE) - 1L}, ",
            "but there are {length(transcripts)} transcripts."
          ),
          i = "Use the {.file quant.sf} of the same Salmon index."
        ),
        class = "isoformic_invalid_treeterminus"
      )
    }
    parsed$label <- transcripts[position]
  } else if (isTRUE(ignore_after_bar)) {
    parsed$label <- sub("\\|.*$", "", parsed$label)
  }
  res <- tibble::tibble(
    tree_id = parsed$tree,
    node_id = parsed$node,
    transcript_id = parsed$label,
    parent_id = parsed$parent
  )
  attr(res, "source") <- "treeterminus"
  res
}

#' @rdname treeterminus
#' @export
cut_transcript_trees <- function(
  trees,
  x,
  cutoff = NULL,
  unresolved = c("transcript", "root")
) {
  rlang::check_required(trees)
  rlang::check_required(x)
  unresolved <- rlang::arg_match(unresolved)
  cutoff <- cutoff %||% uncertainty_settings(x)$inf_rv_cutoff
  required <- c("tree_id", "node_id", "transcript_id", "parent_id")
  if (!is.data.frame(trees) || !all(required %in% names(trees))) {
    cli::cli_abort(
      message = c(x = "{.arg trees} must come from {.fn read_treeterminus}."),
      class = "isoformic_invalid_treeterminus"
    )
  }
  arr <- draws_array(x)
  n_sample <- dim(arr)[2L]
  n_draw <- dim(arr)[3L]
  # Transcripts by (sample within draw) columns, as the replicates are summed
  flat <- matrix(arr, nrow = dim(arr)[1L], dimnames = list(rownames(arr), NULL))
  leaves <- trees[!is.na(trees$transcript_id), ]
  leaves <- leaves[leaves$transcript_id %in% rownames(flat), ]
  nodes <- trees[!is.na(trees$node_id), ]
  if (!nrow(leaves)) {
    cli::cli_abort(
      message = c(
        x = "None of the transcripts of the trees have replicates in the object.",
        i = "Use the trees of the samples quantified in {.arg x}."
      ),
      class = "isoformic_invalid_treeterminus"
    )
  }
  defaults <- uncertainty_defaults()
  rv_of <- function(m) {
    rv <- matrix(0, nrow(m), n_sample)
    for (s in seq_len(n_sample)) {
      cols <- s + n_sample * (seq_len(n_draw) - 1L)
      v <- m[, cols, drop = FALSE]
      mu <- rowMeans(v)
      va <- if (n_draw > 1L) rowSums((v - mu)^2) / (n_draw - 1L) else 0
      rv[, s] <- pmax(defaults$shift, va - mu) / (mu + defaults$pseudocount)
    }
    rowMeans(rv)
  }

  # Node sums: leaves into their parents, then every node into its parent.
  # Node ids are given in pre-order, so children have larger ids.
  node_ids <- sort(unique(nodes$node_id))
  node_index <- match(node_ids, node_ids)
  parent_of <- nodes$parent_id[match(node_ids, nodes$node_id)]
  direct <- rowsum(
    flat[leaves$transcript_id, , drop = FALSE],
    leaves$parent_id,
    reorder = FALSE
  )
  sums <- matrix(0, length(node_ids), ncol(flat))
  sums[match(as.integer(rownames(direct)), node_ids), ] <- direct
  for (i in rev(node_index)) {
    p <- parent_of[[i]]
    if (!is.na(p)) {
      j <- match(p, node_ids)
      sums[j, ] <- sums[j, ] + sums[i, ]
    }
  }
  has_leaf <- rowSums(sums) > 0
  node_rv <- rep(NA_real_, length(node_ids))
  node_rv[has_leaf] <- rv_of(sums[has_leaf, , drop = FALSE])
  leaf_rv <- rv_of(flat[leaves$transcript_id, , drop = FALSE])

  # Lowest well-determined node on the way up from every node (top-down)
  good <- !is.na(node_rv) & node_rv < cutoff
  lowest_good <- rep(NA_integer_, length(node_ids))
  for (i in node_index) {
    p <- parent_of[[i]]
    lowest_good[[i]] <- if (good[[i]]) {
      node_ids[[i]]
    } else if (!is.na(p)) {
      lowest_good[[match(p, node_ids)]]
    } else {
      NA_integer_
    }
  }
  parent_idx <- match(leaves$parent_id, node_ids)
  chosen <- ifelse(leaf_rv < cutoff, NA_integer_, lowest_good[parent_idx])
  # Keep the largest of nested chosen nodes: the topmost chosen node on the path
  selected <- node_ids %in% chosen
  top <- rep(NA_integer_, length(node_ids))
  for (i in node_index) {
    p <- parent_of[[i]]
    up <- if (!is.na(p)) top[[match(p, node_ids)]] else NA_integer_
    top[[i]] <- if (!is.na(up)) {
      up
    } else if (selected[[i]]) {
      node_ids[[i]]
    } else {
      NA_integer_
    }
  }
  group_node <- top[parent_idx]
  root_of <- tapply(
    nodes$node_id[is.na(nodes$parent_id)],
    nodes$tree_id[is.na(nodes$parent_id)],
    function(v) v[[1L]]
  )
  if (identical(unresolved, "root")) {
    missing <- is.na(group_node) & leaf_rv >= cutoff
    trees_unresolved <- unique(leaves$tree_id[missing])
    in_tree <- leaves$tree_id %in% trees_unresolved
    group_node[in_tree] <- as.integer(root_of[as.character(leaves$tree_id[
      in_tree
    ])])
  }
  res <- tibble::tibble(
    transcript_id = leaves$transcript_id,
    tree_id = leaves$tree_id,
    node = group_node
  )
  res <- res[!is.na(res$node), ]
  size <- table(res$node)
  res <- res[res$node %in% as.integer(names(size)[size >= 2L]), ]
  res$group_size <- as.integer(size[as.character(res$node)])
  res$group_inf_rv <- node_rv[match(res$node, node_ids)]
  res$group_id <- paste0("tt", res$tree_id, "_n", res$node)
  res <- res[, c(
    "transcript_id",
    "group_id",
    "group_size",
    "group_inf_rv",
    "tree_id"
  )]
  attr(res, "source") <- "treeterminus"
  attr(res, "cutoff") <- cutoff
  res
}

#' Locate the Newick file of a TreeTerminus run
#' @keywords internal
#' @noRd
treeterminus_file <- function(path) {
  if (isTRUE(fs::is_file(path))) {
    return(path)
  }
  assert_dir_exists(path, arg = "path")
  for (name in c("cluster_nwk.txt", "group_nwk.txt")) {
    found <- fs::dir_ls(
      path,
      recurse = 1L,
      glob = paste0("*", name),
      type = "file"
    )
    if (length(found) == 1L) {
      return(as.character(found))
    }
    if (length(found) > 1L) {
      cli::cli_abort(
        message = c(
          x = "Several {.file {name}} files in {.path {path}}.",
          i = "Pass the file of the trees you want."
        ),
        class = "isoformic_invalid_treeterminus"
      )
    }
  }
  cli::cli_abort(
    message = c(
      x = "No TreeTerminus trees found in {.path {path}}.",
      i = "Expected {.file cluster_nwk.txt} or {.file group_nwk.txt}."
    ),
    class = "isoformic_invalid_treeterminus"
  )
}

#' Parse a forest of Newick trees (one per line) into parent pointers
#'
#' Internal nodes get increasing integer ids in pre-order, so a child always
#' has a larger id than its parent. Branch lengths and node labels are ignored.
#' @keywords internal
#' @noRd
parse_newick_forest <- function(lines) {
  tokens <- regmatches(
    lines,
    gregexpr("[(),;]|[^(),;:[:space:]]+(:[^(),;]*)?", lines)
  )
  n_tok <- sum(lengths(tokens))
  tree <- node <- parent <- integer(n_tok)
  label <- character(n_tok)
  is_leaf <- logical(n_tok)
  n <- 0L
  next_id <- 0L
  for (t in seq_along(tokens)) {
    stack <- integer()
    after_close <- FALSE
    for (tok in tokens[[t]]) {
      if (identical(tok, "(")) {
        next_id <- next_id + 1L
        n <- n + 1L
        tree[[n]] <- t
        node[[n]] <- next_id
        parent[[n]] <- if (length(stack)) {
          stack[[length(stack)]]
        } else {
          NA_integer_
        }
        stack <- c(stack, next_id)
        after_close <- FALSE
      } else if (identical(tok, ")")) {
        stack <- stack[-length(stack)]
        after_close <- TRUE
      } else if (tok %in% c(",", ";")) {
        after_close <- FALSE
      } else if (!after_close) {
        name <- sub(":.*$", "", tok)
        if (!length(stack)) {
          cli::cli_abort(
            message = c(x = "Tree {t} is not valid Newick."),
            class = "isoformic_invalid_treeterminus"
          )
        }
        n <- n + 1L
        tree[[n]] <- t
        node[[n]] <- NA_integer_
        parent[[n]] <- stack[[length(stack)]]
        label[[n]] <- name
        is_leaf[[n]] <- TRUE
      }
    }
    if (length(stack)) {
      cli::cli_abort(
        message = c(x = "Tree {t} has unbalanced parentheses."),
        class = "isoformic_invalid_treeterminus"
      )
    }
  }
  keep <- seq_len(n)
  list(
    tree = tree[keep],
    node = node[keep],
    parent = parent[keep],
    label = ifelse(is_leaf[keep], label[keep], NA_character_)
  )
}
