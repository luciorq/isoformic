#' Which transcripts can the equivalence classes tell apart?
#'
#' @description
#' Salmon summarises the reads as equivalence classes: sets of transcripts
#' with the number of fragments compatible with exactly that set. The
#' abundance of a transcript is determined by these counts only when the
#' classes separate it from the other transcripts. Two isoforms that always
#' appear together in the same classes, for example, can exchange fragments
#' without any change to the class counts.
#'
#' `equivalence_identifiability()` reads the classes (`salmon quant
#' --dumpEq`) and reports, for every transcript, whether it is identifiable
#' from them: whether its abundance is fixed once the total of every class is
#' fixed. This is a property of the classes observed in a sample, not of the
#' posterior draws: a transcript can be non-identifiable and still show a
#' small inferential variance (see Details).
#'
#' @details
#' With `M` the class by transcript membership matrix of one sample, a
#' transcript is identifiable when its unit vector lies in the row space of
#' `M`, that is, when no vector of the null space of `M` changes its
#' abundance. The function avoids the decomposition of the full matrix:
#'
#' 1. A class with a single unresolved transcript resolves it (a class of
#'    one transcript does so directly); repeated until nothing changes.
#' 2. Transcripts that are members of exactly the same classes cannot be
#'    told apart: they are marked non-identifiable, replaced by their sum,
#'    and step 1 is repeated.
#' 3. What is left splits into independent blocks. Blocks of up to
#'    `max_dense` transcripts are solved with a singular value
#'    decomposition. Larger blocks need the `Matrix` package: two random
#'    vectors are projected on the null space with conjugate gradients, and
#'    a transcript with a non-zero coordinate is non-identifiable.
#'
#' Limits:
#'
#' * Only class membership is used, not the conditional probabilities of
#'   `--dumpEqWeights`. Isoforms that share all their fragments have
#'   proportional weights, so they stay non-identifiable with the weights;
#'   other cases may be resolved by them.
#' * The annotation is taken as complete, and the classes are those observed
#'   in the sample: a transcript can be identifiable in a deep sample and not
#'   in a shallow one. Transcripts without fragments in a sample are not
#'   assessed in it.
#' * This is not the range of optima of Zheng et al. (2022), which needs
#'   their sufficient statistics and a linear program per transcript.
#' * Human samples have one block of tens of thousands of transcripts
#'   linked by multi-mapping fragments; expect minutes per sample.
#'
#' In the six GSE129540 samples of the quantification stability article,
#' the estimates of non-identifiable transcripts differed more between two
#' halves of the reads of a sample than those of identifiable transcripts of
#' the same expression: a median absolute log2 ratio of 2.4 against 0.7 at 20
#' to 100 reads. The posterior draws call three quarters of both kinds
#' uncertain, so the two measures are close to independent.
#'
#' @inheritParams inferential-io
#' @param max_dense Largest block (number of transcripts) solved with a
#'   dense decomposition. Larger blocks use the iterative solver.
#' @param max_iterations Maximum number of iterations of the iterative
#'   solver per block. The transcripts of a block that does not converge get
#'   `NA`.
#'
#' @returns A tibble with one row per transcript that has fragments in at
#'   least one sample: `transcript_id`, `n_samples` (samples in which it has
#'   fragments and was resolved), `n_identifiable` (in how many of them it is
#'   identifiable), `identifiable_fraction` and `identifiable` (`TRUE` when
#'   identifiable in all of them, `NA` when it was resolved in none).
#'
#' @references
#' Zheng H, Ma C, Kingsford C (2022). Deriving ranges of optimal estimated
#' transcript expression due to nonidentifiability. Journal of Computational
#' Biology 29:121-139. \doi{10.1089/cmb.2021.0444}
#'
#' @seealso [equivalence_support()] for the share of unique fragments,
#'   [add_uncertainty()] for the uncertainty measured by the posterior draws.
#'
#' @examples
#' # Two transcripts that always share their classes, and one with its own
#' dir <- file.path(tempdir(), "identifiability_example")
#' dir.create(file.path(dir, "aux_info"), recursive = TRUE, showWarnings = FALSE)
#' writeLines(
#'   c("3", "2", "tx1", "tx2", "tx3", "2\t0\t1\t10", "1\t2\t5"),
#'   file.path(dir, "aux_info", "eq_classes.txt")
#' )
#' equivalence_identifiability(c(sample1 = dir))
#'
#' @export
equivalence_identifiability <- function(
  paths,
  transcripts = NULL,
  ignore_after_bar = TRUE,
  max_dense = 400L,
  max_iterations = 50000L
) {
  rlang::check_required(paths)
  dirs <- resolve_quant_dirs(paths)
  if (identical(unique(dirs$kind), "oarfish")) {
    cli::cli_abort(
      message = c(
        x = "Identifiability needs the equivalence classes of Salmon.",
        i = "Oarfish output does not hold them."
      ),
      class = "isoformic_no_equivalence_classes"
    )
  }
  ids <- NULL
  n_resolved <- NULL
  n_identifiable <- NULL
  n_read <- 0L
  for (dir in dirs$dir) {
    ec <- read_eq_classes(dir)
    if (is.null(ec)) {
      next
    }
    if (isTRUE(ignore_after_bar)) {
      ec$transcripts <- sub("\\|.*$", "", ec$transcripts)
    }
    if (is.null(ids)) {
      ids <- ec$transcripts
      n_resolved <- integer(length(ids))
      n_identifiable <- integer(length(ids))
    }
    status <- identifiable_from_classes(
      tid = ec$members$tid,
      cls = ec$members$class,
      n_tx = length(ec$transcripts),
      max_dense = max_dense,
      max_iterations = max_iterations
    )
    at <- match(ec$transcripts, ids)
    ok <- !is.na(status) & !is.na(at)
    n_resolved[at[ok]] <- n_resolved[at[ok]] + 1L
    n_identifiable[at[ok]] <- n_identifiable[at[ok]] + as.integer(status[ok])
    n_read <- n_read + 1L
  }
  if (n_read == 0L) {
    cli::cli_abort(
      message = c(
        x = "No equivalence classes found in the Salmon output.",
        i = "Run {.code salmon quant --dumpEq} to write {.file aux_info/eq_classes.txt.gz}."
      ),
      class = "isoformic_no_equivalence_classes"
    )
  }
  res <- tibble::tibble(
    transcript_id = ids,
    n_samples = n_resolved,
    n_identifiable = n_identifiable,
    identifiable_fraction = n_identifiable / n_resolved,
    identifiable = n_identifiable == n_resolved
  )
  res <- res[res$n_samples > 0L, ]
  if (!is.null(transcripts)) {
    res <- res[res$transcript_id %in% transcripts, ]
  }
  res
}

#' Identifiability of every transcript of one sample
#'
#' @param tid,cls Transcript and class index of every class member.
#' @returns A logical vector of length `n_tx`: `NA` for transcripts in no
#'   class and for blocks that the iterative solver did not resolve.
#' @keywords internal
#' @noRd
identifiable_from_classes <- function(
  tid,
  cls,
  n_tx,
  max_dense = 400L,
  max_iterations = 50000L
) {
  seen <- tabulate(tid, n_tx) > 0L
  status <- rep(NA, n_tx)
  if (!length(tid)) {
    return(status)
  }
  n_cls <- max(cls)
  reduced <- reduce_class_system(tid, cls, n_tx, n_cls)
  col <- reduced$col
  cls <- reduced$cls
  # A column is the sum of the transcripts merged into it
  col_status <- rep(TRUE, n_tx)
  open <- !reduced$known[col]
  if (any(open)) {
    left <- sort(unique(col[open]))
    tt <- match(col[open], left)
    cc <- match(cls[open], unique(cls[open]))
    col_status[left] <- solve_class_blocks(
      tt,
      cc,
      length(left),
      max_dense = max_dense,
      max_iterations = max_iterations
    )
  }
  status[seen] <- col_status[reduced$rep_of[seen]] & !reduced$merged[seen]
  status
}

#' Peel resolved transcripts and merge those with identical classes
#' @keywords internal
#' @noRd
reduce_class_system <- function(tid, cls, n_tx, n_cls) {
  hash <- withr::with_seed(
    1L,
    list(stats::runif(n_cls), stats::runif(n_cls))
  )
  col <- tid
  known <- logical(n_tx)
  merged <- logical(n_tx)
  rep_of <- seq_len(n_tx)
  repeat {
    changed <- FALSE
    repeat {
      open <- !known[col]
      n_open <- sum_by_index(as.numeric(open), cls, n_cls)
      resolved <- unique(col[open & n_open[cls] == 1])
      if (!length(resolved)) {
        break
      }
      known[resolved] <- TRUE
      changed <- TRUE
    }
    open <- !known[col]
    h1 <- sum_by_index(hash[[1L]][cls[open]], col[open], n_tx)
    h2 <- sum_by_index(hash[[2L]][cls[open]], col[open], n_tx)
    ids <- which(h1 > 0)
    key <- paste(h1[ids], h2[ids])
    first <- ids[match(key, key)]
    dup <- first != ids
    if (any(dup)) {
      changed <- TRUE
      merged[rep_of %in% c(ids[dup], first[dup])] <- TRUE
      map <- seq_len(n_tx)
      map[ids] <- first
      rep_of <- map[rep_of]
      keep <- !col %in% ids[dup]
      col <- col[keep]
      cls <- cls[keep]
    }
    if (!changed) {
      break
    }
  }
  list(col = col, cls = cls, known = known, merged = merged, rep_of = rep_of)
}

#' Smallest value per group
#' @keywords internal
#' @noRd
min_by_index <- function(x, index, n) {
  o <- order(index, x)
  index <- index[o]
  first <- !duplicated(index)
  out <- rep(NA_integer_, n)
  out[index[first]] <- x[o][first]
  out
}

#' Connected blocks of a sparse 0/1 system and the identifiability of its
#' columns
#' @keywords internal
#' @noRd
solve_class_blocks <- function(tt, cc, n, max_dense, max_iterations) {
  n_cls <- max(cc)
  label <- seq_len(n)
  repeat {
    class_label <- min_by_index(label[tt], cc, n_cls)
    new <- pmin(label, min_by_index(class_label[cc], tt, n))
    new <- new[new]
    if (identical(new, label)) {
      break
    }
    label <- new
  }
  status <- rep(TRUE, n)
  rows <- split(seq_along(tt), label[tt])
  for (r in rows) {
    cols <- unique(tt[r])
    if (length(cols) < 2L) {
      next
    }
    j <- match(tt[r], cols)
    i <- match(cc[r], unique(cc[r]))
    status[cols] <- if (length(cols) <= max_dense) {
      dense_identifiable(i, j, length(cols))
    } else {
      sparse_identifiable(i, j, length(cols), max_iterations)
    }
  }
  status
}

#' @keywords internal
#' @noRd
dense_identifiable <- function(i, j, n) {
  a <- matrix(0, max(i), n)
  a[cbind(i, j)] <- 1
  a <- unique(a)
  dec <- svd(a, nu = 0L)
  rank <- sum(dec$d > max(dim(a)) * dec$d[[1L]] * .Machine$double.eps)
  if (rank == n) {
    return(rep(TRUE, n))
  }
  rowSums(dec$v[, seq_len(rank), drop = FALSE]^2) > 1 - 1e-8
}

#' Projection of random vectors on the null space, by conjugate gradients on
#' the least squares problem (CGLS)
#' @keywords internal
#' @noRd
sparse_identifiable <- function(i, j, n, max_iterations) {
  rlang::check_installed(
    "Matrix",
    reason = "to resolve large blocks of equivalence classes."
  )
  m <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(max(i), n))
  project <- function(z) {
    r <- z
    s <- as.numeric(m %*% r)
    p <- s
    gamma <- sum(s^2)
    target <- 1e-24 * gamma
    for (it in seq_len(max_iterations)) {
      q <- as.numeric(Matrix::crossprod(m, p))
      r <- r - (gamma / sum(q^2)) * q
      s <- as.numeric(m %*% r)
      gamma_new <- sum(s^2)
      if (gamma_new <= target) {
        return(r)
      }
      p <- s + (gamma_new / gamma) * p
      gamma <- gamma_new
    }
    NULL
  }
  z <- withr::with_seed(2L, list(stats::rnorm(n), stats::rnorm(n)))
  r1 <- project(z[[1L]])
  r2 <- if (is.null(r1)) NULL else project(z[[2L]])
  if (is.null(r2)) {
    cli::cli_inform(
      message = c(
        "!" = "A block of {n} transcripts was not resolved in {max_iterations} iterations.",
        i = "Its transcripts get {.val {NA}}; raise {.arg max_iterations}."
      ),
      class = "isoformic_identifiability_not_converged"
    )
    return(rep(NA, n))
  }
  pmax(abs(r1), abs(r2)) < 1e-7
}
