# Shared statistics ----------------------------------------------------------

#' Two-group test of a per-gene proportion from counts
#'
#' Closed-form quasi-binomial generalised linear model with one binary
#' covariate: the fitted proportions are the pooled proportions of each
#' group, the dispersion is the Pearson statistic of the per-sample counts
#' around them (floored at 1, the binomial value) and the log odds ratio
#' between the groups is tested with a t statistic on
#' `n_samples - 2` degrees of freedom. Proportions use a 0.5 continuity
#' correction so that genes at 0 or 1 remain testable.
#'
#' @param num Matrix genes x samples of numerator counts.
#' @param den Matrix genes x samples of totals.
#' @param groups Factor of length `ncol(num)` with exactly two levels; the
#'   effect is level 2 versus level 1.
#'
#' @returns A tibble with `proportion_<level1>`, `proportion_<level2>`,
#'   `delta`, `log_odds_ratio`, `se`, `statistic`, `dispersion`, `df`,
#'   `pvalue`, `padj`, one row per row of `num`.
#'
#' @keywords internal
#' @noRd
two_group_proportion_test <- function(num, den, groups) {
  lv <- levels(groups)
  if (length(lv) != 2L) {
    cli::cli_abort(
      message = c(
        x = "{.arg groups} must have exactly two levels, got {length(lv)}."
      ),
      class = "isoformic_group_levels"
    )
  }
  a <- groups == lv[[1L]]
  b <- groups == lv[[2L]]
  y_a <- rowSums(num[, a, drop = FALSE])
  y_b <- rowSums(num[, b, drop = FALSE])
  n_a <- rowSums(den[, a, drop = FALSE])
  n_b <- rowSums(den[, b, drop = FALSE])
  p_a <- ifelse(n_a > 0, y_a / n_a, NA_real_)
  p_b <- ifelse(n_b > 0, y_b / n_b, NA_real_)
  fitted <- matrix(NA_real_, nrow(num), ncol(num))
  fitted[, a] <- p_a
  fitted[, b] <- p_b
  var <- den * fitted * (1 - fitted)
  resid2 <- (num - den * fitted)^2 / var
  usable <- is.finite(resid2) & den > 0
  resid2[!usable] <- 0
  n_used <- rowSums(den > 0 & !is.na(fitted))
  df <- n_used - 2
  pearson <- rowSums(resid2)
  dispersion <- ifelse(df > 0, pmax(1, pearson / df), NA_real_)
  pc_a <- (y_a + 0.5) / (n_a + 1)
  pc_b <- (y_b + 0.5) / (n_b + 1)
  logit <- function(p) log(p / (1 - p))
  lor <- logit(pc_b) - logit(pc_a)
  se <- sqrt(
    dispersion *
      (1 /
        ((n_a + 1) * pc_a * (1 - pc_a)) +
        1 / ((n_b + 1) * pc_b * (1 - pc_b)))
  )
  statistic <- lor / se
  pvalue <- ifelse(
    df > 0 & is.finite(statistic),
    2 * stats::pt(-abs(statistic), df = df),
    NA_real_
  )
  res <- tibble::tibble(
    p1 = unname(p_a),
    p2 = unname(p_b),
    delta = unname(p_b - p_a),
    log_odds_ratio = unname(lor),
    se = unname(se),
    statistic = unname(statistic),
    dispersion = unname(dispersion),
    df = unname(df),
    pvalue = unname(pvalue),
    padj = stats::p.adjust(unname(pvalue), method = "BH")
  )
  names(res)[1:2] <- paste0("proportion_", lv)
  res
}

#' Summarise a numerator/denominator pair into group means and deltas
#' @keywords internal
#' @noRd
group_fraction_columns <- function(frac, groups, prefix) {
  lv <- levels(groups)
  out <- list()
  for (l in lv) {
    cols <- which(groups == l)
    v <- rowMeans(frac[, cols, drop = FALSE], na.rm = TRUE)
    v[is.nan(v)] <- NA_real_
    out[[paste0(prefix, "_", l)]] <- as.vector(v)
  }
  if (length(lv) >= 2L) {
    out[[paste0("delta_", prefix)]] <- out[[paste0(
      prefix,
      "_",
      lv[[length(lv)]]
    )]] -
      out[[paste0(prefix, "_", lv[[1L]])]]
  }
  out
}
