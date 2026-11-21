# WP-A step 1: the null splits of the placenta controls.
#
# Purpose: draw the label permutations of the 15 control units of the S0
# store (17 control runs; the A/B pairs PHNPR27A/B and PHNPR38A/B are one
# placenta each: S0 keeps the deeper member as the unit and the other as the
# sensitivity member). No preeclampsia run is ever tested. Arms hold units;
# fetal sex counts differ by at most 1 between arms (except in the
# unrestricted set) and no partition repeats.
#   7v7_louo   15 left-out units x 14 splits = 210. Delete-one-unit
#              jackknife: the splits without unit u are the 14 that leave it out.
#   7v7_mc     another 14 per left-out unit (210), used with 7v7_louo for the
#              28-per-unit Monte Carlo check of M1 and M2
#   5v5, 3v3   100 splits each (random units)
#   7v7_unres  30 splits without the sex rule (two per left-out unit)
# Each split records the imbalance between arms in depth (mean log2 processed
# fragments), mapping rate and number of male placentas.
#
# Inputs: <out>/store/qc.parquet and control_units.parquet
# (10_placenta_prepare.R), <out>/prereg.yaml.
# Env vars: ISOFORMIC_UE_OUT (required), ISOFORMIC_UE_PKG (ue_common.R).
# Outputs: <out>/null_placenta/splits/splits.parquet (one row per split) and
# split_samples.parquet (one row per split, member set and run; member set
# "primary" uses the S0 units, "sensitivity" swaps in the other member of
# each A/B pair).
#
# Resources: seconds, under 2 GB (loads the package snapshot).
# Run from the repository root:
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/20_null_splits.R

source("data-raw/uncertainty_evidence/null_common.R")
ue <- ue_setup("20_null_splits")
p <- wpa_paths()
dir.create(p$splits, recursive = TRUE, showWarnings = FALSE)
# Seed rule of prereg.yaml: base + 1000 * package index (A = 1) + unit index
seed <- ue$seed_base + 1000L * ue$prereg$seeds$package_index$A + 1L
set.seed(seed)

qc <- wpa_read(file.path(p$store, "qc.parquet"))
cu <- wpa_read(file.path(p$store, "control_units.parquet"))
units <- sort(cu$sample_id[cu$unit_primary])
stopifnot(length(units) == 15L)
# Sensitivity member of each unit (the other member of an A/B pair)
pairs <- ue$prereg$data$technical_pairs
sens_of <- stats::setNames(units, units)
for (pr in pairs) {
  prim <- intersect(pr, units)
  other <- setdiff(pr, prim)
  stopifnot(length(prim) == 1L, length(other) == 1L)
  sens_of[[prim]] <- other
}
sex <- stats::setNames(qc$fetal_sex, qc$sample_id)
# Checked 2026-11-18: each A/B pair has one female and one male fetus
# (XIST/chrY scores of S0), and their transcript correlation is not above
# that of unrelated control pairs (S0 ab_check), so the pairs look like
# opposite-sex twin placentas rather than technical replicates. The sex rule
# uses the primary members; the sensitivity member set records its own sex
# imbalance (sex_diff_sensitivity) and can break the rule.
sex_mismatch <- units[sex[units] != sex[sens_of[units]]]
if (length(sex_mismatch)) {
  wpa_log(
    "A/B members with a different fetal sex: ",
    paste(sex_mismatch, sens_of[sex_mismatch], sep = "/", collapse = ", ")
  )
}
male <- sex[units] == "male"
names(male) <- units
male_any <- stats::setNames(sex == "male", names(sex))
wpa_log("units: ", length(units), "; male: ", sum(male))

canon <- function(a, b) {
  a <- sort(a)
  b <- sort(b)
  x <- c(paste(a, collapse = ","), paste(b, collapse = ","))
  paste(sort(x), collapse = "|")
}
draw_splits <- function(pool, k, n, sex_rule, seen) {
  out <- list()
  tries <- 0L
  while (length(out) < n) {
    tries <- tries + 1L
    if (tries > 1e5) {
      stop("could not draw ", n, " splits", call. = FALSE)
    }
    pick <- sample(pool, 2L * k)
    a <- sort(pick[seq_len(k)])
    b <- sort(pick[k + seq_len(k)])
    if (sex_rule && abs(sum(male[a]) - sum(male[b])) > 1L) {
      next
    }
    key <- canon(a, b)
    if (key %in% seen) {
      next
    }
    seen <- c(seen, key)
    out[[length(out) + 1L]] <- list(a = a, b = b, key = key)
  }
  out
}

rows <- list()
seen <- character()
add <- function(design, left_out, rank, s) {
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  rows[[length(rows) + 1L]] <<- data.frame(
    design = design,
    size = paste0(length(s$a), "v", length(s$b)),
    left_out = left_out,
    rank_in_unit = rank,
    units_a = paste(s$a, collapse = ","),
    units_b = paste(s$b, collapse = ","),
    partition = s$key
  )
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  seen <<- c(seen, s$key)
}
for (u in units) {
  for (r in seq_len(28L)) {
    s <- draw_splits(setdiff(units, u), 7L, 1L, TRUE, seen)[[1L]]
    add(if (r <= 14L) "7v7_louo" else "7v7_mc", u, r, s)
  }
}
for (s in draw_splits(units, 5L, 100L, TRUE, seen)) {
  add("5v5", NA_character_, NA_integer_, s)
}
for (s in draw_splits(units, 3L, 100L, TRUE, seen)) {
  add("3v3", NA_character_, NA_integer_, s)
}
for (u in units) {
  for (s in draw_splits(setdiff(units, u), 7L, 2L, FALSE, seen)) {
    add("7v7_unres", u, NA_integer_, s)
  }
}
splits <- do.call(rbind, rows)
stopifnot(!anyDuplicated(splits$partition))
splits$split_id <- sprintf("s%04d", seq_len(nrow(splits)))
splits <- splits[, c("split_id", setdiff(colnames(splits), "split_id"))]

# Runs of each split and imbalance ------------------------------------------------
ss <- list()
for (i in seq_len(nrow(splits))) {
  ua <- strsplit(splits$units_a[[i]], ",", fixed = TRUE)[[1L]]
  ub <- strsplit(splits$units_b[[i]], ",", fixed = TRUE)[[1L]]
  for (ms in c("primary", "sensitivity")) {
    ra <- if (ms == "primary") ua else unname(sens_of[ua])
    rb <- if (ms == "primary") ub else unname(sens_of[ub])
    ss[[length(ss) + 1L]] <- data.frame(
      split_id = splits$split_id[[i]],
      member_set = ms,
      sample_id = c(ra, rb),
      unit = c(ua, ub),
      arm = rep(c("a", "b"), c(length(ra), length(rb)))
    )
  }
}
split_samples <- do.call(rbind, ss)
split_samples$swapped <- split_samples$sample_id != split_samples$unit
ldepth <- stats::setNames(log2(qc$num_processed), qc$sample_id)
pmap <- stats::setNames(qc$percent_mapped, qc$sample_id)
sp <- split_samples[split_samples$member_set == "primary", ]
imb <- do.call(
  rbind,
  lapply(split(sp, sp$split_id), function(d) {
    a <- d$sample_id[d$arm == "a"]
    b <- d$sample_id[d$arm == "b"]
    data.frame(
      split_id = d$split_id[[1L]],
      depth_diff = mean(ldepth[b]) - mean(ldepth[a]),
      mapping_diff = mean(pmap[b]) - mean(pmap[a]),
      male_a = sum(male[d$unit[d$arm == "a"]]),
      male_b = sum(male[d$unit[d$arm == "b"]])
    )
  })
)
imb$sex_diff <- imb$male_b - imb$male_a
ssens <- split_samples[split_samples$member_set == "sensitivity", ]
imb_s <- do.call(
  rbind,
  lapply(split(ssens, ssens$split_id), function(d) {
    data.frame(
      split_id = d$split_id[[1L]],
      sex_diff_sensitivity = sum(male_any[d$sample_id[d$arm == "b"]]) -
        sum(male_any[d$sample_id[d$arm == "a"]])
    )
  })
)
imb <- merge(imb, imb_s, by = "split_id")
splits <- merge(splits, imb, by = "split_id")
swapped_splits <- unique(split_samples$split_id[split_samples$swapped])
splits$has_sensitivity_member <- splits$split_id %in% swapped_splits
ue_write(splits, file.path(p$splits, "splits.parquet"), ue, seed)
ue_write(split_samples, file.path(p$splits, "split_samples.parquet"), ue, seed)
ue_write_session(p$splits, ue)
print(table(splits$design))
print(table(splits$design, abs(splits$sex_diff)))
wpa_log("splits written: ", nrow(splits))
