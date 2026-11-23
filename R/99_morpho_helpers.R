###############################################################################
# 99_morpho_helpers.R
# Small helper functions shared by the morphology and life-history scripts.
# Sourced at the top of the other scripts; you don't run this file on its own.
###############################################################################

# -----------------------------------------------------------------------------
# write_sheet(): write a data frame to one sheet of the supplementary workbook.
#
# The old scripts did loadWorkbook() + addWorksheet(), which fails with
# "sheet already exists" the second time you run them. This creates the
# workbook if it is missing and replaces the sheet if it already exists, so
# every script can be re-run as often as you like.
# -----------------------------------------------------------------------------
write_sheet <- function(df, sheet, path = "res/supplementary_table.xlsx") {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  wb <- if (file.exists(path)) openxlsx::loadWorkbook(path) else openxlsx::createWorkbook()
  if (sheet %in% names(wb)) openxlsx::removeWorksheet(wb, sheet)
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, df)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  invisible(path)
}

# -----------------------------------------------------------------------------
# align_to_shapes(): put the metadata rows in the same order as the shapes.
#
# geomorph.data.frame() pairs shapes and metadata purely by position, so if the
# CSV and the TPS file list fish in a different order, every fish silently gets
# someone else's ecotype and sex. Here we match on the specimen ID instead, in
# three steps:
#   1. exact match of the IDs;
#   2. "loose" match: ignore EARLY/LATE, case, spaces and punctuation, so the
#      TPS "NOT_1" matches the CSV "NOT_1EARLY" (only if that is unambiguous);
#   3. otherwise keep row order (as the old scripts did), but check row by row
#      that the loose IDs agree, and list any rows where they don't.
# -----------------------------------------------------------------------------
loose_id <- function(x) {
  toupper(gsub("[^A-Za-z0-9]", "", gsub("EARLY|LATE", "", x, ignore.case = TRUE)))
}

align_to_shapes <- function(data, shape_ids, id_col = "sample_id") {
  ids <- as.character(data[[id_col]])

  # 1. exact
  if (all(shape_ids %in% ids) && !anyDuplicated(ids)) {
    message("Metadata matched to shapes by specimen ID (", length(shape_ids), " fish).")
    return(data[match(shape_ids, ids), , drop = FALSE])
  }

  # 2. loose (only if every loose ID is unique on both sides)
  ls <- loose_id(shape_ids)
  ld <- loose_id(ids)
  if (!anyDuplicated(ls) && !anyDuplicated(ld) && all(ls %in% ld)) {
    message("Metadata matched to shapes by specimen ID, ignoring EARLY/LATE, case and ",
            "punctuation (e.g. TPS '", shape_ids[1], "' = CSV '", ids[match(ls[1], ld)], "').")
    return(data[match(ls, ld), , drop = FALSE])
  }

  # 3. row order, checked
  if (nrow(data) != length(shape_ids)) {
    stop("The CSV has ", nrow(data), " rows but there are ", length(shape_ids),
         " shapes, and the IDs can't be matched (TPS: ", paste(head(shape_ids, 3), collapse = ", "),
         "; CSV: ", paste(head(ids, 3), collapse = ", "), ").", call. = FALSE)
  }
  bad <- which(ls != ld)
  if (length(bad) == 0) {
    message("IDs can't be matched uniquely (e.g. the same number used for an EARLY and a LATE fish), ",
            "so row order is kept. Row-by-row check: every row pairs the same fish number ",
            "in the TPS and the CSV.")
  } else {
    warning(length(bad), " of ", length(shape_ids), " rows pair DIFFERENT fish in the TPS and the CSV. ",
            "First mismatches (TPS vs CSV): ",
            paste0(shape_ids[head(bad, 5)], " vs ", ids[head(bad, 5)], collapse = "; "),
            ". Fix the order or the IDs before trusting the shape results.", call. = FALSE)
  }
  data
}

# -----------------------------------------------------------------------------
# lda_loo(): linear discriminant analysis with leave-one-out cross-validation.
#
# MASS::lda(..., CV = TRUE) predicts each fish from a model fitted to all the
# other fish, so no fish is ever classified by a model that saw it.
#
# Balanced accuracy = mean of the per-class success rates. Plain accuracy
# rewards always guessing the bigger group (there are ~2x more fall-run fish);
# balanced accuracy doesn't, and its chance level is always 0.5 with 2 groups.
#
# prior: by default lda() sets the prior to the group proportions, which pushes
# predictions toward the bigger group. We use equal priors so the result
# reflects the traits, not the sample sizes.
# -----------------------------------------------------------------------------
lda_loo <- function(X, group, prior = NULL) {
  group <- droplevels(factor(group))
  if (is.null(prior)) prior <- rep(1 / nlevels(group), nlevels(group))
  fit <- MASS::lda(x = as.matrix(X), grouping = group, prior = prior, CV = TRUE)
  conf <- table(Observed = group, Predicted = factor(fit$class, levels = levels(group)))
  sens <- diag(prop.table(conf, 1))           # proportion correct within each group
  list(confusion = conf,
       sensitivity = sens,
       balanced_accuracy = mean(sens),
       accuracy = sum(diag(conf)) / sum(conf))
}

# Permutation test for balanced accuracy: shuffle the group labels many times,
# redo the whole leave-one-out LDA each time, and ask how often a shuffled
# dataset does at least as well as the real one. (1 + count) / (n + 1) avoids a
# p-value of exactly 0 (Phipson & Smyth 2010).
lda_perm_test <- function(X, group, observed, n_perm = 1000, seed = 123) {
  set.seed(seed)
  perm <- replicate(n_perm, lda_loo(X, sample(group))$balanced_accuracy)
  (1 + sum(perm >= observed)) / (n_perm + 1)
}

# -----------------------------------------------------------------------------
# broken_stick(): how many principal components to keep.
#
# If total variance were split at random into p pieces, the k-th largest piece
# would on average get (1/p) * sum_{i=k..p} 1/i of the total. We keep the
# leading PCs that explain more than that expectation, stopping at the first
# one that doesn't.
# -----------------------------------------------------------------------------
broken_stick <- function(eigenvalues) {
  p <- length(eigenvalues)
  observed <- eigenvalues / sum(eigenvalues)
  expected <- rev(cumsum(1 / (p:1))) / p
  keep <- observed > expected
  n_keep <- if (all(keep)) p else which(!keep)[1] - 1
  list(n = n_keep,
       table = data.frame(PC = seq_len(p), observed = observed, expected = expected))
}

# -----------------------------------------------------------------------------
# wilks_table(): multivariate tests (Wilks' lambda) from car::Manova().
#
# car::Manova(type = "II") gives Type II hypothesis (H) and error (E)
# sums-of-squares-and-products matrices for each term. Wilks' lambda is
# prod(1 / (1 + eigenvalues of E^-1 H)); the F approximation below is Rao's,
# the same formula base R uses in summary.manova().
# -----------------------------------------------------------------------------
wilks_table <- function(av) {
  E <- av$SSPE
  do.call(rbind, lapply(names(av$SSP), function(term) {
    H   <- av$SSP[[term]]
    eig <- Re(eigen(qr.coef(qr(E), H), symmetric = FALSE, only.values = TRUE)$values)
    q   <- av$df[[term]]
    df_res <- av$error.df
    p   <- length(eig)
    wilks <- prod(1 / (1 + eig))
    tmp1 <- df_res - 0.5 * (p - q + 1)
    tmp2 <- (p * q - 2) / 4
    tmp3 <- p^2 + q^2 - 5
    tmp3 <- if (tmp3 > 0) sqrt(((p * q)^2 - 4) / tmp3) else 1
    df1 <- p * q
    df2 <- tmp1 * tmp3 - 2 * tmp2
    Fval <- ((wilks^(-1 / tmp3) - 1) * df2) / df1
    data.frame(Term = term, Df = q, Wilks = wilks, approx_F = Fval,
               num_Df = df1, den_Df = df2,
               p_value = pf(Fval, df1, df2, lower.tail = FALSE))
  }))
}
