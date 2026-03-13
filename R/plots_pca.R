# R/plots_pca.R
#
# Functions for Principal Component Analysis (PCA) and reverse scree plot
# visualization of py-GC-MS concentration data.
#
# Two complementary PCA modes are provided:
#
#   1. Per-batch PCA (export_pca_by_batch):
#      Runs a separate PCA for each batch using all replicate types
#      (Unknown, QC, Standard, Blank). Useful for detecting within-batch
#      outliers, assessing replicate clustering by type, and identifying
#      analytes that drive within-batch variance.
#
#   2. Pooled Unknown-only PCA (export_pooled_unknown_pca):
#      Runs a single PCA across all batches using only Unknown replicates.
#      Points are colored by Batch_ID to reveal between-batch clustering,
#      systematic batch effects, and project-wide sample groupings.
#
# PCA details:
#   - Input: ug/g concentration values (Norm_ug_Per_Gram or Unnorm_ug_Per_Gram)
#   - Negative values are clipped to 0 before PCA (non-negative constraint)
#   - Data are centered and scaled (prcomp center = TRUE, scale. = TRUE)
#     so that analytes with large absolute concentrations do not dominate
#     the variance decomposition
#   - NA handling: median imputation (default) replaces missing values with
#     the column median; alternatively, rows with any NA are dropped ("drop")
#   - Zero-variance analytes are removed before fitting (prcomp requires
#     at least some variance to scale)
#   - Minimum 2 analyte features and 2 sample replicates required per PCA
#
# Scree plots:
#   A "reverse scree" layout is used — components are ordered from smallest
#   to largest variance explained (left to right). This is a deliberate
#   design choice that makes it easier to identify the point where explained
#   variance begins to level off (the "elbow"), which can be obscured in a
#   traditional descending scree when many small components are compressed
#   on the right side.
#
# Output layout: portrait PDF (8.5 × 11 inches), PCA biplot stacked above
#   reverse scree, one page per batch.
#
# Output subfolder: pca_output/
#
# Contents:
#   Internal helpers (not exported):
#     prepare_pca_by_batch()              — preps and runs per-batch PCA
#     prepare_pooled_unknown_pca()        — preps and runs pooled Unknown PCA
#     make_pca_plot_for_batch()           — builds per-batch PCA biplot
#     make_reverse_scree_plot_for_batch() — builds per-batch scree plot
#     make_pooled_unknown_pca_plot()      — builds pooled PCA biplot
#     make_pooled_unknown_reverse_scree_plot() — builds pooled scree plot
#
#   Exported pipeline functions:
#     export_pca_by_batch()         — per-batch PCA + scree PDFs
#     export_pooled_unknown_pca()   — pooled Unknown PCA + scree PDFs


# ══ Internal helpers: data preparation ═══════════════════════════════════════

# ── prepare_pca_by_batch ──────────────────────────────────────────────────────
#
# Prepares concentration data and runs PCA separately for each batch.
# Returns PCA scores and scree data for all batches that had sufficient
# data to fit a model.
#
# @param df        A data frame from calculate_concentrations().
# @param value_col Character. Concentration column: "Norm_ug_Per_Gram" or
#                  "Unnorm_ug_Per_Gram".
# @param na_method Character. "median" (impute NAs with column median) or
#                  "drop" (remove rows with any NA). Default "median".
# @param min_features Integer. Minimum number of analyte columns required.
# @param min_samples  Integer. Minimum number of replicate rows required.
# @return A named list with elements:
#         $scores — data frame of PC1/PC2 scores with metadata
#         $scree  — data frame of variance explained per component
#
#' @keywords internal
prepare_pca_by_batch <- function(df, value_col,
                                 na_method    = c("median", "drop"),
                                 min_features = 2,
                                 min_samples  = 2) {
  na_method <- match.arg(na_method)
  stopifnot(all(c("Batch_ID", "Replicate_Name", "Replicate_Type",
                  "Molecule_Prefix") %in% names(df)))
  stopifnot(value_col %in% names(df))

  # Clip negatives to 0 — PCA on concentration data should not have
  # negative values; they represent below-blank measurements
  df_sanitized <- df %>%
    dplyr::mutate(
      !!rlang::sym(value_col) := ifelse(!!rlang::sym(value_col) < 0, 0,
                                        !!rlang::sym(value_col))
    )

  # Pivot to wide format: rows = replicates, columns = analytes
  df_wide <- df_sanitized %>%
    dplyr::select(Batch_ID, Replicate_Name, Replicate_Type,
                  Molecule_Prefix, !!rlang::sym(value_col)) %>%
    tidyr::pivot_wider(
      names_from  = Molecule_Prefix,
      values_from = !!rlang::sym(value_col)
    ) %>%
    dplyr::arrange(Batch_ID, Replicate_Name)

  batch_list  <- split(df_wide, df_wide$Batch_ID)
  scores_list <- list()
  scree_list  <- list()

  for (bid in names(batch_list)) {
    batch_df <- batch_list[[bid]]
    meta     <- batch_df[, c("Batch_ID", "Replicate_Name", "Replicate_Type")]
    mat      <- batch_df %>%
      dplyr::select(-Batch_ID, -Replicate_Name, -Replicate_Type)

    # Drop analytes that are entirely NA for this batch
    mat <- mat[, colSums(!is.na(mat)) > 0, drop = FALSE]
    if (ncol(mat) < min_features) next

    # NA handling
    if (na_method == "median") {
      mat <- as.data.frame(lapply(mat, function(col) {
        m <- stats::median(col, na.rm = TRUE)
        if (is.finite(m)) col[is.na(col)] <- m
        col
      }))
      # Re-clip negatives introduced by imputation
      mat <- as.data.frame(lapply(mat, function(col) ifelse(col < 0, 0, col)))
    }

    # Remove rows with remaining NA (those not resolved by median imputation)
    complete_rows <- stats::complete.cases(mat)
    mat_complete  <- mat[complete_rows, , drop = FALSE]
    meta_complete <- meta[complete_rows, , drop = FALSE]

    # Drop zero-variance analytes — prcomp(scale. = TRUE) requires variance > 0
    if (ncol(mat_complete) == 0) next
    var_cols     <- sapply(mat_complete, function(x) stats::var(x, na.rm = TRUE))
    mat_complete <- mat_complete[, var_cols > 0, drop = FALSE]

    if (ncol(mat_complete) < min_features ||
        nrow(mat_complete) < min_samples) next

    # Run PCA: center and scale so all analytes contribute equally
    pca      <- stats::prcomp(mat_complete, center = TRUE, scale. = TRUE)
    var_expl <- pca$sdev^2 / sum(pca$sdev^2)

    scores_list[[bid]] <- tibble::tibble(
      Batch_ID       = meta_complete$Batch_ID,
      Replicate_Name = meta_complete$Replicate_Name,
      Replicate_Type = meta_complete$Replicate_Type,
      PC1            = pca$x[, 1],
      PC2            = if (ncol(pca$x) >= 2) pca$x[, 2] else NA_real_,
      PC1_var_pct    = round(if (length(var_expl) >= 1) var_expl[1] * 100
                             else NA_real_, 1),
      PC2_var_pct    = round(if (length(var_expl) >= 2) var_expl[2] * 100
                             else NA_real_, 1)
    )

    # Reverse scree: arrange smallest → largest for elbow identification
    scree_list[[bid]] <- tibble::tibble(
      Batch_ID     = unique(meta_complete$Batch_ID)[1],
      Component    = paste0("PC", seq_along(var_expl)),
      Variance_Pct = round(var_expl * 100, 2)
    ) %>%
      dplyr::arrange(Variance_Pct) %>%
      dplyr::mutate(Component = factor(Component, levels = Component))
  }

  list(
    scores = dplyr::bind_rows(scores_list),
    scree  = dplyr::bind_rows(scree_list)
  )
}


# ── prepare_pooled_unknown_pca ────────────────────────────────────────────────
#
# Prepares Unknown-only replicate data pooled across all batches and runs
# a single PCA. Returns scores colored by Batch_ID and scree data.
#
# @param df        A data frame from calculate_concentrations().
# @param value_col Character. Concentration column to use.
# @param na_method Character. "median" or "drop".
# @param min_features Integer. Minimum analyte columns required.
# @param min_samples  Integer. Minimum replicate rows required.
# @return A named list with $scores and $scree data frames.
#
#' @keywords internal
prepare_pooled_unknown_pca <- function(df, value_col,
                                       na_method    = c("median", "drop"),
                                       min_features = 2,
                                       min_samples  = 2) {
  na_method <- match.arg(na_method)
  stopifnot(all(c("Batch_ID", "Replicate_Name", "Replicate_Type",
                  "Molecule_Prefix") %in% names(df)))
  stopifnot(value_col %in% names(df))

  empty <- list(scores = tibble::tibble(), scree = tibble::tibble())

  # Filter to Unknown replicates only and clip negatives
  df_unknown <- df %>%
    dplyr::filter(Replicate_Type == "Unknown") %>%
    dplyr::mutate(
      value = !!rlang::sym(value_col),
      value = dplyr::if_else(value < 0, 0, value)
    )

  if (nrow(df_unknown) == 0) {
    warning("No Unknown replicates found. Skipping pooled Unknown PCA.")
    return(empty)
  }

  # Pivot wide: use mean to resolve any duplicate replicate/analyte pairs
  wide <- df_unknown %>%
    dplyr::select(Batch_ID, Replicate_Name, Molecule_Prefix, value) %>%
    tidyr::pivot_wider(
      names_from  = Molecule_Prefix,
      values_from = value,
      values_fn   = list(value = mean)
    ) %>%
    dplyr::arrange(Batch_ID, Replicate_Name)

  meta <- wide %>% dplyr::select(Batch_ID, Replicate_Name)
  mat  <- wide %>% dplyr::select(-Batch_ID, -Replicate_Name)

  # Drop all-NA analytes
  mat <- mat[, colSums(!is.na(mat)) > 0, drop = FALSE]
  if (ncol(mat) < min_features) {
    warning("Insufficient analyte features for pooled Unknown PCA.")
    return(empty)
  }

  # NA handling
  if (na_method == "median") {
    mat <- as.data.frame(lapply(mat, function(col) {
      m <- stats::median(col, na.rm = TRUE)
      if (is.finite(m)) col[is.na(col)] <- m
      col
    }))
    mat <- as.data.frame(lapply(mat, function(col) ifelse(col < 0, 0, col)))
  }

  complete_rows <- stats::complete.cases(mat)
  mat_complete  <- mat[complete_rows, , drop = FALSE]
  meta_complete <- meta[complete_rows, , drop = FALSE]

  if (ncol(mat_complete) == 0) {
    warning("No analytes remain after filtering.")
    return(empty)
  }

  var_cols     <- sapply(mat_complete, function(x) stats::var(x, na.rm = TRUE))
  mat_complete <- mat_complete[, var_cols > 0, drop = FALSE]

  if (ncol(mat_complete) < min_features) {
    warning("Insufficient features after variance filtering.")
    return(empty)
  }
  if (nrow(mat_complete) < min_samples) {
    warning("Insufficient Unknown samples for pooled PCA.")
    return(empty)
  }

  pca      <- stats::prcomp(mat_complete, center = TRUE, scale. = TRUE)
  var_expl <- pca$sdev^2 / sum(pca$sdev^2)

  scores <- tibble::tibble(
    Batch_ID       = meta_complete$Batch_ID,
    Replicate_Name = meta_complete$Replicate_Name,
    PC1            = pca$x[, 1],
    PC2            = if (ncol(pca$x) >= 2) pca$x[, 2] else NA_real_,
    PC1_var_pct    = round(if (length(var_expl) >= 1) var_expl[1] * 100
                           else NA_real_, 1),
    PC2_var_pct    = round(if (length(var_expl) >= 2) var_expl[2] * 100
                           else NA_real_, 1)
  )

  scree <- tibble::tibble(
    Component    = paste0("PC", seq_along(var_expl)),
    Variance_Pct = round(var_expl * 100, 2)
  ) %>%
    dplyr::arrange(Variance_Pct) %>%
    dplyr::mutate(Component = factor(Component, levels = Component))

  list(scores = scores, scree = scree)
}


# ══ Internal helpers: plot builders ══════════════════════════════════════════

# ── make_pca_plot_for_batch ───────────────────────────────────────────────────
#
# Builds a PCA biplot for a single batch. Points are colored by Replicate_Type
# and labeled with Replicate_Name using ggrepel for non-overlapping labels.
#
# @param scores_df_batch A data frame of PCA scores for one batch.
# @param title           Character. Plot title.
# @param base_size       Numeric. Base font size. Default 10.
# @return A ggplot object, or NULL if scores_df_batch is empty.
#
#' @keywords internal
make_pca_plot_for_batch <- function(scores_df_batch, title, base_size = 10) {
  if (nrow(scores_df_batch) == 0) return(NULL)

  pc1_lab <- paste0("PC1 (", unique(scores_df_batch$PC1_var_pct)[1], "%)")
  pc2_lab <- paste0("PC2 (", unique(scores_df_batch$PC2_var_pct)[1], "%)")

  ggplot2::ggplot(
    scores_df_batch,
    ggplot2::aes(x = PC1, y = PC2,
                 color = Replicate_Type, label = Replicate_Name)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    ggplot2::geom_point(size = 1.8, alpha = 0.75) +
    ggrepel::geom_text_repel(
      min.segment.length = 0, size = 1.8, fontface = "bold",
      max.overlaps = 300, box.padding = 0.1,
      point.padding = 0.05, force = 0.3
    ) +
    ggplot2::labs(
      title = title,
      x     = pc1_lab,
      y     = pc2_lab,
      color = "Replicate Type"
    ) +
    theme_minimal_bold(base_size = base_size)
}


# ── make_reverse_scree_plot_for_batch ─────────────────────────────────────────
#
# Builds a reverse scree bar chart for one batch. Components are ordered
# smallest to largest variance explained (left to right) to highlight
# the elbow more clearly than a traditional descending scree.
#
# @param scree_df_batch A data frame of scree data for one batch.
# @param title          Character. Plot title.
# @param base_size      Numeric. Base font size. Default 10.
# @return A ggplot object, or NULL if scree_df_batch is empty.
#
#' @keywords internal
make_reverse_scree_plot_for_batch <- function(scree_df_batch, title,
                                              base_size = 10) {
  if (nrow(scree_df_batch) == 0) return(NULL)

  ggplot2::ggplot(scree_df_batch,
                  ggplot2::aes(x = Component, y = Variance_Pct)) +
    ggplot2::geom_col(fill = "steelblue", color = "black",
                      linewidth = 1, width = 0.7) +
    ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 1) +
    ggplot2::geom_text(
      ggplot2::aes(label = Variance_Pct),
      vjust = -0.25, size = 3, fontface = "bold"
    ) +
    ggplot2::labs(
      title = title,
      x     = "Principal Components (reverse: smallest \u2192 largest)",
      y     = "Variance Explained (%)"
    ) +
    theme_minimal_bold(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::coord_cartesian(clip = "off")
}


# ── make_pooled_unknown_pca_plot ──────────────────────────────────────────────
#
# Builds a PCA biplot for the pooled Unknown-only analysis. Points are
# colored by Batch_ID (not Replicate_Type) to reveal between-batch effects.
# Labels show "Replicate_Name (Batch_ID)" when use_combined_labels = TRUE.
#
# @param scores_df          A data frame of pooled PCA scores.
# @param title              Character. Plot title.
# @param base_size          Numeric. Base font size.
# @param use_combined_labels Logical. If TRUE, labels show
#                            "Replicate_Name (Batch_ID)". Default TRUE.
# @param legend_position    Character. ggplot2 legend position.
# @return A ggplot object, or NULL if scores_df is empty.
#
#' @keywords internal
make_pooled_unknown_pca_plot <- function(scores_df, title,
                                         base_size           = 10,
                                         use_combined_labels = TRUE,
                                         legend_position     = c("right",
                                                                 "bottom",
                                                                 "top",
                                                                 "left",
                                                                 "none")) {
  if (nrow(scores_df) == 0) return(NULL)
  legend_position <- match.arg(legend_position)

  pc1_lab <- paste0("PC1 (", unique(scores_df$PC1_var_pct)[1], "%)")
  pc2_lab <- paste0("PC2 (", unique(scores_df$PC2_var_pct)[1], "%)")

  # Combined labels make it easier to identify which batch each point belongs
  # to without relying solely on color (useful for colorblind accessibility)
  if (use_combined_labels) {
    scores_df <- scores_df %>%
      dplyr::mutate(Label = paste0(Replicate_Name, " (", Batch_ID, ")"))
    label_col <- "Label"
  } else {
    label_col <- "Replicate_Name"
  }

  ggplot2::ggplot(
    scores_df,
    ggplot2::aes(x = PC1, y = PC2,
                 color = Batch_ID, label = !!rlang::sym(label_col))
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    ggplot2::geom_point(size = 1.8, alpha = 0.75) +
    ggrepel::geom_text_repel(
      min.segment.length = 0, size = 1.8, fontface = "bold",
      max.overlaps = 300, box.padding = 0.1,
      point.padding = 0.05, force = 0.3
    ) +
    ggplot2::labs(
      title = title,
      x     = pc1_lab,
      y     = pc2_lab,
      color = "Batch ID"
    ) +
    theme_minimal_bold(base_size = base_size) +
    ggplot2::theme(legend.position = legend_position)
}


# ── make_pooled_unknown_reverse_scree_plot ────────────────────────────────────
#
# Builds a reverse scree plot for the pooled Unknown-only PCA.
#
# @param scree_df       A data frame of pooled scree data.
# @param title          Character. Plot title.
# @param base_size      Numeric. Base font size.
# @param legend_position Character. ggplot2 legend position.
# @return A ggplot object, or NULL if scree_df is empty.
#
#' @keywords internal
make_pooled_unknown_reverse_scree_plot <- function(
    scree_df, title, base_size = 10,
    legend_position = c("right", "bottom", "top", "left", "none")) {

  if (nrow(scree_df) == 0) return(NULL)
  legend_position <- match.arg(legend_position)

  ggplot2::ggplot(scree_df,
                  ggplot2::aes(x = Component, y = Variance_Pct)) +
    ggplot2::geom_col(fill = "steelblue", color = "black",
                      linewidth = 1, width = 0.7) +
    ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 1) +
    ggplot2::geom_text(
      ggplot2::aes(label = Variance_Pct),
      vjust = -0.25, size = 3, fontface = "bold"
    ) +
    ggplot2::labs(
      title = title,
      x     = "Principal Components (reverse: smallest \u2192 largest)",
      y     = "Variance Explained (%)"
    ) +
    theme_minimal_bold(base_size = base_size) +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = legend_position
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::coord_cartesian(clip = "off")
}


# ══ Exported pipeline functions ═══════════════════════════════════════════════

# ── export_pca_by_batch ───────────────────────────────────────────────────────
#
#' Export per-batch PCA and reverse scree PDFs
#'
#' Runs PCA separately for each batch using all replicate types, then
#' exports two portrait-layout PDFs — one for normalized ug/g and one for
#' unnormalized ug/g — each containing one page per batch with the PCA
#' biplot stacked above the reverse scree plot.
#'
#' Output filenames:
#' \itemize{
#'   \item \code{PCA_Scree_by_batch_normalized.pdf}
#'   \item \code{PCA_Scree_by_batch_unnormalized.pdf}
#' }
#'
#' @param all_batches_subtracted A data frame as returned by
#'   \code{calculate_concentrations()}.
#' @param out_dir   Character. Base results directory. PDFs are written to
#'   a \code{pca_output/} subfolder. Default \code{"."}.
#' @param na_method Character. NA handling: \code{"median"} (impute with
#'   column median) or \code{"drop"} (remove rows with any NA).
#'   Default \code{"median"}.
#' @param base_size Numeric. Base font size for plots. Default \code{10}.
#' @param file_norm   Character. Filename for normalized PDF.
#' @param file_unnorm Character. Filename for unnormalized PDF.
#' @param width  Numeric. PDF width in inches. Default \code{8.5}.
#' @param height Numeric. PDF height in inches. Default \code{11}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   writing PDFs.
#'
#' @details
#' \strong{Insufficient data:} Batches with fewer than 2 analyte features
#' or 2 replicate rows after filtering are skipped with a placeholder page
#' noting "insufficient data".
#'
#' \strong{What to modify:} \code{na_method} to switch between median
#' imputation and complete-case analysis; \code{out_dir} for output location.
#'
#' @examples
#' \dontrun{
#' export_pca_by_batch(
#'   all_data,
#'   out_dir   = "results/",
#'   na_method = "median"
#' )
#' }
#'
#' @export
export_pca_by_batch <- function(all_batches_subtracted,
                                out_dir     = ".",
                                na_method   = "median",
                                base_size   = 10,
                                file_norm   = "PCA_Scree_by_batch_normalized.pdf",
                                file_unnorm = "PCA_Scree_by_batch_unnormalized.pdf",
                                width       = 8.5,
                                height      = 11) {

  # Create pca_output subfolder automatically if it doesn't exist
  plot_dir <- file.path(out_dir, "pca_output")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    message("Created output subfolder: ", plot_dir)
  }

  # Placeholder for batches with insufficient data
  empty_plot <- function(msg) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg,
                        size = 4, fontface = "bold") +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", color = NA)
      )
  }

  # ── Normalized ─────────────────────────────────────────────────────────────
  message("Preparing per-batch PCA (normalized)...")
  data_norm <- prepare_pca_by_batch(all_batches_subtracted,
                                    "Norm_ug_Per_Gram",
                                    na_method = na_method)

  norm_path <- file.path(plot_dir, file_norm)
  grDevices::pdf(norm_path, width = width, height = height)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  all_bids_norm <- sort(unique(c(data_norm$scores$Batch_ID,
                                 data_norm$scree$Batch_ID)))

  for (bid in all_bids_norm) {
    p_pca <- make_pca_plot_for_batch(
      dplyr::filter(data_norm$scores, Batch_ID == bid),
      title     = paste0("PCA (Normalized ug/g) \u2014 Batch ", bid),
      base_size = base_size
    )
    p_scree <- make_reverse_scree_plot_for_batch(
      dplyr::filter(data_norm$scree, Batch_ID == bid),
      title     = paste0("Reverse Scree (Normalized ug/g) \u2014 Batch ", bid),
      base_size = base_size
    )
    if (is.null(p_pca))   p_pca   <- empty_plot("Normalized PCA: insufficient data")
    if (is.null(p_scree)) p_scree <- empty_plot("Normalized Scree: insufficient data")

    # Stack PCA above scree using patchwork; annotate with batch title
    combined <- (p_pca / p_scree) +
      patchwork::plot_annotation(
        title = paste0("Batch ", bid,
                       " \u2014 Normalized PCA & Reverse Scree"),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold",
                                             size = base_size + 3)
        )
      )
    print(combined)
  }
  grDevices::dev.off()
  on.exit()  # clear the on.exit handler after explicit dev.off
  message("Written: ", norm_path)

  # ── Unnormalized ───────────────────────────────────────────────────────────
  message("Preparing per-batch PCA (unnormalized)...")
  data_unnorm <- prepare_pca_by_batch(all_batches_subtracted,
                                      "Unnorm_ug_Per_Gram",
                                      na_method = na_method)

  unnorm_path <- file.path(plot_dir, file_unnorm)
  grDevices::pdf(unnorm_path, width = width, height = height)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  all_bids_unnorm <- sort(unique(c(data_unnorm$scores$Batch_ID,
                                   data_unnorm$scree$Batch_ID)))

  for (bid in all_bids_unnorm) {
    p_pca <- make_pca_plot_for_batch(
      dplyr::filter(data_unnorm$scores, Batch_ID == bid),
      title     = paste0("PCA (Unnormalized ug/g) \u2014 Batch ", bid),
      base_size = base_size
    )
    p_scree <- make_reverse_scree_plot_for_batch(
      dplyr::filter(data_unnorm$scree, Batch_ID == bid),
      title     = paste0("Reverse Scree (Unnormalized ug/g) \u2014 Batch ", bid),
      base_size = base_size
    )
    if (is.null(p_pca))   p_pca   <- empty_plot("Unnormalized PCA: insufficient data")
    if (is.null(p_scree)) p_scree <- empty_plot("Unnormalized Scree: insufficient data")

    combined <- (p_pca / p_scree) +
      patchwork::plot_annotation(
        title = paste0("Batch ", bid,
                       " \u2014 Unnormalized PCA & Reverse Scree"),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold",
                                             size = base_size + 3)
        )
      )
    print(combined)
  }
  grDevices::dev.off()
  on.exit()
  message("Written: ", unnorm_path)

  message("Per-batch PCA export complete.")
  invisible(NULL)
}


# ── export_pooled_unknown_pca ─────────────────────────────────────────────────
#
#' Export pooled Unknown-only PCA and reverse scree PDFs
#'
#' Runs a single PCA across all batches using only Unknown replicates,
#' then exports two portrait-layout PDFs — one normalized, one unnormalized
#' — each containing one page with the PCA biplot stacked above the reverse
#' scree plot. Points are colored by Batch_ID to reveal between-batch
#' clustering and systematic batch effects.
#'
#' Output filenames:
#' \itemize{
#'   \item \code{PCA_Scree_pooled_normalized_UNKNOWN.pdf}
#'   \item \code{PCA_Scree_pooled_unnormalized_UNKNOWN.pdf}
#' }
#'
#' @param all_batches_subtracted A data frame as returned by
#'   \code{calculate_concentrations()}.
#' @param out_dir   Character. Base results directory. PDFs written to
#'   \code{pca_output/} subfolder. Default \code{"."}.
#' @param na_method Character. \code{"median"} or \code{"drop"}.
#'   Default \code{"median"}.
#' @param base_size Numeric. Base font size. Default \code{10}.
#' @param file_norm   Character. Filename for normalized PDF.
#' @param file_unnorm Character. Filename for unnormalized PDF.
#' @param width  Numeric. PDF width in inches. Default \code{8.5}.
#' @param height Numeric. PDF height in inches. Default \code{11}.
#' @param use_combined_labels Logical. If TRUE, labels show
#'   "Replicate_Name (Batch_ID)". Default \code{TRUE}.
#' @param legend_position Character. Legend position for both plots.
#'   Default \code{"bottom"}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   writing PDFs.
#'
#' @details
#' \strong{Why Unknown-only:} Including standards and QC replicates in a
#' pooled PCA would dominate the variance decomposition since their
#' concentration profiles differ systematically from environmental samples.
#' Restricting to Unknowns gives a cleaner view of between-sample and
#' between-batch variability in the actual study samples.
#'
#' \strong{What to modify:} \code{use_combined_labels = FALSE} for cleaner
#' plots when batch IDs are already shown in the legend; \code{legend_position}
#' to manage legend crowding with many batches.
#'
#' @examples
#' \dontrun{
#' export_pooled_unknown_pca(
#'   all_data,
#'   out_dir             = "results/",
#'   na_method           = "median",
#'   use_combined_labels = TRUE,
#'   legend_position     = "bottom"
#' )
#' }
#'
#' @export
export_pooled_unknown_pca <- function(all_batches_subtracted,
                                      out_dir             = ".",
                                      na_method           = "median",
                                      base_size           = 10,
                                      file_norm   = "PCA_Scree_pooled_normalized_UNKNOWN.pdf",
                                      file_unnorm = "PCA_Scree_pooled_unnormalized_UNKNOWN.pdf",
                                      width               = 8.5,
                                      height              = 11,
                                      use_combined_labels = TRUE,
                                      legend_position     = c("bottom", "right",
                                                              "top", "left",
                                                              "none")) {

  legend_position <- match.arg(legend_position)

  # Create pca_output subfolder automatically if it doesn't exist
  plot_dir <- file.path(out_dir, "pca_output")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    message("Created output subfolder: ", plot_dir)
  }

  empty_plot <- function(msg) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg,
                        size = 4, fontface = "bold") +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", color = NA)
      )
  }

  # ── Normalized ─────────────────────────────────────────────────────────────
  message("Preparing pooled Unknown-only PCA (normalized)...")
  pooled_norm <- prepare_pooled_unknown_pca(all_batches_subtracted,
                                            "Norm_ug_Per_Gram",
                                            na_method = na_method)

  norm_path <- file.path(plot_dir, file_norm)
  grDevices::pdf(norm_path, width = width, height = height)
  {
    p_pca <- make_pooled_unknown_pca_plot(
      pooled_norm$scores,
      title               = "Pooled Unknown-only PCA (Normalized ug/g) \u2014 Color by Batch ID",
      base_size           = base_size,
      use_combined_labels = use_combined_labels,
      legend_position     = legend_position
    )
    p_scree <- make_pooled_unknown_reverse_scree_plot(
      pooled_norm$scree,
      title           = "Reverse Scree (Pooled Unknown-only, Normalized ug/g)",
      base_size       = base_size,
      legend_position = legend_position
    )
    if (is.null(p_pca))   p_pca   <- empty_plot("Pooled Unknown PCA (Normalized): insufficient data")
    if (is.null(p_scree)) p_scree <- empty_plot("Reverse Scree (Normalized): insufficient data")

    combined <- (p_pca / p_scree) +
      patchwork::plot_annotation(
        title = "Project-wide Unknown-only \u2014 Normalized PCA & Reverse Scree",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold",
                                             size = base_size + 3)
        )
      )
    print(combined)
  }
  grDevices::dev.off()
  message("Written: ", norm_path)

  # ── Unnormalized ───────────────────────────────────────────────────────────
  message("Preparing pooled Unknown-only PCA (unnormalized)...")
  pooled_unnorm <- prepare_pooled_unknown_pca(all_batches_subtracted,
                                              "Unnorm_ug_Per_Gram",
                                              na_method = na_method)

  unnorm_path <- file.path(plot_dir, file_unnorm)
  grDevices::pdf(unnorm_path, width = width, height = height)
  {
    p_pca <- make_pooled_unknown_pca_plot(
      pooled_unnorm$scores,
      title               = "Pooled Unknown-only PCA (Unnormalized ug/g) \u2014 Color by Batch ID",
      base_size           = base_size,
      use_combined_labels = use_combined_labels,
      legend_position     = legend_position
    )
    p_scree <- make_pooled_unknown_reverse_scree_plot(
      pooled_unnorm$scree,
      title           = "Reverse Scree (Pooled Unknown-only, Unnormalized ug/g)",
      base_size       = base_size,
      legend_position = legend_position
    )
    if (is.null(p_pca))   p_pca   <- empty_plot("Pooled Unknown PCA (Unnormalized): insufficient data")
    if (is.null(p_scree)) p_scree <- empty_plot("Reverse Scree (Unnormalized): insufficient data")

    combined <- (p_pca / p_scree) +
      patchwork::plot_annotation(
        title = "Project-wide Unknown-only \u2014 Unnormalized PCA & Reverse Scree",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold",
                                             size = base_size + 3)
        )
      )
    print(combined)
  }
  grDevices::dev.off()
  message("Written: ", unnorm_path)

  message("Pooled Unknown-only PCA export complete.")
  invisible(NULL)
}
