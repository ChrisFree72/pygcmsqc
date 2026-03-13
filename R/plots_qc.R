# R/plots_qc.R
#
# Functions for visualizing Quality Control (QC) replicate signal intensity
# across batches and within individual batches.
#
# QC replicates are prepared reference materials injected alongside unknown
# samples at regular intervals throughout a batch sequence. Their signal
# should remain consistent across injections and batches if the analytical
# method is under control. Significant drift or variability in QC signal
# indicates method performance issues that may affect sample quantification.
#
# Four complementary plot types are produced:
#   1. Cross-batch trend (per analyte): mean QC intensity per batch with
#      error bars — detects systematic batch-to-batch drift per polymer.
#   2. Replicate-level scatter (per analyte): individual QC intensities
#      across all replicates and batches with ±SD bands — detects outliers
#      and between-batch variability at the replicate level.
#   3. All-batches bar chart: mean QC intensity per analyte pooled across
#      all batches — gives a single-page overview of relative polymer
#      response magnitudes.
#   4. Per-batch bar charts: mean QC intensity per analyte within each
#      batch — allows batch-specific QC review.
#
# All plots use Area_BG_Sub (background-subtracted raw signal).
#
# Output subfolder: qc_plots/
#
# Contents:
#   Internal helpers (not exported):
#     .qc_bar_theme()    — shared ggplot2 theme (defined in utils.R)
#     .qc_cleanup()      — shared filter/coerce step (defined in utils.R)
#
#   Exported pipeline functions:
#     plot_qc_across_batches()    — cross-batch trend per analyte
#     plot_qc_replicate_level()   — replicate-level scatter per analyte
#     plot_qc_bar_all_batches()   — all-batches summary bar chart
#     plot_qc_bar_by_batch()      — per-batch bar charts


# ── plot_qc_across_batches ────────────────────────────────────────────────────
#
#' Plot mean QC intensity trend across batches, one page per analyte
#'
#' For each analyte (Molecule_Prefix), computes the mean and SD of QC
#' replicate intensity within each batch and plots a connected
#' point-and-error-bar chart showing how QC intensity varies across the
#' batch sequence. Produces a multi-page PDF with one page per analyte.
#'
#' @param all_batches_subtracted A data frame as returned by
#'   \code{calculate_concentrations()}. Must contain \code{Batch_ID},
#'   \code{Replicate_Name}, \code{Replicate_Type}, \code{Molecule_Prefix},
#'   and \code{Area_BG_Sub}.
#' @param out_dir Character. Base results directory. Output goes to the
#'   \code{qc_plots/} subfolder. Default \code{"."}.
#' @param output_file Character. PDF filename.
#'   Default \code{"QC_Averages_Per_Batch.pdf"}.
#' @param width  Numeric. PDF width in inches. Default \code{8}.
#' @param height Numeric. PDF height in inches. Default \code{6}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   writing a PDF.
#'
#' @details
#' \strong{What to modify:} \code{out_dir}, \code{output_file},
#' \code{width}, \code{height}. Plot styling is controlled by
#' \code{.qc_bar_theme()} in \code{utils.R}.
#'
#' @examples
#' \dontrun{
#' plot_qc_across_batches(all_data, out_dir = "results/")
#' }
#'
#' @export
plot_qc_across_batches <- function(all_batches_subtracted,
                                   out_dir     = ".",
                                   output_file = "QC_Averages_Per_Batch.pdf",
                                   width       = 8,
                                   height      = 6) {

  # Create qc_plots subfolder automatically if it doesn't exist
  plot_dir <- file.path(out_dir, "qc_plots")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    message("Created output subfolder: ", plot_dir)
  }

  # ── Prepare data ───────────────────────────────────────────────────────────
  #
  # .qc_cleanup() (utils.R): filters to QC replicates, coerces column types,
  # and drops rows with missing or non-finite intensity values.
  data <- .qc_cleanup(all_batches_subtracted)

  if (nrow(data) == 0) {
    warning("No QC replicates found after filtering. ",
            "Skipping QC cross-batch plot.")
    return(invisible(NULL))
  }

  # Summarize QC intensity per Batch_ID × Molecule_Prefix.
  # Deduplication with distinct() guards against any residual duplicate rows
  # that could inflate means or produce multiple lines per analyte.
  qc_summary <- data %>%
    dplyr::group_by(Batch_ID, Molecule_Prefix) %>%
    dplyr::summarise(
      Avg_Area = mean(Area_BG_Sub, na.rm = TRUE),
      Std_Area = sd(Area_BG_Sub,   na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    dplyr::distinct(Batch_ID, Molecule_Prefix, .keep_all = TRUE) %>%
    dplyr::mutate(
      # Replace NA/NaN (from single-replicate groups) with 0
      dplyr::across(where(is.numeric), ~ ifelse(is.na(.) | is.nan(.), 0, .))
    )

  # ── Build one plot per analyte ─────────────────────────────────────────────
  out_path <- file.path(plot_dir, output_file)
  grDevices::pdf(out_path, width = width, height = height)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  for (mp in unique(qc_summary$Molecule_Prefix)) {

    dsub <- qc_summary %>%
      dplyr::filter(Molecule_Prefix == mp) %>%
      dplyr::arrange(Batch_ID)  # connect points in batch order

    p <- ggplot2::ggplot(dsub, ggplot2::aes(x = Batch_ID, y = Avg_Area)) +
      ggplot2::geom_line(
        ggplot2::aes(group = 1),
        color = "red", linewidth = 1.1
      ) +
      ggplot2::geom_point(shape = 21, fill = "black", size = 3) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = Avg_Area - Std_Area, ymax = Avg_Area + Std_Area),
        width = 0.1, linewidth = 1, color = "black"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(
        title = paste0("QC Analyte Intensity Across Batches \u2014 ", mp),
        x     = "Batch ID",
        y     = "Average Intensity (Area, BG Subtracted)"
      ) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line.x      = ggplot2::element_line(linewidth = 1.1, color = "black"),
        axis.line.y      = ggplot2::element_line(linewidth = 1.1, color = "black"),
        panel.border     = ggplot2::element_blank(),
        axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
        axis.title       = ggplot2::element_text(face = "bold"),
        plot.title       = ggplot2::element_text(face = "bold"),
        axis.text        = ggplot2::element_text(face = "bold", color = "black")
      )

    print(p)
  }

  message("Written: ", out_path,
          " (", dplyr::n_distinct(qc_summary$Molecule_Prefix), " page(s))")
  invisible(NULL)
}


# ── plot_qc_replicate_level ───────────────────────────────────────────────────
#
#' Plot QC intensity at the replicate level across all batches, per analyte
#'
#' For each analyte, plots every individual QC replicate intensity value
#' across all batches. Points are colored by Batch_ID. A global mean and
#' ±SD band (computed across all replicates of all batches for that analyte)
#' is overlaid to identify replicates falling outside the expected range.
#'
#' Produces a multi-page PDF with one page per analyte.
#'
#' @param all_batches_subtracted A data frame as returned by
#'   \code{calculate_concentrations()}.
#' @param out_dir Character. Base results directory. Output goes to the
#'   \code{qc_plots/} subfolder. Default \code{"."}.
#' @param output_file Character. PDF filename.
#'   Default \code{"QC_Replicate_Plots_Per_Analyte.pdf"}.
#' @param sd_multiplier Numeric. Number of SDs to use for the acceptance
#'   band. Default \code{1} (±1 SD). Increase to \code{2} or \code{3}
#'   for a wider acceptance window.
#' @param batch_order Character vector. Optional custom ordering of Batch_IDs
#'   for the legend color assignment. Default \code{NULL} uses data order.
#' @param width  Numeric. PDF width in inches. Default \code{9}.
#' @param height Numeric. PDF height in inches. Default \code{6}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   writing a PDF.
#'
#' @details
#' \strong{SD band:} The mean and SD are computed across ALL QC replicates
#' for a given analyte, pooling across batches. This gives a project-wide
#' reference range rather than a batch-specific one, making it easier to
#' identify batches that are systematically high or low.
#'
#' \strong{What to modify:} \code{sd_multiplier} to widen or narrow the
#' acceptance band; \code{batch_order} to control legend color assignment.
#'
#' @examples
#' \dontrun{
#' plot_qc_replicate_level(
#'   all_data,
#'   out_dir       = "results/",
#'   sd_multiplier = 2
#' )
#' }
#'
#' @export
plot_qc_replicate_level <- function(all_batches_subtracted,
                                    out_dir       = ".",
                                    output_file   = "QC_Replicate_Plots_Per_Analyte.pdf",
                                    sd_multiplier = 1,
                                    batch_order   = NULL,
                                    width         = 9,
                                    height        = 6) {

  # Create qc_plots subfolder automatically if it doesn't exist
  plot_dir <- file.path(out_dir, "qc_plots")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    message("Created output subfolder: ", plot_dir)
  }

  # .qc_cleanup() filters to QC, coerces types, drops non-finite rows
  data <- .qc_cleanup(all_batches_subtracted)

  if (nrow(data) == 0) {
    warning("No QC replicates found after filtering. ",
            "Skipping QC replicate-level plot.")
    return(invisible(NULL))
  }

  if (is.null(batch_order)) batch_order <- unique(data$Batch_ID)

  # Remove existing file to avoid device-lock issues on Windows
  out_path <- file.path(plot_dir, output_file)
  if (file.exists(out_path)) {
    try(suppressWarnings(file.remove(out_path)), silent = TRUE)
  }

  grDevices::pdf(out_path, width = width, height = height)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  for (mp in unique(data$Molecule_Prefix)) {

    dsub <- data %>% dplyr::filter(Molecule_Prefix == mp)
    if (nrow(dsub) == 0) next

    # ── Global mean and SD across all replicates for this analyte ─────────────
    #
    # Pooled across all batches so the band represents the project-wide
    # reference range. Replicates falling outside are flagged visually.
    mu    <- mean(dsub$Area_BG_Sub, na.rm = TRUE)
    sigma <- sd(dsub$Area_BG_Sub,   na.rm = TRUE)
    mu    <- ifelse(is.finite(mu),    mu,    0)
    sigma <- ifelse(is.finite(sigma), sigma, 0)

    lower <- pmax(0, mu - sd_multiplier * sigma)
    upper <- mu + sd_multiplier * sigma

    # Set x-axis factor levels to preserve replicate appearance order
    repl_levels <- unique(dsub$Replicate_Name)
    dsub <- dsub %>%
      dplyr::mutate(
        Replicate_Name = factor(Replicate_Name, levels = repl_levels),
        Batch_ID       = factor(Batch_ID, levels = batch_order)
      ) %>%
      dplyr::arrange(Replicate_Name)

    # y-axis upper limit with 10% headroom
    top_y <- max(c(dsub$Area_BG_Sub, upper), na.rm = TRUE)
    top_y <- if (is.finite(top_y) && top_y > 0) top_y * 1.10 else 1
    n_x   <- length(levels(dsub$Replicate_Name))

    p <- ggplot2::ggplot(
      dsub,
      ggplot2::aes(x = Replicate_Name, y = Area_BG_Sub)
    ) +
      # Green acceptance band
      ggplot2::annotate("rect",
                        xmin = 0.5, xmax = n_x + 0.5,
                        ymin = lower, ymax = upper,
                        alpha = 0.20, fill = "green") +
      # Red band below lower bound
      ggplot2::annotate("rect",
                        xmin = 0.5, xmax = n_x + 0.5,
                        ymin = 0, ymax = lower,
                        alpha = 0.20, fill = "red") +
      # Red band above upper bound
      ggplot2::annotate("rect",
                        xmin = 0.5, xmax = n_x + 0.5,
                        ymin = upper, ymax = top_y,
                        alpha = 0.20, fill = "red") +
      # Individual replicate points colored by Batch_ID
      ggplot2::geom_point(
        ggplot2::aes(fill = Batch_ID),
        shape = 21, color = "black", size = 3
      ) +
      # Mean line
      ggplot2::geom_hline(
        yintercept = mu, color = "black",
        linetype = "dashed", linewidth = 1.1
      ) +
      # Lower bound line
      ggplot2::geom_hline(
        yintercept = lower, color = "red",
        linetype = "dashed", linewidth = 1.1
      ) +
      # Upper bound line
      ggplot2::geom_hline(
        yintercept = upper, color = "red",
        linetype = "dashed", linewidth = 1.1
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(
        title    = "QC Intensity (Replicate Level)",
        subtitle = paste0(
          "Analyte: ", mp,
          "  |  Mean = ", round(mu, 2),
          "  |  SD = ",   round(sigma, 2),
          if (sd_multiplier != 1)
            paste0("  |  \u00b1", sd_multiplier, " SD") else ""
        ),
        x    = "Replicate Name",
        y    = "Intensity (Area, BG Subtracted)",
        fill = "Batch ID"
      ) +
      ggplot2::scale_fill_viridis_d(option = "B") +
      ggplot2::scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, top_y)
      ) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line.x      = ggplot2::element_line(linewidth = 1.1, color = "black"),
        axis.line.y      = ggplot2::element_line(linewidth = 1.1, color = "black"),
        panel.border     = ggplot2::element_blank(),
        axis.title       = ggplot2::element_text(face = "bold"),
        plot.title       = ggplot2::element_text(face = "bold"),
        plot.subtitle    = ggplot2::element_text(face = "bold"),
        axis.text        = ggplot2::element_text(face = "bold", color = "black"),
        axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
        legend.title     = ggplot2::element_text(face = "bold"),
        legend.text      = ggplot2::element_text(face = "bold"),
        plot.margin      = ggplot2::margin(10, 20, 10, 20)
      )

    print(p)
  }

  message("Written: ", out_path,
          " (", dplyr::n_distinct(data$Molecule_Prefix), " page(s))")
  invisible(NULL)
}


# ── plot_qc_bar_all_batches ───────────────────────────────────────────────────
#
#' Bar chart of mean QC intensity per analyte pooled across all batches
#'
#' Produces a single-page bar chart showing the mean QC replicate intensity
#' for each analyte (Molecule_Prefix), pooled across all batches. Error bars
#' represent ±1 SD of individual replicate intensities. Useful for a rapid
#' overview of relative polymer signal magnitudes across the full project.
#'
#' @param all_batches_subtracted A data frame as returned by
#'   \code{calculate_concentrations()}.
#' @param out_dir Character. Base results directory. Output goes to the
#'   \code{qc_plots/} subfolder. Default \code{"."}.
#' @param output_file Character. PDF filename.
#'   Default \code{"QC_Bar_AllBatches.pdf"}.
#' @param sort_bars Character. Bar ordering: \code{"none"} (alphabetical),
#'   \code{"asc"} (ascending mean), or \code{"desc"} (descending mean).
#'   Default \code{"desc"}.
#' @param fill_palette Named character vector. Optional custom color palette
#'   mapping Molecule_Prefix values to hex colors. Default \code{NULL} uses
#'   viridis D.
#' @param width  Numeric. PDF width in inches. Default \code{11}.
#' @param height Numeric. PDF height in inches. Default \code{7}.
#'
#' @return Invisibly returns the ggplot object. Called primarily for its
#'   side effect of writing a PDF.
#'
#' @details
#' \strong{What to modify:} \code{sort_bars} to control bar ordering;
#' \code{fill_palette} to use custom polymer colors for publication figures.
#'
#' @examples
#' \dontrun{
#' plot_qc_bar_all_batches(
#'   all_data,
#'   out_dir   = "results/",
#'   sort_bars = "desc"
#' )
#' }
#'
#' @export
plot_qc_bar_all_batches <- function(all_batches_subtracted,
                                    out_dir      = ".",
                                    output_file  = "QC_Bar_AllBatches.pdf",
                                    sort_bars    = c("desc", "asc", "none"),
                                    fill_palette = NULL,
                                    width        = 11,
                                    height       = 7) {

  sort_bars <- match.arg(sort_bars)

  # Create qc_plots subfolder automatically if it doesn't exist
  plot_dir <- file.path(out_dir, "qc_plots")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    message("Created output subfolder: ", plot_dir)
  }

  data          <- .qc_cleanup(all_batches_subtracted)
  prefix_levels <- sort(unique(data$Molecule_Prefix))

  # Summarize across ALL batches (replicate-level stats)
  dd <- data %>%
    dplyr::group_by(Molecule_Prefix) %>%
    dplyr::summarise(
      mean_val = mean(Area_BG_Sub, na.rm = TRUE),
      sd_val   = sd(Area_BG_Sub,   na.rm = TRUE),
      n        = dplyr::n(),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      sd_val = ifelse(is.finite(sd_val), sd_val, 0),
      lower  = pmax(0, mean_val - sd_val),
      upper  = mean_val + sd_val
    )

  if (sort_bars == "asc")  dd <- dd %>% dplyr::arrange(mean_val)
  if (sort_bars == "desc") dd <- dd %>% dplyr::arrange(dplyr::desc(mean_val))

  # Lock in factor levels after sorting so ggplot respects the order
  dd$Molecule_Prefix <- factor(dd$Molecule_Prefix, levels = prefix_levels)

  y_top <- max(dd$upper, na.rm = TRUE)
  y_top <- if (is.finite(y_top) && y_top > 0) y_top * 1.05 else 1

  g <- ggplot2::ggplot(
    dd,
    ggplot2::aes(x = Molecule_Prefix, y = mean_val, fill = Molecule_Prefix)
  ) +
    ggplot2::geom_col(width = 1, color = "black", linewidth = 0.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      width = 0.5, linewidth = 0.5
    ) +
    ggplot2::labs(
      title    = "QC Mean Intensity by Analyte (All Batches)",
      subtitle = stringr::str_wrap(
        "Bar = mean of replicate intensities across all batches; error bars = \u00b11 SD",
        width = 100
      ),
      x    = "Molecule Prefix",
      y    = "Average Intensity (Area, BG Subtracted)",
      fill = "Molecule Prefix"
    ) +
    .qc_bar_theme() +
    ggplot2::scale_y_continuous(
      limits = c(0, y_top),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::scale_x_discrete(
      expand = ggplot2::expansion(mult = c(0.05, 0.05))
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::coord_cartesian(clip = "off")

  if (!is.null(fill_palette)) {
    g <- g + ggplot2::scale_fill_manual(
      values = fill_palette, limits = prefix_levels, drop = FALSE
    )
  } else {
    g <- g + ggplot2::scale_fill_viridis_d(
      option = "D", begin = 0.20, end = 0.90,
      limits = prefix_levels, drop = FALSE
    )
  }

  out_path <- file.path(plot_dir, output_file)
  ggplot2::ggsave(filename = out_path, plot = g, width = width, height = height)
  message("Written: ", out_path)
  invisible(g)
}


# ── plot_qc_bar_by_batch ──────────────────────────────────────────────────────
#
#' Bar charts of mean QC intensity per analyte, one page per batch
#'
#' Produces a multi-page PDF with one bar chart per batch. Each chart shows
#' the mean QC replicate intensity per analyte within that batch, with ±1 SD
#' error bars. Analyte colors are consistent across all batch pages to allow
#' easy visual comparison.
#'
#' @param all_batches_subtracted A data frame as returned by
#'   \code{calculate_concentrations()}.
#' @param out_dir Character. Base results directory. Output goes to the
#'   \code{qc_plots/} subfolder. Default \code{"."}.
#' @param output_file Character. PDF filename.
#'   Default \code{"QC_Bar_ByBatch.pdf"}.
#' @param sort_bars Character. Bar ordering within each batch page:
#'   \code{"none"}, \code{"asc"}, or \code{"desc"}. Default \code{"desc"}.
#' @param fill_palette Named character vector. Optional custom color palette.
#'   Default \code{NULL} uses viridis D.
#' @param width  Numeric. PDF width in inches. Default \code{11}.
#' @param height Numeric. PDF height in inches. Default \code{7}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   writing a multi-page PDF.
#'
#' @details
#' \strong{Consistent colors:} Analyte fill colors are assigned from the
#' global set of all Molecule_Prefix values across the full dataset, not
#' just those present in each batch. This ensures colors are consistent
#' across batch pages even if some analytes are absent in certain batches.
#'
#' \strong{What to modify:} \code{sort_bars} and \code{fill_palette}.
#'
#' @examples
#' \dontrun{
#' plot_qc_bar_by_batch(
#'   all_data,
#'   out_dir   = "results/",
#'   sort_bars = "desc"
#' )
#' }
#'
#' @export
plot_qc_bar_by_batch <- function(all_batches_subtracted,
                                 out_dir      = ".",
                                 output_file  = "QC_Bar_ByBatch.pdf",
                                 sort_bars    = c("desc", "asc", "none"),
                                 fill_palette = NULL,
                                 width        = 11,
                                 height       = 7) {

  sort_bars <- match.arg(sort_bars)

  # Create qc_plots subfolder automatically if it doesn't exist
  plot_dir <- file.path(out_dir, "qc_plots")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    message("Created output subfolder: ", plot_dir)
  }

  data          <- .qc_cleanup(all_batches_subtracted)
  prefix_levels <- sort(unique(data$Molecule_Prefix))

  # Summarize within each batch
  dd <- data %>%
    dplyr::group_by(Batch_ID, Molecule_Prefix) %>%
    dplyr::summarise(
      mean_val = mean(Area_BG_Sub, na.rm = TRUE),
      sd_val   = sd(Area_BG_Sub,   na.rm = TRUE),
      n        = dplyr::n(),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(
      sd_val = ifelse(is.finite(sd_val), sd_val, 0),
      lower  = pmax(0, mean_val - sd_val),
      upper  = mean_val + sd_val
    )

  batches <- unique(dd$Batch_ID)
  if (length(batches) == 0) {
    warning("No QC batches found. Skipping per-batch QC bar plot.")
    return(invisible(NULL))
  }

  out_path <- file.path(plot_dir, output_file)
  grDevices::pdf(out_path, width = width, height = height)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  for (b in batches) {
    sub <- dd %>% dplyr::filter(Batch_ID == b)

    if (sort_bars == "asc")  sub <- sub %>% dplyr::arrange(mean_val)
    if (sort_bars == "desc") sub <- sub %>% dplyr::arrange(dplyr::desc(mean_val))

    # Use global prefix_levels so colors are consistent across batch pages
    sub$Molecule_Prefix <- factor(sub$Molecule_Prefix, levels = prefix_levels)

    y_top <- max(sub$upper, na.rm = TRUE)
    y_top <- if (is.finite(y_top) && y_top > 0) y_top * 1.05 else 1

    g <- ggplot2::ggplot(
      sub,
      ggplot2::aes(x = Molecule_Prefix, y = mean_val, fill = Molecule_Prefix)
    ) +
      ggplot2::geom_col(width = 1, color = "black", linewidth = 0.5) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = lower, ymax = upper),
        width = 0.5, linewidth = 0.5
      ) +
      ggplot2::labs(
        title    = paste0("QC Mean Intensity by Analyte \u2014 Batch ", b),
        subtitle = stringr::str_wrap(
          "Bar = mean of replicate intensities within this batch; error bars = \u00b11 SD",
          width = 100
        ),
        x    = "Molecule Prefix",
        y    = "Average Intensity (Area, BG Subtracted)",
        fill = "Molecule Prefix"
      ) +
      .qc_bar_theme() +
      ggplot2::scale_y_continuous(
        limits = c(0, y_top),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::scale_x_discrete(
        expand = ggplot2::expansion(mult = c(0.05, 0.05))
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
      ggplot2::coord_cartesian(clip = "off")

    if (!is.null(fill_palette)) {
      g <- g + ggplot2::scale_fill_manual(
        values = fill_palette, limits = prefix_levels, drop = FALSE
      )
    } else {
      g <- g + ggplot2::scale_fill_viridis_d(
        option = "D", begin = 0.20, end = 0.90,
        limits = prefix_levels, drop = FALSE
      )
    }

    print(g)
  }

  message("Written: ", out_path,
          " (", length(batches), " page(s))")
  invisible(NULL)
}
