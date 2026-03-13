# R/plots_is.R
#
# Functions for visualizing internal standard (IS) signal intensity across
# batches and within individual batches.
#
# Internal standards are co-injected reference compounds added at a fixed
# known amount to every sample. Their signal intensity should remain
# consistent across all replicates and batches if instrument performance
# is stable. Significant drift or variability in IS signal indicates:
#   - Instrument sensitivity changes between injections
#   - Matrix suppression effects varying between samples
#   - Sample preparation inconsistencies
#
# Two complementary plot types are produced:
#   1. Cross-batch trend plot: mean IS intensity per batch with error bars,
#      useful for detecting systematic batch-to-batch instrument drift.
#   2. Per-batch replicate plot: individual IS intensities per replicate
#      within a single batch with CV annotation and ±1 SD bands, useful
#      for detecting within-batch outliers and injection order effects.
#
# Both plots use Area_BG_Sub (background-subtracted raw area) rather than
# Normalized_Area_BG_Sub because the IS itself is the normalization reference
# and cannot be self-normalized.
#
# Output subfolder: internal_standard/
#
# Contents:
#   Exported pipeline functions:
#     plot_is_across_batches()  — cross-batch trend plot
#     plot_is_per_batch()       — per-batch replicate-level plots


# ── plot_is_across_batches ────────────────────────────────────────────────────
#
#' Plot internal standard intensity trend across all batches
#'
#' Computes the mean and standard deviation of IS signal intensity within
#' each batch, then plots a connected point-and-error-bar chart showing how
#' IS intensity varies across the full run sequence. A consistent flat trend
#' indicates stable instrument performance; a declining or rising trend
#' suggests instrument drift that may require normalization or batch
#' correction.
#'
#' @param all_batches_subtracted A data frame as returned by
#'   \code{calculate_concentrations()}. Must contain columns \code{Batch_ID},
#'   \code{Molecule_Prefix}, and \code{Area_BG_Sub}.
#' @param out_dir Character. Base results directory. The plot will be written
#'   to a \code{internal_standard/} subfolder within this directory.
#'   Default \code{"."}.
#' @param output_file Character. Filename for the output PDF.
#'   Default \code{"IS_Averages_Per_Batch.pdf"}.
#' @param width  Numeric. PDF page width in inches. Default \code{8}.
#' @param height Numeric. PDF page height in inches. Default \code{6}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   writing a PDF file.
#'
#' @details
#' \strong{IS identification:} Rows are filtered to \code{Molecule_Prefix == "IS"}.
#' Ensure your IS molecule is named with the prefix "IS" in Skyline (e.g.,
#' "IS_1", "IS_deuterated") for it to be captured correctly.
#'
#' \strong{Error bars:} Represent ±1 SD of replicate IS intensities within
#' each batch. Wide error bars within a batch indicate within-batch injection
#' variability.
#'
#' \strong{What to modify:} The \code{out_dir}, \code{output_file},
#' \code{width}, and \code{height} arguments. To change the IS molecule
#' prefix from "IS", modify the filter inside this function.
#'
#' @examples
#' \dontrun{
#' plot_is_across_batches(all_data, out_dir = "results/")
#' }
#'
#' @export
plot_is_across_batches <- function(all_batches_subtracted,
                                   out_dir     = ".",
                                   output_file = "IS_Averages_Per_Batch.pdf",
                                   width       = 8,
                                   height      = 6) {

  # Create internal_standard subfolder if it doesn't exist
  plot_dir <- file.path(out_dir, "internal_standard")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    message("Created output subfolder: ", plot_dir)
  }

  # ── Filter to IS rows and summarize per batch ──────────────────────────────
  #
  # IS signal is summarized by batch rather than replicate to show the
  # overall batch-level trend. Per-replicate variation is shown separately
  # in plot_is_per_batch().
  is_data <- all_batches_subtracted %>%
    dplyr::filter(Molecule_Prefix == "IS") %>%
    dplyr::select(Batch_ID, Replicate_Name, Molecule_Prefix, Area_BG_Sub)

  if (nrow(is_data) == 0) {
    warning("No IS replicates found (Molecule_Prefix == 'IS'). ",
            "Skipping IS cross-batch plot. Check that your IS molecule ",
            "has prefix 'IS' in Skyline.")
    return(invisible(NULL))
  }

  # Compute mean and SD of IS intensity per batch in a single summarise call
  is_summary <- is_data %>%
    dplyr::group_by(Batch_ID) %>%
    dplyr::summarise(
      Avg_Area = mean(Area_BG_Sub, na.rm = TRUE),
      Std_Area = sd(Area_BG_Sub,   na.rm = TRUE),
      .groups  = "drop"
    )

  # ── Build plot ─────────────────────────────────────────────────────────────
  #
  # geom_line() connects batch means in the order they appear on the x-axis,
  # showing the directional trend across the run sequence.
  # geom_errorbar() shows ±1 SD of within-batch IS variability.
  p <- ggplot2::ggplot(
    is_summary,
    ggplot2::aes(x = Batch_ID, y = Avg_Area, fill = Batch_ID)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(group = 1),
      color     = "red",
      linewidth = 1.1
    ) +
    ggplot2::geom_point(
      stat  = "identity",
      color = "black",
      size  = 3
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = Avg_Area - Std_Area,
        ymax = Avg_Area + Std_Area
      ),
      width     = 0.1,
      linewidth = 1,
      color     = "black"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = "Internal Standard Intensity Across Batches",
      x     = "Batch ID",
      y     = "Average Intensity (Area, BG Subtracted)"
    ) +
    ggplot2::scale_fill_viridis_d(option = "B") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line.x      = ggplot2::element_line(linewidth = 1.1, color = "black"),
      axis.line.y      = ggplot2::element_line(linewidth = 1.1, color = "black"),
      panel.border     = ggplot2::element_blank(),
      axis.title       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold"),
      axis.text        = ggplot2::element_text(face = "bold", color = "black"),
      legend.position  = "none",
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 9)
    )

  # ── Write PDF ──────────────────────────────────────────────────────────────
  out_path <- file.path(plot_dir, output_file)
  grDevices::pdf(out_path, width = width, height = height)
  print(p)
  grDevices::dev.off()

  message("Written: ", out_path)
  invisible(NULL)
}


# ── plot_is_per_batch ─────────────────────────────────────────────────────────
#
#' Plot IS intensity for every replicate within each batch
#'
#' Produces a multi-page PDF with one page per batch. Each page shows
#' individual IS intensity values for every replicate in that batch,
#' annotated with the batch mean (dashed black line), ±1 SD bounds
#' (dashed red lines), color-coded ±1 SD bands (green = within bounds,
#' red = outside bounds), and a subtitle reporting mean, SD, and CV%.
#'
#' This plot is used to identify within-batch outlier injections and to
#' assess the consistency of IS recovery across the injection sequence.
#'
#' @param all_batches_subtracted A data frame as returned by
#'   \code{calculate_concentrations()}. Must contain columns \code{Batch_ID},
#'   \code{Replicate_Name}, \code{Molecule_Prefix}, and \code{Area_BG_Sub}.
#' @param out_dir Character. Base results directory. The plot will be written
#'   to a \code{internal_standard/} subfolder within this directory.
#'   Default \code{"."}.
#' @param output_file Character. Filename for the output PDF.
#'   Default \code{"IS_Batch_Plots.pdf"}.
#' @param width  Numeric. PDF page width in inches. Default \code{8}.
#' @param height Numeric. PDF page height in inches. Default \code{6}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   writing a multi-page PDF file.
#'
#' @details
#' \strong{CV calculation:} Coefficient of variation is computed as
#' \code{(SD / Mean) × 100}. If the mean is zero, CV is reported as 0
#' to avoid division by zero.
#'
#' \strong{SD of zero:} If a batch contains only one IS replicate, SD is
#' undefined (\code{NA}). This is replaced with 0 so the plot still renders;
#' the SD band will collapse to a line at the mean.
#'
#' \strong{Batch ordering:} Batches are plotted in the order they appear in
#' \code{unique(all_batches_subtracted$Batch_ID)} — typically the order they
#' were loaded, which reflects the run sequence.
#'
#' \strong{What to modify:} The \code{out_dir}, \code{output_file},
#' \code{width}, and \code{height} arguments.
#'
#' @examples
#' \dontrun{
#' plot_is_per_batch(all_data, out_dir = "results/")
#' }
#'
#' @export
plot_is_per_batch <- function(all_batches_subtracted,
                              out_dir     = ".",
                              output_file = "IS_Batch_Plots.pdf",
                              width       = 8,
                              height      = 6) {

  # Create internal_standard subfolder if it doesn't exist
  plot_dir <- file.path(out_dir, "internal_standard")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    message("Created output subfolder: ", plot_dir)
  }

  # ── Filter to IS rows ──────────────────────────────────────────────────────
  is_data <- all_batches_subtracted %>%
    dplyr::filter(Molecule_Prefix == "IS") %>%
    dplyr::select(Batch_ID, Replicate_Name, Molecule_Prefix, Area_BG_Sub)

  if (nrow(is_data) == 0) {
    warning("No IS replicates found (Molecule_Prefix == 'IS'). ",
            "Skipping IS per-batch plots.")
    return(invisible(NULL))
  }

  # Derive batch list from the actual IS data so we only loop over batches
  # that have IS replicates — avoids empty pages in the PDF
  batch_list <- unique(is_data$Batch_ID)

  out_path <- file.path(plot_dir, output_file)
  grDevices::pdf(out_path, width = width, height = height)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  for (batch in batch_list) {

    # Filter to this batch's IS replicates
    batch_data <- is_data %>% dplyr::filter(Batch_ID == batch)

    # ── Compute summary statistics ───────────────────────────────────────────
    avg_area <- mean(batch_data$Area_BG_Sub, na.rm = TRUE)
    sd_area  <- sd(batch_data$Area_BG_Sub,   na.rm = TRUE)

    # Replace NA SD (single replicate) with 0 so the plot still renders.
    # The SD bands will collapse to a line at the mean in this case.
    if (!is.finite(sd_area)) sd_area <- 0

    upper_bound <- avg_area + sd_area
    lower_bound <- avg_area - sd_area

    # CV: (SD / Mean) × 100. Guard against mean = 0.
    cv_percent <- if (avg_area != 0) (sd_area / avg_area) * 100 else 0

    # x-axis span for the annotated SD band rectangles
    n_replicates <- length(unique(batch_data$Replicate_Name))
    x_min        <- 0.5
    x_max        <- n_replicates + 0.5

    # y-axis upper limit: max observed value with 5% padding
    top_value <- max(batch_data$Area_BG_Sub, na.rm = TRUE)
    if (!is.finite(top_value) || top_value <= 0) top_value <- 1

    # ── Build plot ───────────────────────────────────────────────────────────
    #
    # Three annotated rectangles divide the plot into:
    #   Green band: within ±1 SD of the batch mean (acceptable range)
    #   Red bands: above upper bound and below lower bound (potential outliers)
    # Dashed horizontal lines mark the mean and ±1 SD bounds explicitly.
    p <- ggplot2::ggplot(
      batch_data,
      ggplot2::aes(x = Replicate_Name, y = Area_BG_Sub, fill = Replicate_Name)
    ) +
      # Green band: ±1 SD acceptance zone
      ggplot2::annotate("rect",
                        xmin = x_min, xmax = x_max,
                        ymin = lower_bound, ymax = upper_bound,
                        alpha = 0.2, fill = "green") +
      # Red band: below lower bound
      ggplot2::annotate("rect",
                        xmin = x_min, xmax = x_max,
                        ymin = 0, ymax = lower_bound,
                        alpha = 0.2, fill = "red") +
      # Red band: above upper bound
      ggplot2::annotate("rect",
                        xmin = x_min, xmax = x_max,
                        ymin = upper_bound, ymax = top_value,
                        alpha = 0.2, fill = "red") +
      ggplot2::geom_point(size = 2) +
      # Mean line
      ggplot2::geom_hline(
        yintercept = avg_area, color = "black",
        linetype = "dashed", linewidth = 1.1
      ) +
      # Upper ±1 SD line
      ggplot2::geom_hline(
        yintercept = upper_bound, color = "red",
        linetype = "dashed", linewidth = 1.1
      ) +
      # Lower ±1 SD line
      ggplot2::geom_hline(
        yintercept = lower_bound, color = "red",
        linetype = "dashed", linewidth = 1.1
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(
        title    = "Internal Standard Intensity Per Replicate",
        subtitle = paste0(
          "Batch: ", batch,
          "  |  Mean = ", round(avg_area, 2),
          "  |  SD = ",   round(sd_area, 2),
          "  |  CV = ",   round(cv_percent, 2), "%"
        ),
        x = paste("Batch ID:", batch),
        y = "Intensity (Area, BG Subtracted)"
      ) +
      ggplot2::scale_fill_viridis_d(option = "B") +
      ggplot2::scale_y_continuous(
        limits = c(0, top_value * 1.05),
        expand = c(0, 0)
      ) +
      ggplot2::scale_x_discrete(
        expand = ggplot2::expansion(mult = c(0.05, 0.05))
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
        legend.position  = "none",
        axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
        plot.margin      = ggplot2::margin(10, 20, 10, 20)
      ) +
      ggplot2::coord_cartesian(clip = "off")

    print(p)
  }

  message("Written: ", out_path,
          " (", length(batch_list), " page(s))")
  invisible(NULL)
}
