# R/plots_calibration.R
#
# Functions for generating calibration curve diagnostic plots for each
# analyte in each batch, exported as per-batch PDFs.
#
# Each plot shows:
#   - Standard replicates plotted at their Scaled_Concentration (x-axis)
#   - Non-standard replicates (Unknowns, QC, Blanks) plotted at their
#     back-calculated Calc_Concentration on the same x-axis, allowing
#     visual assessment of where samples fall on the calibration curve
#   - The fitted weighted 1/x regression line drawn across the standard range
#   - Vertical LOD and LOQ threshold lines (from calculate_lod_loq())
#   - R², weighting type, and removed calibrator annotations in the subtitle
#
# Two PDFs are produced per batch:
#   - Normalized: uses Normalized_Area_BG_Sub (IS-normalized signal)
#   - Unnormalized: uses Area_BG_Sub (raw signal)
#
# Contents:
#   Internal helpers (not exported):
#     prepare_plot_data()   — assigns x and y plotting columns by signal type
#     make_curve_plot()     — builds a single ggplot calibration curve
#
#   Exported pipeline function:
#     export_calibration_pdfs() — loops over batches and analytes, writes PDFs


# ══ Internal helpers ══════════════════════════════════════════════════════════

# ── prepare_plot_data ─────────────────────────────────────────────────────────
#
# Assigns x_to_plot and y_to_plot columns for a calibration curve plot.
#
# For standards, the x-axis is Scaled_Concentration (the known calibrator
# level scaled by Calibration_Factor) — this is what the model was fit on.
# For all other replicates (Unknowns, QC, Blanks), the x-axis is the
# back-calculated concentration from the fitted model, placing them on the
# same axis as the standards for visual comparison.
#
# This dual x-axis assignment is purely cosmetic — it allows the analyst to
# see visually where samples fall relative to the calibration range without
# changing any underlying values.
#
# @param df_group A data frame for a single Batch_ID × Molecule_Prefix group,
#                 as produced by filtering All_Batches_Subtracted.
# @param y_col    Character. Which signal column to plot on the y-axis.
#                 One of "Normalized_Area_BG_Sub" or "Area_BG_Sub".
# @return The input data frame with two new columns: x_to_plot and y_to_plot.
#
#' @keywords internal
prepare_plot_data <- function(df_group,
                              y_col = c("Normalized_Area_BG_Sub",
                                        "Area_BG_Sub")) {
  y_col <- match.arg(y_col)

  if (y_col == "Normalized_Area_BG_Sub") {
    df_group %>%
      dplyr::mutate(
        # Standards: known scaled concentration (model x-axis)
        # Others: back-calculated from normalized signal
        x_to_plot = dplyr::if_else(
          Replicate_Type == "Standard",
          Scaled_Concentration,
          Norm_Calc_Concentration
        ),
        y_to_plot = Normalized_Area_BG_Sub
      )
  } else {
    df_group %>%
      dplyr::mutate(
        x_to_plot = dplyr::if_else(
          Replicate_Type == "Standard",
          Scaled_Concentration,
          Unnorm_Calc_Concentration
        ),
        y_to_plot = Area_BG_Sub
      )
  }
}


# ── make_curve_plot ───────────────────────────────────────────────────────────
#
# Builds a single calibration curve ggplot for one analyte in one batch.
#
# Plot elements:
#   - Points: colored by Replicate_Type using a fixed color map
#   - Regression line: drawn from min to max standard Scaled_Concentration
#     using the fitted slope and intercept; annotated "No fit available"
#     if model parameters are NA
#   - LOD/LOQ: vertical dashed/solid red lines at the back-calculated
#     threshold concentrations; sourced from the lod_loq table
#   - Subtitle: weighting type, R², and removed calibrator names
#
# @param df_group A data frame for a single Batch_ID × Molecule_Prefix group.
# @param lod_loq  The full LOD/LOQ table from calculate_lod_loq(), used to
#                 look up threshold values for this batch and analyte.
# @param y_col    Character. Signal column to plot. One of
#                 "Normalized_Area_BG_Sub" or "Area_BG_Sub".
# @return A ggplot object.
#
#' @keywords internal
make_curve_plot <- function(df_group,
                            lod_loq,
                            y_col = c("Normalized_Area_BG_Sub",
                                      "Area_BG_Sub")) {
  y_col <- match.arg(y_col)

  # Pull scalar metadata from the first row of the group
  # (all rows share the same Batch_ID, Molecule_Prefix, and model params)
  pdat      <- df_group %>% dplyr::slice(1)
  batch     <- pdat$Batch_ID
  analyte   <- pdat$Molecule_Prefix
  weighting <- if (!is.null(pdat$Weighting_Type) &&
                   !is.na(pdat$Weighting_Type)) pdat$Weighting_Type else "N/A"

  # Select signal-type-specific model parameters and labels
  if (y_col == "Normalized_Area_BG_Sub") {
    slope     <- pdat$Norm_Slope
    intercept <- pdat$Norm_Y_Int
    r2        <- pdat$New_R2_Norm
    title     <- paste0("Calibration: ", analyte,
                        " (Batch ", batch, ") \u2014 Normalized")
    ylab      <- "Normalized Area (BG Subtracted)"
  } else {
    slope     <- pdat$Unnorm_Slope
    intercept <- pdat$Unnorm_Y_Int
    r2        <- pdat$New_R2_Unnorm
    title     <- paste0("Calibration: ", analyte,
                        " (Batch ", batch, ") \u2014 Unnormalized")
    ylab      <- "Area (BG Subtracted)"
  }

  # Build subtitle: weighting, R², and any removed calibrators
  removed_txt <- if (!is.null(pdat$Removed_Calibrators) &&
                     !is.na(pdat$Removed_Calibrators) &&
                     pdat$Removed_Calibrators != "None") {
    paste0(" | Removed Calibrators: ", pdat$Removed_Calibrators)
  } else ""

  subtitle_text <- paste0(
    "Weighting: ", weighting,
    if (!is.na(r2)) paste0("   |   R\u00b2 = ", formatC(r2, format = "f",
                                                        digits = 4)) else "",
    removed_txt
  )

  # Prepare plotting data: assigns x_to_plot and y_to_plot
  dd <- prepare_plot_data(df_group, y_col) %>%
    dplyr::filter(!is.na(x_to_plot), !is.na(y_to_plot))

  # Remove the excluded calibrator replicates from the plot so they don't
  # appear as points even though their rows exist in the data
  if (!is.null(pdat$Removed_Calibrators) &&
      !is.na(pdat$Removed_Calibrators) &&
      pdat$Removed_Calibrators != "None") {
    removed_names <- stringr::str_split(
      pdat$Removed_Calibrators, ","
    )[[1]] %>% trimws()
    dd <- dd %>% dplyr::filter(!(Replicate_Name %in% removed_names))
  }

  # Fixed color map and legend order for Replicate_Type.
  # limits = replicate_categories enforces a consistent legend order
  # across all plots regardless of which types are present in this group.
  replicate_categories <- c("Unknown", "QC", "Blank", "Cup_Blank", "Standard")
  color_map <- c(
    "Unknown"   = "maroon",
    "QC"        = "darkgreen",
    "Standard"  = "darkgray",
    "Blank"     = "blue",
    "Cup_Blank" = "orange"
  )

  # ── Base plot ──────────────────────────────────────────────────────────────
  g <- ggplot2::ggplot(
    dd,
    ggplot2::aes(x = x_to_plot, y = y_to_plot, color = Replicate_Type)
  ) +
    ggplot2::geom_point(alpha = 0.9, size = 2) +
    ggplot2::scale_color_manual(
      values = color_map,
      limits = replicate_categories,
      name   = "Replicate Type"
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = stringr::str_wrap(subtitle_text, width = 80),
      x        = "Concentration (Standards: Scaled; Samples: Back-Calculated)",
      y        = ylab
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      axis.title       = ggplot2::element_text(face = "bold"),
      axis.line        = ggplot2::element_line(linewidth = 1.1, color = "black"),
      axis.text        = ggplot2::element_text(face = "bold", color = "black"),
      plot.title       = ggplot2::element_text(face = "bold"),
      legend.text      = ggplot2::element_text(face = "bold")
    )

  # ── Regression line ────────────────────────────────────────────────────────
  #
  # The line is drawn only across the standard concentration range
  # (min to max Scaled_Concentration of standards), not extrapolated
  # beyond it. This correctly represents the validated calibration range.
  if (!is.na(slope) && !is.na(intercept)) {
    x_std <- df_group %>%
      dplyr::filter(Replicate_Type == "Standard") %>%
      dplyr::pull(Scaled_Concentration)

    x_min <- if (length(x_std) > 0 && any(!is.na(x_std))) {
      min(x_std, na.rm = TRUE)
    } else {
      min(dd$x_to_plot, na.rm = TRUE)
    }
    x_max <- if (length(x_std) > 0 && any(!is.na(x_std))) {
      max(x_std, na.rm = TRUE)
    } else {
      max(dd$x_to_plot, na.rm = TRUE)
    }

    line_df   <- data.frame(x = c(x_min, x_max))
    line_df$y <- slope * line_df$x + intercept

    g <- g + ggplot2::geom_line(
      data         = line_df,
      ggplot2::aes(x = x, y = y),
      inherit.aes  = FALSE,
      color        = "black",
      linewidth    = 0.8
    )
  } else {
    # No fit available — annotate rather than leaving the plot blank
    g <- g + ggplot2::annotate(
      "text", x = Inf, y = Inf,
      label = "No fit available",
      hjust = 1.1, vjust = 1.5, size = 3.5
    )
  }

  # ── LOD / LOQ vertical lines ───────────────────────────────────────────────
  #
  # LOD and LOQ are looked up from the lod_loq table by Batch_ID and
  # Molecule_Prefix. The x-axis of calibration curve plots is concentration
  # (ug/g space), so the LOD/LOQ values from calculate_lod_loq() are
  # plotted directly as vertical lines without any further conversion.
  #
  # Lines are only drawn when the threshold is finite and > 0. A threshold
  # of 0 (from zero blank SD) or NA (from missing blanks) is silently skipped.
  lod_row <- lod_loq %>%
    dplyr::filter(Batch_ID == batch, Molecule_Prefix == analyte) %>%
    dplyr::slice(1)

  if (nrow(lod_row) == 1) {
    if (y_col == "Normalized_Area_BG_Sub") {
      lod_x <- lod_row$Norm_LOD
      loq_x <- lod_row$Norm_LOQ
    } else {
      lod_x <- lod_row$Unnorm_LOD
      loq_x <- lod_row$Unnorm_LOQ
    }

    thresh_df <- tibble::tibble(
      type = c("LOD", "LOQ"),
      x    = c(lod_x, loq_x)
    ) %>%
      dplyr::filter(is.finite(x), x > 0)

    if (nrow(thresh_df) > 0) {
      g <- g +
        ggplot2::geom_vline(
          data        = thresh_df,
          ggplot2::aes(xintercept = x, linetype = type),
          color       = "red",
          linewidth   = 0.8,
          show.legend = TRUE
        ) +
        ggplot2::scale_linetype_manual(
          name   = "Thresholds",
          breaks = c("LOD", "LOQ"),
          values = c("LOD" = "solid", "LOQ" = "dotted")
        ) +
        ggplot2::guides(
          # Prevent red threshold lines from appearing in the color legend keys
          color    = ggplot2::guide_legend(
            order        = 1,
            override.aes = list(linetype = "blank")
          ),
          # Separate compact legend for LOD/LOQ
          linetype = ggplot2::guide_legend(
            order     = 2,
            keywidth  = grid::unit(1.4, "lines"),
            keyheight = grid::unit(0.5, "lines"),
            override.aes = list(color = "red", linewidth = 0.8)
          )
        )
    }
  }

  g
}


# ══ Exported pipeline function ════════════════════════════════════════════════

# ── export_calibration_pdfs ───────────────────────────────────────────────────
#
#' Export calibration curve PDFs for all batches and analytes
#'
#' Iterates over every batch and analyte in the dataset and generates two
#' PDF files per batch — one for normalized signal and one for unnormalized
#' signal — each containing one calibration curve page per analyte.
#'
#' Each plot shows standards at their known scaled concentrations, non-standard
#' replicates at their back-calculated concentrations, the fitted weighted 1/x
#' regression line, and LOD/LOQ vertical threshold lines.
#'
#' Output filenames follow the pattern:
#' \itemize{
#'   \item \code{Cal_Curves_Normalized_Batch_<BatchID>.pdf}
#'   \item \code{Cal_Curves_Unnormalized_Batch_<BatchID>.pdf}
#' }
#'
#' @param all_batches_subtracted A data frame as returned by
#'   \code{calculate_concentrations()}.
#' @param lod_loq A data frame as returned by \code{calculate_lod_loq()}.
#' @param out_dir Character. Directory where PDF files will be written.
#'   Defaults to the current working directory (\code{"."}).
#'   The directory must already exist.
#' @param width  Numeric. PDF page width in inches. Default \code{9.25}.
#' @param height Numeric. PDF page height in inches. Default \code{6}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   writing PDF files to \code{out_dir}.
#'
#' @details
#' \strong{File naming:} Batch IDs are sanitized for use in filenames —
#' any character that is not alphanumeric, a dot, underscore, or hyphen
#' is replaced with an underscore.
#'
#' \strong{Page order:} Within each PDF, analytes are ordered
#' alphabetically by \code{Molecule_Prefix}.
#'
#' \strong{What to modify:} The \code{out_dir}, \code{width}, and
#' \code{height} arguments. Plot colors, theme, and layout are controlled
#' by \code{make_curve_plot()} and are not exposed as arguments here to
#' keep the user interface simple. To change plot styling, modify
#' \code{make_curve_plot()} directly.
#'
#' @examples
#' \dontrun{
#' export_calibration_pdfs(
#'   all_data,
#'   lod_loq,
#'   out_dir = "results/calibration_plots"
#' )
#' }
#'
#' @export
export_calibration_pdfs <- function(all_batches_subtracted,
                                    lod_loq,
                                    out_dir = ".",
                                    width   = 9.25,
                                    height  = 6) {

  # Validate base output directory
  if (!dir.exists(out_dir)) {
    stop("Output directory does not exist: ", out_dir,
         "\nCreate it first with dir.create('", out_dir, "')")
  }

  # Create calibration_curves subfolder automatically if it doesn't exist.
  # All calibration curve PDFs are written here to keep the results folder clean.
  plot_dir <- file.path(out_dir, "calibration_curves")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    message("Created output subfolder: ", plot_dir)
  }

  # Sanitize batch IDs for use in filenames
  sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

  batches <- sort(unique(all_batches_subtracted$Batch_ID))
  message("Exporting calibration curve PDFs for ", length(batches),
          " batch(es) to: ", plot_dir)

  for (b in batches) {
    df_b     <- all_batches_subtracted %>%
      dplyr::filter(Batch_ID == b)
    analytes <- sort(unique(df_b$Molecule_Prefix))

    # ── Normalized PDF ────────────────────────────────────────────────────
    norm_path <- file.path(
      plot_dir,
      paste0("Cal_Curves_Normalized_Batch_", sanitize(b), ".pdf")
    )
    grDevices::pdf(file = norm_path, width = width, height = height)
    for (m in analytes) {
      df_g <- df_b %>% dplyr::filter(Molecule_Prefix == m)
      print(make_curve_plot(df_g, lod_loq,
                            y_col = "Normalized_Area_BG_Sub"))
    }
    grDevices::dev.off()
    message("  Written: ", basename(norm_path))

    # ── Unnormalized PDF ──────────────────────────────────────────────────
    unnorm_path <- file.path(
      plot_dir,
      paste0("Cal_Curves_Unnormalized_Batch_", sanitize(b), ".pdf")
    )
    grDevices::pdf(file = unnorm_path, width = width, height = height)
    for (m in analytes) {
      df_g <- df_b %>% dplyr::filter(Molecule_Prefix == m)
      print(make_curve_plot(df_g, lod_loq,
                            y_col = "Area_BG_Sub"))
    }
    grDevices::dev.off()
    message("  Written: ", basename(unnorm_path))
  }

  message("Calibration curve export complete.")
  invisible(NULL)
}
