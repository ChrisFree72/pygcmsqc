# R/plots_concentrations.R
#
# Functions for generating bar plots of calculated analyte concentrations
# (ug/g) for all replicates in each batch, exported as per-batch PDFs.
#
# Each plot shows:
#   - One bar per replicate colored by Replicate_Type (Unknown, QC, Blank,
#     Cup_Blank, Standard), with bar height equal to the calculated ug/g
#     concentration (negative values and NAs are treated as 0 for display)
#   - Value labels printed above each bar for rapid data review
#   - Horizontal LOD and LOQ threshold lines from calculate_lod_loq(),
#     allowing immediate visual assessment of which samples exceed detection
#     and quantification limits
#   - Removed calibrator annotation in the subtitle
#
# Two PDFs are produced per batch:
#   - Normalized: uses Norm_ug_Per_Gram (IS-normalized signal)
#   - Unnormalized: uses Unnorm_ug_Per_Gram (raw signal)
#
# Note on negative concentrations:
#   Negative ug/g values (samples below the blank level) are valid analytical
#   results preserved in All_Batches_Subtracted. For display purposes only,
#   they are plotted as 0 in these bar charts (a bar cannot have negative
#   height). The underlying data is unchanged.
#
# Output subfolder: concentration_bars/
#
# Contents:
#   Internal helpers (not exported):
#     prepare_bar_data()   — assigns display value and x-axis label columns
#     make_bar_plot()      — builds a single ggplot concentration bar chart
#
#   Exported pipeline function:
#     export_concentration_bar_pdfs() — loops over batches/analytes, writes PDFs


# ══ Internal helpers ══════════════════════════════════════════════════════════

# ── prepare_bar_data ──────────────────────────────────────────────────────────
#
# Prepares a data frame for bar plotting by:
#   1. Replacing NA and negative concentration values with 0 for display
#      (negative values are analytically valid but cannot be shown as bars)
#   2. Assigning an x-axis label column from the specified label column
#      (typically Replicate_Name), preserving input row order via fct_inorder()
#
# @param df_group  A data frame for a single Batch_ID × Molecule_Prefix group.
# @param conc_col  Character. Which concentration column to plot.
#                  One of "Norm_ug_Per_Gram" or "Unnorm_ug_Per_Gram".
# @param label_col Character or NULL. Column to use for x-axis labels.
#                  If NULL, attempts to use "Replicate_Name"; falls back to
#                  row numbers if that column is absent.
# @return The input data frame with two new columns:
#         value_to_plot (numeric, >= 0) and x_label (ordered factor).
#
#' @keywords internal
prepare_bar_data <- function(df_group,
                             conc_col  = c("Norm_ug_Per_Gram",
                                           "Unnorm_ug_Per_Gram"),
                             label_col = NULL) {
  conc_col <- match.arg(conc_col)

  # Auto-detect label column if not specified
  if (is.null(label_col)) {
    candidates <- c("Replicate_Name")
    label_col  <- candidates[candidates %in% names(df_group)][1]
  }

  df_group %>%
    dplyr::mutate(
      # Clip NA and negative values to 0 for bar display only.
      # The underlying concentration columns are not modified.
      value_to_plot = ifelse(
        is.na(.data[[conc_col]]) | .data[[conc_col]] < 0,
        0,
        .data[[conc_col]]
      ),
      # x-axis label: use the specified column if present, else row number.
      # fct_inorder() preserves the input row sequence so bars appear in
      # the same order as the data (injection sequence order).
      x_label = forcats::fct_inorder(
        as.factor(
          if (!is.null(label_col) && label_col %in% names(df_group)) {
            .data[[label_col]]
          } else {
            as.character(dplyr::row_number())
          }
        )
      )
    )
}


# ── make_bar_plot ─────────────────────────────────────────────────────────────
#
# Builds a single concentration bar chart ggplot for one analyte in one batch.
#
# Plot elements:
#   - Bars: filled by Replicate_Type, height = ug/g (clipped to 0 if negative)
#   - Value labels: printed above each bar rotated 90° for readability
#   - LOD line: solid red horizontal line at the LOD threshold (ug/g)
#   - LOQ line: dotted red horizontal line at the LOQ threshold (ug/g)
#   - Subtitle: removed calibrators annotation if any were dropped
#
# LOD/LOQ thresholds come from the lod_loq table (calculate_lod_loq()).
# They are in ug/g units and are plotted directly as horizontal lines
# on the ug/g y-axis without any further conversion.
#
# @param df_group  A data frame for one Batch_ID × Molecule_Prefix group.
# @param lod_loq   The full LOD/LOQ table from calculate_lod_loq().
# @param conc_col  Character. Concentration column to plot.
# @param label_col Character or NULL. Column for x-axis labels.
# @return A ggplot object.
#
#' @keywords internal
make_bar_plot <- function(df_group,
                          lod_loq,
                          conc_col  = c("Norm_ug_Per_Gram",
                                        "Unnorm_ug_Per_Gram"),
                          label_col = NULL) {
  conc_col <- match.arg(conc_col)

  # Pull scalar metadata from the first row
  pdat    <- df_group %>% dplyr::slice(1)
  batch   <- pdat$Batch_ID
  analyte <- pdat$Molecule_Prefix

  # Signal-type-specific titles and y-axis labels
  is_norm <- conc_col == "Norm_ug_Per_Gram"
  title   <- paste0(
    "Concentrations ug/g: ", analyte,
    " (Batch ", batch, ") \u2014 ",
    if (is_norm) "Normalized" else "Unnormalized"
  )
  ylab <- if (is_norm) "Normalized ug/g" else "Unnormalized ug/g"

  # Subtitle: show removed calibrators if any were dropped during optimization
  removed_txt <- if (!is.null(pdat$Removed_Calibrators) &&
                     !is.na(pdat$Removed_Calibrators) &&
                     pdat$Removed_Calibrators != "None") {
    paste0("Removed Calibrators: ", pdat$Removed_Calibrators)
  } else ""

  # Prepare display data (clips negatives/NAs to 0, assigns x_label)
  dd <- prepare_bar_data(df_group, conc_col = conc_col, label_col = label_col)

  # Fixed color map and legend order for Replicate_Type
  replicate_categories <- c("Unknown", "QC", "Blank", "Cup_Blank", "Standard")
  color_map <- c(
    "Unknown"   = "maroon",
    "QC"        = "darkgreen",
    "Standard"  = "darkgray",
    "Blank"     = "blue",
    "Cup_Blank" = "orange"
  )

  # ── Base bar plot ──────────────────────────────────────────────────────────
  g <- ggplot2::ggplot(
    dd,
    ggplot2::aes(x = x_label, y = value_to_plot, fill = Replicate_Type)
  ) +
    ggplot2::geom_col(width = 1, color = "black", linewidth = 0.2) +
    # Value labels rotated 90°, positioned just above each bar
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", value_to_plot)),
      vjust = 0.5, hjust = -0.2, angle = 90, size = 3
    ) +
    ggplot2::scale_fill_manual(
      values = color_map,
      limits = replicate_categories,
      drop   = FALSE
    ) +
    # scale_x_discrete MUST be outside theme() — it is a scale layer, not
    # a theme element. In the original script this was incorrectly placed
    # inside theme() where it was silently ignored.
    ggplot2::scale_x_discrete(
      expand = ggplot2::expansion(mult = c(0.05, 0.05))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.2))
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = stringr::str_wrap(removed_txt, width = 80),
      x        = "Samples",
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
      legend.text      = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      plot.margin      = ggplot2::margin(10, 20, 10, 20)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(nrow = 1)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  # ── LOD / LOQ horizontal lines ─────────────────────────────────────────────
  #
  # LOD and LOQ are in ug/g units (from calculate_lod_loq()) and are plotted
  # directly as horizontal lines on the ug/g y-axis.
  #
  # Calibration_Factor is NOT applied here — it is already embedded in the
  # fitted slope and therefore already reflected in the ug/g values from
  # which the blank SD (and hence LOD/LOQ) were derived.
  #
  # Lines are drawn only when the threshold is finite, non-NA, and > 0.
  # A threshold of 0 (zero blank SD) or NA (missing blanks) is skipped.
  # If a threshold exactly equals 0, a tiny epsilon nudge makes it visible
  # just above the x-axis rather than sitting on it invisibly.
  lod_row <- lod_loq %>%
    dplyr::filter(Batch_ID == batch, Molecule_Prefix == analyte) %>%
    dplyr::slice(1)

  if (nrow(lod_row) == 1) {
    lod_y <- if (is_norm) lod_row$Norm_LOD   else lod_row$Unnorm_LOD
    loq_y <- if (is_norm) lod_row$Norm_LOQ   else lod_row$Unnorm_LOQ

    # Compute a small epsilon for zero-threshold nudging
    y_max <- suppressWarnings(max(dd$value_to_plot, na.rm = TRUE))
    y_min <- suppressWarnings(min(dd$value_to_plot, na.rm = TRUE))
    eps   <- (y_max - y_min) * 0.01
    if (!is.finite(eps) || eps <= 0) eps <- 1e-6

    # Nudge zero thresholds to eps so they remain visible on the plot
    lod_y_plot <- if (is.finite(lod_y) && !is.na(lod_y)) {
      if (lod_y <= 0) eps else lod_y
    } else NA_real_

    loq_y_plot <- if (is.finite(loq_y) && !is.na(loq_y)) {
      if (loq_y <= 0) eps else loq_y
    } else NA_real_

    thresh_df <- data.frame(
      type = c("LOD", "LOQ"),
      y    = c(lod_y_plot, loq_y_plot)
    )
    thresh_df <- thresh_df[
      is.finite(thresh_df$y) & !is.na(thresh_df$y), ,
      drop = FALSE
    ]

    if (nrow(thresh_df) > 0) {
      g <- g +
        ggplot2::geom_hline(
          data        = thresh_df,
          ggplot2::aes(yintercept = y, linetype = type),
          color       = "red",
          linewidth   = 0.8,
          show.legend = TRUE
        ) +
        ggplot2::scale_linetype_manual(
          name   = "Thresholds",
          breaks = c("LOD", "LOQ"),
          values = c(LOD = "solid", LOQ = "dotted")
        ) +
        ggplot2::guides(
          # Prevent red lines from appearing inside the fill legend keys
          fill     = ggplot2::guide_legend(
            nrow         = 1,
            override.aes = list(linetype = "blank")
          ),
          # Separate compact legend for LOD/LOQ lines
          linetype = ggplot2::guide_legend(
            order        = 2,
            keywidth     = grid::unit(1.4, "lines"),
            keyheight    = grid::unit(0.6, "lines"),
            override.aes = list(color = "red", linewidth = 0.8)
          ),
          # Suppress any auto-generated color legend for the red lines
          color    = "none"
        )
    }
  }

  g
}


# ══ Exported pipeline function ════════════════════════════════════════════════

# ── export_concentration_bar_pdfs ─────────────────────────────────────────────
#
#' Export concentration bar chart PDFs for all batches and analytes
#'
#' Iterates over every batch and analyte and generates two PDF files per
#' batch — one for normalized ug/g concentrations and one for unnormalized
#' — each containing one bar chart page per analyte. Bars are colored by
#' Replicate_Type and annotated with LOD/LOQ threshold lines.
#'
#' Output filenames follow the pattern:
#' \itemize{
#'   \item \code{Bar_CalcConcs_Normalized_Batch_<BatchID>.pdf}
#'   \item \code{Bar_CalcConcs_Unnormalized_Batch_<BatchID>.pdf}
#' }
#'
#' @param all_batches_subtracted A data frame as returned by
#'   \code{calculate_concentrations()}.
#' @param lod_loq A data frame as returned by \code{calculate_lod_loq()}.
#' @param out_dir Character. Base results directory. PDFs are written to a
#'   \code{concentration_bars/} subfolder. Default \code{"."}.
#' @param width     Numeric. PDF page width in inches. Default \code{9.25}.
#' @param height    Numeric. PDF page height in inches. Default \code{6}.
#' @param label_col Character or NULL. Column name to use for x-axis bar
#'   labels. Default \code{"Replicate_Name"}. Pass \code{NULL} to use
#'   row numbers.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   writing PDF files.
#'
#' @details
#' \strong{Negative concentrations:} Values below zero are displayed as 0
#' in the bar charts (a bar cannot have negative height). The underlying
#' data in \code{all_batches_subtracted} is not modified.
#'
#' \strong{Page order:} Within each PDF, analytes are ordered
#' alphabetically by \code{Molecule_Prefix}.
#'
#' \strong{What to modify:} \code{out_dir}, \code{width}, \code{height},
#' and \code{label_col}. To change plot styling, modify \code{make_bar_plot()}
#' directly.
#'
#' @examples
#' \dontrun{
#' export_concentration_bar_pdfs(
#'   all_data,
#'   lod_loq,
#'   out_dir   = "results/",
#'   label_col = "Replicate_Name"
#' )
#' }
#'
#' @export
export_concentration_bar_pdfs <- function(all_batches_subtracted,
                                          lod_loq,
                                          out_dir   = ".",
                                          width     = 9.25,
                                          height    = 6,
                                          label_col = "Replicate_Name") {

  # Validate base output directory
  if (!dir.exists(out_dir)) {
    stop("Output directory does not exist: ", out_dir,
         "\nCreate it first with dir.create('", out_dir, "')")
  }

  # Create concentration_bars subfolder automatically if it doesn't exist
  plot_dir <- file.path(out_dir, "concentration_bars")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    message("Created output subfolder: ", plot_dir)
  }

  sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

  batches <- sort(unique(all_batches_subtracted$Batch_ID))
  message("Exporting concentration bar PDFs for ", length(batches),
          " batch(es) to: ", plot_dir)

  for (b in batches) {
    df_b     <- all_batches_subtracted %>% dplyr::filter(Batch_ID == b)
    analytes <- sort(unique(df_b$Molecule_Prefix))

    # ── Normalized PDF ────────────────────────────────────────────────────────
    norm_path <- file.path(
      plot_dir,
      paste0("Bar_CalcConcs_Normalized_Batch_", sanitize(b), ".pdf")
    )
    grDevices::pdf(file = norm_path, width = width, height = height)
    for (m in analytes) {
      df_g <- df_b %>% dplyr::filter(Molecule_Prefix == m)
      print(make_bar_plot(df_g, lod_loq,
                          conc_col  = "Norm_ug_Per_Gram",
                          label_col = label_col))
    }
    grDevices::dev.off()
    message("  Written: ", basename(norm_path))

    # ── Unnormalized PDF ──────────────────────────────────────────────────────
    unnorm_path <- file.path(
      plot_dir,
      paste0("Bar_CalcConcs_Unnormalized_Batch_", sanitize(b), ".pdf")
    )
    grDevices::pdf(file = unnorm_path, width = width, height = height)
    for (m in analytes) {
      df_g <- df_b %>% dplyr::filter(Molecule_Prefix == m)
      print(make_bar_plot(df_g, lod_loq,
                          conc_col  = "Unnorm_ug_Per_Gram",
                          label_col = label_col))
    }
    grDevices::dev.off()
    message("  Written: ", basename(unnorm_path))
  }

  message("Concentration bar export complete.")
  invisible(NULL)
}
