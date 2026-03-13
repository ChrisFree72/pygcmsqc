# R/export.R
#
# Functions for exporting pipeline results to styled Excel workbooks.
#
# Three workbooks are produced:
#
#   1. R2_Comparison.xlsx (export_r2_table):
#      One sheet containing the calibration model parameter table from
#      run_calibration() — R² values, slopes, intercepts, removed
#      calibrators, and improvement flags for every Batch_ID ×
#      Molecule_Prefix combination. Used to review and audit calibration
#      quality across batches.
#
#   2. Total_Batch_Info.xlsx (export_batch_info):
#      One sheet containing the full All_Batches_Subtracted data frame —
#      every row and column including raw signals, background-subtracted
#      signals, calibration metadata, calculated concentrations, and ug/g
#      values. This is the master output table for downstream analysis.
#
#   3. Final_Sample_Data_Output.xlsx (export_wide_concentrations):
#      Two sheets — "Normalized" and "Unnormalized" — each containing a
#      wide-format pivot of Unknown replicate ug/g concentrations with one
#      row per replicate and one column per analyte. Rows are banded by
#      Batch_ID with a light pastel color for readability. Negative values
#      and NAs are replaced with 0. Values are rounded to 3 significant
#      figures. This is the primary deliverable table for reporting.
#
# All workbooks use consistent styling:
#   - Bold, gray-filled header row
#   - Frozen header row (scroll-safe)
#   - Auto-fitted column widths
#   - Batch-banded row colors (Final_Sample_Data_Output only)
#
# Output subfolder: All three files are written directly to out_dir
#   (no subfolder) since they are primary data deliverables, not
#   diagnostic plots.
#
# Contents:
#   Internal helpers (not exported):
#     build_wide()               — pivots Unknown ug/g to wide format
#     .make_light_palette()      — generates pastel batch band colors
#     .style_sheet_by_batch()    — applies header and batch band styles
#
#   Exported pipeline functions:
#     export_r2_table()           — writes R2_Comparison.xlsx
#     export_batch_info()         — writes Total_Batch_Info.xlsx
#     export_wide_concentrations() — writes Final_Sample_Data_Output.xlsx


# ══ Internal helpers ══════════════════════════════════════════════════════════

# ── build_wide ────────────────────────────────────────────────────────────────
#
# Pivots Unknown replicate ug/g concentration data from long to wide format.
# Each row is one replicate; each column is one analyte (Molecule_Prefix).
# Only Unknown replicates are included — Standards, QC, and Blanks are
# excluded from the reportable concentration table.
#
# Post-pivot processing:
#   - NA and negative values are replaced with 0 (below-blank or missing
#     measurements reported as non-detect)
#   - Values are rounded to 3 significant figures for reporting clarity
#   - Replicate_Type column is dropped (all rows are Unknown)
#   - Analyte columns are sorted alphabetically
#
# @param df       A data frame from calculate_concentrations().
# @param calc_col Character. Which concentration column to pivot.
#                 "Norm_ug_Per_Gram" or "Unnorm_ug_Per_Gram".
# @return A wide data frame with columns: Batch_ID, Replicate_Name,
#         then one column per Molecule_Prefix.
#
#' @keywords internal
build_wide <- function(df, calc_col) {

  # Filter to Unknown replicates only
  df_unknown <- df %>% dplyr::filter(Replicate_Type == "Unknown")

  if (nrow(df_unknown) == 0) {
    warning("No Unknown replicates found. Wide export will be empty.")
    return(data.frame())
  }

  # Prepare long form for pivoting.
  # Convert calc_col to character before pivoting to avoid type conflicts
  # when some batches have numeric and others have NA for the same analyte.
  combined_long <- df_unknown %>%
    dplyr::transmute(
      Batch_ID, Replicate_Name, Replicate_Type,
      Molecule_Prefix,
      value = as.character(.data[[calc_col]])
    )

  # Pivot wide: one column per Molecule_Prefix.
  # first_non_na() (utils.R) resolves any duplicate replicate/analyte pairs
  # without arbitrary aggregation — takes the first non-NA value encountered.
  wide <- combined_long %>%
    tidyr::pivot_wider(
      names_from  = Molecule_Prefix,
      values_from = value,
      values_fn   = list(value = first_non_na)
    ) %>%
    dplyr::arrange(Batch_ID, Replicate_Name)

  # Column ordering: ID columns first, then analytes alphabetically
  id_cols  <- c("Batch_ID", "Replicate_Name", "Replicate_Type")
  mol_cols <- sort(setdiff(names(wide), id_cols))
  wide     <- wide[, c(id_cols, mol_cols)]

  # Post-pivot cleaning:
  #   1. Parse character back to numeric (pivot coerces to character)
  #   2. Replace NA and negative values with 0 (non-detect reporting)
  #   3. Round to 3 significant figures for reporting
  #   4. Drop Replicate_Type (all rows are Unknown, redundant in wide format)
  wide %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(mol_cols),
        ~ {
          x           <- suppressWarnings(as.numeric(.))
          x[is.na(x) | x < 0] <- 0
          x
        }
      )
    ) %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(mol_cols), ~ signif(.x, 3))
    ) %>%
    dplyr::select(-Replicate_Type)
}


# ── .make_light_palette ───────────────────────────────────────────────────────
#
# Generates a vector of light pastel hex colors for batch row banding.
# Uses a predefined set of 10 distinctive pastels; for more than 10 batches,
# additional very-light HCL tints are generated programmatically.
#
# Excel does not support alpha transparency in cell fills, so pre-tinted
# light colors are used directly rather than applying alpha to base colors.
#
# @param n Integer. Number of distinct colors needed (one per batch).
# @return A character vector of n hex color codes.
#
#' @keywords internal
.make_light_palette <- function(n) {
  base <- c("#FDEBD0", "#D6EAF8", "#E8DAEF", "#D5F5E3", "#FADBD8",
            "#E5E7E9", "#FCF3CF", "#E8F8F5", "#F5EEF8", "#EBDEF0")
  if (n <= length(base)) {
    base[seq_len(n)]
  } else {
    # Generate additional very-light tints using HCL color space.
    # Low chroma (c = 20) and high luminance (l = 95) ensure all generated
    # colors are pale enough to not interfere with text readability.
    grDevices::hcl(h = seq(0, 1, length.out = n) * 360, c = 20, l = 95)
  }
}


# ── .style_sheet_by_batch ─────────────────────────────────────────────────────
#
# Applies consistent Excel styling to a worksheet:
#   - Bold, light gray (#BFBFBF) header row (row 1)
#   - Frozen header row so it stays visible while scrolling
#   - Light pastel row banding by Batch_ID (different color per batch block)
#   - Auto-fitted column widths for readability
#
# Row indices are offset by +1 to account for the header row occupying row 1
# in the Excel sheet (data rows start at row 2).
#
# @param wb    An openxlsx workbook object.
# @param sheet Character. Sheet name to style.
# @param df    The data frame written to this sheet (used to compute row
#              ranges per Batch_ID).
# @return Invisibly returns NULL. Called for its side effects on wb.
#
#' @keywords internal
.style_sheet_by_batch <- function(wb, sheet, df) {
  n_rows <- nrow(df)
  n_cols <- ncol(df)

  if (n_rows == 0 || n_cols == 0) {
    message("Sheet '", sheet, "': no rows or columns to style.")
    return(invisible(NULL))
  }

  # ── Header row: bold text, light gray fill ─────────────────────────────────
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill         = "#BFBFBF"
  )
  openxlsx::addStyle(
    wb, sheet = sheet, style = header_style,
    rows = 1, cols = 1:n_cols, gridExpand = TRUE
  )

  # Freeze the header row so column names stay visible when scrolling down
  openxlsx::freezePane(wb, sheet = sheet, firstRow = TRUE)

  # ── Batch row banding ──────────────────────────────────────────────────────
  #
  # Add a row_idx column (+1 offset for header) to map each data row to its
  # Excel row number, then group by Batch_ID to get contiguous row ranges.
  batches_df   <- df %>% dplyr::mutate(row_idx = dplyr::row_number() + 1)
  uniq_batches <- unique(batches_df$Batch_ID)
  pal          <- .make_light_palette(length(uniq_batches))

  # Create one style object per batch and name them by Batch_ID
  batch_styles <- lapply(pal, function(col) {
    openxlsx::createStyle(fgFill = col)
  })
  names(batch_styles) <- as.character(uniq_batches)

  # Apply each batch's style to its rows.
  # stack = TRUE preserves the bold header style applied above.
  for (b in uniq_batches) {
    rows_b <- batches_df %>%
      dplyr::filter(Batch_ID == b) %>%
      dplyr::pull(row_idx)

    if (length(rows_b) > 0) {
      openxlsx::addStyle(
        wb, sheet = sheet,
        style     = batch_styles[[as.character(b)]],
        rows      = rows_b, cols = 1:n_cols,
        gridExpand = TRUE, stack = TRUE
      )
    }
  }

  # Auto-fit column widths based on content
  openxlsx::setColWidths(wb, sheet = sheet,
                         cols = 1:n_cols, widths = "auto")

  message("  Sheet '", sheet, "': styled ", length(uniq_batches),
          " batch group(s).")
  invisible(NULL)
}


# ══ Exported pipeline functions ═══════════════════════════════════════════════

# ── export_r2_table ───────────────────────────────────────────────────────────
#
#' Export calibration R² comparison table to Excel
#'
#' Writes the calibration model parameter table from \code{run_calibration()}
#' to a single-sheet Excel workbook with a bold, gray header row and
#' auto-fitted column widths. Used to audit calibration quality and review
#' which calibrators were removed by the optimizer for each batch and analyte.
#'
#' @param model_params A data frame as returned by \code{run_calibration()}.
#' @param out_dir Character. Directory where the file will be written.
#'   Default \code{"."}.
#' @param output_file Character. Filename. Default \code{"R2_Comparison.xlsx"}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   writing an Excel file.
#'
#' @details
#' \strong{Columns included:} Batch_ID, Molecule_Prefix, Original_R2_Norm,
#' Original_R2_Unnorm, New_R2_Norm, New_R2_Unnorm, Avg_R2_Gain_Batch,
#' Avg_R2_Gain_Rounded, Improvement_Flag, Removed_Calibrators, Norm_Slope,
#' Unnorm_Slope, Norm_Y_Int, Unnorm_Y_Int.
#'
#' \strong{What to modify:} \code{out_dir} and \code{output_file} only.
#'
#' @examples
#' \dontrun{
#' models <- run_calibration(all_data, max_remove = 2)
#' export_r2_table(models, out_dir = "results/")
#' }
#'
#' @export
export_r2_table <- function(model_params,
                            out_dir     = ".",
                            output_file = "R2_Comparison.xlsx") {

  if (!dir.exists(out_dir)) {
    stop("Output directory does not exist: ", out_dir)
  }

  out_path <- file.path(out_dir, output_file)

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "R2_Comparison")
  openxlsx::writeData(wb, sheet = "R2_Comparison",
                      x = model_params, withFilter = FALSE)

  # Apply header styling (bold + gray) and auto-fit columns
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill         = "#BFBFBF"
  )
  openxlsx::addStyle(
    wb, sheet = "R2_Comparison", style = header_style,
    rows = 1, cols = 1:ncol(model_params), gridExpand = TRUE
  )
  openxlsx::freezePane(wb, sheet = "R2_Comparison", firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet = "R2_Comparison",
                         cols = 1:ncol(model_params), widths = "auto")

  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
  message("Written: ", out_path)
  invisible(NULL)
}


# ── export_batch_info ─────────────────────────────────────────────────────────
#
#' Export the full batch data table to Excel
#'
#' Writes the complete \code{All_Batches_Subtracted} data frame — every row
#' and column — to a single-sheet Excel workbook with a bold, gray header
#' row and auto-fitted column widths. This is the master output table
#' containing all raw signals, background-subtracted signals, calibration
#' metadata, calculated concentrations, and ug/g values.
#'
#' @param all_batches_subtracted A data frame as returned by
#'   \code{calculate_concentrations()}.
#' @param out_dir Character. Directory where the file will be written.
#'   Default \code{"."}.
#' @param output_file Character. Filename.
#'   Default \code{"Total_Batch_Info.xlsx"}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   writing an Excel file.
#'
#' @details
#' \strong{File size note:} This table contains one row per replicate per
#' analyte per batch and can be large for multi-batch projects. If Excel
#' performance is slow, consider exporting as CSV instead using
#' \code{write.csv(all_data, "Total_Batch_Info.csv")}.
#'
#' \strong{What to modify:} \code{out_dir} and \code{output_file} only.
#'
#' @examples
#' \dontrun{
#' export_batch_info(all_data, out_dir = "results/")
#' }
#'
#' @export
export_batch_info <- function(all_batches_subtracted,
                              out_dir     = ".",
                              output_file = "Total_Batch_Info.xlsx") {

  if (!dir.exists(out_dir)) {
    stop("Output directory does not exist: ", out_dir)
  }

  out_path <- file.path(out_dir, output_file)

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "All_Batches")
  openxlsx::writeData(wb, sheet = "All_Batches",
                      x = all_batches_subtracted, withFilter = FALSE)

  # Bold gray header + freeze + auto-width
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill         = "#BFBFBF"
  )
  openxlsx::addStyle(
    wb, sheet = "All_Batches", style = header_style,
    rows = 1, cols = 1:ncol(all_batches_subtracted), gridExpand = TRUE
  )
  openxlsx::freezePane(wb, sheet = "All_Batches", firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet = "All_Batches",
                         cols    = 1:ncol(all_batches_subtracted),
                         widths  = "auto")

  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
  message("Written: ", out_path)
  invisible(NULL)
}


# ── export_wide_concentrations ────────────────────────────────────────────────
#
#' Export wide-format Unknown-only concentration table to Excel
#'
#' Builds wide-format pivot tables of Unknown replicate ug/g concentrations
#' for both normalized and unnormalized signal, then writes them to a
#' two-sheet Excel workbook with full batch-banded styling. This is the
#' primary deliverable table for reporting sample concentrations.
#'
#' Sheet structure:
#' \describe{
#'   \item{Normalized}{Norm_ug_Per_Gram values pivoted wide.}
#'   \item{Unnormalized}{Unnorm_ug_Per_Gram values pivoted wide.}
#' }
#'
#' Each sheet has:
#' \itemize{
#'   \item One row per Unknown replicate
#'   \item One column per analyte (Molecule_Prefix), sorted alphabetically
#'   \item Negative values and NAs replaced with 0
#'   \item Values rounded to 3 significant figures
#'   \item Bold gray header row, frozen for scrolling
#'   \item Row bands colored by Batch_ID (different pastel per batch)
#'   \item Auto-fitted column widths
#' }
#'
#' @param all_batches_subtracted A data frame as returned by
#'   \code{calculate_concentrations()}.
#' @param out_dir Character. Directory where the file will be written.
#'   Default \code{"."}.
#' @param output_file Character. Filename.
#'   Default \code{"Final_Sample_Data_Output.xlsx"}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   writing an Excel file.
#'
#' @details
#' \strong{Unknown-only filter:} Only rows where \code{Replicate_Type ==
#' "Unknown"} are included. Standards, QC replicates, and Blanks are
#' excluded from the reportable concentration table.
#'
#' \strong{Negative concentrations:} Below-blank measurements (negative
#' ug/g) are reported as 0 in this table. The full data including negatives
#' is preserved in \code{Total_Batch_Info.xlsx}.
#'
#' \strong{What to modify:} \code{out_dir} and \code{output_file} only.
#' To change the significant figures, modify the \code{signif(.x, 3)} call
#' inside \code{build_wide()}.
#'
#' @examples
#' \dontrun{
#' export_wide_concentrations(all_data, out_dir = "results/")
#' }
#'
#' @export
export_wide_concentrations <- function(all_batches_subtracted,
                                       out_dir     = ".",
                                       output_file = "Final_Sample_Data_Output.xlsx") {

  if (!dir.exists(out_dir)) {
    stop("Output directory does not exist: ", out_dir)
  }

  out_path <- file.path(out_dir, output_file)

  # Build wide tables for both signal types
  message("Building wide concentration tables...")
  normalized_wide   <- build_wide(all_batches_subtracted,
                                  calc_col = "Norm_ug_Per_Gram")
  unnormalized_wide <- build_wide(all_batches_subtracted,
                                  calc_col = "Unnorm_ug_Per_Gram")

  # Create workbook and add both sheets
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Normalized")
  openxlsx::addWorksheet(wb, "Unnormalized")

  # Write data without Excel auto-filters (cleaner for deliverable tables)
  openxlsx::writeData(wb, sheet = "Normalized",
                      x = normalized_wide,   withFilter = FALSE)
  openxlsx::writeData(wb, sheet = "Unnormalized",
                      x = unnormalized_wide, withFilter = FALSE)

  # Apply batch-banded styling to both sheets
  message("Applying styling...")
  .style_sheet_by_batch(wb, sheet = "Normalized",   df = normalized_wide)
  .style_sheet_by_batch(wb, sheet = "Unnormalized", df = unnormalized_wide)

  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)

  message("Written: ", out_path,
          " (Normalized: ", nrow(normalized_wide), " rows",
          " | Unnormalized: ", nrow(unnormalized_wide), " rows)")
  invisible(NULL)
}
