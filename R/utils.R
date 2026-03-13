# R/utils.R
#
# Internal helper functions shared across the pygcmsqc pipeline.
# None of these are exported to the user — they are workhorses called
# internally by the pipeline functions in other R/ files.
#
# Contents:
#   1. to_numeric_clean()       — robust character-to-numeric coercion
#   2. should_be_numeric()      — column-type decision rule
#   3. coerce_mixed_columns()   — applies type coercion across a data frame
#   4. first_non_na()           — first non-NA value from a vector
#   5. theme_minimal_bold()     — shared ggplot2 theme for PCA plots
#   6. .qc_bar_theme()          — shared ggplot2 theme for QC bar plots
#   7. .qc_cleanup()            — filters and coerces QC data before plotting


# ── 1. to_numeric_clean ───────────────────────────────────────────────────────
#
# Converts a character or factor vector to numeric, handling common data
# artifacts found in Skyline CSV exports:
#   - Factor levels are dropped first to avoid integer encoding
#   - Whitespace is stripped
#   - Common NA-like tokens (#N/A, N/A, NA, empty string) become NA
#   - Infinity symbols (∞, Inf, +Inf, -Inf) become NA rather than Inf,
#     which would silently corrupt downstream sum/mean calculations
#   - Any numeric Inf that survives conversion is also clamped to NA
#
# @param x A character or factor vector.
# @return A numeric vector of the same length; non-parsable entries are NA.
#
#' @keywords internal
to_numeric_clean <- function(x) {
  # Drop factor levels to avoid R encoding factor as integer indices
  if (is.factor(x)) x <- as.character(x)

  # Normalize whitespace (tabs, multiple spaces, leading/trailing)
  x <- stringr::str_squish(as.character(x))

  # Treat common NA-like tokens as true NA before numeric conversion.
  # Case-insensitive match via toupper().
  na_tokens <- c("#N/A", "N/A", "NA", "NA", "")
  x[trimws(toupper(x)) %in% toupper(na_tokens)] <- NA_character_

  # Infinity representations become NA — Inf would silently propagate through
  # sum() and mean() calls and corrupt concentration calculations.
  x[x %in% c("\u221e", "Inf", "+Inf", "-Inf")] <- NA_character_

  # Suppress warnings: non-numeric strings that survived the above filters
  # will produce NA with a warning; we expect and accept that.
  out <- suppressWarnings(as.numeric(x))

  # Final safety net: catch any numeric Inf that slipped through conversion
  out[is.infinite(out)] <- NA_real_
  out
}


# ── 2. should_be_numeric ─────────────────────────────────────────────────────
#
# Decides whether a character or factor column contains predominantly numeric
# data and should be coerced. Only character and factor columns are evaluated;
# already-numeric columns are ignored.
#
# The threshold (default 0.5) is the minimum proportion of values that must
# parse successfully as numeric. A value of 0.5 means: if more than half the
# non-NA entries look like numbers, treat the whole column as numeric.
# This is intentionally conservative to avoid coercing mixed ID columns.
#
# @param x         A vector (any type).
# @param threshold Numeric scalar in (0, 1). Minimum success rate to coerce.
# @return Logical scalar.
#
#' @keywords internal
should_be_numeric <- function(x, threshold = 0.5) {
  # Only consider character and factor columns
  if (!(is.character(x) || is.factor(x))) return(FALSE)

  x_chr <- as.character(x)

  # Quick screen: if there are no digits at all, skip the full parse
  if (!any(grepl("\\d", x_chr))) return(FALSE)

  # Attempt parse and compute the proportion that succeeded
  num          <- suppressWarnings(as.numeric(trimws(gsub(",", "", x_chr))))
  success_rate <- mean(!is.na(num))
  success_rate > threshold
}


# ── 3. coerce_mixed_columns ───────────────────────────────────────────────────
#
# Applies column-type coercion across an entire data frame:
#   - Factors are first converted to character to prevent accidental integer encoding
#   - Numeric columns are left untouched
#   - Character columns are coerced to numeric if should_be_numeric() returns TRUE
#   - Remaining character columns stay as character
#
# This is applied to each batch CSV after loading, before binding rows, because
# Skyline may export the same column as character in one batch and numeric in
# another depending on the presence of error strings like "#N/A".
#
# @param df        A data frame.
# @param threshold Passed to should_be_numeric(); default 0.5.
# @return A data frame with corrected column types.
#
#' @keywords internal
coerce_mixed_columns <- function(df, threshold = 0.5) {
  # Convert all factors to character first
  df <- df %>% dplyr::mutate(dplyr::across(where(is.factor), as.character))

  for (nm in names(df)) {
    col <- df[[nm]]

    # Skip columns that are already numeric or non-character types (dates, logicals)
    if (is.numeric(col)) next

    if (is.character(col) || is.factor(col)) {
      if (should_be_numeric(col, threshold = threshold)) {
        # Coerce: strip commas and whitespace, then parse
        df[[nm]] <- suppressWarnings(
          as.numeric(trimws(gsub(",", "", as.character(col))))
        )
      } else {
        df[[nm]] <- as.character(col)
      }
    }
  }
  df
}


# ── 4. first_non_na ───────────────────────────────────────────────────────────
#
# Returns the first non-NA element of a vector, or NA if all elements are NA.
# Used in pivot_wider() to resolve duplicate (replicate, analyte) pairs cleanly
# without introducing arbitrary aggregation functions like mean().
#
# @param x A vector of any type.
# @return A scalar of the same type as x, or NA.
#
#' @keywords internal
first_non_na <- function(x) {
  out <- x[!is.na(x)]
  if (length(out) == 0) NA else out[1]
}


# ── 5. theme_minimal_bold ─────────────────────────────────────────────────────
#
# A ggplot2 theme used for PCA and scree plots throughout the package.
# Built on theme_minimal() with bold axis/strip/legend text, no grid lines,
# and explicit black axis lines to match the style of other plots in the pipeline.
#
# @param base_size Numeric. Base font size passed to theme_minimal(). Default 10.
# @return A ggplot2 theme object.
#
#' @keywords internal
theme_minimal_bold <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = base_size + 2),
      strip.text       = ggplot2::element_text(face = "bold", size = base_size + 1),
      axis.title       = ggplot2::element_text(face = "bold", size = base_size + 1),
      axis.text        = ggplot2::element_text(size = base_size, face = "bold", color = "black"),
      legend.title     = ggplot2::element_text(face = "bold", size = base_size + 1),
      legend.text      = ggplot2::element_text(size = base_size),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line.x      = ggplot2::element_line(color = "black", linewidth = 1),
      axis.line.y      = ggplot2::element_line(color = "black", linewidth = 1)
    )
}


# ── 6. .qc_bar_theme ─────────────────────────────────────────────────────────
#
# A ggplot2 theme used for all QC bar plots (both all-batches and per-batch).
# Centralizes the shared theme so changes propagate to all QC bar plots at once.
# Named with a leading dot to signal that it is an internal convention, not
# a user-facing function.
#
# @return A ggplot2 theme object.
#
#' @keywords internal
.qc_bar_theme <- function() {
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
    )
}


# ── 7. .qc_cleanup ───────────────────────────────────────────────────────────
#
# Filters a data frame to QC replicates only (if Replicate_Type is present),
# coerces key columns to the expected types, and drops rows with missing or
# non-finite intensity values. Called at the top of every QC plotting function
# to guarantee a clean, consistent input regardless of what the user passes in.
#
# @param df A data frame. Expected columns: Batch_ID, Molecule_Prefix,
#           Replicate_Name, Area_BG_Sub. Optionally: Replicate_Type.
# @return A filtered and type-corrected data frame.
#
#' @keywords internal
.qc_cleanup <- function(df) {
  # Filter to QC replicates if the column is present
  if ("Replicate_Type" %in% names(df)) {
    df <- df %>% dplyr::filter(Replicate_Type == "QC")
  }

  df %>%
    dplyr::mutate(
      Batch_ID        = as.character(Batch_ID),
      Molecule_Prefix = as.character(Molecule_Prefix),
      Replicate_Name  = as.character(Replicate_Name),
      # Suppress warnings: non-finite values are intentionally filtered below
      Area_BG_Sub     = suppressWarnings(as.numeric(Area_BG_Sub))
    ) %>%
    # Drop rows with missing keys or non-finite intensity.
    # is.finite() catches NA, NaN, and Inf in a single check.
    dplyr::filter(
      !is.na(Batch_ID),
      !is.na(Molecule_Prefix),
      !is.na(Replicate_Name),
      is.finite(Area_BG_Sub)
    )
}
