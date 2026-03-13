# R/lod_loq.R
#
# Computes Limit of Detection (LOD) and Limit of Quantification (LOQ)
# for each analyte in each batch from the distribution of blank replicates.
#
# Method basis:
#   LOD and LOQ are derived from the variability of back-calculated blank
#   concentrations in ug/g units. This approach is preferred over signal-
#   based derivation for py-GC-MS data because:
#
#     1. The output plots (bar charts and calibration curves) display values
#        in ug/g, so LOD/LOQ lines must be in the same units to be
#        interpretable on those axes.
#
#     2. Back-calculating blank signal through the calibration equation
#        (Signal - Intercept) / (Weight x Slope) converts both the blank
#        measurement and its variability into ug/g space, giving thresholds
#        that are directly comparable to sample concentrations.
#
#     3. The SD of back-calculated blank ug/g values correctly propagates
#        the uncertainty from the calibration model (slope, intercept) into
#        the detection limit estimate.
#
#   The intercept IS subtracted during back-calculation of each individual
#   blank replicate because this is a point estimate (position on the
#   calibration curve), not a spread calculation. The SD of those point
#   estimates is then taken, and because SD is invariant to constant shifts,
#   the intercept cancels out in the final LOD/LOQ values - but it is still
#   correct to include it in the per-replicate back-calculation step.
#
#   Standard multipliers (ICH Q2(R1), EPA):
#     LOD = lod_multiplier x SD(blank ug/g)   [default: 3]
#     LOQ = loq_multiplier x SD(blank ug/g)   [default: 10]
#
#   Calibration_Factor is NOT applied to LOD/LOQ after calculation. It is
#   already embedded in the fitted Slope (via Scaled_Concentration), so the
#   back-calculated ug/g values are already in the correct concentration
#   space without any further factor multiplication.
#
# Contents:
#   1. calculate_lod_loq()  - exported pipeline function


# -- calculate_lod_loq ---------------------------------------------------------
#
#' Compute LOD and LOQ from blank replicate distributions
#'
#' Back-calculates ug/g concentrations for all blank replicates using the
#' calibration model parameters from \code{calculate_concentrations()}, then
#' computes the standard deviation of those values within each
#' \code{Batch_ID} x \code{Molecule_Prefix} group. LOD and LOQ are derived
#' as multiples of that SD (default: 3x and 10x respectively).
#'
#' The returned table is used by the plotting functions
#' (\code{export_calibration_pdfs()} and \code{export_concentration_bar_pdfs()})
#' to overlay LOD and LOQ threshold lines on calibration curve and
#' concentration bar plots.
#'
#' @param all_batches_subtracted A data frame as returned by
#'   \code{calculate_concentrations()}, containing columns \code{Batch_ID},
#'   \code{Replicate_Name}, \code{Molecule_Prefix}, \code{Replicate_Type},
#'   \code{Unknown_Weight}, \code{Area}, \code{Normalized_Area},
#'   \code{Norm_Slope}, \code{Unnorm_Slope}, \code{Norm_Y_Int},
#'   \code{Unnorm_Y_Int}, and \code{Calibration_Factor}.
#' @param lod_multiplier Numeric. Multiplier applied to the blank SD to
#'   compute the LOD. Default is \code{3}, per ICH Q2(R1) and EPA guidance.
#'   Change only if your laboratory SOP specifies a different multiplier.
#' @param loq_multiplier Numeric. Multiplier applied to the blank SD to
#'   compute the LOQ. Default is \code{10}, per ICH Q2(R1) and EPA guidance.
#'   Change only if your laboratory SOP specifies a different multiplier.
#'
#' @return A data frame with one row per \code{Batch_ID} x
#'   \code{Molecule_Prefix} combination, containing:
#' \describe{
#'   \item{Batch_ID}{Batch identifier.}
#'   \item{Molecule_Prefix}{Analyte polymer prefix.}
#'   \item{Blank_SD_Unnorm_ug_per_g}{SD of unnormalized blank ug/g values
#'         across blank replicates. Used to derive Unnorm_LOD and Unnorm_LOQ.}
#'   \item{Blank_SD_Norm_ug_per_g}{SD of normalized blank ug/g values.
#'         Used to derive Norm_LOD and Norm_LOQ.}
#'   \item{Norm_Slope / Unnorm_Slope}{Calibration slopes, carried forward
#'         for reference and diagnostic use.}
#'   \item{Norm_Y_Int / Unnorm_Y_Int}{Calibration intercepts, carried
#'         forward for reference.}
#'   \item{Calibration_Factor}{Carried forward for reference.}
#'   \item{Unnorm_LOD}{LOQ threshold in ug/g for unnormalized signal:
#'         \code{lod_multiplier x Blank_SD_Unnorm_ug_per_g}.}
#'   \item{Unnorm_LOQ}{LOQ threshold in ug/g for unnormalized signal:
#'         \code{loq_multiplier x Blank_SD_Unnorm_ug_per_g}.}
#'   \item{Norm_LOD}{LOD threshold in ug/g for normalized signal.}
#'   \item{Norm_LOQ}{LOQ threshold in ug/g for normalized signal.}
#' }
#'
#' @details
#' \strong{Missing blanks:} If no blank replicates exist for a given
#' \code{Batch_ID} x \code{Molecule_Prefix} combination, the SD is
#' \code{NA}, replaced with \code{0} via \code{replace_na()}. This produces
#' LOD = LOQ = 0, which means no threshold line will be drawn for that
#' analyte (the plotting functions filter for \code{x > 0} and
#' \code{is.finite(x)} before drawing lines).
#'
#' \strong{Single blank replicates:} \code{sd()} returns \code{NA} for
#' groups with only one blank replicate (SD is undefined for n = 1). This
#' is also replaced with 0. Having at least 3 blank replicates per batch
#' is recommended for reliable LOD/LOQ estimation.
#'
#' \strong{What to modify:} The \code{lod_multiplier} and
#' \code{loq_multiplier} arguments if your laboratory SOP requires different
#' multipliers. The default values of 3 and 10 are the international
#' standard (ICH Q2(R1)) and should not be changed without documented
#' scientific justification.
#'
#' @examples
#' \dontrun{
#' lod_loq <- calculate_lod_loq(all_data)
#'
#' # Inspect LOD/LOQ values for a specific batch
#' lod_loq %>%
#'   filter(Batch_ID == "Batch_230615_XX") %>%
#'   select(Molecule_Prefix, Unnorm_LOD, Unnorm_LOQ, Norm_LOD, Norm_LOQ)
#'
#' # Use non-standard multipliers (e.g., some EPA methods use 3.3 and 10)
#' lod_loq <- calculate_lod_loq(all_data, lod_multiplier = 3.3, loq_multiplier = 10)
#' }
#'
#' @export
calculate_lod_loq <- function(all_batches_subtracted,
                              lod_multiplier = 3,
                              loq_multiplier = 10) {

  # Validate multipliers - must be positive numbers, and LOQ must exceed LOD
  if (!is.numeric(lod_multiplier) || lod_multiplier <= 0) {
    stop("lod_multiplier must be a positive number. Received: ", lod_multiplier)
  }
  if (!is.numeric(loq_multiplier) || loq_multiplier <= 0) {
    stop("loq_multiplier must be a positive number. Received: ", loq_multiplier)
  }
  if (loq_multiplier <= lod_multiplier) {
    warning("loq_multiplier (", loq_multiplier, ") is not greater than ",
            "lod_multiplier (", lod_multiplier, "). ",
            "LOQ should always exceed LOD by convention.")
  }

  # -- Step A: Select blank replicates and back-calculate ug/g ---------------
  #
  # The calibration equation is applied to each blank replicate individually:
  #   ug/g = (Signal - Intercept) / (Unknown_Weight x Slope)
  #
  # The intercept IS subtracted here because this is a point estimate of
  # where each blank falls on the calibration curve. The SD of these
  # point estimates is taken in Step B - and because SD is invariant to
  # constant additive shifts, the intercept cancels in the final LOD/LOQ.
  # It is nonetheless correct to include it in the per-replicate step.
  #
  # Rows where Slope is NA/0 or Unknown_Weight is NA/0 receive NA for the
  # ug/g back-calculation, which propagates to NA SD, then to 0 after
  # replace_na() - producing no threshold line in plots (by design).
  blank_concs <- all_batches_subtracted %>%
    dplyr::select(
      Batch_ID, Replicate_Name, Molecule_Prefix,
      Calibration_Factor, Replicate_Type, Unknown_Weight,
      Area, Normalized_Area,
      Norm_Slope, Unnorm_Slope, Norm_Y_Int, Unnorm_Y_Int
    ) %>%
    dplyr::filter(Replicate_Type == "Blank") %>%
    dplyr::mutate(

      # Unnormalized blank ug/g: uses raw Area signal
      Blank_Unnorm_ug_per_g = dplyr::if_else(
        is.na(Unnorm_Slope)    | Unnorm_Slope    == 0 |
          is.na(Unknown_Weight)  | Unknown_Weight  == 0,
        NA_real_,
        (Area - Unnorm_Y_Int) / (Unknown_Weight * Unnorm_Slope)
      ),

      # Normalized blank ug/g: uses IS-normalized Area signal
      Blank_Norm_ug_per_g = dplyr::if_else(
        is.na(Norm_Slope)      | Norm_Slope      == 0 |
          is.na(Unknown_Weight)  | Unknown_Weight  == 0,
        NA_real_,
        (Normalized_Area - Norm_Y_Int) / (Unknown_Weight * Norm_Slope)
      )
    )

  # Warn if no blank replicates are found at all
  if (nrow(blank_concs) == 0) {
    warning("No Blank replicates found in the data. ",
            "LOD and LOQ will be NA for all analytes. ",
            "Check Replicate_Type assignments in Batch_Info.")
    return(
      all_batches_subtracted %>%
        dplyr::distinct(Batch_ID, Molecule_Prefix) %>%
        dplyr::mutate(
          Blank_SD_Unnorm_ug_per_g = NA_real_,
          Blank_SD_Norm_ug_per_g   = NA_real_,
          Unnorm_LOD = NA_real_, Unnorm_LOQ = NA_real_,
          Norm_LOD   = NA_real_, Norm_LOQ   = NA_real_
        )
    )
  }

  # -- Step B: Compute SD of blank ug/g values per Batch_ID x Molecule_Prefix
  #
  # SD is the spread of blank concentrations - the basis for LOD and LOQ.
  # na.rm = TRUE: blanks with missing signal (e.g., below instrument noise)
  # are excluded from the SD calculation rather than propagating NA.
  # If all blanks are NA for a group, sd() returns NA -> replaced with 0.
  #
  # Also carry forward the calibration model parameters using first() since
  # they are identical for all rows within a Batch_ID x Molecule_Prefix group.
  lod_loq_table <- blank_concs %>%
    dplyr::group_by(Batch_ID, Molecule_Prefix) %>%
    dplyr::summarise(
      Blank_SD_Unnorm_ug_per_g = sd(Blank_Unnorm_ug_per_g, na.rm = TRUE),
      Blank_SD_Norm_ug_per_g   = sd(Blank_Norm_ug_per_g,   na.rm = TRUE),
      # Carry forward for reference/diagnostics
      Norm_Slope         = dplyr::first(Norm_Slope),
      Unnorm_Slope       = dplyr::first(Unnorm_Slope),
      Norm_Y_Int         = dplyr::first(Norm_Y_Int),
      Unnorm_Y_Int       = dplyr::first(Unnorm_Y_Int),
      Calibration_Factor = dplyr::first(Calibration_Factor),
      .groups = "drop"
    ) %>%
    # Replace NA SD (from missing or single-replicate blanks) with 0.
    # A SD of 0 produces LOD = LOQ = 0, and the plotting functions will
    # not draw threshold lines at 0 (filtered by x > 0).
    dplyr::mutate(
      Blank_SD_Unnorm_ug_per_g = tidyr::replace_na(Blank_SD_Unnorm_ug_per_g, 0),
      Blank_SD_Norm_ug_per_g   = tidyr::replace_na(Blank_SD_Norm_ug_per_g,   0)
    )

  # -- Step C: Warn about groups with zero or near-zero SD -------------------
  #
  # Zero SD means all blanks had identical (or all-missing) ug/g values.
  # This is unusual and may indicate a data issue worth investigating.
  zero_sd <- lod_loq_table %>%
    dplyr::filter(Blank_SD_Unnorm_ug_per_g == 0 | Blank_SD_Norm_ug_per_g == 0)

  if (nrow(zero_sd) > 0) {
    message("Note: ", nrow(zero_sd), " Batch_ID x Molecule_Prefix group(s) ",
            "have zero blank SD (possibly only one blank replicate or all ",
            "blanks below detection). No LOD/LOQ line will be drawn for these:")
    print(zero_sd %>% dplyr::select(Batch_ID, Molecule_Prefix,
                                    Blank_SD_Unnorm_ug_per_g,
                                    Blank_SD_Norm_ug_per_g))
  }

  # -- Step D: Apply LOD and LOQ multipliers ---------------------------------
  #
  # LOD = lod_multiplier x SD(blank ug/g)
  # LOQ = loq_multiplier x SD(blank ug/g)
  #
  # These are already in ug/g units - no further conversion or factor
  # multiplication is needed. The values are used directly as horizontal
  # (bar plots) or vertical (calibration curves) threshold lines.
  lod_loq_table <- lod_loq_table %>%
    dplyr::mutate(
      Unnorm_LOD = lod_multiplier * Blank_SD_Unnorm_ug_per_g,
      Unnorm_LOQ = loq_multiplier * Blank_SD_Unnorm_ug_per_g,
      Norm_LOD   = lod_multiplier * Blank_SD_Norm_ug_per_g,
      Norm_LOQ   = loq_multiplier * Blank_SD_Norm_ug_per_g
    )

  message(
    "LOD/LOQ calculated for ", nrow(lod_loq_table), " Batch_ID x ",
    "Molecule_Prefix combination(s). ",
    "Multipliers used - LOD: ", lod_multiplier, "x, LOQ: ", loq_multiplier, "x."
  )

  lod_loq_table
}
