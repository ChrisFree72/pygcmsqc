# R/calibration.R
#
# Functions for weighted 1/x linear calibration curve fitting, automated
# calibrator optimization by R² gain, and back-calculation of analyte
# concentrations and ug/g values.
#
# Background on the calibration approach:
#   py-GC-MS data is heteroscedastic — signal variance increases with
#   concentration. A standard unweighted OLS regression would give excessive
#   influence to high-concentration calibrators and underfit the low end of
#   the calibration range where environmental samples typically fall.
#   Weighted 1/x regression (weights = 1 / Scaled_Concentration) compensates
#   by giving proportionally more weight to low-concentration calibrators,
#   producing a better fit across the full dynamic range.
#
#   Reference: US EPA Method guidance and ICH Q2(R1) both endorse 1/x
#   weighting for analytical calibration with heteroscedastic signal response.
#
# Scaled_Concentration:
#   Each calibrator level (mass injected or known concentration) is multiplied
#   by a per-polymer, per-batch Calibration_Factor before fitting. This factor
#   accounts for differences in detector response between polymers (e.g., PE
#   vs. PET) and batch-to-batch instrument variation. Because the factor is
#   embedded in Scaled_Concentration, the fitted slope already reflects it —
#   no further Calibration_Factor multiplication is needed in back-calculation
#   or LOD/LOQ derivation.
#
# Contents:
#   Internal helpers (not exported):
#     scale_fun()                — computes Scaled_Concentration
#     fit_lm_w1x()               — fits a weighted 1/x linear model
#     safe_coef()                — safely extracts model coefficients
#     safe_r2()                  — safely extracts R²
#     optimize_batch()           — finds the optimal calibrator subset per batch
#
#   Exported pipeline functions:
#     run_calibration()          — runs optimization across all batches
#     calculate_concentrations() — back-calculates concentrations and ug/g


# ══ Internal helpers ══════════════════════════════════════════════════════════

# ── scale_fun ─────────────────────────────────────────────────────────────────
#
# Computes Scaled_Concentration = Calibrator_Level × Calibration_Factor.
#
# Calibration_Factor is a per-polymer, per-batch multiplier that converts
# the nominal calibrator level (e.g., ng injected) into the concentration
# space used for curve fitting. It accounts for polymer-specific detector
# response differences and is sourced from the Calibrator_Factors sheet.
#
# @param level  Numeric. The raw calibrator level (e.g., ng injected).
# @param factor Numeric. The Calibration_Factor for this polymer and batch.
# @return Numeric. The scaled concentration used as the x-axis in calibration.
#
#' @keywords internal
scale_fun <- function(level, factor) level * factor


# ── fit_lm_w1x ────────────────────────────────────────────────────────────────
#
# Fits a weighted 1/x linear regression of signal on Scaled_Concentration.
#
# Model: signal = slope × Scaled_Concentration + intercept
# Weights: 1 / Scaled_Concentration (gives more influence to low-concentration
#   calibrators where most environmental samples fall, and compensates for the
#   proportional increase in signal variance at higher concentrations).
#
# Returns NULL (rather than erroring) if the data are insufficient to fit a
# model — this allows the optimizer to skip degenerate subsets gracefully.
#
# Minimum requirements:
#   - At least 2 non-NA data points
#   - At least 2 distinct Scaled_Concentration values (otherwise slope is
#     undefined)
#   - Scaled_Concentration > 0 (required for 1/x weights to be finite)
#
# @param data  A data frame containing at least Scaled_Concentration and the
#              column named by y_col.
# @param y_col Character. Name of the signal column to model (either
#              "Normalized_Area_BG_Sub" or "Area_BG_Sub").
# @return An lm object, or NULL if fitting requirements are not met.
#
#' @keywords internal
fit_lm_w1x <- function(data, y_col) {
  # Filter to rows where both x and y are finite and x is strictly positive
  # (x <= 0 would produce non-finite weights and cause lm() to error)
  dd <- data[
    !is.na(data[[y_col]]) &
      is.finite(data$Scaled_Concentration) &
      data$Scaled_Concentration > 0,
  ]

  # Need at least 2 points with at least 2 distinct x values to define a line
  if (nrow(dd) < 2 || dplyr::n_distinct(dd$Scaled_Concentration) < 2) {
    return(NULL)
  }

  lm(
    as.formula(paste(y_col, "~ Scaled_Concentration")),
    data    = dd,
    weights = 1 / dd$Scaled_Concentration
  )
}


# ── safe_coef ─────────────────────────────────────────────────────────────────
#
# Extracts a named coefficient from a fitted lm object without erroring if
# the model is NULL or the coefficient name is absent. Returns NA instead.
# Used to safely extract slope ("Scaled_Concentration") and intercept
# ("(Intercept)") from potentially NULL model objects.
#
# @param model An lm object or NULL.
# @param name  Character. The coefficient name to extract.
# @return Numeric scalar, or NA_real_ if model is NULL or name not found.
#
#' @keywords internal
safe_coef <- function(model, name) {
  if (is.null(model)) return(NA_real_)
  co <- coef(model)
  if (!name %in% names(co)) return(NA_real_)
  unname(co[[name]])
}


# ── safe_r2 ───────────────────────────────────────────────────────────────────
#
# Extracts R² from a fitted lm object without erroring if the model is NULL.
# R² is used as the optimization criterion — higher R² indicates a better
# fit of the calibration curve to the selected calibrator subset.
#
# @param model An lm object or NULL.
# @return Numeric scalar R² in [0, 1], or NA_real_ if model is NULL.
#
#' @keywords internal
safe_r2 <- function(model) {
  if (is.null(model)) NA_real_ else summary(model)$r.squared
}


# ── optimize_batch ────────────────────────────────────────────────────────────
#
# Finds the optimal subset of calibrators for a single batch by testing all
# valid removal combinations and selecting the one that maximizes average R²
# gain across all analytes and both signal types (normalized and unnormalized).
#
# Algorithm:
#   1. The lowest-concentration calibrator is always retained as an anchor
#      point. Removing the lowest calibrator would leave the regression
#      unconstrained near zero and risk extrapolating into negative
#      concentration space.
#   2. All combinations of up to max_remove calibrators (from the non-anchor
#      set) are tested. Only subsets leaving >= 5 calibrators total are
#      considered (fewer than 5 points is insufficient for a reliable weighted
#      regression).
#   3. For each removal set, the R² gain vs. the full calibrator set is
#      computed for both normalized and unnormalized fits, averaged across
#      all analytes in the batch.
#   4. The removal set with the highest average R² gain is selected.
#      If no removal improves R² (gain <= 0), the full calibrator set is kept.
#
# @param batch_df   A data frame of standards for a single batch, containing
#                   Scaled_Concentration, Molecule_Prefix, Calibrator_Level,
#                   Replicate_Name, Replicate_Type, Normalized_Area_BG_Sub,
#                   Area_BG_Sub, and Batch_ID.
# @param max_remove Integer. Maximum number of non-anchor calibrators the
#                   optimizer may consider removing. Default 2. Runtime grows
#                   combinatorially with this value — increase cautiously.
# @return A data frame (one row per analyte) containing model parameters
#         (slopes, intercepts), R² values, and metadata about which
#         calibrators were removed and the R² gain achieved.
#
#' @keywords internal
optimize_batch <- function(batch_df, max_remove = 2) {

  # Compute Scaled_Concentration if not already present
  # (it should be pre-computed by run_calibration, but this guards against
  # direct calls to optimize_batch() in isolation)
  if (!"Scaled_Concentration" %in% names(batch_df)) {
    batch_df <- batch_df %>%
      dplyr::mutate(
        Scaled_Concentration = scale_fun(Calibrator_Level, Calibration_Factor)
      )
  }

  # Sort by Scaled_Concentration ascending so the minimum is always first
  sorted <- batch_df %>% dplyr::arrange(Scaled_Concentration)

  # Identify the lowest-concentration calibrator — this is the anchor point
  # that is always retained regardless of the removal set being tested
  lowest_level <- sorted$Calibrator_Level[
    which.min(sorted$Scaled_Concentration)
  ]

  # All other unique calibrator levels are candidates for removal testing
  candidates <- sorted$Calibrator_Level[
    sorted$Calibrator_Level != lowest_level
  ]

  # Build the list of removal sets to test:
  #   NULL = no removal (full calibrator set, always tested as baseline)
  #   Then all combinations of 1..max_remove candidates, subject to the
  #   constraint that at least 5 calibrators remain after removal.
  removal_sets <- list(NULL)
  for (k in seq_len(min(max_remove, length(candidates)))) {
    combs <- combn(candidates, k, simplify = FALSE)
    for (set in combs) {
      if ((nrow(sorted) - length(set)) >= 5) {
        removal_sets <- append(removal_sets, list(set))
      }
    }
  }

  # ── Evaluate each removal set ──────────────────────────────────────────────
  results <- purrr::map(removal_sets, function(rem_set) {

    # Build the test calibrator data frame for this removal set
    test_df <- if (is.null(rem_set)) {
      sorted
    } else {
      sorted %>% dplyr::filter(!Calibrator_Level %in% rem_set)
    }

    # Split by analyte and fit models for each
    per_prefix <- batch_df %>%
      dplyr::group_by(Molecule_Prefix) %>%
      dplyr::group_split()

    metrics <- purrr::map_dfr(per_prefix, function(df) {
      # Original fit: full calibrator set (baseline R²)
      orig_norm   <- fit_lm_w1x(df, "Normalized_Area_BG_Sub")
      orig_unnorm <- fit_lm_w1x(df, "Area_BG_Sub")

      # New fit: restricted to calibrators in test_df for this analyte
      new_df      <- df %>%
        dplyr::filter(Calibrator_Level %in% test_df$Calibrator_Level)
      new_norm    <- fit_lm_w1x(new_df, "Normalized_Area_BG_Sub")
      new_unnorm  <- fit_lm_w1x(new_df, "Area_BG_Sub")

      tibble::tibble(
        Molecule_Prefix    = df$Molecule_Prefix[1],
        Original_R2_Norm   = safe_r2(orig_norm),
        Original_R2_Unnorm = safe_r2(orig_unnorm),
        New_R2_Norm        = safe_r2(new_norm),
        New_R2_Unnorm      = safe_r2(new_unnorm),
        # Slope and intercept from the NEW (potentially reduced) fit.
        # These are the parameters used for concentration back-calculation.
        Norm_Slope         = safe_coef(new_norm,   "Scaled_Concentration"),
        Unnorm_Slope       = safe_coef(new_unnorm, "Scaled_Concentration"),
        Norm_Y_Int         = safe_coef(new_norm,   "(Intercept)"),
        Unnorm_Y_Int       = safe_coef(new_unnorm, "(Intercept)")
      )
    })

    # Average R² gain across all analytes and both signal types.
    # This single scalar is the optimization criterion — the removal set
    # with the highest avg_gain is selected as the best calibration.
    avg_gain <- mean(
      ((metrics$New_R2_Norm + metrics$New_R2_Unnorm) / 2) -
        ((metrics$Original_R2_Norm + metrics$Original_R2_Unnorm) / 2),
      na.rm = TRUE
    )

    list(removal = rem_set, avg_gain = avg_gain, metrics = metrics)
  })

  # ── Select best removal set ────────────────────────────────────────────────
  #
  # If no removal improves R² (all gains <= 0), fall back to the full
  # calibrator set (results[[1]], which is always the NULL/no-removal case).
  best <- results[[which.max(purrr::map_dbl(results, "avg_gain"))]]
  if (best$avg_gain <= 0) best <- results[[1]]

  # Resolve the names of removed calibrator replicates for reporting
  removed_names <- character(0)
  if (!is.null(best$removal) && length(best$removal) > 0) {
    removed_names <- batch_df %>%
      dplyr::filter(
        Replicate_Type == "Standard",
        Calibrator_Level %in% best$removal
      ) %>%
      dplyr::distinct(Replicate_Name) %>%
      dplyr::arrange(Replicate_Name) %>%
      dplyr::pull(Replicate_Name)
  }

  # Annotate the metrics table with batch-level metadata
  best$metrics %>%
    dplyr::mutate(
      Batch_ID            = batch_df$Batch_ID[1],
      Removed_Calibrators = ifelse(
        is.null(best$removal), "None",
        paste(removed_names, collapse = ",")
      ),
      Avg_R2_Gain_Batch   = best$avg_gain,
      Avg_R2_Gain_Rounded = round(best$avg_gain, 3),
      Improvement_Flag    = best$avg_gain > 0
    )
}


# ══ Exported pipeline functions ═══════════════════════════════════════════════

# ── run_calibration ───────────────────────────────────────────────────────────
#
#' Fit calibration curves and optimize calibrator selection across all batches
#'
#' For each batch, fits weighted 1/x linear calibration curves for every
#' analyte (Molecule_Prefix) using both normalized and unnormalized
#' background-subtracted signal. Runs the calibrator optimization algorithm
#' to find the subset of calibrators that maximizes average R² across all
#' analytes, subject to the constraint that the lowest-concentration
#' calibrator is always retained and at least 5 calibrators remain.
#'
#' Prints an estimated execution time before running and reports actual
#' elapsed time on completion.
#'
#' @param all_batches_subtracted A data frame as returned by
#'   \code{subtract_blanks()}, containing columns \code{Batch_ID},
#'   \code{Molecule_Prefix}, \code{Replicate_Type}, \code{Calibrator_Level},
#'   \code{Calibration_Factor}, \code{Normalized_Area_BG_Sub}, and
#'   \code{Area_BG_Sub}.
#' @param max_remove Integer. Maximum number of non-anchor calibrators the
#'   optimizer may remove per batch. Default 2. Runtime grows combinatorially
#'   — increasing to 3 or more is not recommended without first checking the
#'   estimated execution time.
#'
#' @return A data frame of calibration model parameters, one row per
#'   \code{Batch_ID} × \code{Molecule_Prefix} combination, containing:
#' \describe{
#'   \item{Batch_ID}{Batch identifier.}
#'   \item{Molecule_Prefix}{Analyte polymer prefix.}
#'   \item{Original_R2_Norm / Original_R2_Unnorm}{R² of the full-set fit.}
#'   \item{New_R2_Norm / New_R2_Unnorm}{R² of the optimized fit.}
#'   \item{Avg_R2_Gain_Batch}{Mean R² improvement from optimization.}
#'   \item{Improvement_Flag}{TRUE if optimization improved R².}
#'   \item{Removed_Calibrators}{Comma-separated names of removed replicates,
#'         or "None".}
#'   \item{Norm_Slope / Unnorm_Slope}{Fitted slopes (signal / concentration).}
#'   \item{Norm_Y_Int / Unnorm_Y_Int}{Fitted intercepts.}
#' }
#'
#' @details
#' \strong{What to modify:} Only the \code{max_remove} argument. The
#' optimization algorithm, weighting scheme, and minimum calibrator count
#' (5) are analytically fixed.
#'
#' \strong{Execution time:} The function times one batch as a sample and
#' extrapolates. Actual time may differ based on the number of analytes,
#' calibrator levels, and system speed.
#'
#' @examples
#' \dontrun{
#' models <- run_calibration(all_data, max_remove = 2)
#' # Inspect R² results
#' models %>% select(Batch_ID, Molecule_Prefix, New_R2_Norm, New_R2_Unnorm)
#' }
#'
#' @export
run_calibration <- function(all_batches_subtracted, max_remove = 2) {

  # Pre-compute Scaled_Concentration for all rows.
  # This is done once here rather than inside optimize_batch() to avoid
  # redundant computation across the many optimize_batch() calls.
  all_batches_subtracted <- all_batches_subtracted %>%
    dplyr::mutate(
      Scaled_Concentration = scale_fun(Calibrator_Level, Calibration_Factor)
    )

  # Filter to standard replicates only — only standards have known
  # Scaled_Concentration values and are used for curve fitting.
  standards <- all_batches_subtracted %>%
    dplyr::filter(Replicate_Type == "Standard")

  if (nrow(standards) == 0) {
    stop("No Standard replicates found in the data. ",
         "Check that Replicate_Type is correctly assigned in Batch_Info.")
  }

  batch_splits <- standards %>%
    dplyr::group_by(Batch_ID) %>%
    dplyr::group_split()

  # ── Estimate execution time before running ─────────────────────────────────
  #
  # Times a single batch as a representative sample and multiplies by the
  # total number of batches to give an upfront estimate. This lets the user
  # abort early if the runtime is unexpectedly long before committing to the
  # full run.
  sample_batch <- batch_splits[[1]]
  t0           <- Sys.time()
  optimize_batch(sample_batch, max_remove = max_remove)
  t1           <- Sys.time()

  sample_secs  <- as.numeric(t1 - t0, units = "secs")
  est_secs     <- sample_secs * length(batch_splits)

  cat("Estimated execution time:", round(est_secs, 2), "seconds (",
      round(est_secs / 60, 2), "minutes)\n")
  cat("Running calibration optimization across", length(batch_splits),
      "batch(es)...\n")

  if (est_secs > 600) {
    warning("Estimated execution time exceeds 10 minutes. ",
            "Consider reducing max_remove or dataset size.")
  }

  start_time <- Sys.time()

  # ── Run optimization across all batches ────────────────────────────────────
  final_results <- purrr::map_dfr(batch_splits, function(batch_df) {
    cat("  Processing Batch_ID:", batch_df$Batch_ID[1], "\n")
    optimize_batch(batch_df, max_remove = max_remove)
  })

  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
  cat("Calibration complete. Elapsed:", round(elapsed, 2), "seconds\n")

  # ── Reorder columns for readability in exported Excel ─────────────────────
  final_results %>%
    dplyr::select(
      Batch_ID, Molecule_Prefix,
      Original_R2_Norm, Original_R2_Unnorm,
      New_R2_Norm,      New_R2_Unnorm,
      Avg_R2_Gain_Batch, Avg_R2_Gain_Rounded, Improvement_Flag,
      Removed_Calibrators,
      Norm_Slope, Unnorm_Slope, Norm_Y_Int, Unnorm_Y_Int
    )
}


# ── calculate_concentrations ──────────────────────────────────────────────────
#
#' Back-calculate analyte concentrations and ug/g values from calibration models
#'
#' Joins the calibration model parameters (slopes and intercepts) from
#' \code{run_calibration()} back onto every row of the full dataset, then
#' applies the inverse calibration equation to compute concentrations and
#' sample-weight-normalized ug/g values for all replicates.
#'
#' Also joins R² values and removed calibrator information for use in
#' diagnostic plots, and adds a \code{Weighting_Type} column for plot
#' annotations.
#'
#' @param all_batches_subtracted A data frame as returned by
#'   \code{subtract_blanks()}.
#' @param model_params A data frame as returned by \code{run_calibration()},
#'   containing \code{Batch_ID}, \code{Molecule_Prefix}, \code{Norm_Slope},
#'   \code{Unnorm_Slope}, \code{Norm_Y_Int}, \code{Unnorm_Y_Int},
#'   \code{New_R2_Norm}, \code{New_R2_Unnorm}, and
#'   \code{Removed_Calibrators}.
#'
#' @return The input data frame with the following new columns:
#' \describe{
#'   \item{Scaled_Concentration}{Calibrator_Level × Calibration_Factor.
#'         Non-NA only for Standard replicates.}
#'   \item{Norm_Calc_Concentration}{Back-calculated concentration from
#'         normalized signal: (Normalized_Area_BG_Sub - Norm_Y_Int) /
#'         Norm_Slope.}
#'   \item{Unnorm_Calc_Concentration}{Back-calculated concentration from
#'         unnormalized signal: (Area_BG_Sub - Unnorm_Y_Int) / Unnorm_Slope.}
#'   \item{Norm_ug_Per_Gram}{Normalized concentration normalized to sample
#'         weight: (Normalized_Area_BG_Sub - Norm_Y_Int) /
#'         (Unknown_Weight × Norm_Slope).}
#'   \item{Unnorm_ug_Per_Gram}{Unnormalized concentration normalized to
#'         sample weight: (Area_BG_Sub - Unnorm_Y_Int) /
#'         (Unknown_Weight × Unnorm_Slope).}
#'   \item{Weighting_Type}{Character. Always "1/x", used for plot subtitles.}
#'   \item{New_R2_Norm / New_R2_Unnorm}{Joined from model_params for plots.}
#'   \item{Removed_Calibrators}{Joined from model_params for plots.}
#' }
#'
#' @details
#' \strong{Concentration back-calculation:}
#' The calibration model is: Signal = Slope × Scaled_Concentration + Intercept
#' Rearranged: Concentration = (Signal - Intercept) / Slope
#' Normalized to sample weight: ug/g = (Signal - Intercept) /
#'   (Unknown_Weight × Slope)
#'
#' \strong{Calibration_Factor and Scaled_Concentration:}
#' Calibration_Factor is already embedded in Scaled_Concentration, which was
#' used as the x-axis when fitting the model. Therefore the fitted Slope
#' already accounts for Calibration_Factor — no further multiplication is
#' needed here. This is a common source of error and is intentionally
#' documented.
#'
#' \strong{NA propagation:} Rows where Slope is NA or zero (analytes with
#' no calibration data) receive NA for all concentration columns rather
#' than producing Inf or NaN.
#'
#' \strong{Negative concentrations:} Values below the blank level will
#' produce negative concentrations. These are retained as analytically
#' valid results indicating the sample is below the blank level. They are
#' clipped to zero only in export and PCA functions.
#'
#' \strong{What to modify:} Nothing. The back-calculation equations are
#' fixed by the calibration model. To change which replicates get
#' concentrations calculated, modify the upstream Replicate_Type assignments
#' in the Batch_Info sheet.
#'
#' @examples
#' \dontrun{
#' models   <- run_calibration(all_data, max_remove = 2)
#' all_data <- calculate_concentrations(all_data, models)
#'
#' # Check Unknown sample concentrations
#' all_data %>%
#'   filter(Replicate_Type == "Unknown") %>%
#'   select(Batch_ID, Replicate_Name, Molecule_Prefix,
#'          Norm_ug_Per_Gram, Unnorm_ug_Per_Gram)
#' }
#'
#' @export
calculate_concentrations <- function(all_batches_subtracted, model_params) {

  # ── Step A: Pre-compute Scaled_Concentration ───────────────────────────────
  #
  # Scaled_Concentration = Calibrator_Level × Calibration_Factor.
  # This is the x-axis value used during calibration curve fitting.
  # For non-standard replicates, Calibrator_Level is NA, so
  # Scaled_Concentration will also be NA — this is expected and correct.
  # The resulting column is used by plotting functions to position standards
  # on the calibration curve x-axis.
  all_batches_subtracted <- all_batches_subtracted %>%
    dplyr::mutate(
      Scaled_Concentration = scale_fun(Calibrator_Level, Calibration_Factor)
    )

  # ── Step B: Join model parameters ─────────────────────────────────────────
  #
  # Join slopes, intercepts, R² values, and removed calibrator info by
  # (Batch_ID, Molecule_Prefix). Every row for a given analyte in a given
  # batch receives the same model parameters — this is correct because the
  # calibration model is defined at the batch × analyte level.
  all_batches_subtracted <- all_batches_subtracted %>%
    dplyr::left_join(
      model_params %>%
        dplyr::select(
          Batch_ID, Molecule_Prefix,
          Norm_Slope, Unnorm_Slope,
          Norm_Y_Int, Unnorm_Y_Int,
          New_R2_Norm, New_R2_Unnorm,
          Removed_Calibrators
        ),
      by = c("Batch_ID", "Molecule_Prefix")
    )

  # ── Step C: Back-calculate concentrations and ug/g ────────────────────────
  #
  # Inverse of the calibration model:
  #   Model:  Signal = Slope × Concentration + Intercept
  #   Inverse: Concentration = (Signal - Intercept) / Slope
  #   ug/g:    ug/g = (Signal - Intercept) / (Unknown_Weight × Slope)
  #
  # IMPORTANT: Calibration_Factor is NOT applied here. It is already
  # embedded in the fitted Slope via Scaled_Concentration. Applying it
  # again would double-count the factor and inflate concentrations.
  #
  # Rows where Slope is NA or zero receive NA to avoid Inf/NaN values.
  # Unknown_Weight = 0 also produces NA (can't divide by zero weight).
  all_batches_subtracted <- all_batches_subtracted %>%
    dplyr::mutate(

      # Concentration in the same units as Scaled_Concentration
      # (useful for plotting against the calibration curve x-axis)
      Norm_Calc_Concentration = dplyr::if_else(
        is.na(Norm_Slope) | Norm_Slope == 0,
        NA_real_,
        (Normalized_Area_BG_Sub - Norm_Y_Int) / Norm_Slope
      ),

      Unnorm_Calc_Concentration = dplyr::if_else(
        is.na(Unnorm_Slope) | Unnorm_Slope == 0,
        NA_real_,
        (Area_BG_Sub - Unnorm_Y_Int) / Unnorm_Slope
      ),

      # ug per gram of sample — the primary reportable quantity.
      # Dividing by Unknown_Weight converts from total concentration to
      # concentration per unit sample mass.
      Norm_ug_Per_Gram = dplyr::if_else(
        is.na(Norm_Slope) | Norm_Slope == 0 |
          is.na(Unknown_Weight) | Unknown_Weight == 0,
        NA_real_,
        (Normalized_Area_BG_Sub - Norm_Y_Int) / (Unknown_Weight * Norm_Slope)
      ),

      Unnorm_ug_Per_Gram = dplyr::if_else(
        is.na(Unnorm_Slope) | Unnorm_Slope == 0 |
          is.na(Unknown_Weight) | Unknown_Weight == 0,
        NA_real_,
        (Area_BG_Sub - Unnorm_Y_Int) / (Unknown_Weight * Unnorm_Slope)
      ),

      # Weighting type label for calibration curve plot subtitles
      Weighting_Type = "1/x"
    )

  message("Concentrations calculated. New columns: ",
          "Scaled_Concentration, Norm_Calc_Concentration, ",
          "Unnorm_Calc_Concentration, Norm_ug_Per_Gram, Unnorm_ug_Per_Gram.")

  all_batches_subtracted
}
