# R/data_processing.R
#
# Functions for joining calibration metadata, summarizing signal areas by
# molecule prefix, and performing blank background subtraction.
#
# This file covers everything between raw harmonized data (output of
# data_loading.R) and the calibration step (calibration.R). The output
# of subtract_blanks() is the primary data frame that flows into all
# downstream calibration, LOD/LOQ, plotting, and export functions.
#
# Pipeline order:
#   harmonize_batches()           [data_loading.R]
#        v
#   join_calibration_metadata()   [this file]
#        v
#   summarize_batches()           [this file]
#        v
#   subtract_blanks()             [this file]
#        v
#   run_calibration()             [calibration.R]
#
# Contents:
#   1. join_calibration_metadata()  - joins calibration factors, weights,
#                                     replicate types, and sample weights
#   2. summarize_batches()          - sums Area/Normalized_Area by molecule
#                                     prefix; drops excluded molecules
#   3. subtract_blanks()            - computes and subtracts blank averages


# -- 1. join_calibration_metadata ---------------------------------------------
#
#' Join all calibration reference metadata onto the harmonized batch data
#'
#' Performs four sequential join/annotation operations on the harmonized
#' batch data:
#' \enumerate{
#'   \item Extracts a plastic polymer token from \code{Molecule_Name} and
#'         joins the per-batch, per-polymer \code{Multiplication_Factor}
#'         from \code{Calibrator_Factors} as \code{Calibration_Factor}.
#'   \item Matches \code{Replicate_Name} against \code{Calibrator_Name}
#'         entries to assign a \code{Calibrator_Level} (weight) to each
#'         standard replicate.
#'   \item Joins \code{Replicate_Type} (Unknown, QC, Standard, Blank,
#'         Cup_Blank) from \code{Batch_Info}.
#' }
#'
#' @param all_batches A data frame as returned by \code{harmonize_batches()}.
#' @param calibrator_factors A data frame. The \code{Calibrator_Factors}
#'   sheet from \code{load_reference_data()}, containing columns
#'   \code{Batch_ID}, \code{Plastic}, and \code{Multiplication_Factor}.
#' @param calibrator_values A data frame. The \code{Calibrator_Values}
#'   sheet, containing columns \code{Batch_ID}, \code{Calibrator_Name},
#'   and \code{Calibrator_Weight}.
#' @param batch_info A data frame. The \code{Batch_Info} sheet, containing
#'   columns \code{Batch_ID}, \code{Replicates}, \code{Replicate_Type},
#'   and \code{Unknown_Weight}.
#'
#' @return The input data frame with four new columns appended:
#'   \code{Calibration_Factor}, \code{Calibrator_Level}, and
#'   \code{Replicate_Type}.
#'
#' @details
#' \strong{Plastic token matching:} Polymer names are extracted from
#' \code{Molecule_Name} using a longest-match-first regex built from the
#' unique \code{Plastic} values in \code{Calibrator_Factors}. Longest-first
#' ordering prevents shorter tokens (e.g., \code{N6}) from stealing matches
#' intended for longer ones (e.g., \code{N66}). Matching is case-insensitive
#' and word-boundary-anchored to avoid partial matches within longer strings.
#'
#' \strong{Calibrator level matching:} Within each batch, replicate names
#' are matched against calibrator names using a case-insensitive substring
#' search (\code{str_detect(..., fixed(...))}). When a replicate name matches
#' multiple calibrator names, the longest (most specific) match is used to
#' avoid ambiguous assignments.
#'
#' \strong{Diagnostics:} Duplicate entries in \code{Calibrator_Factors} and
#' unmatched rows are reported via \code{message()} so issues can be
#' identified and corrected in the reference workbook without stopping
#' execution.
#'
#' \strong{What to modify:} Nothing in this function - all configuration
#' comes from the reference Excel workbook. To add a polymer alias
#' (e.g., map "NYLON6" -> "N6"), add the alias directly to the
#' \code{Calibrator_Factors} sheet or uncomment and extend the alias
#' recode block inside this function.
#'
#' @examples
#' \dontrun{
#' ref      <- load_reference_data("data/reference.xlsx")
#' batches  <- load_batch_data("data/Batches/")
#' all_data <- harmonize_batches(batches, ref)
#' all_data <- join_calibration_metadata(
#'   all_data,
#'   ref$Calibrator_Factors,
#'   ref$Calibrator_Values,
#'   ref$Batch_Info
#' )
#' }
#'
#' @export
join_calibration_metadata <- function(all_batches,
                                      calibrator_factors,
                                      calibrator_values,
                                      batch_info) {

  # -- Step A: Build polymer token regex --------------------------------------
  #
  # Extract unique polymer names from Calibrator_Factors and sort them
  # longest-first. This prevents shorter tokens (N6) from matching before
  # longer ones (N66) when both appear in the regex alternation.
  plastics <- calibrator_factors %>%
    dplyr::transmute(Plastic = stringr::str_to_upper(stringr::str_squish(Plastic))) %>%
    dplyr::filter(!is.na(Plastic), Plastic != "") %>%
    dplyr::distinct() %>%
    dplyr::pull(Plastic)

  plastics_sorted <- plastics[order(nchar(plastics), decreasing = TRUE)]

  # str_escape() handles regex metacharacters in polymer names (e.g., +, -)
  plastics_esc <- stringr::str_escape(plastics_sorted)

  # Word-boundary anchors ((?<![A-Za-z0-9]) and (?![A-Za-z0-9])) prevent
  # "PE" from matching inside "HDPE" or "PET".
  token_pattern <- stringr::regex(
    stringr::str_c(
      "(?<![A-Za-z0-9])(",
      stringr::str_c(plastics_esc, collapse = "|"),
      ")(?![A-Za-z0-9])"
    ),
    ignore_case = TRUE
  )

  # -- Step B: Extract Plastic_Join from Molecule_Name ------------------------
  #
  # Plastic_Join is a temporary join key used to match rows in all_batches
  # to their corresponding Multiplication_Factor in Calibrator_Factors.
  # It is dropped after the join.
  all_batches_clean <- all_batches %>%
    dplyr::mutate(
      Batch_ID      = stringr::str_squish(Batch_ID),
      Molecule_Name = stringr::str_squish(Molecule_Name),
      Plastic_Join  = stringr::str_to_upper(
        stringr::str_extract(Molecule_Name, token_pattern)
      )
    )

  # Optional alias mapping - uncomment and extend if your Molecule_Name
  # strings use non-standard polymer abbreviations that don't match the
  # tokens in Calibrator_Factors exactly.
  # all_batches_clean <- all_batches_clean %>%
  #   dplyr::mutate(
  #     Plastic_Join = dplyr::recode(
  #       Plastic_Join,
  #       "NYLON6"   = "N6",
  #       "NYLON 6"  = "N6",
  #       "NYLON66"  = "N66",
  #       "NYLON 66" = "N66",
  #       .default   = Plastic_Join
  #     )
  #   )

  # -- Step C: Clean Calibrator_Factors and check for duplicates --------------
  #
  # Duplicates in (Batch_ID, Plastic_Join) would cause a many-to-many join
  # and silently multiply rows. We detect and report them before joining.
  cal_factors_clean <- calibrator_factors %>%
    dplyr::mutate(
      Batch_ID     = stringr::str_squish(Batch_ID),
      Plastic_Join = stringr::str_to_upper(stringr::str_squish(Plastic))
    ) %>%
    dplyr::select(Batch_ID, Plastic_Join, Multiplication_Factor)

  dups <- cal_factors_clean %>%
    dplyr::count(Batch_ID, Plastic_Join) %>%
    dplyr::filter(n > 1)

  if (nrow(dups) > 0) {
    message("Warning: duplicate Calibrator_Factors entries detected for ",
            "some (Batch_ID, Plastic_Join) combinations. ",
            "Resolve in the reference workbook to ensure correct joining:")
    print(dups)
  }

  # -- Step D: Join Calibration_Factor ----------------------------------------
  #
  # Left join preserves all rows in all_batches_clean; rows without a match
  # receive NA for Calibration_Factor, which will propagate through
  # Scaled_Concentration and flag those analytes as uncalibrated.
  tmp <- all_batches_clean %>%
    dplyr::left_join(cal_factors_clean, by = c("Batch_ID", "Plastic_Join"))

  # Ensure Calibration_Factor column exists even if join produced nothing
  if (!"Calibration_Factor" %in% names(tmp)) tmp[["Calibration_Factor"]] <- NA_real_

  # Coerce Multiplication_Factor to numeric in case it was read as character
  if (is.character(tmp$Multiplication_Factor)) {
    suppressWarnings({
      tmp$Multiplication_Factor <- as.numeric(tmp$Multiplication_Factor)
    })
  }

  # Fill Calibration_Factor from the joined factor; preserve any pre-existing
  # values where the join produced no match (coalesce = first non-NA wins)
  tmp[["Calibration_Factor"]] <- dplyr::coalesce(
    tmp$Multiplication_Factor,
    tmp[["Calibration_Factor"]]
  )

  # Drop temporary join columns
  all_batches <- tmp %>%
    dplyr::select(-Multiplication_Factor, -Plastic_Join)

  # Diagnostic: report rows with no calibration factor match so the user
  # can correct spelling/case issues in the reference workbook
  no_match <- all_batches_clean %>%
    dplyr::left_join(cal_factors_clean, by = c("Batch_ID", "Plastic_Join")) %>%
    dplyr::filter(is.na(Multiplication_Factor))

  if (nrow(no_match) > 0) {
    message("Rows with no Calibration_Factor match ",
            "(check Plastic spelling/case in Calibrator_Factors sheet):")
    print(
      no_match %>%
        dplyr::count(Batch_ID, Plastic_Join, name = "n") %>%
        dplyr::arrange(dplyr::desc(n))
    )
  }

  # -- Step E: Join Calibrator_Level from Calibrator_Values -------------------
  #
  # Each standard replicate has a known concentration level (weight).
  # We match Replicate_Name against Calibrator_Name using substring detection.
  # When multiple calibrators match a single replicate name, the longest
  # (most specific) calibrator name is chosen to avoid ambiguous assignments.
  # Non-standard replicates (Unknowns, QC, Blanks) receive NA.
  cal_values_clean <- calibrator_values %>%
    dplyr::mutate(
      Batch_ID        = stringr::str_squish(Batch_ID),
      Calibrator_Name = stringr::str_squish(Calibrator_Name)
    ) %>%
    dplyr::select(Batch_ID, Calibrator_Name, Calibrator_Weight)

  # Add a temporary row ID to allow safe re-joining after the candidate filter
  all_batches_tmp <- all_batches %>%
    dplyr::mutate(
      Batch_ID       = stringr::str_squish(Batch_ID),
      Replicate_Name = stringr::str_squish(Replicate_Name),
      .rowid         = dplyr::row_number()
    )

  # Find candidate matches: for each row, find all Calibrator_Names whose
  # string appears inside Replicate_Name (case-insensitive, fixed string
  # to avoid regex interpretation of calibrator name characters).
  # Then keep only the longest matching calibrator name per row.
  candidates <- all_batches_tmp %>%
    dplyr::select(.rowid, Batch_ID, Replicate_Name) %>%
    dplyr::inner_join(cal_values_clean, by = "Batch_ID") %>%
    dplyr::filter(
      stringr::str_detect(
        Replicate_Name,
        stringr::fixed(Calibrator_Name, ignore_case = TRUE)
      )
    ) %>%
    dplyr::mutate(name_len = nchar(Calibrator_Name)) %>%
    dplyr::arrange(.rowid, dplyr::desc(name_len)) %>%
    dplyr::group_by(.rowid) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(.rowid, Calibrator_Weight)

  # Re-join the winning calibrator weight back onto every row by .rowid.
  # Rows with no matching calibrator (Unknowns, QCs, Blanks) receive NA.
  all_batches <- all_batches_tmp %>%
    dplyr::left_join(candidates, by = ".rowid") %>%
    dplyr::mutate(Calibrator_Level = Calibrator_Weight) %>%
    dplyr::select(-.rowid, -Calibrator_Weight)

  # -- Step F: Join Replicate_Type and Unknown_Weight from Batch_Info ----------
  #
  # Replicate_Type classifies each replicate as one of:
  #   Unknown, QC, Standard, Blank, Cup_Blank, Internal Standard
  # Unknown_Weight is the sample mass (grams) used to compute ug/g.
  # Both are joined by (Batch_ID, Replicate_Name).
  batch_info_clean <- batch_info %>%
    dplyr::mutate(
      Batch_ID       = stringr::str_squish(Batch_ID),
      Replicates     = stringr::str_squish(Replicates),
      Replicate_Type = stringr::str_squish(Replicate_Type)
    ) %>%
    dplyr::select(Batch_ID, Replicates, Replicate_Type, Unknown_Weight) %>%
    dplyr::distinct()

  all_batches <- all_batches %>%
    dplyr::mutate(
      Batch_ID       = stringr::str_squish(Batch_ID),
      Replicate_Name = stringr::str_squish(Replicate_Name)
    ) %>%
    dplyr::left_join(
      batch_info_clean,
      by = c("Batch_ID", "Replicate_Name" = "Replicates")
    )

  message("Calibration metadata joined. Columns added: ",
          "Calibration_Factor, Calibrator_Level, Replicate_Type, Unknown_Weight.")

  all_batches
}


# -- 2. summarize_batches ------------------------------------------------------
#
#' Summarize signal areas by molecule prefix across fragment ions
#'
#' Extracts the polymer prefix from each \code{Molecule_Name} (the portion
#' before the first underscore, e.g., \code{"PE"} from \code{"PE_1"}),
#' optionally removes excluded molecules, then sums \code{Area} and
#' \code{Normalized_Area} within each unique combination of batch, replicate,
#' and molecule prefix. This collapses multiple fragment-level rows per
#' molecule into a single summarized row per analyte per replicate.
#'
#' Sample weights from \code{Batch_Info} are also joined here for use in
#' ug/g concentration calculations downstream.
#'
#' @param all_batches A data frame as returned by
#'   \code{join_calibration_metadata()}.
#' @param excluded_molecules A character vector of \code{Molecule_Name}
#'   values to remove before summarizing. Useful for dropping known
#'   problematic or redundant molecules (e.g., duplicate standards).
#'   Default is \code{c("N66_2", "PC_1", "PET_2", "PVC_1")}.
#'
#' @return A summarized data frame with one row per unique combination of
#'   \code{Batch_ID}, \code{Replicate_Name}, and \code{Molecule_Prefix},
#'   containing summed \code{Area} and \code{Normalized_Area} columns
#'   and all metadata columns preserved.
#'
#' @details
#' \strong{Prefix extraction:} If \code{Molecule_Name} contains an underscore,
#' the prefix is everything before the first underscore. If no underscore is
#' present, the full \code{Molecule_Name} is used as the prefix (handles
#' molecules like \code{"IS"} or \code{"PE"} without numeric suffixes).
#'
#' \strong{Robust summation:} The sum is computed with \code{na.rm = TRUE}.
#' If all values for a group are \code{NA}, the result is \code{NA} rather
#' than 0, preserving the distinction between a true zero signal and a
#' missing measurement.
#'
#' \strong{What to modify:} The \code{excluded_molecules} argument. Pass a
#' character vector of any \code{Molecule_Name} strings you want to drop.
#' This is the primary user configuration point for this function.
#'
#' @examples
#' \dontrun{
#' all_data <- summarize_batches(
#'   all_data,
#'   excluded_molecules = c("N66_2", "PC_1", "PET_2", "PVC_1")
#' )
#' }
#'
#' @export
summarize_batches <- function(all_batches,
                              excluded_molecules = c("N66_2", "PC_1",
                                                     "PET_2", "PVC_1")) {

  # -- Step A: Coerce Area columns to numeric ---------------------------------
  #
  # to_numeric_clean() (defined in utils.R) handles "#N/A", infinity, and other
  # Skyline export artifacts that would otherwise produce NA with a warning
  # or silently corrupt sums. Applied here rather than at load time because
  # type_convert() in harmonize_batches() may not catch all edge cases.
  all_batches <- all_batches %>%
    dplyr::mutate(
      Area            = to_numeric_clean(Area),
      Normalized_Area = to_numeric_clean(Normalized_Area)
    )

  # -- Step B: Remove excluded molecules -------------------------------------
  #
  # Excluded molecules are dropped before summarization so they do not
  # contribute to prefix-level sums. Common use case: redundant calibrators
  # or known interference peaks that should not be quantified.
  if (length(excluded_molecules) > 0) {
    n_before <- nrow(all_batches)
    all_batches <- all_batches %>%
      dplyr::filter(!Molecule_Name %in% excluded_molecules)
    n_dropped <- n_before - nrow(all_batches)
    if (n_dropped > 0) {
      message("Excluded ", n_dropped, " rows matching: ",
              paste(excluded_molecules, collapse = ", "))
    }
  }

  # -- Step C: Extract Molecule_Prefix ---------------------------------------
  #
  # The prefix groups all fragment-level measurements for a polymer together.
  # For example, "PE_1", "PE_2", "PE_3" all have prefix "PE" and will be
  # summed into a single row. Molecules without an underscore (e.g., "IS")
  # use the full name as the prefix.
  all_batches <- all_batches %>%
    dplyr::mutate(
      Molecule_Prefix = dplyr::case_when(
        stringr::str_detect(Molecule_Name, "_") ~
          stringr::str_extract(Molecule_Name, "^[^_]+"),
        TRUE ~ Molecule_Name
      )
    )

  # Define which columns to carry through vs. which to sum.
  # columns_to_keep: metadata that should be identical within a group
  #   (if they differ, the first value is used - see group_by below)
  # columns_to_sum: signal columns that are additive across fragment ions
  columns_to_keep <- c("Replicate_Name", "Calibration_Factor",
                       "Calibrator_Level", "Replicate_Type")
  columns_to_sum  <- c("Area", "Normalized_Area")

  # Select only the columns needed for summarization
  valid_columns <- c("Batch_ID", "Replicate_Name", "Molecule_Prefix",
                     columns_to_keep, columns_to_sum)

  all_batches_fil <- all_batches %>%
    dplyr::select(dplyr::any_of(valid_columns))

  # -- Step D: Sum Area and Normalized_Area within each group -----------------
  #
  # Groups: Batch_ID x Replicate_Name x Molecule_Prefix x metadata columns.
  # The robust sum returns NA only if ALL values in a group are NA;
  # otherwise it sums the non-NA values. This preserves the difference
  # between a true zero measurement and a completely missing one.
  all_batches_sum <- all_batches_fil %>%
    dplyr::group_by(
      dplyr::across(c("Batch_ID", "Replicate_Name", "Molecule_Prefix",
                      dplyr::all_of(columns_to_keep)))
    ) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(columns_to_sum),
        ~ {
          x        <- suppressWarnings(as.numeric(.x))
          n_non_na <- sum(!is.na(x))
          if (n_non_na == 0) NA_real_ else sum(x, na.rm = TRUE)
        }
      ),
      .groups = "drop"
    )

  # -- Step E: Join Unknown_Weight --------------------------------------------
  #
  # Unknown_Weight (sample mass in grams) is needed for ug/g calculations
  # in calibration.R. Joined here by (Batch_ID, Replicate_Name) from the
  # same Batch_Info used in join_calibration_metadata().
  # Note: this requires Batch_Info to still be accessible; it is passed
  # through all_batches as a joined column from Step F of
  # join_calibration_metadata(). If Unknown_Weight is already present,
  # this step is a no-op.
  if (!"Unknown_Weight" %in% names(all_batches_sum) &&
      "Unknown_Weight" %in% names(all_batches)) {

    sample_weights <- all_batches %>%
      dplyr::select(Batch_ID, Replicate_Name, Unknown_Weight) %>%
      dplyr::distinct()

    all_batches_sum <- all_batches_sum %>%
      dplyr::left_join(sample_weights,
                       by = c("Batch_ID", "Replicate_Name"))
  }

  message("Summarized to ", nrow(all_batches_sum), " rows across ",
          dplyr::n_distinct(all_batches_sum$Molecule_Prefix),
          " molecule prefix(es).")

  all_batches_sum
}


# -- 3. subtract_blanks --------------------------------------------------------
#
#' Subtract average blank signal from all non-blank, non-standard replicates
#'
#' Computes the mean \code{Area} and mean \code{Normalized_Area} of blank
#' replicates within each \code{Batch_ID} and \code{Molecule_Prefix} group,
#' then subtracts those averages from all other replicate types to produce
#' background-corrected signal columns (\code{Area_BG_Sub} and
#' \code{Normalized_Area_BG_Sub}).
#'
#' Standards, Cup_Blanks, and Internal Standard molecules are excluded from
#' background subtraction because:
#' \itemize{
#'   \item Standards are used to build the calibration curve and must retain
#'         their original signal values.
#'   \item Cup_Blanks represent procedural blanks and are not environmental
#'         samples.
#'   \item Internal Standards (IS prefix) are used for normalization and
#'         their absolute signal should not be blank-corrected.
#' }
#'
#' @param all_batches_sum A data frame as returned by
#'   \code{summarize_batches()}, containing columns \code{Batch_ID},
#'   \code{Molecule_Prefix}, \code{Replicate_Type}, \code{Area}, and
#'   \code{Normalized_Area}.
#'
#' @return The input data frame with four new columns:
#' \describe{
#'   \item{Area_Blank}{Mean blank Area for this Batch_ID x Molecule_Prefix.}
#'   \item{Normalized_Area_Blank}{Mean blank Normalized_Area for this group.}
#'   \item{Area_BG_Sub}{Background-subtracted Area. Equal to \code{Area} for
#'         Standards, Cup_Blanks, and IS molecules.}
#'   \item{Normalized_Area_BG_Sub}{Background-subtracted Normalized_Area.
#'         Equal to \code{Normalized_Area} for the same excluded types.}
#' }
#'
#' @details
#' \strong{Missing blank handling:} If no blank exists for a given
#' (Batch_ID, Molecule_Prefix) combination, the blank average is treated as
#' 0 (via \code{coalesce(..., 0)}), so the background-subtracted value equals
#' the raw signal. A message is printed listing any such combinations.
#'
#' \strong{Negative values:} Background subtraction can produce negative
#' values when a sample signal falls below the blank average. These are
#' analytically valid (indicating the sample is below the blank level) and
#' are preserved in \code{Area_BG_Sub} and \code{Normalized_Area_BG_Sub}.
#' They are clipped to zero only in export and PCA functions where
#' non-negative values are required.
#'
#' \strong{What to modify:} Nothing in this function. The excluded replicate
#' types (Standard, Cup_Blank) and the IS prefix exclusion are analytically
#' fixed and should not be changed without careful consideration of the
#' calibration model.
#'
#' @examples
#' \dontrun{
#' all_data <- subtract_blanks(all_data)
#' # Check subtracted values
#' all_data %>%
#'   filter(Replicate_Type == "Unknown") %>%
#'   select(Batch_ID, Replicate_Name, Molecule_Prefix,
#'          Area, Area_Blank, Area_BG_Sub)
#' }
#'
#' @export
subtract_blanks <- function(all_batches_sum) {

  # -- Step A: Compute average blank signal per Batch_ID x Molecule_Prefix ---
  #
  # Only rows classified as "Blank" contribute to the blank average.
  # na.rm = TRUE: if some blank replicates have NA signal (e.g., below
  # detection), the average is computed from the non-NA values only.
  blank_avg <- all_batches_sum %>%
    dplyr::filter(Replicate_Type == "Blank") %>%
    dplyr::group_by(Batch_ID, Molecule_Prefix) %>%
    dplyr::summarise(
      Area            = mean(Area,            na.rm = TRUE),
      Normalized_Area = mean(Normalized_Area, na.rm = TRUE),
      .groups         = "drop"
    )

  # -- Step B: Warn about analyte x batch combinations with no blank ----------
  all_combos <- all_batches_sum %>%
    dplyr::distinct(Batch_ID, Molecule_Prefix)

  missing_blanks <- dplyr::anti_join(all_combos, blank_avg,
                                     by = c("Batch_ID", "Molecule_Prefix"))

  if (nrow(missing_blanks) > 0) {
    message("No blank found for ", nrow(missing_blanks),
            " Batch_ID x Molecule_Prefix combination(s). ",
            "Blank will be treated as 0 for these groups:")
    print(missing_blanks)
  }

  # -- Step C: Join blank averages and compute background-subtracted columns --
  #
  # Suffix "_Blank" is added to the joined blank columns to distinguish them
  # from the original signal columns. coalesce(..., 0) handles the case where
  # no blank exists for a group (sets blank to 0 rather than NA).
  all_batches_subtracted <- all_batches_sum %>%
    dplyr::left_join(
      blank_avg,
      by     = c("Batch_ID", "Molecule_Prefix"),
      suffix = c("", "_Blank")
    ) %>%
    dplyr::mutate(
      # Replace missing blank averages with 0
      Area_Blank            = dplyr::coalesce(Area_Blank,            0),
      Normalized_Area_Blank = dplyr::coalesce(Normalized_Area_Blank, 0),

      # Background subtraction:
      # Applied to all replicate types EXCEPT Standards, Cup_Blanks, and IS.
      # Standards: must retain raw signal for calibration curve fitting.
      # Cup_Blanks: procedural blanks, not environmental samples.
      # IS (Internal Standard): normalization reference; absolute signal
      #   must not be blank-corrected.
      Area_BG_Sub = dplyr::if_else(
        !(Replicate_Type %in% c("Standard", "Cup_Blank")) &
          tolower(Molecule_Prefix) != "is",
        Area - Area_Blank,
        Area    # pass through unchanged for excluded types
      ),

      Normalized_Area_BG_Sub = dplyr::if_else(
        !(Replicate_Type %in% c("Standard", "Cup_Blank")) &
          tolower(Molecule_Prefix) != "is",
        Normalized_Area - Normalized_Area_Blank,
        Normalized_Area
      )
    )

  message("Blank subtraction complete. New columns: ",
          "Area_BG_Sub, Normalized_Area_BG_Sub.")

  all_batches_subtracted
}
