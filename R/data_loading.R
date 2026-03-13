# R/data_loading.R
#
# Functions for reading input files and preparing the raw data structure
# for the pygcmsqc pipeline. This is always the first step in the pipeline.
#
# The three functions here replace the monolithic file-reading block at the
# top of the original script. Crucially, they do NOT use list2env() — instead
# data is returned as named lists and passed explicitly between pipeline steps.
# This keeps the user's global environment clean and makes the pipeline
# reproducible and testable.
#
# Contents:
#   1. load_reference_data()   — reads the reference Excel workbook
#   2. load_batch_data()       — reads all batch CSVs from a directory
#   3. harmonize_batches()     — coerces types, binds rows, filters to targets


# ── 1. load_reference_data ────────────────────────────────────────────────────
#
#' Load the QC reference Excel workbook
#'
#' Reads all sheets from the reference Excel file into a named list. The
#' pipeline expects three specific sheets to be present:
#' \describe{
#'   \item{Batch_Info}{One row per replicate. Contains Batch_ID, Replicates,
#'         Replicate_Type, and Unknown_Weight.}
#'   \item{Calibrator_Values}{Calibrator names and their associated weights
#'         (concentration levels) per batch.}
#'   \item{Calibrator_Factors}{Per-batch, per-polymer multiplication factors
#'         that scale Calibrator_Level to Scaled_Concentration.}
#' }
#'
#' @param reference_path Character. Full path to the reference Excel file
#'   (e.g., \code{"path/to/Environage_QC_Sample_Information.xlsx"}).
#'
#' @return A named list with elements \code{Batch_Info},
#'   \code{Calibrator_Values}, and \code{Calibrator_Factors}, each a
#'   data frame read directly from the corresponding sheet.
#'
#' @details
#' All sheets in the workbook are read. If additional sheets are present
#' beyond the three required ones, they are included in the returned list
#' but ignored by downstream pipeline functions.
#'
#' \strong{What to modify:} Only the \code{reference_path} argument. The
#' sheet names (\code{Batch_Info}, \code{Calibrator_Values},
#' \code{Calibrator_Factors}) must match exactly — if your workbook uses
#' different sheet names, rename them in the Excel file to match.
#'
#' @examples
#' \dontrun{
#' ref <- load_reference_data("data/Environage_QC_Sample_Information.xlsx")
#' head(ref$Batch_Info)
#' head(ref$Calibrator_Factors)
#' }
#'
#' @export
load_reference_data <- function(reference_path) {
  # Validate that the file exists before attempting to read it.
  # Gives a clear error message instead of a cryptic readxl error.
  if (!file.exists(reference_path)) {
    stop("Reference file not found: ", reference_path,
         "\nCheck that the path is correct and the file is not open in Excel.")
  }

  # Read all sheet names from the workbook
  sheet_names <- readxl::excel_sheets(reference_path)

  # Read each sheet into a named list element.
  # setNames() ensures the list element names match the sheet names exactly.
  ref_list <- setNames(
    lapply(sheet_names, function(sheet) {
      readxl::read_excel(reference_path, sheet = sheet)
    }),
    sheet_names
  )

  # Verify that the three required sheets are present.
  # Stops early with a descriptive message rather than failing silently later.
  required_sheets <- c("Batch_Info", "Calibrator_Values", "Calibrator_Factors")
  missing <- setdiff(required_sheets, sheet_names)
  if (length(missing) > 0) {
    stop("The following required sheets are missing from the reference workbook:\n  ",
         paste(missing, collapse = ", "),
         "\nFound sheets: ", paste(sheet_names, collapse = ", "))
  }

  message("Loaded reference workbook: ", basename(reference_path))
  message("  Sheets found: ", paste(sheet_names, collapse = ", "))

  ref_list
}


# ── 2. load_batch_data ────────────────────────────────────────────────────────
#
#' Load all batch CSV files from a directory
#'
#' Reads every \code{.csv} file in the specified directory into a named list.
#' Each list element is a raw data frame corresponding to one batch file.
#' List element names are derived from the file basenames, sanitized with
#' \code{make.names()} to ensure they are valid R names (e.g., spaces
#' replaced with dots).
#'
#' @param data_dir Character. Path to the directory containing batch CSV files.
#'   All \code{.csv} files in this directory will be read.
#'
#' @return A named list of data frames, one per CSV file. Names are the
#'   sanitized file basenames without extensions.
#'
#' @details
#' \strong{Expected CSV format:} Each CSV should be a Skyline export containing
#' at minimum the columns \code{Replicate_Name}, \code{Molecule_Name},
#' \code{Fragment_Ion}, \code{Area}, and \code{Normalized_Area}. Column names
#' with spaces are handled downstream by \code{harmonize_batches()}.
#'
#' \strong{What to modify:} Only the \code{data_dir} argument. File naming
#' does not need to follow a specific convention beyond containing the 6-digit
#' batch number that matches entries in \code{Batch_Info$Batch_ID}.
#'
#' @examples
#' \dontrun{
#' batches <- load_batch_data("data/Batches/")
#' names(batches)   # shows all loaded batch file names
#' }
#'
#' @export
load_batch_data <- function(data_dir) {
  # Validate directory exists
  if (!dir.exists(data_dir)) {
    stop("Data directory not found: ", data_dir,
         "\nCheck that the path is correct.")
  }

  # Find all CSV files in the directory (non-recursive: only top-level files)
  file_list <- list.files(path = data_dir, pattern = "\\.csv$",
                          full.names = TRUE)

  if (length(file_list) == 0) {
    stop("No CSV files found in: ", data_dir)
  }

  message("Found ", length(file_list), " CSV file(s) in: ", data_dir)

  # Read each file. Names are sanitized basenames (spaces -> dots, etc.)
  # so they can be used as valid R list element names.
  batch_list <- setNames(
    lapply(file_list, function(f) {
      readr::read_csv(f, show_col_types = FALSE)
    }),
    make.names(tools::file_path_sans_ext(basename(file_list)))
  )

  message("Loaded batches: ", paste(names(batch_list), collapse = ", "))
  batch_list
}


# ── 3. harmonize_batches ─────────────────────────────────────────────────────
#
#' Match, harmonize, and bind batch data to the reference Batch_Info
#'
#' Takes the raw named list of batch data frames and the reference list,
#' matches each batch CSV to its Batch_ID by extracting the 6-digit batch
#' number from the Batch_ID string, coerces mixed-type columns consistently
#' across all batches, binds all batches into a single data frame, and
#' filters rows to Target fragment ions only.
#'
#' This function replaces the combination of \code{Batch_Data_List},
#' \code{list2env()}, \code{coerce_mixed_columns()}, \code{bind_rows()},
#' and the \code{filter(str_detect(Fragment_Ion, "Target"))} call from the
#' original script — all in one step, with no global environment side effects.
#'
#' @param batch_list A named list of raw batch data frames, as returned by
#'   \code{load_batch_data()}.
#' @param ref A named list as returned by \code{load_reference_data()},
#'   containing at minimum a \code{Batch_Info} element.
#' @param type_threshold Numeric scalar in (0, 1). Passed to
#'   \code{coerce_mixed_columns()}. Columns where more than this proportion
#'   of values parse as numeric will be coerced to numeric. Default 0.5.
#'
#' @return A single data frame with all matched batches bound together,
#'   column names with spaces replaced by underscores, \code{#N/A} strings
#'   converted to \code{NA}, and rows filtered to Target fragment ions only.
#'   A \code{Batch_ID} column is prepended identifying the source batch.
#'
#' @details
#' \strong{Batch matching:} The 6-digit numeric batch number is extracted
#' from each \code{Batch_ID} in \code{Batch_Info} (e.g., \code{"Batch_230615_XX"}
#' → \code{"230615"}) and matched against the names of \code{batch_list}.
#' Batch IDs with no matching CSV are silently dropped with a warning.
#'
#' \strong{What to modify:} The \code{type_threshold} argument if your data
#' has columns that are borderline numeric/character. In practice the default
#' of 0.5 works well for Skyline exports.
#'
#' @examples
#' \dontrun{
#' ref     <- load_reference_data("data/reference.xlsx")
#' batches <- load_batch_data("data/Batches/")
#' all_data <- harmonize_batches(batches, ref)
#' nrow(all_data)
#' unique(all_data$Batch_ID)
#' }
#'
#' @export
harmonize_batches <- function(batch_list, ref, type_threshold = 0.5) {
  batch_info <- ref[["Batch_Info"]]

  # Extract unique Batch_IDs from the reference sheet
  batch_ids <- unique(batch_info$Batch_ID)

  # For each Batch_ID, extract the 6-digit number and find the matching
  # CSV name in batch_list. This tolerates flexible file naming as long as
  # the 6-digit batch number appears somewhere in the filename.
  matched <- data.frame(
    Batch_ID           = batch_ids,
    Batch_Num          = stringr::str_extract(batch_ids, "\\d{6}"),
    stringsAsFactors   = FALSE
  )

  matched$Matched_Name <- vapply(matched$Batch_Num, function(num) {
    if (is.na(num)) return(NA_character_)
    hits <- names(batch_list)[stringr::str_detect(names(batch_list), num)]
    if (length(hits) > 0) hits[1] else NA_character_
  }, character(1))

  # Warn about any Batch_IDs that could not be matched to a CSV
  unmatched <- matched[is.na(matched$Matched_Name), "Batch_ID"]
  if (length(unmatched) > 0) {
    warning("The following Batch_IDs had no matching CSV and will be excluded:\n  ",
            paste(unmatched, collapse = ", "))
  }

  matched <- matched[!is.na(matched$Matched_Name), ]

  if (nrow(matched) == 0) {
    stop("No batch CSVs could be matched to Batch_Info entries. ",
         "Check that batch file names contain the 6-digit batch number ",
         "present in Batch_Info$Batch_ID.")
  }

  # Coerce column types within each matched batch, then bind all into one
  # data frame with a Batch_ID column prepended (.id argument of bind_rows).
  #
  # coerce_mixed_columns() is defined in utils.R. It handles the common
  # Skyline export issue where "#N/A" causes numeric columns to be read
  # as character in some batches but numeric in others.
  batch_data_named <- setNames(
    lapply(seq_len(nrow(matched)), function(i) {
      df <- batch_list[[ matched$Matched_Name[i] ]]
      coerce_mixed_columns(df, threshold = type_threshold)
    }),
    matched$Batch_ID
  )

  # Bind all batches; Batch_ID becomes a column from the list names
  all_batches <- dplyr::bind_rows(batch_data_named, .id = "Batch_ID")

  # Standardize column names: replace spaces with underscores.
  # Skyline sometimes exports headers with spaces (e.g., "Molecule Name").
  names(all_batches) <- gsub(" ", "_", names(all_batches))

  # Filter to Target fragment ions only.
  # Skyline exports both target and decoy/other fragment ions; only rows
  # where Fragment_Ion starts with "Target" are used for quantification.
  if (!"Fragment_Ion" %in% names(all_batches)) {
    warning("Column 'Fragment_Ion' not found. Skipping Target filter. ",
            "Check that your Skyline export includes a Fragment Ion column.")
  } else {
    all_batches <- all_batches %>%
      dplyr::filter(stringr::str_detect(Fragment_Ion, "^Target"))
  }

  # Convert "#N/A" strings (a common Skyline export artifact) to true NA,
  # then attempt to parse any remaining character columns that look numeric.
  all_batches <- all_batches %>%
    dplyr::mutate(
      dplyr::across(where(is.character), ~ dplyr::na_if(.x, "#N/A"))
    ) %>%
    readr::type_convert(col_types = readr::cols())

  message("Harmonized ", length(batch_data_named), " batch(es) into ",
          nrow(all_batches), " rows.")

  all_batches
}
