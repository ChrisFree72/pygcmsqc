# example_config.R
#
# pygcmsqc pipeline configuration and run script.
#
# ── HOW TO USE ────────────────────────────────────────────────────────────────
#
#   1. Copy this file to your project folder (outside the package directory).
#   2. Fill in the three path variables and two analytical settings below.
#   3. Run the entire script with source("example_config.R") or step through
#      it section by section in RStudio.
#
# That is the ONLY file you need to modify to run the full pipeline.
# Do not edit files inside the pygcmsqc/R/ folder — all user configuration
# is controlled from here.
#
# ── WHAT THIS SCRIPT DOES ────────────────────────────────────────────────────
#
#   Step 1  Load reference metadata and batch CSV files
#   Step 2  Harmonize and bind all batches into one data frame
#   Step 3  Join calibration factors, calibrator levels, replicate types,
#           and sample weights
#   Step 4  Summarize signal areas by molecule prefix; remove excluded molecules
#   Step 5  Background-subtract blank signal from all non-standard replicates
#   Step 6  Fit weighted 1/x calibration curves; optimize calibrator selection
#   Step 7  Back-calculate concentrations and ug/g values
#   Step 8  Compute LOD and LOQ from blank replicate distributions
#   Step 9  Export calibration curve PDFs
#   Step 10 Export internal standard plots
#   Step 11 Export QC plots
#   Step 12 Export concentration bar chart PDFs
#   Step 13 Export PCA and scree plots
#   Step 14 Export Excel workbooks


# ══ 0. Install / load package ═════════════════════════════════════════════════
#
# Install from GitHub (run once, then comment out):
# devtools::install_github("yourusername/pygcmsqc")

library(pygcmsqc)


# ══ 1. CONFIGURE PATHS ════════════════════════════════════════════════════════
#
# reference_path: Full path to the QC reference Excel workbook.
#   Required sheets: Batch_Info, Calibrator_Values, Calibrator_Factors
#   (sheet names must match exactly)
#
reference_path <- "path/to/Environage_QC_Sample_Information.xlsx"

# data_dir: Full path to the folder containing batch CSV files.
#   All .csv files in this folder will be loaded.
#   File names must contain the 6-digit batch number that matches
#   Batch_Info$Batch_ID (e.g., a file named "Batch_230615_skyline.csv"
#   will match Batch_ID entries containing "230615").
#
data_dir <- "path/to/Data/"

# output_dir: Full path to the folder where all results will be written.
#   Subfolders are created automatically:
#     output_dir/calibration_curves/
#     output_dir/internal_standard/
#     output_dir/qc_plots/
#     output_dir/concentration_bars/
#     output_dir/pca_output/
#   Excel workbooks are written directly to output_dir/.
#
output_dir <- "path/to/Results/"


# ══ 2. CONFIGURE ANALYTICAL SETTINGS ══════════════════════════════════════════

# excluded_molecules: Molecule_Name values to remove before summarization.
#   Add or remove entries as needed for your dataset.
#   These are full Molecule_Name strings (not prefixes), e.g. "PE_2" not "PE".
#   Set to character(0) to exclude nothing.
#
excluded_molecules <- c("N66_2", "PC_1", "PET_2", "PVC_1")

# max_remove: Maximum number of non-anchor calibrators the optimizer may
#   remove per batch to improve calibration curve R².
#   Default: 2. The lowest-concentration calibrator is always retained.
#   Runtime grows combinatorially — increasing beyond 3 is not recommended
#   without first checking the estimated execution time printed at runtime.
#
max_remove <- 2

# lod_multiplier / loq_multiplier: Multipliers applied to blank SD for LOD/LOQ.
#   Default values (3 and 10) are the international standard (ICH Q2(R1)).
#   Change only if your laboratory SOP specifies different multipliers.
#
lod_multiplier <- 3
loq_multiplier <- 10


# ══ 3. CREATE OUTPUT DIRECTORY ════════════════════════════════════════════════
#
# Create the top-level output directory if it does not already exist.
# Subfolders (calibration_curves/, pca_output/, etc.) are created
# automatically by each plotting function.
#
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}


# ══ 4. LOAD DATA ══════════════════════════════════════════════════════════════

message("\n── Step 1: Loading reference data ──────────────────────────────────")
ref <- load_reference_data(reference_path)

message("\n── Step 2: Loading batch CSV files ─────────────────────────────────")
batches <- load_batch_data(data_dir)

message("\n── Step 3: Harmonizing batches ─────────────────────────────────────")
all_data <- harmonize_batches(batches, ref)


# ══ 5. PROCESS DATA ═══════════════════════════════════════════════════════════

message("\n── Step 4: Joining calibration metadata ─────────────────────────────")
all_data <- join_calibration_metadata(
  all_data,
  calibrator_factors = ref$Calibrator_Factors,
  calibrator_values  = ref$Calibrator_Values,
  batch_info         = ref$Batch_Info
)

message("\n── Step 5: Summarizing by molecule prefix ───────────────────────────")
all_data <- summarize_batches(all_data, excluded_molecules = excluded_molecules)

message("\n── Step 6: Subtracting blank background ─────────────────────────────")
all_data <- subtract_blanks(all_data)


# ══ 6. CALIBRATION ════════════════════════════════════════════════════════════

message("\n── Step 7: Running calibration optimization ─────────────────────────")
# run_calibration() prints an estimated execution time before starting and
# reports actual elapsed time on completion. Adjust max_remove above if the
# estimate is unexpectedly long.
models <- run_calibration(all_data, max_remove = max_remove)

message("\n── Step 8: Calculating concentrations ───────────────────────────────")
all_data <- calculate_concentrations(all_data, models)


# ══ 7. LOD / LOQ ══════════════════════════════════════════════════════════════

message("\n── Step 9: Computing LOD and LOQ ────────────────────────────────────")
lod_loq <- calculate_lod_loq(
  all_data,
  lod_multiplier = lod_multiplier,
  loq_multiplier = loq_multiplier
)


# ══ 8. PLOTS ══════════════════════════════════════════════════════════════════

message("\n── Step 10: Exporting calibration curve PDFs ────────────────────────")
# Output: output_dir/calibration_curves/
# Files:  Cal_Curves_Normalized_Batch_<ID>.pdf
#         Cal_Curves_Unnormalized_Batch_<ID>.pdf  (one per batch)
export_calibration_pdfs(all_data, lod_loq, out_dir = output_dir)

message("\n── Step 11: Exporting internal standard plots ───────────────────────")
# Output: output_dir/internal_standard/
# Files:  IS_Averages_Per_Batch.pdf  (cross-batch trend)
#         IS_Batch_Plots.pdf          (per-batch replicate level)
plot_is_across_batches(
  all_data,
  out_dir     = output_dir,
  output_file = "IS_Averages_Per_Batch.pdf"
)
plot_is_per_batch(
  all_data,
  out_dir     = output_dir,
  output_file = "IS_Batch_Plots.pdf"
)

message("\n── Step 12: Exporting QC plots ──────────────────────────────────────")
# Output: output_dir/qc_plots/
# Files:  QC_Averages_Per_Batch.pdf              (batch trend per analyte)
#         QC_Replicate_Plots_Per_Analyte.pdf     (replicate level per analyte)
#         QC_Bar_AllBatches.pdf                  (all-batches bar summary)
#         QC_Bar_ByBatch.pdf                     (per-batch bar charts)
plot_qc_across_batches(
  all_data,
  out_dir     = output_dir,
  output_file = "QC_Averages_Per_Batch.pdf"
)
plot_qc_replicate_level(
  all_data,
  out_dir       = output_dir,
  output_file   = "QC_Replicate_Plots_Per_Analyte.pdf",
  sd_multiplier = 1,                               # change to 2 for wider band
  batch_order   = sort(unique(all_data$Batch_ID))
)
plot_qc_bar_all_batches(
  all_data,
  out_dir     = output_dir,
  output_file = "QC_Bar_AllBatches.pdf",
  sort_bars   = "desc"
)
plot_qc_bar_by_batch(
  all_data,
  out_dir     = output_dir,
  output_file = "QC_Bar_ByBatch.pdf",
  sort_bars   = "desc"
)

message("\n── Step 13: Exporting concentration bar PDFs ────────────────────────")
# Output: output_dir/concentration_bars/
# Files:  Bar_CalcConcs_Normalized_Batch_<ID>.pdf
#         Bar_CalcConcs_Unnormalized_Batch_<ID>.pdf  (one per batch)
export_concentration_bar_pdfs(
  all_data,
  lod_loq,
  out_dir   = output_dir,
  label_col = "Replicate_Name"
)

message("\n── Step 14: Exporting PCA and scree plots ───────────────────────────")
# Output: output_dir/pca_output/
# Files:  PCA_Scree_by_batch_normalized.pdf        (per-batch, all rep types)
#         PCA_Scree_by_batch_unnormalized.pdf
#         PCA_Scree_pooled_normalized_UNKNOWN.pdf   (pooled Unknowns only)
#         PCA_Scree_pooled_unnormalized_UNKNOWN.pdf
export_pca_by_batch(
  all_data,
  out_dir   = output_dir,
  na_method = "median"      # or "drop" to use complete cases only
)
export_pooled_unknown_pca(
  all_data,
  out_dir             = output_dir,
  na_method           = "median",
  use_combined_labels = TRUE,   # FALSE for cleaner labels with many batches
  legend_position     = "bottom"
)


# ══ 9. EXPORT EXCEL WORKBOOKS ═════════════════════════════════════════════════

message("\n── Step 15: Exporting Excel workbooks ───────────────────────────────")

# R2_Comparison.xlsx: calibration model parameters and R² audit table
# Written to: output_dir/R2_Comparison.xlsx
export_r2_table(models, out_dir = output_dir)

# Total_Batch_Info.xlsx: full master data table (all replicates, all columns)
# Written to: output_dir/Total_Batch_Info.xlsx
export_batch_info(all_data, out_dir = output_dir)

# Final_Sample_Data_Output.xlsx: wide-format Unknown-only ug/g concentrations
#   Sheet "Normalized"   — Norm_ug_Per_Gram pivoted wide
#   Sheet "Unnormalized" — Unnorm_ug_Per_Gram pivoted wide
# Written to: output_dir/Final_Sample_Data_Output.xlsx
export_wide_concentrations(all_data, out_dir = output_dir)

message("\n── Pipeline complete ─────────────────────────────────────────────────")
message("All outputs written to: ", normalizePath(output_dir))
