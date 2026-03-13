# pygcmsqc

**Quality Control Pipeline for Pyrolysis GC-MS Data**

`pygcmsqc` is an R package for processing, calibrating, and visualizing pyrolysis gas chromatography-mass spectrometry (py-GC-MS) data across multiple analytical batches. It takes raw Skyline CSV exports and a reference Excel workbook as input, and produces publication-ready PDFs and styled Excel workbooks as output.

---

## What the pipeline does

| Step | Function | Output |
|------|----------|--------|
| Load reference metadata | `load_reference_data()` | Named list of reference sheets |
| Load batch CSVs | `load_batch_data()` | Named list of raw data frames |
| Harmonize & bind batches | `harmonize_batches()` | Single combined data frame |
| Join calibration metadata | `join_calibration_metadata()` | Calibration factors, replicate types, sample weights |
| Summarize by molecule prefix | `summarize_batches()` | Summed Area / Normalized_Area per analyte |
| Blank background subtraction | `subtract_blanks()` | `Area_BG_Sub`, `Normalized_Area_BG_Sub` |
| Weighted 1/x calibration | `run_calibration()` | Slopes, intercepts, R² per batch × analyte |
| Back-calculate concentrations | `calculate_concentrations()` | `Norm_μg_Per_Gram`, `Unnorm_μg_Per_Gram` |
| LOD / LOQ | `calculate_lod_loq()` | Detection & quantification limits in μg/g |
| Calibration curve PDFs | `export_calibration_pdfs()` | Per-batch PDFs with regression lines + LOD/LOQ |
| IS diagnostic plots | `plot_is_across_batches()`, `plot_is_per_batch()` | Internal standard QC PDFs |
| QC diagnostic plots | `plot_qc_across_batches()`, `plot_qc_replicate_level()`, `plot_qc_bar_all_batches()`, `plot_qc_bar_by_batch()` | QC replicate PDFs |
| Concentration bar charts | `export_concentration_bar_pdfs()` | Per-batch bar chart PDFs |
| PCA & scree plots | `export_pca_by_batch()`, `export_pooled_unknown_pca()` | PCA diagnostic PDFs |
| Excel exports | `export_r2_table()`, `export_batch_info()`, `export_wide_concentrations()` | Styled Excel workbooks |

---

## Installation

Install the development version from GitHub using `devtools`:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install pygcmsqc
devtools::install_github("yourusername/pygcmsqc")
```

---

## Input files required

### 1. Reference Excel workbook
A single `.xlsx` file with **three required sheets** (names must match exactly):

| Sheet | Contents |
|-------|----------|
| `Batch_Info` | One row per replicate. Columns: `Batch_ID`, `Replicates`, `Replicate_Type`, `Unknown_Weight` |
| `Calibrator_Values` | Calibrator names and weights per batch. Columns: `Batch_ID`, `Calibrator_Name`, `Calibrator_Weight` |
| `Calibrator_Factors` | Per-batch, per-polymer multiplication factors. Columns: `Batch_ID`, `Plastic`, `Multiplication_Factor` |

### 2. Batch CSV files (Skyline exports)
- One `.csv` file per analytical batch, exported from Skyline
- All CSV files placed in a single folder
- File names must contain the **6-digit batch number** that matches entries in `Batch_Info$Batch_ID`
  - Example: a file named `Batch_230615_skyline.csv` will match any `Batch_ID` containing `230615`
- Required columns: `Replicate Name`, `Molecule Name`, `Fragment Ion`, `Area`, `Normalized Area`

---

## Quick start

Copy `example_config.R` to your project folder, fill in the three paths, and run it:

```r
# ── Configure these three paths ───────────────────────────────────────────────
reference_path <- "path/to/Environage_QC_Sample_Information.xlsx"
data_dir       <- "path/to/Data/"
output_dir     <- "path/to/Results/"

# ── Analytical settings (defaults are fine for most projects) ─────────────────
excluded_molecules <- c("N66_2", "PC_1", "PET_2", "PVC_1")
max_remove         <- 2    # max calibrators optimizer may remove per batch
lod_multiplier     <- 3    # ICH Q2(R1) standard
loq_multiplier     <- 10   # ICH Q2(R1) standard
```

Then run the full pipeline with:

```r
source("example_config.R")
```

Or run it step by step in RStudio — each step is clearly labeled in the script.

---

## Output folder structure

All outputs are written to `output_dir/`. Subfolders are created automatically.

```
output_dir/
├── calibration_curves/
│   ├── Cal_Curves_Normalized_Batch_<ID>.pdf
│   └── Cal_Curves_Unnormalized_Batch_<ID>.pdf
├── internal_standard/
│   ├── IS_Averages_Per_Batch.pdf
│   └── IS_Batch_Plots.pdf
├── qc_plots/
│   ├── QC_Averages_Per_Batch.pdf
│   ├── QC_Replicate_Plots_Per_Analyte.pdf
│   ├── QC_Bar_AllBatches.pdf
│   └── QC_Bar_ByBatch.pdf
├── concentration_bars/
│   ├── Bar_CalcConcs_Normalized_Batch_<ID>.pdf
│   └── Bar_CalcConcs_Unnormalized_Batch_<ID>.pdf
├── pca_output/
│   ├── PCA_Scree_by_batch_normalized.pdf
│   ├── PCA_Scree_by_batch_unnormalized.pdf
│   ├── PCA_Scree_pooled_normalized_UNKNOWN.pdf
│   └── PCA_Scree_pooled_unnormalized_UNKNOWN.pdf
├── R2_Comparison.xlsx
├── Total_Batch_Info.xlsx
└── Final_Sample_Data_Output.xlsx
```

---

## Calibration method

The pipeline uses **weighted 1/x linear regression** for calibration curve fitting. This is appropriate for py-GC-MS data because signal variance increases with concentration (heteroscedastic response). Weighting by `1/x` gives proportionally more influence to low-concentration calibrators, improving accuracy across the full dynamic range — particularly at the low end where most environmental samples fall.

Reference: ICH Q2(R1), US EPA Method guidance on weighted calibration.

### Calibrator optimization

For each batch, the pipeline tests all valid combinations of removing up to `max_remove` calibrators (excluding the lowest-concentration anchor point) and selects the subset that maximizes the average R² across all analytes. This is useful when one or two calibrators are outliers due to instrument noise or sample preparation issues.

---

## LOD and LOQ

LOD and LOQ are calculated from the standard deviation of back-calculated blank concentrations (in μg/g):

```
LOD = lod_multiplier × SD(blank μg/g)    [default: 3×]
LOQ = loq_multiplier × SD(blank μg/g)    [default: 10×]
```

These are in μg/g units and are overlaid directly on calibration curve and concentration bar chart plots.

---

## Dependencies

```r
tidyverse, readxl, openxlsx, ggrepel, patchwork,
forcats, purrr, rlang, stringr, dplyr, tidyr,
ggplot2, tibble, stats, grDevices, grid
```

All dependencies are installed automatically when you run `devtools::install_github()`.

---

## Package structure

```
pygcmsqc/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── utils.R                  # Internal helpers
│   ├── data_loading.R           # load_reference_data, load_batch_data, harmonize_batches
│   ├── data_processing.R        # join_calibration_metadata, summarize_batches, subtract_blanks
│   ├── calibration.R            # run_calibration, calculate_concentrations
│   ├── lod_loq.R                # calculate_lod_loq
│   ├── plots_calibration.R      # export_calibration_pdfs
│   ├── plots_concentrations.R   # export_concentration_bar_pdfs
│   ├── plots_is.R               # plot_is_across_batches, plot_is_per_batch
│   ├── plots_qc.R               # plot_qc_* functions
│   ├── plots_pca.R              # export_pca_by_batch, export_pooled_unknown_pca
│   └── export.R                 # export_r2_table, export_batch_info, export_wide_concentrations
└── example_config.R             # User-facing run script — edit this to run the pipeline
```

---

## License

MIT
