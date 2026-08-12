# snowgum_temperature_response

Thermal performance and heat-tolerance data/analysis for *Eucalyptus pauciflora*
(snow gum) populations across an elevation gradient (MSFSS project, Chapter 7).
Combines gas-exchange photosynthesis-temperature curves (FASTER), PSII thermal
death time (TDT), and Fv/Fm thermal-decay (Arrhenius) data into one per-curve
dataset, applies data-quality checks, and produces the manuscript's figures/tables.

## Pipeline

Scripts are numbered in run order. Open `snowgum_temperature_response.Rproj`
in RStudio first -- it sets the working directory to the project root, which
every script relies on (no absolute paths).

1. **`00_faster_fitting.R`** -- reads raw LI-6800 gas-exchange files from
   `licor_data/`, reconciles individual IDs against `data/Faster_key_msfss_v2.csv`,
   trims each curve, and fits a Schoolfield thermal-performance curve
   (`fit_mod_SS.R`) per curve.
   - Writes `data/at.clean.trimmed.csv` (trimmed per-curve A-Tleaf gas-exchange
     data) and `data/faster_parameters_SS.csv` (per-curve fit parameters:
     `T_opt`, `Amax`, `breadth_90`, etc.)
   - Also writes per-curve diagnostic plots to `./SS_fits/` (not included in
     this repo -- regenerated on rerun).

2. **`00_tdt_fitting.R`** -- reads `exp_all_wide.csv` (Fv/Fm vs. heat-exposure
   time/temperature), fits a logistic T50 model per individual x temperature,
   then a thermal death time (TDT) model per individual.
   - Writes `data/tdt_clean.csv` (`z`, `T50_prime`, `T50` per individual).

3. **`01_merge_faster_tdt.R`** -- merges the two outputs above (bridged via
   `data/at.clean.trimmed.csv`, which carries both `curveID` and individual
   id), refits an Arrhenius model (`Ea_kJmol`, `lnA`) from `exp_all_wide.csv`,
   attaches `gsw_mean_above30` and source-site metadata
   (`data/snowgum_metadata_source_pops.csv`), and applies data-quality checks:
   1. delete curves where `T_opt` is fit outside the measured Tleaf range
   2. set `breadth_90` to NA if its 90%-of-Amax breadth extends beyond the
      measured Tleaf range
   3. set `T50` to NA outside 34-54°C
   4. null out `z`/`T50_prime`/`Ea_kJmol`/`lnA` (all four) if any is negative
   - Writes `data/faster_tdt_merged.csv` (final per-curve dataset),
     `data/at_clean_trimmed_with_site.csv` (raw curve data + `site`/
     `population_category`), and `data/qc_exclusion_summary.csv` (exclusion
     tally by criterion).

4. **`01_analyses.R`** -- reads `data/faster_tdt_merged.csv` and
   `data/at_clean_trimmed_with_site.csv`, produces the manuscript's figures
   and tables into `figures_updated/`.

## Supporting files

- `correct_6800_functions.R` -- LI-6800 file-reading/correction helpers,
  sourced by `00_faster_fitting.R`.
- `fit_mod_SS.R` -- Schoolfield curve function and per-curve fitting routine,
  sourced by `00_faster_fitting.R` and `01_merge_faster_tdt.R`.

## Data

| File | Description |
|---|---|
| `licor_data/` | Raw LI-6800 export files, one per gas-exchange curve. |
| `data/Faster_key_msfss_v2.csv` | Field metadata / curveID-to-barcode key for the gas-exchange curves. |
| `exp_all_wide.csv` | Fv/Fm heat-exposure time-course data (input to TDT and Arrhenius fits). |
| `data/snowgum_metadata_source_pops.csv` | Source-population site/tree/climate (BIOx) metadata, ground truth for individual IDs. |
| `data/at.clean.trimmed.csv` | Trimmed per-curve gas-exchange data (output of step 1, input to step 3). |
| `data/faster_parameters_SS.csv` | Per-curve Schoolfield fit parameters (output of step 1). |
| `data/tdt_clean.csv` | Per-individual TDT fit parameters (output of step 2). |
| `data/faster_tdt_merged.csv` | Final merged, quality-checked per-curve dataset (output of step 3). |
| `data/at_clean_trimmed_with_site.csv` | `at.clean.trimmed.csv` + `site`/`population_category` (output of step 3). |
| `data/qc_exclusion_summary.csv` | Row/curve exclusion counts by QC criterion (output of step 3). |

## Requirements

R (4.x) with: `tidyverse`, `nls.multstart`, `nlstools`, `lubridate`,
`segmented`, `tidylog`, `readxl`, `rTPC`, `minpack.lm`, `car`, `patchwork`,
`gam`, `mgcv`, `lme4`, `lmerTest`, `Hmisc`, `broom`, `maps`, `ggrepel`.
