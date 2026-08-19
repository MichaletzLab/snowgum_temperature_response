# ============================================================
# 01_merge_faster_tdt.R
#
# Merges the two curve-fit outputs:
#   - data/faster_parameters_SS.csv  (Schoolfield A-Tleaf curve fits, per curveID)
#   - data/tdt_clean.csv             (TDT logistic fits: z, CTmax_1min, T50, per individual id)
#
# and attaches an Arrhenius fit (Ea_PSII, lnA_PSII) of the Fv/Fm thermal-decay
# time course (from exp_all_wide.csv) plus gsw_mean_above30 (mean stomatal
# conductance for Tleaf > 30C, per curve, from at.clean.trimmed.csv), then
# applies data-quality checks and finally joins in source-site metadata
# (data/snowgum_metadata_source_pops.csv, on code == id):
#
#   1) delete curves where Topt is fit outside the measured Tleaf range
#   2) set Tbreadth to NA if the 95%-of-Amax breadth extends beyond the
#      measured Tleaf range
#   3) set T50 to NA where the estimate falls outside 34-54 C
#   4) set z, CTmax_1min, Ea_PSII, lnA_PSII to NA (all four, row-wise) if any
#      of the four is <0
#
# Individuals whose Schoolfield fit is entirely missing (Amax and Topt
# both NA -- i.e. no curve ever existed, or it was deleted by check 1) are
# dropped from the final output altogether: if there's no surviving
# faster curve, the associated TDT/Arrhenius data for that individual is
# excluded too, not kept as a TDT-only row.
#
# data/at.clean.trimmed.csv is used only as a bridge: it's the one file
# that carries both curveID and individual id, and it's the source of
# measured Tleaf per curve for checks 1 & 2. It's also the exact data
# fit_mod_SS.R fit the Schoolfield curves against (see at.clean.trimmed
# in 00_faster_fitting.R), so its Tleaf values are the authoritative
# "measured range" for every one of the 63 curves -- unlike
# gsw_by_tleaf_clean, which was missing Tleaf data for 7 curveIDs.
# ============================================================

library(tidyverse)
library(minpack.lm)   # nlsLM, for the Arrhenius decay-rate fits

# Assumes the working directory is the project root (true by default when
# this project's .Rproj is opened in RStudio).
source("fit_mod_SS.R")   # reuse the schoolfield() curve function for check 2

# ------------------------------------------------------------
# 1. Load source data
#
# faster_parameters_SS.csv and tdt_clean.csv are raw fitting output --
# their on-disk column names (T_opt, E, E_D, breadth_90, T_opt_SE, T_opt_var,
# Topt_ED_cov, Topt_ED_cor, T50_prime) predate this naming cleanup and
# aren't touched by re-running 00_faster_fitting.R/00_tdt_fitting.R, so
# they're renamed here to the clearer names used everywhere downstream:
#   T_opt -> Topt, E -> Ea_A (apparent activation energy of the AT
#     response below Topt; Table 1), E_D -> Ed_A, breadth_90 -> Tbreadth
#     (this one was also misleading on its own terms -- fit_mod_SS.R
#     actually uses a 95% threshold, not 90%, despite the old column
#     name), T50_prime (already CTmax,1min in meaning; see
#     00_tdt_fitting.R) -> CTmax_1min
# ------------------------------------------------------------
faster <- read.csv("data/faster_parameters_SS.csv", row.names = 1, stringsAsFactors = FALSE) %>%
  rename(
    Topt = T_opt, Topt_SE = T_opt_SE, Topt_var = T_opt_var,
    Ea_A = E, Ea_A_SE = E_SE, Ea_A_var = E_var,
    Ed_A = E_D, Ed_A_SE = E_D_SE, Ed_A_var = E_D_var,
    Tbreadth = breadth_90,
    Topt_EdA_cov = Topt_ED_cov, Topt_EdA_cor = Topt_ED_cor
  )
tdt <- read.csv("data/tdt_clean.csv", stringsAsFactors = FALSE) %>%
  rename(CTmax_1min = T50_prime)
at_clean <- read.csv("data/at.clean.trimmed.csv", stringsAsFactors = FALSE)

# ------------------------------------------------------------
# at.clean.trimmed, with site + population_category added, saved as its
# own file (data/at_clean_trimmed_with_site.csv) rather than overwriting
# at.clean.trimmed.csv. Site number is parsed from individualID.x (the
# same fully-populated id column used to build the curveID<->id crosswalk
# below); all 8 sites present in this file (5,6,7,9,10,12,13,14) are
# covered by the site->population_category mapping, so nothing falls
# through to NA. This doesn't feed into the merge/QC pipeline below --
# it's just saved as a standalone output.
# ------------------------------------------------------------
at_clean_with_site <- at_clean %>%
  mutate(
    site = as.numeric(str_extract(individualID.x, "(?<=S)\\d+")),

    population_category = case_when(
      site == 13 ~ "tableland",
      site %in% c(14, 12, 7) ~ "montane",
      site %in% c(5, 6, 9, 10) ~ "sub alpine",
      TRUE ~ NA_character_
    ),

    population_category = factor(
      population_category,
      levels = c("tableland", "montane", "sub alpine")
    )
  )

write.csv(at_clean_with_site, "data/at_clean_trimmed_with_site.csv", row.names = FALSE)

# Exclusion tally, built up as the pipeline runs -- see the summary table
# printed/written at the end of the script.
exclusion_log <- list()

n_curves_start <- nrow(faster)

# curveID 8 is the base-run S05T01P01 curve; curveID 9 (S05T01P01R2) is the
# repeat measurement for the same plant. Per the same keep-the-repeat
# pattern already applied to S13T07P01/S13T02P07/S13T10P03 in
# 00_faster_fitting.R, drop the base run and keep only the R2 curve.
faster <- faster %>% filter(curveID != 8)
exclusion_log[["Duplicate base curve (S05T01P01, keep R2 only)"]] <- n_curves_start - nrow(faster)

# ------------------------------------------------------------
# 2. curveID <-> individual-id crosswalk (bridge table)
#
# faster_parameters_SS is indexed by curveID; tdt_clean is indexed by
# individual id. at.clean.trimmed.csv carries both (as individualID.x --
# the fully-populated copy from the at.df/at.meta join in
# 00_faster_fitting.R; individualID.y has 712 NAs from unresolved
# metadata matches, so it's not used here), so it's used to build the
# crosswalk. The S09T08P01 -> S09T08T01 fix is already baked into this
# file (it's applied upstream, before at.clean.trimmed is built -- see
# 00_faster_fitting.R). What's left is the same S09T02P12 -> S09T02P012
# leading-zero fix from 00_tdt_fitting.R, and stripping the R2/V2
# repeat-measurement suffixes that tdt_clean$id doesn't carry.
# ------------------------------------------------------------
id_fixes <- c(
  "S09T02P12" = "S09T02P012"
)

crosswalk <- at_clean %>%
  distinct(curveID, individualID = individualID.x) %>%
  mutate(
    id = dplyr::recode(individualID, !!!id_fixes),
    id = str_remove(id, "(R2|V2)$")   # repeat-measurement suffixes not used in tdt_clean$id
  ) %>%
  select(curveID, id)

# tdt_clean ids with no corresponding curve in at.clean.trimmed at all --
# i.e. never fit as a Schoolfield curve to begin with (as opposed to fit
# and then deleted by QC1 below). Tracked separately in the exclusion
# table so "no surviving faster curve" doesn't read as a mismatch against
# QC1's count -- QC1 only counts curves that existed and got deleted.
never_had_curve <- setdiff(tdt$id, crosswalk$id)
exclusion_log[["Individual never had a faster curve (e.g. excluded upstream as Bad, row dropped)"]] <-
  length(never_had_curve)
if (length(never_had_curve) > 0) {
  message(
    "NOTE: ", length(never_had_curve), " tdt_clean id(s) never had a faster ",
    "curve fit at all: ", paste(never_had_curve, collapse = ", ")
  )
}

still_unmatched <- setdiff(crosswalk$id, tdt$id)
if (length(still_unmatched) > 0) {
  message(
    "NOTE: ", length(still_unmatched), " curveID(s) have no tdt_clean match ",
    "after ID reconciliation: ", paste(still_unmatched, collapse = ", ")
  )
}

no_curve_match <- setdiff(faster$curveID, crosswalk$curveID)
if (length(no_curve_match) > 0) {
  message(
    "NOTE: curveID(s) not found in at.clean.trimmed.csv (no measured Tleaf range, ",
    "so checks 1 & 2 cannot be evaluated for them): ",
    paste(no_curve_match, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 3. Measured Tleaf range per curve (for checks 1 & 2), and mean gsw for
# leaf temperatures above 30C per curve
# ------------------------------------------------------------
tleaf_range <- at_clean %>%
  filter(!is.na(Tleaf)) %>%
  group_by(curveID) %>%
  summarise(
    Tleaf_min = min(Tleaf),
    Tleaf_max = max(Tleaf),
    .groups = "drop"
  )

gsw_above30 <- at_clean %>%
  filter(!is.na(Tleaf), Tleaf > 30, !is.na(gsw)) %>%
  group_by(curveID) %>%
  summarise(
    gsw_mean_above30 = mean(gsw),
    .groups = "drop"
  )

faster <- faster %>%
  left_join(crosswalk, by = "curveID") %>%
  left_join(tleaf_range, by = "curveID") %>%
  left_join(gsw_above30, by = "curveID")

# ------------------------------------------------------------
# QC 1: delete curves where Topt is fit outside the measured Tleaf range
# (curves with no measured range at all are left as-is -- there's nothing
# to check them against -- and were already flagged above)
# ------------------------------------------------------------
topt_out_of_range <- faster %>%
  filter(!is.na(Tleaf_min), Topt < Tleaf_min | Topt > Tleaf_max) %>%
  pull(curveID)

if (length(topt_out_of_range) > 0) {
  message(
    "QC1: dropping ", length(topt_out_of_range), " curve(s) with Topt outside ",
    "measured Tleaf range: ", paste(topt_out_of_range, collapse = ", ")
  )
}

faster <- faster %>%
  filter(is.na(Tleaf_min) | (Topt >= Tleaf_min & Topt <= Tleaf_max))
exclusion_log[["QC1: Topt outside measured Tleaf range (curve deleted)"]] <- length(topt_out_of_range)

# ------------------------------------------------------------
# QC 2: set Tbreadth to NA if the 95%-of-Amax breadth extends beyond
# the measured Tleaf range.
#
# faster_parameters_SS only stores the breadth WIDTH, not its lower/upper
# bounds, so the bounds are reconstructed here using the same method
# fit_mod_SS.R used to originally compute Tbreadth: rebuild the fitted
# schoolfield() curve from the saved coefficients, evaluate it over
# [Tleaf_min - 5, Tleaf_max + 5] (the same extrapolation window used at
# fitting time), and take the range of temperatures where predicted A is
# >= 95% of Amax. That reproduces Tbreadth exactly; it's then compared
# against the *strictly measured* Tleaf_min/Tleaf_max (no +/-5 padding)
# to decide whether the breadth estimate leans on extrapolated, not
# measured, temperatures.
# ------------------------------------------------------------
recompute_breadth_bounds <- function(J_ref, Ea_A, Ed_A, Topt, Amax, Tleaf_min, Tleaf_max) {
  if (any(is.na(c(J_ref, Ea_A, Ed_A, Topt, Amax, Tleaf_min, Tleaf_max)))) {
    return(c(NA_real_, NA_real_))
  }
  T_extrap <- seq(Tleaf_min - 5, Tleaf_max + 5, by = 0.1)
  # schoolfield() (fit_mod_SS.R) still uses its own original parameter
  # names (E, E_D, T_opt) -- only this script's column/variable names changed.
  pred <- schoolfield(temp = T_extrap, J_ref = J_ref, E = Ea_A, E_D = Ed_A, T_opt = Topt)
  above <- T_extrap[pred >= 0.95 * Amax]
  if (length(above) > 1) range(above) else c(NA_real_, NA_real_)
}

n_breadth_before <- sum(!is.na(faster$Tbreadth))

faster <- faster %>%
  rowwise() %>%
  mutate(
    breadth_bounds = list(recompute_breadth_bounds(J_ref, Ea_A, Ed_A, Topt, Amax, Tleaf_min, Tleaf_max)),
    breadth_low  = breadth_bounds[1],
    breadth_high = breadth_bounds[2]
  ) %>%
  ungroup() %>%
  mutate(
    Tbreadth = if_else(
      !is.na(Tleaf_min) & !is.na(breadth_low) &
        (breadth_low < Tleaf_min | breadth_high > Tleaf_max),
      NA_real_,
      Tbreadth
    )
  ) %>%
  select(-breadth_bounds, -breadth_low, -breadth_high)

exclusion_log[["QC2: Tbreadth extends beyond measured Tleaf range (set NA)"]] <-
  n_breadth_before - sum(!is.na(faster$Tbreadth))

# ------------------------------------------------------------
# 4. Arrhenius fit (Ea_PSII, lnA_PSII) from the Fv/Fm thermal-decay time
# course in exp_all_wide.csv -- the same underlying dataset
# 00_tdt_fitting.R uses to fit T50/z, refit here as first-order decay
# rate k per individual x temperature, then Ea/lnA_PSII from the Arrhenius
# regression of ln(k) on 1/Temp_K across each individual's temperatures.
#
# exp_all_wide.csv$id still carries the RAW (uncorrected) individual
# codes, so the same code_fixes table from 00_tdt_fitting.R is reapplied
# here to line the ids up with tdt_clean$id / the crosswalk above.
# ------------------------------------------------------------
code_fixes <- c(
  "S03T03P06" = "S05T03P06",
  "S05P02T23" = "S05T02P23",
  "S09T02P12" = "S09T02P012",
  "S09T05P09" = "S07T05P09",
  "S09T08P01" = "S09T08T01",
  "S12T03P07" = "S12T03P04",
  "S14T03P04" = "S14T03P09",
  "S05T02P09" = "S05T02P07",
  "S9T07P01"  = "S09T07P01"
)

df <- read.csv("exp_all_wide.csv", stringsAsFactors = FALSE) %>%
  mutate(
    plant_code = dplyr::recode(id, !!!code_fixes),
    Temp_K     = Temp + 273.15,
    fvfm_norm  = final / initial
  )

k_estimates <- df %>%
  group_by(plant_code, Temp_K) %>%
  nest() %>%
  mutate(
    fit = purrr::map(data, ~{
      tryCatch(
        nlsLM(
          fvfm_norm ~ fvfm0 * exp(-k * Time),
          data = .x,
          start = list(fvfm0 = max(.x$fvfm_norm, na.rm = TRUE), k = 0.1),
          lower = c(fvfm0 = 0, k = 0),
          upper = c(fvfm0 = 1, k = Inf),
          control = nls.lm.control(maxiter = 200)
        ),
        error = function(e) NULL
      )
    }),
    k = purrr::map_dbl(fit, ~{
      if (!is.null(.x)) coef(.x)["k"] else NA_real_
    })
  ) %>%
  ungroup() %>%   # group_by(plant_code, Temp_K) %>% nest() leaves the tibble grouped by
                   # both vars; without ungroup() the count(plant_code) below tallies
                   # within those 1-row leftover groups instead of across each plant's
                   # temperatures, so every plant would fail the n >= 3 threshold below.
  select(plant_code, Temp_K, k) %>%
  mutate(
    inv_T = 1 / Temp_K,
    ln_k  = log(k)
  )

k_estimates_clean <- k_estimates %>%
  filter(is.finite(k), k > 0)

valid_plants <- k_estimates_clean %>%
  count(plant_code) %>%
  filter(n >= 3) %>%
  pull(plant_code)

k_estimates_clean <- k_estimates_clean %>%
  filter(plant_code %in% valid_plants)

arrhenius_fits <- k_estimates_clean %>%
  group_by(plant_code) %>%
  summarise(fit = list(lm(ln_k ~ inv_T)), .groups = "drop") %>%
  mutate(
    Ea_Jmol  = -purrr::map_dbl(fit, ~ coef(.x)[["inv_T"]] * 8.314),
    Ea_PSII = Ea_Jmol / 1000,
    lnA_PSII      = purrr::map_dbl(fit, ~ coef(.x)[["(Intercept)"]])
  ) %>%
  select(id = plant_code, Ea_PSII, lnA_PSII)

# ------------------------------------------------------------
# 5. Merge faster (post QC1/QC2) + tdt_clean + Arrhenius fit
#
# full_join on id so that individuals with TDT/Arrhenius data but no
# surviving Schoolfield curve (dropped in QC1, or never had a curve at
# all) are still kept, just with the Schoolfield-derived columns NA.
# ------------------------------------------------------------
merged <- faster %>%
  full_join(tdt, by = "id") %>%
  left_join(arrhenius_fits, by = "id")

# ------------------------------------------------------------
# QC 3: T50 outside 34-54 C -> NA
# ------------------------------------------------------------
n_t50_before <- sum(!is.na(merged$T50))

merged <- merged %>%
  mutate(
    T50 = if_else(!is.na(T50) & (T50 < 34 | T50 > 54), NA_real_, T50)
  )

exclusion_log[["QC3: T50 outside 34-54C (set NA)"]] <- n_t50_before - sum(!is.na(merged$T50))

# ------------------------------------------------------------
# QC 4: if any of z, CTmax_1min, Ea_PSII, lnA_PSII is <0 for a row, null out
# all four (not just the negative one)
# ------------------------------------------------------------
merged <- merged %>%
  mutate(
    .any_negative = purrr::pmap_lgl(
      list(z, CTmax_1min, Ea_PSII, lnA_PSII),
      ~ any(c(...) < 0, na.rm = TRUE)
    ),
    z         = if_else(.any_negative, NA_real_, z),
    CTmax_1min = if_else(.any_negative, NA_real_, CTmax_1min),
    Ea_PSII  = if_else(.any_negative, NA_real_, Ea_PSII),
    lnA_PSII       = if_else(.any_negative, NA_real_, lnA_PSII)
  )

exclusion_log[["QC4: negative z/CTmax_1min/Ea_PSII/lnA_PSII (row nulled)"]] <- sum(merged$.any_negative)

merged <- merged %>% select(-.any_negative)

# ------------------------------------------------------------
# Drop rows with no surviving faster curve: if Amax and Topt are both
# missing (no curve fit at all, or QC1 deleted it), the individual's
# TDT/Arrhenius data is excluded from the final output too, rather than
# being kept as a TDT-only row. This total should equal (and is checked
# against) "QC1: Topt outside range" + "never had a faster curve" above
# -- those are the only two ways a row can end up missing Amax/Topt.
# ------------------------------------------------------------
n_before_curve_filter <- nrow(merged)

merged <- merged %>%
  filter(!is.na(Amax) & !is.na(Topt))

n_dropped_no_curve <- n_before_curve_filter - nrow(merged)
exclusion_log[["Total rows dropped for missing Amax/Topt (= QC1 + never-had-curve, check)"]] <-
  n_dropped_no_curve

stopifnot(n_dropped_no_curve == length(topt_out_of_range) + length(never_had_curve))

# ------------------------------------------------------------
# 6. Attach source-site metadata (snowgum_metadata_source_pops.csv), joined
# on code == id. This is the ground-truth site/tree/climate (BIOx) table
# the earlier ID reconciliations in 00_faster_fitting.R / 00_tdt_fitting.R
# were checked against. Nine codes in it are duplicated (different
# site/tree/block combos sharing a code, e.g. two "S04T03P04" rows with
# tree 3 vs tree 5) -- none of them are among our 57 remaining ids, so a
# left_join here is safe and won't fan out any rows. If a future rerun
# adds one of those codes to the merged data, this join would silently
# duplicate that row, so it's worth re-checking the intersect if the set
# of surviving curves changes.
# ------------------------------------------------------------
site_meta <- read.csv("data/snowgum_metadata_source_pops.csv", stringsAsFactors = FALSE)

dup_codes_in_use <- intersect(merged$id, site_meta$code[duplicated(site_meta$code)])
if (length(dup_codes_in_use) > 0) {
  stop(
    "snowgum_metadata_source_pops.csv has duplicate 'code' rows for id(s) ",
    "still in the merged data (", paste(dup_codes_in_use, collapse = ", "),
    ") -- a left_join would fan out these rows. Resolve the duplicates ",
    "before merging."
  )
}

merged <- merged %>%
  left_join(site_meta, by = c("id" = "code"))

# ------------------------------------------------------------
# 7. Exclusion summary table
# ------------------------------------------------------------
exclusion_summary <- data.frame(
  Criterion = names(exclusion_log),
  N_excluded = unlist(exclusion_log),
  row.names = NULL
)

print(exclusion_summary)
write.csv(exclusion_summary, "data/qc_exclusion_summary.csv", row.names = FALSE)

# ------------------------------------------------------------
# 8. Write output
# ------------------------------------------------------------
write.csv(merged, "data/faster_tdt_merged.csv", row.names = FALSE)

message("Done. ", nrow(merged), " rows written to data/faster_tdt_merged.csv")
