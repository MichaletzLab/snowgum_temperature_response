library(tidyverse)

# ============================================================
# 1. Read and prepare data
# ============================================================

dat_wide_clean <- read.csv("exp_all_wide.csv") %>%
  mutate(
    fvfm_prop_if = final / initial
  )

temp_min <- min(dat_wide_clean$Temp, na.rm = TRUE)
temp_max <- max(dat_wide_clean$Temp, na.rm = TRUE)


# ============================================================
# 2. Logistic model
# ============================================================

logistic_model <- function(df) {
  tryCatch(
    nls(
      fvfm_prop_if ~ 1 / (1 + exp(b * (Temp - T50))),
      data = df,
      start = list(
        b = 0.5,
        T50 = mean(df$Temp)
      ),
      control = nls.control(maxiter = 200)
    ),
    error = function(e) NULL
  )
}


# ============================================================
# 3. Bootstrap T50
# ============================================================

boot_T50 <- function(df, nboot = 100) {
  
  m <- logistic_model(df)
  
  if (is.null(m)) {
    return(tibble(
      T50 = NA_real_,
      T50_lower = NA_real_,
      T50_upper = NA_real_
    ))
  }
  
  T50_hat <- coef(m)[["T50"]]
  
  boot_vals <- map_dbl(seq_len(nboot), function(i) {
    
    boot_df <- df %>%
      slice_sample(prop = 1, replace = TRUE)
    
    m_boot <- logistic_model(boot_df)
    
    if (is.null(m_boot)) return(NA_real_)
    
    coef(m_boot)[["T50"]]
  })
  
  tibble(
    T50 = T50_hat,
    T50_lower = quantile(boot_vals, 0.025, na.rm = TRUE),
    T50_upper = quantile(boot_vals, 0.975, na.rm = TRUE)
  )
}


# ============================================================
# 4. Fit T50 for each individual × AOI × exposure duration
# ============================================================

T50_estimates <- dat_wide_clean %>%
  group_by(id, AOI, Time) %>%
  nest() %>%
  mutate(T50_info = map(data, boot_T50)) %>%
  unnest(T50_info) %>%
  ungroup() %>%
  filter(
    T50 >= temp_min,
    T50 <= temp_max
  ) %>%
  mutate(logTime = log10(Time))


# ============================================================
# 5. Average AOI replicates for TDT analysis
# ============================================================

T50_for_TDT <- T50_estimates %>%
  group_by(id, Time) %>%
  summarise(
    T50 = mean(T50, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(logTime = log10(Time))


# ============================================================
# 6. Fit TDT for each individual
# ============================================================

tdt_clean <- T50_for_TDT %>%
  group_by(id) %>%
  group_modify(~ {
    
    mod <- lm(T50 ~ logTime, data = .x)
    
    tibble(
      T50_prime = coef(mod)[1],
      z = -coef(mod)[2],
      r2 = summary(mod)$r.squared
    )
    
  }) %>%
  ungroup()


# ============================================================
# 7. Add classic 30-min T50
# ============================================================

T50_30 <- T50_estimates %>%
  filter(Time == 30) %>%
  group_by(id) %>%
  summarise(
    T50 = mean(T50, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# 8. Final dataset
# ============================================================

tdt_clean <- tdt_clean %>%
  left_join(T50_30, by = "id") %>%
  mutate(
    
    # Negative z = biologically implausible TDT
    T50_prime = if_else(
      z < 0,
      NA_real_,
      T50_prime
    ),
    
    z = if_else(
      z < 0,
      NA_real_,
      z
    ),
    
    # 30-min T50 must be within measured temperature range
    T50 = if_else(
      T50 < temp_min | T50 > temp_max,
      NA_real_,
      T50
    )
  ) %>%
  select(
    id,
    z,
    T50_prime,
    T50
  )
library(dplyr)

# ---------------------------------------------------------------------------
# Sample code corrections identified during data cleaning of the snowgum
# faster/arrhenius/TDT/T50 datasets. Each mapping was confirmed against
# snowgum_metadata_source_pops.csv (ground truth) and/or field_notes.xlsx.
# See conversation notes for full justification of each fix.
# ---------------------------------------------------------------------------
code_fixes <- c(
  "S03T03P06" = "S05T03P06",   # site digit typo (03 -> 05)
  "S05P02T23" = "S05T02P23",   # P/T segment transposition
  "S09T02P12" = "S09T02P012",  # missing leading zero (this tree uses 3-digit padding for individuals 10-14)
  "S09T05P09" = "S07T05P09",   # site digit typo (09 -> 07)
  "S09T08P01" = "S09T08T01",   # metadata uses "T0x" not "P0x" for this tree (site 9, tree 8)
  "S12T03P07" = "S12T03P04",   # individual-number typo (P07 -> P04)
  "S14T03P04" = "S14T03P09",   # individual-number typo (P04 -> P09)
  "S05T02P09" = "S05T02P07",   # individual-number typo (P09 -> P07)
  "S9T07P01"  = "S09T07P01"    # missing leading zero on site number
)

tdt <- tdt_clean

# Keep the original id for an audit trail, then apply corrections
tdt <- tdt %>%
  mutate(
    id_original = id,
    id = recode(id, !!!code_fixes),
    id_was_corrected = id_original != id
  )

# Quick sanity check: report which fixes were actually applied
applied <- tdt %>% filter(id_was_corrected) %>% select(id_original, id)
print(applied)

write.csv(tdt, "./data/tdt_clean.csv", row.names = FALSE)

#write.csv(tdt_clean,"./data/tdt_clean.csv")


