# Code to read LI-6800 files
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(nls.multstart)
library(nlstools)
library(lubridate)
library(segmented)
library(tidylog)
library(readxl)
library(rTPC)

# Assumes the working directory is the project root (true by default when
# this project's .Rproj is opened in RStudio).
source("correct_6800_functions.R")

# load list of file paths for Licor data
at.files = list.files("licor_data", full.names = TRUE)

# extract list of IDs from the file names
at.ids <- lapply(at.files, function(x)strsplit(x, split = "[_]")[[1]][2]) #save the individual IDs, extracted from file names
at.ids <- sub("\\.xlsx$", "", at.ids)

# LOAD LICOR FILES
file.x <- list()
for (i in seq_along(at.files)) {
  x <- at.files[[i]]  # Take a single file name
  
  # Print the file name before processing
  print(paste("Processing file:", x))
  
  # Check if the file size is greater than 2000 bytes
  if (file.info(x)$size > 2000) {
    y <- tryCatch({
      read_6800_with_BLC(x)  # Read file
    }, error = function(e) {
      message(paste("Error reading file:", x))
      message(e)
      return(NULL)
    })
    
    # Only proceed if the file was read successfully
    if (!is.null(y)) {
      at.id <- at.ids[[i]]  # Take corresponding individual ID
      y$individualID <- rep(at.id, nrow(y))  # Label the individual file ID
      y$obs <- as.integer(y$obs)
      y$time <- as.double(y$time)
      
      file.x[[i]] <- y  # Add file (now a dataframe) to the list of files/dataframes
      print(paste("Successfully processed file:", x))
    } else {
      print(paste("Failed to process file:", x))
    }
  } else {
    print(paste("File size is less than or equal to 2000 bytes, skipping:", x))
  }
}

file.x <- file.x[lengths(file.x) != 0] # remove empty placeholders
at.df <- do.call(bind_rows, file.x) # combine files into 1 dataframe
at.df$individualID <- toupper(at.df$individualID )
at.df$curveID <- cumsum(at.df$obs == 1 & at.df$elapsed == 0)

# Add in the key for curveID and barcodes from photosynthesis group
at.meta = read.csv("./data/Faster_key_msfss_v2.csv")

at.meta$individualID <- toupper(paste0(at.meta$Population, at.meta$Individual))
library(dplyr)
library(dplyr)

# ---------------------------------------------------------------------------
# Final reconciliation of at.df$individualID vs at.meta$individualID
# ahead of merging faster (LICOR) data with field metadata.
#
# Run this top-to-bottom on a FRESH load of at.df / at.meta.
# ---------------------------------------------------------------------------

# --- 1. S02T06P01 (metadata only): delete ------------------------------------
# Field notes flagged this measurement for deletion -- drop it from at.meta.
at.meta <- at.meta %>% filter(individualID != "S02T06P01")

# --- 2. S05T05P12, S12T02P01 (LICOR only): keep as-is ------------------------
# No action needed -- these stay in at.df unchanged. They simply won't have
# a matching field metadata note, which is expected and fine. They will end
# up with NA metadata columns after the join.

# --- 3. S13T04P04 (metadata) vs S13T04P02 (LICOR): keep P02 ------------------
# Correct the metadata entry to match the LICOR-entered ID.
at.meta$individualID <- recode(at.meta$individualID,
                               "S13T04P04" = "S13T04P02"
)

# --- 4. S13T07P01: base + R2 + V2 variants, keep V2 only ---------------------
# at.df actually contains THREE variants: "S13T07P01" (base), "S13T07P01V2",
# and "S13T07P01R2". The R2 variant wasn't visible in the original setdiff()
# because it happened to already match at.meta's "R2" entry exactly -- a
# coincidental match, not evidence it's the right one to keep. Per decision:
# drop both the base and the R2 variant, keep only V2, and repoint the
# metadata note (currently "R2") to "V2" to match.
at.df <- at.df %>% filter(!individualID %in% c("S13T07P01", "S13T07P01R2"))
at.meta$individualID <- recode(at.meta$individualID,
                               "S13T07P01R2" = "S13T07P01V2"
)

# --- 5. S13T02P07 / S13T10P03: base + R2 variants, keep R2 only --------------
# Same pattern as S13T07P01 above. at.meta ONLY has an "R2" row for each --
# no base row exists in metadata at all. For S13T10P03 the base run in
# at.df is also clearly truncated (72 rows vs. R2's full 361), consistent
# with an aborted first attempt that was redone. Drop the base run from
# at.df for both, keep only R2 (matches what metadata actually documents;
# no recode needed since at.meta already uses "R2" naming for these).
at.df <- at.df %>% filter(!individualID %in% c("S13T02P07", "S13T10P03"))

# --- 6. Incomplete metadata entries: exclude ---------------------------------
# "S12", "S09", "", "-" have no tree/individual component -- drop from at.meta.
blank_ids <- c("S12", "S09", "", "-")
at.meta <- at.meta %>% filter(!individualID %in% blank_ids)

# --- Recheck remaining ID discrepancies after fixes --------------------------
setdiff(at.meta$individualID, at.df$individualID)
setdiff(at.df$individualID,  at.meta$individualID)
# Expect only "S05T05P12" and "S12T02P01" left in the at.df-only diff --
# these are the two you decided to keep unmatched. Everything else should
# now be empty.

# ---------------------------------------------------------------------------
# Date standardization (required before building a composite join key)
# ---------------------------------------------------------------------------
#
# at.meta$Date is a single YYYYMMDD value per row.
# at.df$date is a full "YYYYMMDD HH:MM:SS" timestamp, logged every ~5 sec.
#
# Confirmed: at.meta$Date has the same 2023-vs-2025 issue found earlier in
# field_notes.xlsx -- 20250130/20250131 are fine as-is, but 20230201-20230206
# are simply a year typo and should be 20250201-20250206 (month/day unchanged).

# Extract date-only from the raw timestamp
at.df <- at.df %>%
  mutate(date_only = substr(date, 1, 8))   # "20250130 11:02:03" -> "20250130"

at.meta <- at.meta %>%
  mutate(Date = case_when(
    Date == 20230201 ~ 20250201,
    Date == 20230202 ~ 20250202,
    Date == 20230203 ~ 20250203,
    Date == 20230204 ~ 20250204,
    Date == 20230205 ~ 20250205,
    Date == 20230206 ~ 20250206,
    TRUE ~ Date
  ))

# Quick sanity check: confirm the corrected at.meta dates actually show up
# in at.df for a few individuals from that batch.
at.df %>%
  filter(individualID %in% c("S14T03P09", "S06T02P29", "S13T10P03R2")) %>%
  count(individualID, date_only)

# --- 7. Duplicate individualIDs in at.meta -----------------------------------
# NOTE: at.df will have hundreds of rows per individualID -- that's expected
# (raw high-frequency LICOR logging, many points per curve), not an error.
# The concern is specifically at.meta having >1 row per individualID, since
# metadata should be one row per curve/measurement.
at.df   %>% count(individualID) %>% filter(n > 1)
at.meta %>% count(individualID) %>% filter(n > 1)

# Inspect what differs between the duplicated metadata rows -- turns out to
# be two different situations:
#  - S14T03P09 (20250130, Bad=1) and S13T08P02 (20250130, Bad=1) are each
#    duplicated by a single BAD-flagged row on 2025-01-30 -- consistent with
#    the same Jan-30 batch already known (from field_notes) to have several
#    measurements flagged for deletion. Drop these.
#  - S09T08P01 (20250202 & 20250203) has NO Bad flag on either row -- these
#    are two genuinely separate, valid measurements on different days. The
#    join_key (individualID + date) below disambiguates this correctly, so
#    no row needs to be dropped for this one.
at.meta %>% filter(individualID %in% c("S09T08P01", "S13T08P02", "S14T03P09"))

at.meta <- at.meta %>%
  filter(!(individualID %in% c("S14T03P09", "S13T08P02") & Date == 20250130 & Bad == 1))

# The metadata rows above are gone, but the corresponding raw curves are
# still sitting in at.df -- drop those too, or they'll just show up with
# blank metadata after the join instead of being excluded outright.
at.df <- at.df %>%
  filter(!(individualID %in% c("S14T03P09", "S13T08P02") & date_only == "20250130"))

# Conform S09T08P01 to metadata's canonical naming for this tree (site 9 /
# tree 8 uses "T0x" not "P0x" in snowgum_metadata_source_pops.csv -- same
# fix applied to check2_corrected.csv earlier).
at.meta$individualID <- recode(at.meta$individualID, "S09T08P01" = "S09T08T01")
at.df$individualID   <- recode(at.df$individualID,   "S09T08P01" = "S09T08T01")

# --- 8. S09T02P04: keep only 2025-02-06 --------------------------------------
# at.meta only has a metadata row for the 2025-02-06 measurement. The
# 2025-02-03 run (curve 46 in check2_corrected) has no corresponding field
# metadata at all -- decision: drop the 2025-02-03 raw rows, keep only the
# 2025-02-06 ones, which already have a metadata match.
at.df <- at.df %>%
  filter(!(individualID == "S09T02P04" & date_only == "20250203"))

# ---------------------------------------------------------------------------
# Build composite join key and join
# ---------------------------------------------------------------------------
#
# individualID alone is NOT a safe join key -- S09T08T01 was deliberately
# left with 2 legitimate metadata rows (2 real measurement dates). Joining
# on individualID alone would fan out every S09T08T01 raw row against BOTH
# metadata rows, silently doubling that individual's data (no NAs produced,
# so it wouldn't be caught by an is.na() check -- only a row-count mismatch
# reveals it). Always join on join_key, never plain individualID.
at.df   <- at.df   %>% mutate(join_key = paste(individualID, date_only, sep = "_"))
at.meta <- at.meta %>% mutate(join_key = paste(individualID, Date, sep = "_"))

at.df   %>% count(join_key) %>% filter(n > 1) %>% nrow()   # per-curve row counts still vary, fine
at.meta %>% count(join_key) %>% filter(n > 1)                # should now be empty

faster_merged <- at.df %>%
  left_join(at.meta, by = "join_key")

# Verify row count didn't change (i.e. no accidental many-to-many fan-out)
nrow(at.df)
nrow(faster_merged)   # should be IDENTICAL to nrow(at.df) -- any increase means fan-out

# Confirm which rows ended up with no metadata match -- should be exactly
# the two expected ones (S05T05P12 / S12T02P01), 361 raw rows each = 722 total
faster_merged %>%
  filter(is.na(Population)) %>%   # swap "Population" for whatever unique
  distinct(individualID.x)         # non-key column at.meta actually has

at.all <- faster_merged

at.all$`Î”Pcham` <- at.all$Î.Pcham

#at.all$Bad <- as.numeric(at.all$Bad)
at.all <- at.all %>% filter(is.na(Bad))

at.clean.trimmed <- at.all %>%
  group_by(curveID) %>%
  arrange(elapsed, .by_group = TRUE) %>%
  filter(row_number() <= n() - 5) %>%
  ungroup()

write.csv(at.clean.trimmed,"data/at.clean.trimmed.csv")

#source("code/schoolfield_with_vcov.R")
#source("code/fit_mod_schoolfield_4.R")
source("fit_mod_SS.R")
#source("code/test_aug11.R")

#colnames(at.clean.trimmed)
# fit_mod_schoolfield2() creates its own "SS_fits" folder (per-curve
# diagnostic plots) inside whatever the working directory is when it's
# called -- since we're staying in the project root throughout this
# script, that lands at ./SS_fits.
#source("code/fit_mod_schoolfield_updated_3.R")
Sfield.results.discard.hooks <- fit_mod_schoolfield2(
  Data = at.clean.trimmed,
  curve_id_col = "curveID",  # <- make sure this matches your data
  x = "Tleaf",
  y = "A",
  T_ref = 25
)
write.csv(Sfield.results.discard.hooks,"data/faster_parameters_SS.csv")
