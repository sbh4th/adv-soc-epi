# =============================================================================
# Create a trimmed ESS Round 11 dataset for course use
#
# Input:  data/ESS11e04_1.csv  (~80 MB, 600+ variables, all countries)
#         — this file is gitignored; download from https://ess-search.nsd.no
#           (ESS Round 11, Integrated file, Edition 4.1, CSV format)
#
# Output: data/ESS11_slim.csv  (~2–3 MB, ~37 variables, all countries)
#         — committed to git; suitable for course examples
#
# Variable selection rationale:
#   Admin/weights  — needed for survey-correct analyses
#   Demographics   — core sociodemographic covariates
#   Education      — ISCED harmonised scale (eisced)
#   Employment     — activity, occupation (ISCO-08), hours, contract
#   Income         — household decile (hinctnta), source, perceived adequacy
#   Health         — subjective health, activity limitation, happiness, life sat.
#   Behaviours     — smoking, alcohol (freq/quantity/binge), height, weight
#   Geography      — NUTS region
# =============================================================================

library(tidyverse)

in_path  <- here::here("data", "ESS11e04_1.csv")
out_path <- here::here("data", "ESS11_slim.csv")

keep_vars <- c(
  # --- admin / identifiers ---
  "name", "essround", "edition", "idno", "cntry",

  # --- survey weights ---
  "dweight",   # design weight
  "pspwght",   # post-stratification weight (use for most analyses)
  "pweight",   # population size weight (use when pooling countries)
  "anweight",  # analysis weight = pspwght × pweight

  # --- demographics ---
  "gndr",      # gender (1=male, 2=female)
  "agea",      # age in years
  "maritalb",  # legal marital status
  "domicil",   # urbanicity (big city → rural)
  "hhmmb",     # household size
  "brncntr",   # born in country of interview
  "ctzcntr",   # citizen of country of interview

  # --- education ---
  "eisced",    # ISCED harmonised education (1–7)

  # --- employment / occupation ---
  "mnactic",   # main activity last 7 days
  "pdwrk",     # paid work last 7 days (binary)
  "emplrel",   # employment relation (employee / self-employed)
  "wrkctra",   # work contract type
  "wkhtot",    # total hours worked per week
  "isco08",    # occupation (ISCO-08 major group)

  # --- income ---
  "hinctnta",  # household total net income decile (1=lowest, 10=highest)
  "hincsrca",  # main source of household income
  "hincfel",   # feeling about household income

  # --- health outcomes ---
  "health",    # subjective general health (1=very good … 5=very bad)
  "hlthhmp",   # hampered in daily activities by health (1=yes a lot … 3=no)
  "happy",     # how happy are you (0–10)
  "stflife",   # satisfaction with life as a whole (0–10)

  # --- health behaviours ---
  "cgtsmok",   # cigarette smoking frequency (1=every day … 5=never)
  "alcfreq",   # how often drink alcohol
  "alcwkdy",   # units on a weekday
  "alcwknd",   # units on a weekend day
  "alcbnge",   # binge drinking frequency
  "height",    # height in cm
  "weighta",   # weight in kg

  # --- geography ---
  "region"     # NUTS region code
)

raw <- read_csv(in_path, show_col_types = FALSE)

# Confirm all requested variables exist (warn about any that don't)
missing_vars <- setdiff(keep_vars, names(raw))
if (length(missing_vars) > 0) {
  warning("Variables not found in source file (skipped): ",
          paste(missing_vars, collapse = ", "))
}

slim <- raw |> select(any_of(keep_vars))

write_csv(slim, out_path)

# Report size reduction
in_mb  <- file.size(in_path)  / 1e6
out_mb <- file.size(out_path) / 1e6

message(sprintf(
  "Done.\n  Input : %s  (%.1f MB, %d vars, %d rows)\n  Output: %s  (%.1f MB, %d vars)\n  Reduction: %.0f%%",
  basename(in_path),  in_mb,  ncol(raw),  nrow(raw),
  basename(out_path), out_mb, ncol(slim),
  100 * (1 - out_mb / in_mb)
))
