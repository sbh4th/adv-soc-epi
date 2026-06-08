# =============================================================================
# Kitagawa-Blinder-Oaxaca Decomposition using the `oaxaca` package
# Data: European Social Survey Round 11, Italy (IT)
#
# This script replicates the hand-coded decomposition in ess11-kbo-decomp.R
# using the oaxaca package (Hlavac 2022). It is intended to show students
# how the package works and how its output maps to the manual results.
#
# Reference:
#   Hlavac, M. (2022). oaxaca: Blinder-Oaxaca Decomposition in R.
#   CRAN: https://cran.r-project.org/package=oaxaca
#
# Key formula syntax:
#   outcome ~ covariates | group_indicator
#
#   The group indicator must be binary (0/1). Observations where the
#   indicator = 1 form "Group 1" (the numerator of the gap).
#   Gap = mean(Y | group=1) - mean(Y | group=0)
#
# In our case: low_ed = 1 (low education, Group 1) vs low_ed = 0 (high ed)
#   => gap = mean(BMI | low_ed) - mean(BMI | high_ed) > 0
# =============================================================================

library(tidyverse)
library(oaxaca)
library(here)
library(glue)

# ── 0. Settings ───────────────────────────────────────────────────────────────

COUNTRY <- "IT"
DATA_PATH <- here("data", "ESS11_slim.csv")

# ── 1. Load & clean (identical to ess11-kbo-decomp.R) ────────────────────────

raw <- read_csv(DATA_PATH, show_col_types = FALSE)

ess_kbo <- raw |>
  filter(cntry == COUNTRY) |>
  transmute(
    weight_kg = if_else(weighta %in% c(777, 888, 999), NA_real_, as.numeric(weighta)),
    height_cm = if_else(height  %in% c(777, 888, 999), NA_real_, as.numeric(height)),
    bmi       = weight_kg / (height_cm / 100)^2,
    female    = if_else(gndr == 2, 1L, if_else(gndr == 1, 0L, NA_integer_)),
    age       = if_else(agea < 777, as.numeric(agea), NA_real_),
    educ      = if_else(eisced %in% 1:7, as.integer(eisced), NA_integer_),
    low_ed    = if_else(educ %in% 1:4, 1L, if_else(educ %in% 5:7, 0L, NA_integer_)),
    inc_dec   = if_else(hinctnta %in% 1:10, as.integer(hinctnta), NA_integer_),
    smoke     = case_when(
      cgtsmok %in% 1:3 ~ 1L,
      cgtsmok %in% 4:6 ~ 0L,
      TRUE             ~ NA_integer_
    ),
    binge_mo  = case_when(
      alcbnge %in% 1:3       ~ 1L,
      alcbnge %in% c(4:5, 6) ~ 0L,
      TRUE                   ~ NA_integer_
    ),
    married   = if_else(maritalb == 1, 1L, if_else(maritalb %in% 2:6, 0L, NA_integer_)),
    wt        = pspwght
  ) |>
  filter(bmi > 10, bmi < 60) |>
  drop_na(bmi, low_ed, age, female, smoke, binge_mo, married, inc_dec, wt)

message(glue("N = {nrow(ess_kbo)} | Low ed: {sum(ess_kbo$low_ed)} | High ed: {sum(1 - ess_kbo$low_ed)}"))

# ── 2. Run oaxaca decomposition ───────────────────────────────────────────────
#
# The `weight` argument controls which group's betas anchor the endowments
# component. You can pass a vector to get multiple decompositions at once:
#
#   weight = 1     → Group 1 (low_ed) betas for endowments  = our Col1
#   weight = 0     → Group 0 (high_ed) betas for endowments = our Col2
#   weight = 0.5   → Reimers (1983): simple average of betas
#   weight = "Neumark" → Neumark (1988): pooled-regression betas
#   weight = "Cotton"  → Cotton (1988): group-size-weighted betas
#
# R = number of bootstrap replicates for standard errors.

set.seed(42)

result <- oaxaca(
  formula = bmi ~ age + female + smoke + binge_mo + married + inc_dec | low_ed,
  data    = ess_kbo,
  weights = ess_kbo$wt,
  R       = 500
)

# ── 3. Inspect the output structure ──────────────────────────────────────────
#
# result$y            : group means and gap with bootstrap SEs
# result$threefold    : unique threefold decomposition (E, C, I) — same
#                       regardless of reference group choice
# result$twofold      : twofold for each weight option (E + C, no I)
# result$reg          : regression objects for each group and pooled

cat("\n── Object structure ─────────────────────────────────────────────────────\n")
str(result, max.level = 2)

# ── 4. Overall results ────────────────────────────────────────────────────────

cat("\n── Group means and gap ──────────────────────────────────────────────────\n")
print(result$y)

# The threefold decomposition is unique — E, C, I have one value each
cat("\n── Threefold decomposition (unique) ─────────────────────────────────────\n")
print(result$threefold$overall)

# ── 5. Compare to hand-coded results ─────────────────────────────────────────
#
# Our hand-coded Col2 (high_ed / Group 0 betas):
#   E2 = (Xbar_2 - Xbar_1)' * beta_0
#   C2 = Xbar_0' * (beta_1 - beta_0)     [sign: gap = ybar_1 - ybar_0]
#   I2 = (Xbar_2 - Xbar_1)' * (beta_1 - beta_0)
#
# The threefold in oaxaca should match these exactly.

cat("\n── Per-variable threefold contributions ─────────────────────────────────\n")
print(result$threefold$variables)

# ── 6. Twofold decompositions (multiple weight options) ───────────────────────
#
# The twofold collapses I into either E or C depending on the weight:
#   weight = 1  (Group 1 betas):  gap = E_twofold + C_twofold  [no I]
#   weight = 0  (Group 0 betas):  gap = E_twofold + C_twofold  [no I]
#
# These correspond to our Col1 (w=1) and Col2 (w=0) E+C combined (I=0).

result2 <- oaxaca(
  formula = bmi ~ age + female + smoke + binge_mo + married + inc_dec | low_ed,
  data    = ess_kbo,
  weights = ess_kbo$wt,
  R       = 500,
  weight  = c(1, 0, "Neumark", "Cotton")
)

cat("\n── Twofold overall results (four reference choices) ─────────────────────\n")
print(result2$twofold$overall)

cat("\n── Twofold per-variable contributions ───────────────────────────────────\n")
print(result2$twofold$variables)

# ── 7. Built-in plot ──────────────────────────────────────────────────────────
#
# oaxaca has a plot method that shows the decomposition graphically.
# The `components` argument selects which parts to display.

plot(result,
     components     = c("endowments", "coefficients", "interaction"),
     component.labels = c("Endowments\n(characteristics)",
                          "Coefficients\n(returns)",
                          "Interaction"),
     group.labels   = c("High ed", "Low ed"),
     variable.labels = c(
       age      = "Age",
       female   = "Female",
       smoke    = "Smoker",
       binge_mo = "Binge drinker",
       married  = "Married",
       inc_dec  = "Income decile"
     ))

# ── 8. Extract to a tidy tibble ───────────────────────────────────────────────
#
# If you want to use oaxaca results in ggplot or tinytable, you need to
# reshape the output. The per-variable contributions live in a matrix:
#   result$threefold$variables  — rows = variables, cols = E / C / I (each
#                                  with estimate and SE)

threefold_tidy <- as_tibble(result$threefold$variables, rownames = "variable") |>
  rename_with(\(x) str_replace_all(x, " ", "_")) |>
  select(variable,
         E_est = endowments_coef,  E_se = endowments_se,
         C_est = coefficients_coef, C_se = coefficients_se,
         I_est = interaction_coef,  I_se = interaction_se)

cat("\n── Tidy threefold contributions ─────────────────────────────────────────\n")
print(threefold_tidy)

# ── 9. Key takeaways for students ────────────────────────────────────────────
#
# 1. oaxaca() does in one call what took ~100 lines by hand:
#    bootstrap SEs, multiple reference groups, threefold + twofold.
#
# 2. The THREEFOLD (E, C, I) is unique — it does not depend on which
#    group's betas you use as a reference. This is a property of the algebra.
#
# 3. The TWOFOLD folds the interaction term into either E or C depending
#    on the reference group — which is why endowment estimates differ
#    across our three table columns.
#
# 4. The `weight` argument only affects the twofold decomposition. For
#    the threefold, all weight choices give the same E, C, I.
#
# 5. There is no universally "correct" reference group. Showing multiple
#    columns (as in our table) is good practice — it reveals how sensitive
#    the decomposition is to that choice.
