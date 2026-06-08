# =============================================================================
# Concentration Index Decomposition of Smoking Inequality
# Data: European Social Survey Round 11 (ESS11e04_1), Finland
# Method: Wagstaff-van Doorslaer-Watanabe decomposition via rineq
#
# NOTE: Denmark (DK) is not included in ESS Round 11. Finland (FI) is used
# as it has the most complete cases among Nordic countries (n ≈ 957 after
# listwise deletion). Change COUNTRY below to switch (e.g. "NO", "SE").
#
# Outcome    : current smoking (cgtsmok 1–3 = smoker, 4–5 = non-smoker)
# Rank       : household income decile (hinctnta, 1 = lowest, 10 = highest)
# Covariates : age (agea), education (eisced 1–7), gender (gndr)
# =============================================================================

library(tidyverse)
library(rineq)
library(here)

# ── 0. Settings ───────────────────────────────────────────────────────────────

COUNTRY  <- "SE"     # Change to "NO", "SE", or "IS" as needed
DATA_PATH <- here("data", "ESS11_slim.csv")

# ── 1. Load & clean ───────────────────────────────────────────────────────────

raw <- read_csv(DATA_PATH, show_col_types = FALSE)

ess <- raw |>
  filter(cntry == COUNTRY) |>
  transmute(
    # Binary current smoking: 1 = smoker (every day / most days / less often)
    #   cgtsmok: 1=every day, 2=most days, 3=less often, 4=not now, 5=never
    #   6=refused, 7=don't know → NA (excluded by drop_na below)
    smoke    = case_when(
      cgtsmok %in% 1:3 ~ 1L,
      cgtsmok %in% 4:6 ~ 0L,
      TRUE             ~ NA_integer_
    ),
    cfsmoke  = case_when(
      cgtsmok %in% 1:4 ~ 1L,
      cgtsmok %in% 5:6 ~ 0L,
      TRUE             ~ NA_integer_
    ),
    # Income decile (1–10); 77/88/99 = refused / DK / N/A → set to NA
    inc_dec  = if_else(hinctnta %in% 1:10, 
      as.integer(hinctnta), NA_integer_),
    # Age in years (777 = refused → NA)
    age      = if_else(agea < 777, as.numeric(agea), NA_real_),
    # Education: ISCED 1–7 (55 = other, 77/88/99 → NA)
    educ     = if_else(eisced %in% 1:7, as.integer(eisced), NA_integer_),
    low_ed   = if_else(educ %in% 1:4, 1L, 
                if_else(educ %in% 5:7, 0L, NA_integer_)),
    # Gender: 1 = male, 2 = female → binary female indicator
    female   = if_else(gndr == 2, 1L, if_else(gndr == 1, 0L, NA_integer_)),
    alc = case_when(
      alcfreq %in% 1:2 ~ 1L,
      alcfreq %in% 3:6 ~ 2L,
      alcfreq %in% 7 ~ 3L,
      alcfreq %in% 77 ~ NA_integer_),
    no_alc = if_else(alc == 3, 1, 0),
    binge_mo = case_when(
      alcbnge %in% 1:3    ~ 1L,
      alcbnge %in% c(4:5, 6) ~ 0L,
      TRUE                ~ NA_integer_),
    bmi = weighta / (height/100)^2,
    hi_bmi = if_else(weighta / (height/100)^2 >=30, 1L, 0L),
    married = if_else(maritalb == 1, 1L,
                if_else(maritalb %in% 2:6, 0L, NA_integer_)),
    # Survey weight (post-stratification)
    wt       = pspwght
  )  |>
 drop_na()                    # listwise deletion for complete cases

message(glue::glue(
  "Country: {COUNTRY}  |  N = {nrow(ess)}  |  Smokers: {sum(ess$smoke)}  ",
  "({round(100*mean(ess$smoke),1)}%)"
))

# ── 2. Concentration curve ────────────────────────────────────────────────────

# ci() ranks by inc_dec; outcome = smoke; CIw = Wagstaff correction for binary
ci_obj <- ci(
  ineqvar = ess$inc_dec,
  outcome  = ess$smoke,
  weights  = ess$wt,
  type     = "CIw"
)

ci_val  <- ci_obj$concentration_index   # scalar CI estimate
ci_bounds <- confint(ci_obj)
ci_lo   <- ci_bounds[1]                 # 95% CI lower
ci_hi   <- ci_bounds[2]                 # 95% CI upper

# Build the concentration curve data manually:
#   rank respondents by income decile (weighted fractional rank)
#   then compute cumulative share of smoking vs. cumulative share of population
curve_data <- ess |>
  arrange(inc_dec) |>
  mutate(
    # Weighted fractional rank: cumulative weight up to midpoint of observation
    r = (cumsum(wt) - wt / 2) / sum(wt)
  ) |>
  arrange(r) |>
  mutate(
    cum_pop   = cumsum(wt) / sum(wt),
    cum_smoke = cumsum(smoke * wt) / sum(smoke * wt)
  ) |>
  bind_rows(tibble(cum_pop = 0, cum_smoke = 0, r = NA_real_)) |>
  arrange(cum_pop)

# CI annotation label
ci_label <- sprintf(
  "CIw = %.3f\n(95%% CI: %.3f, %.3f)",
  ci_val, ci_lo, ci_hi
)

plot_curve <- ggplot(curve_data, aes(x = cum_pop, y = cum_smoke)) +
  # Line of equality
  geom_abline(slope = 1, intercept = 0,
              colour = "grey50", linetype = "dashed", linewidth = 0.6) +
  # Concentration curve
  geom_line(colour = "#2166ac", linewidth = 1.1) +
  # Shade area between curve and diagonal
  geom_ribbon(aes(ymin = cum_smoke, ymax = cum_pop),
              fill = "#2166ac", alpha = 0.12) +
  # CI text annotation
  annotate("text", x = 0.2, y = 0.85,
           label = ci_label, hjust = 0, size = 4.2,
           colour = "#2166ac", fontface = "bold", family = "sans") +
  scale_x_continuous(labels = scales::percent_format(),
                     breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(
    title    = "Concentration Curve: Current Smoking by Income Decile",
    subtitle = glue::glue("ESS Round 11 · {COUNTRY} · n = {nrow(ess)}"),
    x        = "Cumulative share of population (ranked by income decile)",
    y        = "Cumulative share of smokers"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    aspect.ratio  = 1
  )

print(plot_curve)

# ── 2b. Smoking prevalence by income decile ───────────────────────────────────

smoke_by_dec <- ess |>
  group_by(inc_dec) |>
  summarise(
    pct_smoke = weighted.mean(smoke, wt) * 100,
    n         = n(),
    .groups   = "drop"
  ) |>
  mutate(totn = sum(n))

saveRDS(smoke_by_dec, here("data", "smoke-inc-se.rds"))

plot_smoke_dec <- ggplot(smoke_by_dec,
                         aes(x = factor(inc_dec), y = pct_smoke)) +
  geom_col(fill = "#2166ac", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct_smoke)),
            vjust = -0.4, size = 3.8) +
  scale_y_continuous(
    labels = scales::label_number(suffix = "%"),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title    = "Smoking Prevalence by Income Decile",
    subtitle = glue::glue("ESS Round 11 · {COUNTRY} · n = {nrow(ess)}"),
    x        = "Income decile (1 = lowest, 10 = highest)",
    y        = "% current smokers"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

print(plot_smoke_dec)

# ── 2c. Concentration curve for low education (ranked by income) ──────────────────
# Shows how education is distributed across the income distribution —
# a key input to the decomposition (the CI_k for education).

educ_ci_obj <- ci(
  ineqvar = ess$inc_dec,
  outcome  = ess$low_ed,
  weights  = ess$wt,
  type     = "CIw"
)

educ_ci_val    <- educ_ci_obj$concentration_index
educ_ci_bounds <- confint(educ_ci_obj)
educ_ci_lo     <- educ_ci_bounds[1]
educ_ci_hi     <- educ_ci_bounds[2]

educ_curve_data <- ess |>
  arrange(inc_dec) |>
  mutate(
    r = (cumsum(wt) - wt / 2) / sum(wt)
  ) |>
  arrange(r) |>
  mutate(
    cum_pop  = cumsum(wt) / sum(wt),
    cum_educ = cumsum(low_ed * wt) / sum(low_ed * wt)
  ) |>
  bind_rows(tibble(cum_pop = 0, cum_educ = 0, r = NA_real_)) |>
  arrange(cum_pop)

saveRDS(educ_curve_data, here("data", "rci-educ.rds"))

educ_ci_label <- sprintf(
  "CIw = %.3f\n(95%% CI: %.3f, %.3f)",
  educ_ci_val, educ_ci_lo, educ_ci_hi
)

plot_educ_curve <- ggplot(educ_curve_data, 
  aes(x = cum_pop, y = cum_educ)) +
  geom_abline(slope = 1, intercept = 0,
              colour = "grey50", linetype = "dashed", linewidth = 0.6) +
  geom_line(colour = "#d6604d", linewidth = 1.1) +
  geom_ribbon(aes(ymin = cum_pop, ymax = cum_educ),
              fill = "#d6604d", alpha = 0.12) +
  annotate("text", x = 0.05, y = 0.80,
           label = educ_ci_label, hjust = 0, size = 4.2,
           colour = "#d6604d", fontface = "bold", family = "sans") +
  scale_x_continuous(labels = scales::percent_format(),
                     breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(
    title    = "Concentration Curve: Education by Income Decile",
    subtitle = glue::glue("ESS Round 11 · {COUNTRY} · n = {nrow(ess)}"),
    x        = "Cumulative share of population (ranked by income decile)",
    y        = "Cumulative share of low education"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    aspect.ratio     = 1
  )

print(plot_educ_curve)

# ── 3. WDW decomposition ──────────────────────────────────────────────────────

# Linear probability model is conventional for WDW decomposition;
# rineq uses predicted values to compute the CI when type = "CIw".
# We include education as a factor to capture non-linearity across ISCED levels.

fit <- lm(
  smoke ~ age + female + low_ed + married + binge_mo + hi_bmi,
  data = ess,
  weights = wt
)

decomp <- contribution(
  object     = fit,
  ranker     = ess$inc_dec,
  correction = FALSE,       # no sign-correction: keep raw estimates
  type       = "CIw"
)

decomp_summary <- summary(decomp)

# ── 4. Build tidy decomposition table ─────────────────────────────────────────

# summary(decomp) returns a matrix; actual column names are:
#   "Contribution (%)", "Contribution (Abs)", "Elasticity",
#   "Concentration Index", "lower 5%", "upper 5%", "Corrected"
decomp_mat <- as.data.frame(decomp_summary) |>
  rownames_to_column("term") |>
  rename(
    pct_contrib  = `Contribution (%)`,
    abs_contrib  = `Contribution (Abs)`,
    elasticity   = Elasticity,
    ci_k         = `Concentration Index`,
    ci_k_lo      = `lower 5%`,
    ci_k_hi      = `upper 5%`
  )

# Collapse education factor levels into one "Education (ISCED)" row
decomp_tidy <- decomp_mat |>
  mutate(
    group = case_when(
      term == "age"              ~ "Age",
      term == "low_ed"           ~ "Low education",
      term == "married"           ~ "Married",
      term == "female"           ~ "Female",
      term == "binge_mo"         ~ "Binge drinker",
      term == "hi_bmi"           ~ "Obese",
      term == "residual"         ~ "Residual",
      TRUE                       ~ term
    )
  ) |>
  group_by(group) |>
  summarise(
    elasticity   = sum(elasticity,  na.rm = TRUE),
    # CI_k for collapsed factor groups is reported as the contribution-weighted mean;
    # for single-variable groups (age, female) this equals the variable's own CI.
    ci_k = if (sum(!is.na(ci_k)) == 0) NA_real_
      else sum(abs_contrib * replace_na(ci_k, 0), na.rm = TRUE) /
                        (sum(abs_contrib, na.rm = TRUE) + 1e-15),
    abs_contrib  = sum(abs_contrib, na.rm = TRUE),
    pct_contrib  = sum(pct_contrib, na.rm = TRUE),
    .groups = "drop"
  ) 

# ── 5. Plot decomposition bar chart ───────────────────────────────────────────

decomp_plot_data <- decomp_tidy |>
  filter(group != "Residual") |>
  mutate(
    group    = fct_reorder(group, pct_contrib),
    fill_col = if_else(pct_contrib >= 0, "pro-inequality", "pro-equality")
  )

plot_decomp <- ggplot(decomp_plot_data,
                      aes(x = pct_contrib, y = group, fill = fill_col)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", pct_contrib),
                hjust = if_else(pct_contrib >= 0, -0.15, 1.15)),
            size = 4.2) +
  scale_fill_manual(values = c("pro-inequality" = "#d6604d",
                                "pro-equality"   = "#4393c3")) +
  scale_x_continuous(labels = scales::label_number(suffix = "%"),
                     expand = expansion(mult = 0.2)) +
  labs(
    title    = "Decomposition of Smoking Inequality by Income",
    subtitle = glue::glue(
      "CIw = {round(ci_val, 3)}  |  ESS Round 11 · {COUNTRY} · n = {nrow(ess)}"
    ),
    x = "% contribution to CIw",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

print(plot_decomp)

# ── 6. Print decomposition table ──────────────────────────────────────────────

cat("\n── WDW Decomposition: Smoking Inequality by Income ──────────────────────\n")
cat(sprintf("Country: %s   N = %d   Smokers: %d (%.1f%%)\n",
            COUNTRY, nrow(ess), sum(ess$smoke), 100*mean(ess$smoke)))
cat(sprintf("Overall CIw = %.4f  (95%% CI: %.4f, %.4f)\n\n", ci_val, ci_lo, ci_hi))

decomp_tidy |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  rename(
    "Variable"       = group,
    "Elasticity"     = elasticity,
    "CI_k"           = ci_k,
    "Contribution"   = abs_contrib,
    "% of total"     = pct_contrib
  ) |>
  print(n = Inf)

saveRDS(decomp_tidy, here("data", "ess-decomp.rds"))

# ── 7. Save outputs ───────────────────────────────────────────────────────────

out_dir <- here::here("slides", "03-decomp", "figs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "ess11_concentration_curve.png"),
       plot_curve,      width = 6, height = 6, dpi = 200, bg = "white")

ggsave(file.path(out_dir, "ess11_decomp_bar.png"),
       plot_decomp,     width = 7, height = 4, dpi = 200, bg = "white")

ggsave(file.path(out_dir, "ess11_smoke_by_decile.png"),
       plot_smoke_dec,  width = 8, height = 5, dpi = 200, bg = "white")

ggsave(file.path(out_dir, "ess11_educ_curve.png"),
       plot_educ_curve, width = 6, height = 6, dpi = 200, bg = "white")

message("Figures saved to: ", out_dir)











library(broom)

predictors <- c("age", "female", "low_ed", "no_alc", "binge_mo", "hi_bmi", "married")  # etc.

results <- ess |>
  pivot_longer(cols = all_of(predictors), 
               names_to = "variable", 
               values_to = "value") |>
  group_by(variable) |>
  nest() |>
  mutate(
    # logistic for binary smoke
    mod_smoke = map(data, \(d) glm(smoke ~ value, data = d, family = binomial)),
    # linear for income decile
    mod_inc   = map(data, \(d) lm(inc_dec ~ value, data = d)),
    # tidy up
    tidy_smoke = map(mod_smoke, \(m) tidy(m, exponentiate = TRUE, conf.int = TRUE)),
    tidy_inc   = map(mod_inc,   \(m) tidy(m, conf.int = TRUE))
  ) |>
  select(variable, tidy_smoke, tidy_inc)

# Pull out just the predictor row (not intercept)
smoke_results <- results |>
  unnest(tidy_smoke) |>
  filter(term == "value")

inc_results <- results |>
  unnest(tidy_inc) |>
  filter(term == "value")
