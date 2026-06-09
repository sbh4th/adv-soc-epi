# =============================================================================
# Kitagawa-Blinder-Oaxaca Decomposition — General Template
# Data: European Social Survey Round 11
# Method: Twofold (Blinder-Oaxaca) decomposition, manual implementation
#
# Edit the SETTINGS block below to change country, outcome, group, covariates.
# Everything else runs automatically.
# =============================================================================

library(tidyverse)
library(here)
library(glue)

# ── 0. SETTINGS — edit here ───────────────────────────────────────────────────

COUNTRY    <- "IT"          # ESS country code

OUTCOME    <- "bmi"         # continuous outcome variable name (must exist in ess_kbo)
OUTCOME_LABEL <- "BMI (kg/m²)"

GROUP_VAR  <- "low_ed"      # binary grouping variable (0 = group_1, 1 = group_2)
GROUP_1    <- "High ed"     # label for group coded 0
GROUP_2    <- "Low ed"      # label for group coded 1
GROUP_COLORS <- c("#2166ac", "#d6604d")  # colours for group 1, group 2

COVARIATES <- c(            # covariates for the decomposition model
  "age", "female", "smoke", "binge_mo", "married", "inc_dec"
)
COVARIATE_LABELS <- c(      # display labels (must match COVARIATES order)
  "Age (years)",
  "Female",
  "Current smoker",
  "Binge drinking (monthly+)",
  "Married",
  "Income (decile)"
)

# ── 1. Load & clean ───────────────────────────────────────────────────────────

DATA_PATH <- here("data", "ESS11_slim.csv")
out_dir   <- here("slides", "03-decomp", "figs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
  drop_na(all_of(c(OUTCOME, GROUP_VAR, COVARIATES, "wt")))

# Named colour/label vectors keyed on group labels
grp_colors <- setNames(GROUP_COLORS, c(GROUP_1, GROUP_2))
grp_shapes <- setNames(c(16, 17),    c(GROUP_1, GROUP_2))
cov_labels <- setNames(COVARIATE_LABELS, COVARIATES)

# Helper: factor-ise the grouping variable for plots
as_group <- function(x) factor(x, levels = c(0, 1), labels = c(GROUP_1, GROUP_2))

message(glue(
  "Country: {COUNTRY}  |  N = {nrow(ess_kbo)}  |  ",
  "{GROUP_1}: {sum(ess_kbo[[GROUP_VAR]] == 0)}  ",
  "{GROUP_2}: {sum(ess_kbo[[GROUP_VAR]] == 1)}"
))
message(glue(
  "Mean {OUTCOME}  {GROUP_1}: ",
  "{round(weighted.mean(ess_kbo[[OUTCOME]][ess_kbo[[GROUP_VAR]]==0], ess_kbo$wt[ess_kbo[[GROUP_VAR]]==0]), 2)}  |  ",
  "{GROUP_2}: ",
  "{round(weighted.mean(ess_kbo[[OUTCOME]][ess_kbo[[GROUP_VAR]]==1], ess_kbo$wt[ess_kbo[[GROUP_VAR]]==1]), 2)}"
))

# ── 2. Outcome distribution by group ─────────────────────────────────────────

mean_df <- ess_kbo |>
  mutate(grp = as_group(.data[[GROUP_VAR]])) |>
  group_by(grp) |>
  summarise(mean_out = weighted.mean(.data[[OUTCOME]], wt), .groups = "drop")

plot_density <- ess_kbo |>
  mutate(grp = as_group(.data[[GROUP_VAR]])) |>
  ggplot(aes(x = .data[[OUTCOME]], fill = grp, colour = grp)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  geom_vline(data = mean_df,
             aes(xintercept = mean_out, colour = grp),
             linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values = grp_colors) +
  scale_colour_manual(values = grp_colors) +
  labs(
    title    = glue("Distribution of {OUTCOME_LABEL} by {GROUP_VAR}"),
    subtitle = glue("ESS Round 11 · {COUNTRY} · n = {nrow(ess_kbo)}"),
    x        = OUTCOME_LABEL, y = "Density",
    fill = NULL, colour = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "top")

print(plot_density)

# ── 3. Covariate means by group ───────────────────────────────────────────────

means_df <- ess_kbo |>
  group_by(grp = as_group(.data[[GROUP_VAR]])) |>
  summarise(across(all_of(COVARIATES), \(x) weighted.mean(x, wt)),
            .groups = "drop") |>
  pivot_longer(-grp, names_to = "variable", values_to = "mean") |>
  mutate(variable = factor(variable, levels = COVARIATES, labels = cov_labels))

pct_vars <- c("Female", "Current smoker", "Binge drinking (monthly+)", "Married")

plot_df <- means_df |>
    mutate(
      mean     = if_else(variable %in% pct_vars, mean * 100, mean),
      variable = if_else(variable %in% pct_vars,
                         paste0(variable, " (%)"), variable)
    )

plot_means <- ggplot(plot_df,
                     aes(x = mean, y = variable, colour = grp, shape = grp)) +
  geom_line(aes(group = variable), colour = "grey70", linewidth = 0.6) +
  geom_point(size = 4) +
  scale_colour_manual(values = grp_colors) +
  scale_shape_manual(values  = grp_shapes) +
  labs(
    title    = glue("Covariate Means by {GROUP_VAR}"),
    subtitle = glue("ESS Round 11 · {COUNTRY}"),
    x = "Weighted mean", y = NULL, colour = NULL, shape = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "top")

print(plot_means)

# ── 4. Regression coefficients by group ──────────────────────────────────────

grp1 <- ess_kbo |> filter(.data[[GROUP_VAR]] == 0)
grp2 <- ess_kbo |> filter(.data[[GROUP_VAR]] == 1)

fmla <- reformulate(COVARIATES, response = OUTCOME)

fit_1 <- lm(fmla, data = grp1, weights = wt)
fit_2 <- lm(fmla, data = grp2, weights = wt)

coef_df <- bind_rows(
  broom::tidy(fit_1, conf.int = TRUE) |> mutate(grp = GROUP_1),
  broom::tidy(fit_2, conf.int = TRUE) |> mutate(grp = GROUP_2)
) |>
  filter(term != "(Intercept)") |>
  mutate(
    grp  = factor(grp, levels = c(GROUP_1, GROUP_2)),
    term = factor(term, levels = COVARIATES, labels = cov_labels)
  )

plot_coefs <- ggplot(coef_df,
                     aes(x = estimate, y = term, colour = grp, shape = grp)) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed",
             linewidth = 0.5) +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high),
                 position = position_dodge(width = 0.5), linewidth = 0.6) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = grp_colors) +
  scale_shape_manual(values  = grp_shapes) +
  labs(
    title    = glue("Regression Coefficients by {GROUP_VAR}"),
    subtitle = glue("Outcome: {OUTCOME_LABEL} · ESS Round 11 · {COUNTRY}"),
    x = glue("Coefficient ({OUTCOME_LABEL})"), y = NULL,
    colour = NULL, shape = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "top")

print(plot_coefs)

# ── 5. Three-reference-group KBO decomposition + bootstrap SEs ───────────────
#
# Three versions, each a valid decomposition of the same gap:
#
#   Col1 (GROUP_2 betas for E, GROUP_2 means for C):
#     E1 = (Xbar_2 - Xbar_1)' * beta_2
#     C1 = Xbar_2' * (beta_2 - beta_1)
#     I1 = -(Xbar_2 - Xbar_1)' * (beta_2 - beta_1)   [= -I_threefold]
#
#   Col2 (GROUP_1 betas for E, GROUP_1 means for C):
#     E2 = (Xbar_2 - Xbar_1)' * beta_1
#     C2 = Xbar_1' * (beta_2 - beta_1)
#     I2 = +(Xbar_2 - Xbar_1)' * (beta_2 - beta_1)   [= +I_threefold]
#
#   Col3 (pooled betas for E, Cotton (1988) C — no interaction):
#     E3 = (Xbar_2 - Xbar_1)' * beta_pooled
#     C3 = Xbar_2'*(beta_2 - beta_pooled) + Xbar_1'*(beta_pooled - beta_1)
#
# Intercept contributes only to C (db_int = beta_2_int - beta_1_int) and
# is identical across all three columns.

# Pooled regression (full sample, no group indicator)
fit_pooled <- lm(fmla, data = ess_kbo, weights = wt)

# ── Core decomposition function ───────────────────────────────────────────────
decompose_kbo <- function(dat, group_var, outcome, covariates, wt_var = "wt") {
  g1 <- dat[dat[[group_var]] == 0, ]
  g2 <- dat[dat[[group_var]] == 1, ]

  fm <- reformulate(covariates, response = outcome)
  b1 <- suppressWarnings(coef(lm(fm, data = g1,  weights = g1[[wt_var]])))
  b2 <- suppressWarnings(coef(lm(fm, data = g2,  weights = g2[[wt_var]])))
  bp <- suppressWarnings(coef(lm(fm, data = dat, weights = dat[[wt_var]])))

  wm <- function(d, v) weighted.mean(d[[v]], d[[wt_var]])

  x1_no <- sapply(covariates, \(v) wm(g1, v))   # group1 covariate means
  x2_no <- sapply(covariates, \(v) wm(g2, v))   # group2 covariate means

  ybar1 <- wm(g1, outcome)
  ybar2 <- wm(g2, outcome)
  gap   <- ybar2 - ybar1

  dX  <- x2_no - x1_no                # Xbar_2 - Xbar_1 (no intercept)
  db_no  <- (b2 - b1)[-1]             # beta_2 - beta_1 (no intercept)
  db_int <- (b2 - b1)[["(Intercept)"]] # intercept difference (scalar)

  # Col1: GROUP_2 betas for E, GROUP_2 means for C
  E1_v   <- dX * b2[-1]
  C1_v   <- x2_no * db_no
  C1_int <- db_int
  I1     <- -sum(dX * db_no)          # = -I_threefold

  # Col2: GROUP_1 betas for E, GROUP_1 means for C
  E2_v   <- dX * b1[-1]
  C2_v   <- x1_no * db_no
  C2_int <- db_int                    # same intercept term
  I2     <- sum(dX * db_no)           # = +I_threefold

  # Col3: pooled betas for E, Cotton C (no interaction by construction)
  E3_v   <- dX * bp[-1]
  C3_v   <- x2_no * (b2[-1] - bp[-1]) + x1_no * (bp[-1] - b1[-1])
  C3_int <- db_int                    # b2_int - bp_int + bp_int - b1_int = db_int

  list(
    ybar1 = ybar1, ybar2 = ybar2, gap = gap,
    # Totals
    E1 = sum(E1_v), C1 = sum(C1_v) + C1_int, I1 = I1,
    E2 = sum(E2_v), C2 = sum(C2_v) + C2_int, I2 = I2,
    E3 = sum(E3_v), C3 = sum(C3_v) + C3_int,
    # Per-variable (named vectors; intercept appended to C)
    E1_v = E1_v, C1_v = c(C1_v, `(Intercept)` = C1_int),
    E2_v = E2_v, C2_v = c(C2_v, `(Intercept)` = C2_int),
    E3_v = E3_v, C3_v = c(C3_v, `(Intercept)` = C3_int)
  )
}

# Point estimates
est <- decompose_kbo(ess_kbo, GROUP_VAR, OUTCOME, COVARIATES)

cat(glue("
── KBO Decomposition (three reference groups) ───────────────────────────────
  Country: {COUNTRY}  N = {nrow(ess_kbo)}
  Mean {GROUP_1}: {round(est$ybar1, 3)}   Mean {GROUP_2}: {round(est$ybar2, 3)}
  Gap: {round(est$gap, 3)}

  Version        E         C         I         Sum
  {GROUP_2} betas  {round(est$E1,3)}  {round(est$C1,3)}  {round(est$I1,3)}  {round(est$E1+est$C1+est$I1,3)}
  {GROUP_1} betas  {round(est$E2,3)}  {round(est$C2,3)}  {round(est$I2,3)}  {round(est$E2+est$C2+est$I2,3)}
  Pooled       {round(est$E3,3)}  {round(est$C3,3)}  NA        {round(est$E3+est$C3,3)}
─────────────────────────────────────────────────────────────────────────────
"))

# ── Bootstrap SEs (B = 500) ───────────────────────────────────────────────────
set.seed(42)
B <- 500
message(glue("Running {B} bootstrap replicates..."))

boot_list <- lapply(seq_len(B), \(i) {
  idx <- sample(nrow(ess_kbo), replace = TRUE)
  tryCatch(
    decompose_kbo(ess_kbo[idx, ], GROUP_VAR, OUTCOME, COVARIATES),
    error = \(e) NULL
  )
}) |> Filter(Negate(is.null), x = _)

message(glue("Bootstrap complete: {length(boot_list)} valid replicates"))

bse <- function(key) sd(sapply(boot_list, \(r) r[[key]]), na.rm = TRUE)
bse_v <- function(key, nm) sd(sapply(boot_list, \(r) r[[key]][[nm]]), na.rm = TRUE)

# Scalars
se <- list(
  ybar1 = bse("ybar1"), ybar2 = bse("ybar2"), gap = bse("gap"),
  E1 = bse("E1"), C1 = bse("C1"), I1 = bse("I1"),
  E2 = bse("E2"), C2 = bse("C2"), I2 = bse("I2"),
  E3 = bse("E3"), C3 = bse("C3")
)

# Per-variable SEs
cov_and_int <- c(COVARIATES, "(Intercept)")
se_E1_v <- sapply(COVARIATES,   \(v) bse_v("E1_v", v))
se_C1_v <- sapply(cov_and_int,  \(v) bse_v("C1_v", v))
se_E2_v <- sapply(COVARIATES,   \(v) bse_v("E2_v", v))
se_C2_v <- sapply(cov_and_int,  \(v) bse_v("C2_v", v))
se_E3_v <- sapply(COVARIATES,   \(v) bse_v("E3_v", v))
se_C3_v <- sapply(cov_and_int,  \(v) bse_v("C3_v", v))

# ── Build the display table ───────────────────────────────────────────────────
# Row layout mirrors the image: overall | endowments | coefficients | interaction
# Columns: label | G2_est G2_se | G1_est G1_se | Pooled_est Pooled_se

cov_display  <- cov_labels                         # COVARIATES → display names
int_display  <- c("(Intercept)" = "Intercept")
all_display  <- c(cov_display, int_display)

make_rows <- function(label, section, bold,
                      e1, se1, e2, se2, ep, sep_) {
  tibble(label = label, section = section, bold = bold,
         g2_est = e1, g2_se = se1,
         g1_est = e2, g1_se = se2,
         gp_est = ep, gp_se = sep_)
}

overall_rows <- make_rows(
  label   = c(GROUP_1, GROUP_2, "Difference"),
  section = "overall", bold = FALSE,
  e1  = c(est$ybar1, est$ybar2, est$gap),
  se1 = c(se$ybar1,  se$ybar2,  se$gap),
  e2  = c(est$ybar1, est$ybar2, est$gap),
  se2 = c(se$ybar1,  se$ybar2,  se$gap),
  ep  = c(est$ybar1, est$ybar2, est$gap),
  sep_ = c(se$ybar1, se$ybar2, se$gap)
)

endow_rows <- bind_rows(
  make_rows("Endowments", "endow_tot", TRUE,
            est$E1, se$E1, est$E2, se$E2, est$E3, se$E3),
  make_rows(unname(cov_display), "endow_var", FALSE,
            est$E1_v[COVARIATES], se_E1_v,
            est$E2_v[COVARIATES], se_E2_v,
            est$E3_v[COVARIATES], se_E3_v)
)

coef_rows <- bind_rows(
  make_rows("Coefficients", "coef_tot", TRUE,
            est$C1, se$C1, est$C2, se$C2, est$C3, se$C3),
  make_rows(unname(all_display), "coef_var", FALSE,
            est$C1_v[names(all_display)], se_C1_v,
            est$C2_v[names(all_display)], se_C2_v,
            est$C3_v[names(all_display)], se_C3_v)
)

inter_row <- make_rows(
  "Interaction", "inter", TRUE,
  est$I1, se$I1, est$I2, se$I2,
  NA_real_, NA_real_          # no interaction for pooled
)

decomp_big_tbl <- bind_rows(overall_rows, 
  endow_rows, coef_rows, inter_row)

# Convenience scalars retained for QMD inline text
ybar_1      <- est$ybar1;  
ybar_2     <- est$ybar2;  
gap        <- est$gap
endow_total <- est$E2;     
coef_total <- est$C2;     
inter_total <- est$I2

# Per-variable endowments (col2 = group1 betas) for bar chart
endow_detail <- tibble(
  variable     = COVARIATES,
  contribution = est$E2_v[COVARIATES]
) |>
  mutate(
    variable = factor(variable, levels = COVARIATES, 
    labels = cov_labels),
    pct      = 100 * contribution / gap
  )

# ── 6. Decomposition plots ────────────────────────────────────────────────────

# Summary: total gap → endowments + coefficients
summary_df <- tibble(
  component = factor(
    c("Total gap", "Endowments\n(characteristics)", "Coefficients\n(returns)"),
    levels = c("Total gap", "Endowments\n(characteristics)", "Coefficients\n(returns)")
  ),
  value = c(gap, endow_total, coef_total),
  fill  = c("gap", "endow", "coef")
)

plot_decomp_summary <- ggplot(summary_df,
  aes(x = component, y = value, fill = fill)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_text(aes(
    label = sprintf("%+.2f\n(%+.1f%%)", value, 100 * value / gap),
    vjust = if_else(value >= 0, -0.3, 1.2)
  ), size = 4.2) +
  scale_fill_manual(values = c(gap = "grey40", endow = GROUP_COLORS[1],
                               coef = GROUP_COLORS[2])) +
  scale_y_continuous(expand = expansion(mult = 0.25)) +
  labs(
    title    = glue("KBO Decomposition: {OUTCOME_LABEL} gap ({GROUP_2} − {GROUP_1})"),
    subtitle = glue("ESS Round 11 · {COUNTRY} · Gap = {round(gap, 2)}"),
    x = NULL, y = OUTCOME_LABEL
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank())

print(plot_decomp_summary)

# Detail: per-variable endowments contributions
plot_decomp_detail <- ggplot(endow_detail,
  aes(x = contribution, y = fct_reorder(variable, contribution),
    fill = contribution >= 0)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_text(aes(
    label = sprintf("%+.3f", contribution),
    hjust = if_else(contribution >= 0, -0.15, 1.15)
  ), size = 4) +
  scale_fill_manual(values = c("TRUE" = GROUP_COLORS[1], "FALSE" = GROUP_COLORS[2])) +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  labs(
    title    = "Endowments Component: Per-Variable Contributions",
    subtitle = glue(
      "Total endowments = {round(endow_total, 2)}  ",
      "({round(100*endow_total/gap, 1)}% of gap)"
    ),
    x = glue("Contribution to {OUTCOME_LABEL} gap"), y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank())

print(plot_decomp_detail)

# ── 7. Save figures ───────────────────────────────────────────────────────────

slug <- glue("{COUNTRY}_{OUTCOME}_{GROUP_VAR}")

ggsave(file.path(out_dir, glue("ess11_{slug}_density.png")),
       plot_density,       width = 8, height = 5, dpi = 200, bg = "white")
ggsave(file.path(out_dir, glue("ess11_{slug}_means.png")),
       plot_means,         width = 7, height = 4, dpi = 200, bg = "white")
ggsave(file.path(out_dir, glue("ess11_{slug}_coefs.png")),
       plot_coefs,         width = 7, height = 4, dpi = 200, bg = "white")
ggsave(file.path(out_dir, glue("ess11_{slug}_decomp_summary.png")),
       plot_decomp_summary, width = 7, height = 5, dpi = 200, bg = "white")
ggsave(file.path(out_dir, glue("ess11_{slug}_decomp_detail.png")),
       plot_decomp_detail,  width = 7, height = 4, dpi = 200, bg = "white")

message("Figures saved with prefix: ess11_", slug)

# ── 8. Save objects for QMD presentation ─────────────────────────────────────

rds_dir <- here("slides", "03-decomp", "data")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# Summary scalars
kbo_summary <- tibble(
  country      = COUNTRY,
  outcome      = OUTCOME,
  outcome_label = OUTCOME_LABEL,
  group_var    = GROUP_VAR,
  group_1      = GROUP_1,
  group_2      = GROUP_2,
  ybar_1       = ybar_1,
  ybar_2       = ybar_2,
  gap          = gap,
  endow_total  = endow_total,
  coef_total   = coef_total,
  inter_total  = inter_total
)

saveRDS(kbo_summary,    file.path(rds_dir, "kbo_summary.rds"))
saveRDS(decomp_big_tbl, file.path(rds_dir, "kbo_decomp_big_tbl.rds"))
saveRDS(means_df,       file.path(rds_dir, "kbo_means_df.rds"))
saveRDS(coef_df,        file.path(rds_dir, "kbo_coef_df.rds"))
saveRDS(endow_detail,   file.path(rds_dir, "kbo_endow_detail.rds"))
saveRDS(ess_kbo,        file.path(rds_dir, "kbo_ess_kbo.rds"))

# Save settings list so QMD can rebuild plots with correct labels/colours
kbo_settings <- list(
  COUNTRY        = COUNTRY,
  OUTCOME        = OUTCOME,
  OUTCOME_LABEL  = OUTCOME_LABEL,
  GROUP_VAR      = GROUP_VAR,
  GROUP_1        = GROUP_1,
  GROUP_2        = GROUP_2,
  GROUP_COLORS   = GROUP_COLORS,
  COVARIATES     = COVARIATES,
  COVARIATE_LABELS = COVARIATE_LABELS,
  grp_colors     = grp_colors,
  grp_shapes     = grp_shapes,
  cov_labels     = cov_labels
)
saveRDS(kbo_settings, file.path(rds_dir, "kbo_settings.rds"))

message("RDS objects saved to: ", rds_dir)



