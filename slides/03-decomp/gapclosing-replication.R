## -------------------------------------------------------
## Manual replication of Lundberg (2021) gap-closing estimand
## Using the same DGP as gapclosing::generate_simulated_data()
## -------------------------------------------------------

library(tidyverse)

set.seed(08544)
n <- 1000

# -------------------------------------------------------
# 1. DATA GENERATING PROCESS
# -------------------------------------------------------
# category  X  ~ Uniform({A,B,C})
# confounder L  ~ N(mu_X, 1),  mu_A=-1.5, mu_B=0, mu_C=1.5
# treatment  T  ~ Bernoulli(plogis(L))   [confounded by L]
# outcome    Y  ~ N(L + T*tau_X, 1),     tau_A=1.3, tau_B=0, tau_C=-1.2
#
# Key feature: treatment is positively associated with L,
# and L's mean differs by category, so naively comparing
# treated vs. untreated mixes confounding with the effect.

dat <- data.frame(category = sample(c("A","B","C"), n, replace = TRUE)) |>
  mutate(
    confounder = rnorm(n,
      mean = case_when(category == "A" ~ -1.5,
                       category == "B" ~  0,
                       category == "C" ~  1.5),
      sd = 1),
    treatment = rbinom(n, 1, plogis(confounder)),
    tau       = case_when(category == "A" ~  1.3,
                          category == "B" ~  0,
                          category == "C" ~ -1.2),
    outcome   = rnorm(n, mean = confounder + treatment * tau, sd = 1)
  ) |>
  select(-tau)


# -------------------------------------------------------
# 2. FACTUAL GAP  (observed means by category)
# -------------------------------------------------------
factual <- dat |>
  group_by(category) |>
  summarise(mean_y = mean(outcome), .groups = "drop")

factual
# Gap A-B, A-C, B-C in observed data
factual_gap_AB <- factual$mean_y[factual$category=="A"] -
                  factual$mean_y[factual$category=="B"]
cat("Factual gap A - B:", round(factual_gap_AB, 3), "\n")


# -------------------------------------------------------
# 3. COUNTERFACTUAL: set T = 1 for everyone
#
# Identification assumption: no unmeasured confounders of
# T -> Y within levels of (category, confounder), i.e.
# Y(t) _||_ T | {category, confounder}.
# -------------------------------------------------------


## -------------------------------------------------------
## 3a. OUTCOME MODEL (g-computation / standardization)
## -------------------------------------------------------
# Fit an outcome model on observed data
m_out <- lm(outcome ~ confounder + category * treatment, data = dat)

# Predict each person's potential outcome under T=1
dat_t1 <- dat |> mutate(treatment = 1)
dat$y_hat_t1 <- predict(m_out, newdata = dat_t1)

# Counterfactual mean per category = average predicted Y(T=1)
cf_om <- dat |>
  group_by(category) |>
  summarise(cf_mean = mean(y_hat_t1), .groups = "drop")

cf_om
cf_gap_AB_om <- cf_om$cf_mean[cf_om$category=="A"] -
                cf_om$cf_mean[cf_om$category=="B"]
cat("Counterfactual gap A - B (outcome model):", round(cf_gap_AB_om, 3), "\n")
cat("Gap closed (outcome model):",
    round(1 - cf_gap_AB_om / factual_gap_AB, 3), "\n\n")


## -------------------------------------------------------
## 3b. IPW (inverse probability weighting)
## -------------------------------------------------------
# Fit treatment propensity model
m_ps <- glm(treatment ~ confounder + category,
            data = dat, family = binomial)
dat$ps <- predict(m_ps, type = "response")   # P(T=1 | L, X)

# IPW estimator for E[Y(T=1)] within each category:
# E[Y(1) | X=x] = E[ T*Y / P(T=1|L,X)  |  X=x ]
#                 ---------------------------------
#                 E[ T / P(T=1|L,X)      |  X=x ]
# (Hajek / stabilized form)
cf_ipw <- dat |>
  group_by(category) |>
  summarise(
    num   = sum(treatment * outcome / ps),
    denom = sum(treatment / ps),
    cf_mean = num / denom,
    .groups = "drop"
  ) |>
  select(category, cf_mean)

cf_ipw
cf_gap_AB_ipw <- cf_ipw$cf_mean[cf_ipw$category=="A"] -
                 cf_ipw$cf_mean[cf_ipw$category=="B"]
cat("Counterfactual gap A - B (IPW):", round(cf_gap_AB_ipw, 3), "\n")
cat("Gap closed (IPW):",
    round(1 - cf_gap_AB_ipw / factual_gap_AB, 3), "\n\n")


## -------------------------------------------------------
## 3c. DOUBLY-ROBUST (AIPW) estimator
## -------------------------------------------------------
# Combines outcome model predictions with IPW correction.
# Consistent if EITHER the outcome model OR the propensity
# model is correctly specified.
#
# DR pseudo-outcome for person i targeting E[Y(1)]:
#   psi_i = y_hat_t1_i  +  (T_i / ps_i) * (Y_i - y_hat_t1_i)
#
# Then E[Y(1) | X=x] = mean(psi_i | X=x)

dat <- dat |>
  mutate(
    psi = y_hat_t1 + (treatment / ps) * (outcome - y_hat_t1)
  )

cf_dr <- dat |>
  group_by(category) |>
  summarise(cf_mean = mean(psi), .groups = "drop")

cf_dr
cf_gap_AB_dr <- cf_dr$cf_mean[cf_dr$category=="A"] -
                cf_dr$cf_mean[cf_dr$category=="B"]
cat("Counterfactual gap A - B (DR):", round(cf_gap_AB_dr, 3), "\n")
cat("Gap closed (DR):",
    round(1 - cf_gap_AB_dr / factual_gap_AB, 3), "\n\n")


# -------------------------------------------------------
# 4. SUMMARY TABLE
# -------------------------------------------------------
factual |>
  rename(factual_mean = mean_y) |>
  left_join(cf_om,  by = "category") |> rename(om  = cf_mean) |>
  left_join(cf_ipw, by = "category") |> rename(ipw = cf_mean) |>
  left_join(cf_dr,  by = "category") |> rename(dr  = cf_mean) |>
  print()


# -------------------------------------------------------
# 5. VISUALIZE: factual vs. counterfactual means
# -------------------------------------------------------
plot_df <- bind_rows(
  factual  |> mutate(estimator = "Factual",        mean = mean_y) |> select(category, estimator, mean),
  cf_om    |> mutate(estimator = "Outcome model")  |> rename(mean = cf_mean),
  cf_ipw   |> mutate(estimator = "IPW")            |> rename(mean = cf_mean),
  cf_dr    |> mutate(estimator = "Doubly-robust")  |> rename(mean = cf_mean)
) |>
  mutate(estimator = factor(estimator,
    levels = c("Factual", "Outcome model", "IPW", "Doubly-robust")))

ggplot(plot_df, aes(x = category, y = mean, color = estimator, group = estimator)) +
  geom_point(size = 3) +
  geom_line() +
  labs(
    title = "Factual vs. counterfactual (T=1) outcome means by category",
    subtitle = "Replicating Lundberg (2021) gap-closing estimand",
    x = "Category", y = "Mean outcome", color = NULL
  ) +
  theme_minimal()
