# Gap-closing explainer — session notes

Working file: `slides/03-decomp/gapclosing-explainer.qmd`
Companion R script: `slides/03-decomp/gapclosing-replication.R`

---

## What the document does

Replicates Lundberg (2024) `gapclosing` package simulation manually, walking through:

1. DGP (`set.seed(08544)`, n=1000) — exact clone of `generate_simulated_data()`
2. Factual gap by category (observed means)
3. Estimand clarification (gap-closing ≠ ATE within group; it's E[Y(1)|X=A] − E[Y(1)|X=B])
4. `gapclosing()` call (`se = FALSE`; use `cache: true` if `se = TRUE` to prevent bootstrap re-run on every render)
5. By-hand estimators:
   - **Outcome model** (standardization): predict Y(T=1) for everyone, average within category
   - **IPW** (Hajek): `sum(T*Y/p) / sum(T/p)` within category — treated obs only, 1/p reweighting makes them represent the full population
   - **Doubly-robust (AIPW)**: `psi_i = yhat_t1 + (T/p)*(Y - yhat_t1)`, then `mean(psi)` within category
6. Comparison table + plot

---

## Pending edits

- **"g-computation" label** — user agreed it should be updated. The Overview list says `**Outcome model** (g-computation / standardization)` and the section header says `### 3a. Outcome model (g-computation)`. Should change to just "standardization" or "outcome model (standardization)" — g-computation is the time-varying generalization (Robins 1986); the single time-point version is classical standardization, which predates Robins.

---

## Key conceptual points settled in session

### Gap-closing estimand vs ATE
The estimand is E[Y(1)|X=A] − E[Y(1)|X=B]: a cross-group comparison under a shared counterfactual intervention, NOT a within-group treatment effect (ATE or ATT). Step 2 of the document explains this.

### Why IPW uses treated obs only
Y(1) is only observable for T=1 (consistency assumption). The 1/p reweighting makes the treated subsample represent the full population within each category — so this IS estimating E[Y(1)|X=x] for everyone, not the ATT.

### Stabilized (Hajek) weights
Standard form: `w_i = 1 / ps_i` for treated obs.  
Stabilized form: `w_i = P(T=1) / ps_i` — marginal treatment probability in numerator. Keeps weights centered around 1, reduces variance from extreme propensities.

For weighted regression on treated obs:
```r
p_marginal <- mean(dat$treatment)
dat <- dat |> mutate(sw = if_else(treatment == 1, p_marginal / ps, NA_real_))

# Weighted regression
m_ipw_stab <- lm(outcome ~ category, data = dat[dat$treatment == 1, ], weights = sw)

# Or equivalently, weighted mean within category
cf_ipw_stab <- dat |>
  filter(treatment == 1) |>
  group_by(category) |>
  summarise(cf_mean = weighted.mean(outcome, sw), .groups = "drop")
```
The two are numerically equivalent (P(T=1) cancels in the ratio).

### Hajek vs Horvitz-Thompson in finite samples
`gapclosing` uses Hajek normalization (`sum(T/p*(Y-Ŷ)) / sum(T/p)`) for the DR correction; our code uses `sum(T/p*(Y-Ŷ)) / n`. They converge asymptotically (E[T/p] = 1) but differ in finite samples — explains small discrepancy between our DR and `gapclosing`'s `cf_mean_dr`.

### Terminology history
- **Standardization / outcome model**: classical epidemiology, predates Robins. Robins (1986) "g-computation" is the generalization to time-varying treatments.
- **Hajek estimator**: from survey sampling (Jaroslav Hájek, 1971), not epidemiology. Ratio estimator that divides by Σ(1/π_i) rather than N. Entered causal inference through Robins/Hernán IPW literature. In epidemiology you see "stabilized weights" for the same idea.

---

## marginaleffects notes

`avg_predictions(model, by = "category", newdata = dat_t1)` is the correct call for counterfactual means with SEs. Common mistakes:
- `var = "category"` sets all obs to each level (gives treatment effects, not counterfactual means)
- `datagrid(treatment = 1)` creates one synthetic row, not the full dataset
- Lambda `newdata = \(d) transform(...)` not supported in older marginaleffects versions — pass pre-made `dat_t1` directly
