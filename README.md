# Statistical Modeling Project — Pacific-Arctic Shipping GAMs

Generalized Additive Model (GAM) analysis of monthly Pacific-Arctic vessel traffic
(2015–2020, 6,409 hex cells × 72 months = 461,448 hex-months) as a function of sea-ice
concentration, air temperature, wind speed, distance to port, season, trend and space.
Modelling follows Wood's penalised-likelihood / information-criterion / residual-diagnostic
approach (`mgcv`); no cross-validation or black-box methods.

## Repository layout

```
On Data/
  Code/                  data-construction + EDA scripts (build_dataset_and_eda.R, …)
  Roadmaps/              project planning docs
  Data/                  *** NOT in repo *** read-only third-party inputs (see below)
On experiments/
  00_setup.R             data load, derived responses, shared helpers, smooth/engine conventions
  01_buildup_tweedie.R   nested Tweedie mean build-up (M0–M5) + spatial bake-off
  02_distribution_comparison.R   tw vs nb vs gaussian at the final predictor
  03_dispersion_twlss.R  Tweedie location-scale: does dispersion depend on ice?
  04_zero_process_ziplss.R       zero-inflated count model (did not converge — see below)
  experiments.csv        one row per fitted model (the comparison table)
  experiments_terms.csv  per-term EDF + k.check, keyed by exp_id
  experiments_interpretation.md  per-experiment write-up (so the report writes itself)
  plots_manifest.md      numbered index of the 49 diagnostic figures
  run.log                full console log of the modelling run
  plots/                 49 PNG diagnostics (smooth effects, QQ, residual maps, ACF, …)
  models/                *** NOT in repo *** 2.5 GB of fitted .rds models (regenerable)
```

## Not included (and why)

- **`On experiments/models/`** — fitted GAM objects total ~2.5 GB (the hex random-effect
  baseline alone is ~2 GB, far over GitHub's 100 MB per-file limit). Regenerate by running
  the `On experiments` scripts in order.
- **`On Data/Data/`** — the analytic dataset and raw inputs (NSIDC sea-ice CDR v5, ERA5
  reanalysis, Kapsar et al. 2022 shipping). These are read-only, license-restricted, and
  large, so they are not redistributed here; obtain them from their original sources and
  place them under `On Data/Data/`.

## Reproduce

```r
# from inside "On experiments/", with the analytic dataset present in ../On Data/Data/
Rscript 01_buildup_tweedie.R
Rscript 02_distribution_comparison.R
Rscript 03_dispersion_twlss.R   # twlss (gam + efs)
Rscript 04_zero_process_ziplss.R
```
Requires R with `mgcv`, `DHARMa`, `gratia`, `ggplot2`, `ape`.

## Conventions (cited in the report)

- **Engine:** `bam(method="fREML", discrete=TRUE)` for single-parameter families
  (gaussian, `tw()`, `nb()`); `gam(method="REML"|"efs")` only where `bam` can't fit the
  family (`twlss`, `ziplss`).
- **Smooths:** every continuous covariate is `s(x, bs="ps", m=c(2,2), k=20)` (cubic P-spline,
  2nd-order penalty); season is the one exception — `s(month, bs="cc", k=12)` (cyclic);
  space is `te(centroid_x, centroid_y, bs="ps")`.
- **Responses (CARGO + OTHER only):** `nOprD_CO` (vessel-days, Tweedie) and `nMMSI_CO`
  (unique vessels, count). Tweedie power `p` is estimated, never fixed.

## Results at a glance

- **Mean structure:** every environmental smooth lowered AIC; the full predictor (M5) is best
  in the build-up; ice and temperature effects are clearly non-linear. The tensor `te()` spatial
  smooth was frozen over the hex random effect (parsimony, interpretability, cross-script consistency).
- **Distribution:** at the final predictor, **negative binomial** has the lowest AIC/BIC
  (Tweedie second; Gaussian far worse). DHARMa tests are over-powered at n≈461k, so AIC gaps and
  effect sizes are weighed over p-values.
- **Dispersion (twlss):** the Tweedie scale declines with ice — traffic *variability* is also
  ice-modulated, not just the mean.
- **Zero process (ziplss): did not converge** at 461k rows under the mandated 2nd-order penalty
  with ~43% structural zeros (mgcv recommends a low-order penalty for such data); reported as a
  non-converged model rather than substituted.

_Claims are associational (conditional association, not causal); significance is de-emphasised given the large n._
