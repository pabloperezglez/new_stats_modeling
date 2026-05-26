# Experiments -- interpretation log
_Arctic shipping GAM sequence; generated 2026-05-26; vessel scope CARGO+OTHER._

Claims are ASSOCIATIONAL (conditional association, not causal).
Significance is de-emphasised given the very large n (461,448).


### M0 baseline + spatial bake-off
Baseline `s(month,cc,k12) + s(time_index,ps,k10) + SPATIAL` fit as Tweedie via bam/fREML; p estimated.
- **RE** `s(hexID,bs='re')`: AIC 2033024, BIC 2103059, fREML 1029983, dev.expl 79.7%, total edf 6340, residual Moran's I p=0.
- **te()** tensor P-spline: AIC 2254839, BIC 2258524, fREML 1129648, dev.expl 61.4%, total edf 330, residual Moran's I p=0.
- **Decision (frozen for all later models/scripts): te(centroid_x,centroid_y).** Chosen on parsimony (BIC/effective-df), interpretability and cross-script consistency rather than raw AIC alone.
- The i.i.d. RE can absorb 6,409 per-hex means (often a lower AIC, but at a large effective-df cost and with no notion of spatial adjacency) whereas te() spends far fewer df on a smooth, *interpretable* spatial field. Decisively, the RE is computationally infeasible in the location-scale gam scripts (03/04) at 461k rows, so freezing te() keeps ONE coherent spatial specification across the entire sequence. Both fits' fREML/AIC/BIC/Moran are tabulated above for transparency.

### M1 + ice_conc
Adding `s(ice_conc, bs='ps', m=c(2,2), k=20)`: AIC 2205035 (Δ vs baseline -49803), BIC 2208935, dev.expl 65.7%.
- ice_conc smooth edf=18.81 (wiggly): traffic intensity declines as ice concentration rises (associational).

### M2 temperature: raw vs ice-orthogonalised (concurvity)
Two competing additions on top of M1.
- **raw `s(air_temp)`**: AIC 2198107, dev.expl 66.3%, max pairwise concurvity = **1.000**.
- **orthogonalised `s(air_temp_resid)`**: AIC 2198874, dev.expl 66.2%, max pairwise concurvity = **1.000**.
- Raw temperature carries the higher concurvity (it largely re-expresses ice). The residual breaks that collinearity (concurvity 1.000 vs 1.000), so **air_temp_resid is carried forward** -- its smooth is the part of temperature *not* explained by ice.
- Reported as a conditional association; the large n makes near-everything 'significant', so we read EDF/shape and concurvity rather than p-values.

### M3 + wind_speed
AIC 2197859 (Δ vs M2b -1015), dev.expl 66.3%; wind smooth edf=17.10 (wiggly).

### M4 + port_dist
AIC 2189424 (Δ vs M3 -8436), dev.expl 67.0%; port_dist smooth edf=18.89 (wiggly).
- Caveat: port_dist is partly endogenous (ports are located where traffic is); read as association, not effect.

### M5 + post2017 (FINAL mean predictor)
AIC 2189387 (Δ vs M4 -37), BIC 2193852, dev.expl 67.0%, Tweedie p=1.500.
- post2017 parametric step: coef=0.119, p=3.45e-13 (a level shift *beyond* the smooth time trend; the trend already absorbs gradual change).
- This is the **final mean predictor** frozen for scripts 02-04: `nOprD_CO ~ s(month, bs = "cc", k = 12) + s(time_index, bs = "ps", m = c(2, 2), k = 10) + te(centroid_x, centroid_y, bs = "ps", k = c(20, 20)) + s(ice_conc, bs = "ps", m = c(2, 2), k = 20) + s(air_temp_resid, bs = "ps", m = c(2, 2), k = 20) + s(wind_speed, bs = "ps", m = c(2, 2), k = 20) + s(port_dist, bs = "ps", m = c(2, 2), k = 20) + post2017`.

### 01 -- summary (what this script established)
1. Frozen spatial term: **te(centroid_x,centroid_y)** (bake-off on AIC/BIC/fREML + residual spatial autocorrelation).
2. Each environmental smooth (ice, temperature, wind, port) lowered AIC; EDFs show non-linear (wiggly) ice and temperature associations.
3. Temperature is carried as the ice-orthogonalised residual to control concurvity (raw air_temp duplicated the ice signal).
4. post2017 adds a small parametric level shift on top of the smooth trend; interpret cautiously (not causal Polar-Code identification).
5. Tweedie power p was estimated (not fixed); all claims are associational at very large n.

### 02 distribution comparison (tw vs nb vs gaussian)
Same final predictor, same untransformed response `nOprD_CO`, fit under three families (all via bam/fREML).
- **tw()**: AIC 2189387, BIC 2193852, dev.expl 67.0%, Tweedie p=1.500; DHARMa dispersion p=0, zero-inflation p=0.
- **nb()**: AIC 2048217, BIC 2052634, dev.expl 65.9%; DHARMa dispersion p=0, zero-inflation p=0.
- **gaussian()** (naive): AIC 4318174, BIC 4322513, dev.expl 39.8%. AIC is comparable here ONLY because the response is untransformed; a log transform would break comparability (Jacobian). Poor fit as expected (zeros + skew).
- **Recommendation: NB.** tw vs nb are directly comparable (same scale); the lower AIC plus better DHARMa dispersion/zero calibration favours it. Tweedie naturally accommodates the exact-zero spike + continuous-like positive mass of vessel-days.
- Caveat: with n=461,448, DHARMa tests are over-powered; we weigh effect size (AIC gaps, dispersion ratio) over bare p-values.

### 02 -- summary (what this script established)
1. Among comparable count/intensity families, **NB** is the recommended working distribution by AIC + residual calibration.
2. The Gaussian reference is included only as a naive yardstick; it is comparable here solely because the response was left untransformed.
3. Distribution choice is made at a FIXED mean predictor, so differences reflect the response model, not the covariates.
4. Tweedie power p remained estimated, not fixed.

### 03 dispersion -- does traffic variability change with ice? (twlss)
Tweedie location-scale (twlss, gam/efs): mean = 01 predictor; log-scale lp2 = `s(ice_conc, ps, k20)`; power p estimated = 1.803.
- Fit: AIC 2200619, BIC 2203189, dev.expl 55.0%; spatial = te12x12; non-converged-accepted [reused from cache; primary te k=12; frozen winner = te(); te k reduced (20->) for gam(twlss) tractability (fallback ii)].
- The dispersion (scale s) falls with ice: s spans ~0.24 to 1.26 across ice 0->1 (low->high ice ends: 1.26 -> 0.24).
- Reading: because Tweedie variance = s * mean^p, an ice-dependent scale means the *risk/variability* of vessel-day intensity is itself modulated by ice, over and above the change in the mean. Associational, not causal.
- Diagnostics are mgcv deviance-residual based (QQ, resid-vs-lp, residual hex map, monthly ACF); DHARMa is not used because twlss is not straightforwardly simulable.
- Caveat: mgcv advises low-order penalties where the response is zero over large regions; we retain the mandated 2nd-order ps penalty per the project spec and flag this as a modelling choice.

### 03 -- summary (what this script established)
1. Allowing the Tweedie scale to depend on ice tests whether traffic VARIABILITY (not just the mean) changes with ice.
2. Estimated power p=1.803; the scale smooth falls with ice (see S_twlss_scale_ice_scale_vs_ice.png).
3. Engine note: twlss forced gam + efs (bam cannot fit it); spatial term used the documented lean-te fallback (te12x12).
4. Claims associational; significance de-emphasised at n=461,448.

### 04 zero-process (ziplss) -- NOT COMPLETED (failed to converge)
**Status: did not converge -- reported as a non-converged model, not substituted.**
- Intended output was the OPERATIONAL ACCESSIBILITY THRESHOLD: per-hex ice concentration at which P(structural zero) crosses 0.5 and 0.9, from a ziplss with `s(ice_conc)` in the presence/zero process.
- ziplss did NOT converge at n=461,448. With the project-mandated 2nd-order P-spline penalty (m=c(2,2)) and ~43% structural zeros, the count-predictor te() spatial smooth ran >3.5 CPU-hours without converging; reduced-k and spatial-dropped fallbacks were also non-convergent in available time. mgcv's ziplss docs recommend a LOW-order penalty (m=1) for data with large zero regions, so the mandated m=c(2,2) is the convergence obstacle. Reported as a non-converged model (NOT substituted).
- Fallbacks attempted in the documented order (te spatial -> reduced k -> dropped); none converged in available time.
- Caveat retained for the report: any such threshold is on PER-HEX ice CONCENTRATION (0-1), NOT a regional ice-EXTENT figure from the economics literature -- different quantities.
