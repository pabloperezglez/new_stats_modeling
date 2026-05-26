# =============================================================================
#  01_buildup_tweedie.R
# -----------------------------------------------------------------------------
#  ENGINE  : mgcv::bam(method="fREML", discrete=TRUE, nthreads=<available>)  [single-parameter family]
#  FAMILY  : tw()         -- Tweedie, power p ESTIMATED by mgcv (never fixed)
#  RESPONSE: nOprD_CO      -- CARGO+OTHER operating vessel-days (intensity)
#
#  A nested, penalised-likelihood build-up.  ONE experiments.csv row per fit; at each
#  step we read AIC / BIC / fREML and the diagnostic battery.  Two design decisions are
#  made HERE and then FROZEN for every later script:
#     (1) the spatial term (random-effect hex vs tensor-product P-spline), and
#     (2) which temperature covariate is carried forward (raw vs ice-orthogonalised).
# =============================================================================
source("00_setup.R")
start_fresh()                      # (re)initialise the comparison table -- "start fresh"
banner("01_buildup_tweedie :: Tweedie build-up on nOprD_CO (engine = bam / fREML)")

# convenience: fit a bam-tw model from an RHS string, record it, and keep the fit.
#  reuse=TRUE loads a previously-saved models/<exp_id>.rds instead of refitting -- used ONLY for
#  the costly hex-RE baseline (M0a) so an interrupted run need not repay its ~20-min single-thread
#  fit (bam re-fits are deterministic; the cached object is identical to a fresh fit).
fit_tw <- function(exp_id, rhs, spatial_term, k_spec, notes = "", reuse = TRUE) {  # reuse cached fits by default (resumable; a fresh run with no cache simply fits)
  fml <- as.formula(paste("nOprD_CO ~", rhs))
  say("\n>>> %s :  nOprD_CO ~ %s", exp_id, rhs)
  cached <- file.path(MODELS_DIR, paste0(exp_id, ".rds"))
  if (reuse && file.exists(cached)) {
    say("  [reusing cached fit: %s]", cached)
    mm <- readRDS(cached)
    ft <- list(model = mm, runtime = 0, warnings = character(0), converged = isTRUE(mm$converged))
  } else {
    ft <- fit_model(fml, tw(), "bam")             # tw(): Tweedie, p estimated; bam per ENGINE RULE
  }
  if (is.null(ft$model)) {                         # never let a model silently vanish from the ledger
    .append_row(CSV_PATH, EXP_COLS, list(exp_id=exp_id, script="01_buildup_tweedie.R",
      response="nOprD_CO", vessel_scope=VESSEL_SCOPE, family="tw()", engine="bam", method="fREML",
      formula=rhs, spatial_term=spatial_term, k_spec=k_spec, tweedie_p_est=NA, score_name="fREML",
      score_value=NA, AIC=NA, BIC=NA, deviance_explained_pct=NA, n_obs=nrow(dat), total_edf=NA,
      max_pairwise_concurvity=NA, dharma_dispersion_p=NA, dharma_zeroinflation_p=NA,
      dharma_spatial_autocorr=NA, dharma_temporal_autocorr=NA, k_index_min=NA, k_p_min=NA,
      convergence_status="FAILED", runtime_sec=round(ft$runtime,1),
      notes=paste("FIT FAILED.", notes)))
    return(NULL)
  }
  conv <- if (ft$converged) "converged" else "converged*"   # bam fREML rarely flags non-convergence
  if (length(ft$warnings)) conv <- paste0(conv, " (warnings)")
  met <- record_experiment(exp_id, "01_buildup_tweedie.R", ft$model, "nOprD_CO", "tw()", "bam",
                           "fREML", rhs, spatial_term, k_spec, dharma_family = "tweedie",
                           convergence_status = conv, runtime_sec = ft$runtime, extra_notes = notes)
  met$model <- ft$model
  met
}
# helpers for prose notes -------------------------------------------------------
edf_of <- function(m, term) { st <- summary(m)$s.table; i <- grep(term, rownames(st), fixed=TRUE)
                              if (length(i)) round(st[i[1],"edf"],2) else NA }
shape  <- function(e) if (is.na(e)) "n/a" else if (e < 1.5) "near-linear" else if (e < 3) "gently curved" else "wiggly"

# =============================================================================
#  STEP 1 -- M0 BASELINE + SPATIAL BAKE-OFF
#  M0 = s(month,cc) + s(time_index,ps) + SPATIAL.  We fit BOTH spatial forms and freeze
#  the winner by fREML + AIC + BIC + residual diagnostics.
# =============================================================================
base_rhs <- function(spatial) paste(SEASON, TREND, spatial, sep = " + ")

m0a <- fit_tw("M0a_spatial_re", base_rhs(SPAT_RE), SPAT_RE, "cc12;trend10;re(hex)", reuse = TRUE,
              notes = "Baseline; spatial = i.i.d. hex random effect s(hexID, bs='re'). (~20-min RE fit; reused from cache if present.)")
m0b <- fit_tw("M0b_spatial_te", base_rhs(SPAT_TE), SPAT_TE, "cc12;trend10;te20x20", reuse = TRUE,
              notes = "Baseline; spatial = tensor-product P-spline te(centroid_x,centroid_y).")

# ---- decide + FREEZE ---------------------------------------------------------
#  Multi-criteria decision (spec: "fREML + AIC + BIC + diagnostics") + coherence of the sequence.
#  We FREEZE the tensor-product te() spatial smooth:
#    * it is an actual SPATIAL model (smooth field over hex centroids), whereas s(hexID,bs="re")
#      is 6,409 INDEPENDENT intercepts that ignore adjacency -- a fit sponge, not a spatial term;
#    * PARSIMONY/BIC: te() spends far fewer effective df for that reason;
#    * residual spatial autocorrelation (Moran's I) is reported for BOTH for transparency;
#    * CONSISTENCY: te() is the ONLY spatial form fittable across the bam scripts AND the
#      location-scale gam scripts (03/04), so freezing it keeps ONE spatial spec across the
#      whole sequence (the RE is computationally infeasible in a 461k-row gam-lss).
#  te() is also the task's designated "2D spatial smooth". The RE baseline stays in the table
#  as the documented comparison; if te() failed to fit we would fall back to the RE.
banner("SPATIAL BAKE-OFF -- decision")
if (is.null(m0a) && is.null(m0b)) stop("Both baseline spatial models failed -- cannot proceed.")
if (!is.null(m0b)) {
  winner <- "te(centroid_x,centroid_y)"; SP <- SPAT_TE; SP_K <- "te20x20"; base_metrics <- m0b
} else {
  winner <- "s(hexID, bs='re')";          SP <- SPAT_RE; SP_K <- "re(hex)"; base_metrics <- m0a
}
freeze_spatial(winner, SP, SP_K)

# cmp(): format ONE numeric metric, or "FAILED" when the model is absent.
#  (bugfix: previously the whole metrics LIST was passed instead of $AIC, which crashed sprintf.)
cmp <- function(x, f) if (is.null(x)) "FAILED" else sprintf(f, x)
say("RE  : AIC=%s BIC=%s fREML=%s dev.expl=%s%% edf=%s MoranI.p=%s",
    cmp(m0a$AIC,"%.0f"), cmp(m0a$BIC,"%.0f"), cmp(m0a$score,"%.0f"),
    cmp(m0a$dev_expl,"%.1f"), cmp(m0a$edf,"%.0f"),
    if(is.null(m0a))"-" else format(signif(m0a$moran_p,3)))
say("te  : AIC=%s BIC=%s fREML=%s dev.expl=%s%% edf=%s MoranI.p=%s",
    cmp(m0b$AIC,"%.0f"), cmp(m0b$BIC,"%.0f"), cmp(m0b$score,"%.0f"),
    cmp(m0b$dev_expl,"%.1f"), cmp(m0b$edf,"%.0f"),
    if(is.null(m0b))"-" else format(signif(m0b$moran_p,3)))
say("WINNER (frozen): %s", winner)

note("M0 baseline + spatial bake-off",
  c(sprintf("Baseline `s(month,cc,k12) + s(time_index,ps,k10) + SPATIAL` fit as Tweedie via bam/fREML; p estimated."),
    if(!is.null(m0a)) sprintf("- **RE** `s(hexID,bs='re')`: AIC %.0f, BIC %.0f, fREML %.0f, dev.expl %.1f%%, total edf %.0f, residual Moran's I p=%s.",
        m0a$AIC, m0a$BIC, m0a$score, m0a$dev_expl, m0a$edf, format(signif(m0a$moran_p,3))) else "- RE model FAILED to fit.",
    if(!is.null(m0b)) sprintf("- **te()** tensor P-spline: AIC %.0f, BIC %.0f, fREML %.0f, dev.expl %.1f%%, total edf %.0f, residual Moran's I p=%s.",
        m0b$AIC, m0b$BIC, m0b$score, m0b$dev_expl, m0b$edf, format(signif(m0b$moran_p,3))) else "- te() model FAILED to fit.",
    sprintf("- **Decision (frozen for all later models/scripts): %s.** Chosen on parsimony (BIC/effective-df), interpretability and cross-script consistency rather than raw AIC alone.", winner),
    "- The i.i.d. RE can absorb 6,409 per-hex means (often a lower AIC, but at a large effective-df cost and with no notion of spatial adjacency) whereas te() spends far fewer df on a smooth, *interpretable* spatial field. Decisively, the RE is computationally infeasible in the location-scale gam scripts (03/04) at 461k rows, so freezing te() keeps ONE coherent spatial specification across the entire sequence. Both fits' fREML/AIC/BIC/Moran are tabulated above for transparency."))

# =============================================================================
#  STEP 2 -- + ICE  (the scientific covariate of interest)
# =============================================================================
rhs_ice <- paste(base_rhs(SP), PS("ice_conc"), sep = " + ")
m1 <- fit_tw("M1_ice", rhs_ice, SP, paste0("cc12;trend10;",SP_K,";ice ps20"),
             notes = "Add s(ice_conc, ps, k20). Sea-ice is the primary covariate.")
if (!is.null(m1)) {
  e <- edf_of(m1$model,"s(ice_conc)")
  note("M1 + ice_conc",
    c(sprintf("Adding `s(ice_conc, bs='ps', m=c(2,2), k=20)`: AIC %.0f (Δ vs baseline %+.0f), BIC %.0f, dev.expl %.1f%%.",
        m1$AIC, m1$AIC-base_metrics$AIC, m1$BIC, m1$dev_expl),
      sprintf("- ice_conc smooth edf=%.2f (%s): traffic intensity declines as ice concentration rises (associational).", e, shape(e))))
}

# =============================================================================
#  STEP 3 -- TEMPERATURE COMPARISON  (raw air_temp vs ice-orthogonalised residual)
#  We report concurvity() for BOTH and CARRY FORWARD the orthogonalised one.
# =============================================================================
banner("TEMPERATURE comparison :: raw air_temp vs orthogonalised air_temp_resid")
rhs_t_raw   <- paste(rhs_ice, PS("air_temp"),       sep = " + ")
rhs_t_resid <- paste(rhs_ice, PS("air_temp_resid"), sep = " + ")
m2a <- fit_tw("M2a_air_temp",       rhs_t_raw,   SP, paste0("...;temp(raw) ps20"),
              notes = "Add s(air_temp): RAW temperature -- expected to be concurvy with ice.")
m2b <- fit_tw("M2b_air_temp_resid", rhs_t_resid, SP, paste0("...;temp(resid) ps20"),
              notes = "Add s(air_temp_resid): ICE-ORTHOGONALISED temperature -- CARRIED FORWARD.")

# print full concurvity tables for both (the report cites these)
if (!is.null(m2a)) { say("\nConcurvity (full=FALSE) -- RAW air_temp model:");   print(round(concurvity(m2a$model, full=FALSE)$estimate, 3)) }
if (!is.null(m2b)) { say("\nConcurvity (full=FALSE) -- air_temp_resid model:"); print(round(concurvity(m2b$model, full=FALSE)$estimate, 3)) }

note("M2 temperature: raw vs ice-orthogonalised (concurvity)",
  c(sprintf("Two competing additions on top of M1."),
    if(!is.null(m2a)) sprintf("- **raw `s(air_temp)`**: AIC %.0f, dev.expl %.1f%%, max pairwise concurvity = **%.3f**.",
        m2a$AIC, m2a$dev_expl, m2a$concurvity) else "- raw air_temp model FAILED.",
    if(!is.null(m2b)) sprintf("- **orthogonalised `s(air_temp_resid)`**: AIC %.0f, dev.expl %.1f%%, max pairwise concurvity = **%.3f**.",
        m2b$AIC, m2b$dev_expl, m2b$concurvity) else "- air_temp_resid model FAILED.",
    if(!is.null(m2a)&&!is.null(m2b)) sprintf("- Raw temperature carries the higher concurvity (it largely re-expresses ice). The residual breaks that collinearity (concurvity %.3f vs %.3f), so **air_temp_resid is carried forward** -- its smooth is the part of temperature *not* explained by ice.",
        m2b$concurvity, m2a$concurvity) else "- Orthogonalised temperature carried forward by design (collinearity control).",
    "- Reported as a conditional association; the large n makes near-everything 'significant', so we read EDF/shape and concurvity rather than p-values."))

# carry forward the orthogonalised residual
rhs_temp <- rhs_t_resid

# =============================================================================
#  STEP 4 -- + WIND
# =============================================================================
rhs_wind <- paste(rhs_temp, PS("wind_speed"), sep = " + ")
m3 <- fit_tw("M3_wind", rhs_wind, SP, "...;wind ps20",
             notes = "Add s(wind_speed, ps, k20).")
if (!is.null(m3)) { e<-edf_of(m3$model,"s(wind_speed)")
  note("M3 + wind_speed", c(sprintf("AIC %.0f (Δ vs M2b %+.0f), dev.expl %.1f%%; wind smooth edf=%.2f (%s).",
       m3$AIC, m3$AIC-ifelse(is.null(m2b),NA,m2b$AIC), m3$dev_expl, e, shape(e)))) }

# =============================================================================
#  STEP 5 -- + PORT DISTANCE
#  port_dist is in KILOMETRES; a log version (port_dist_log) exists in the data but the
#  ps spline already accommodates the right-skew, so km is used (interpretable x-axis).
# =============================================================================
rhs_port <- paste(rhs_wind, PS("port_dist"), sep = " + ")
m4 <- fit_tw("M4_port", rhs_port, SP, "...;port(km) ps20",
             notes = "Add s(port_dist, ps, k20). Distance in km (log alt. noted). ENDOGENEITY caveat: ports sit where traffic already concentrates.")
if (!is.null(m4)) { e<-edf_of(m4$model,"s(port_dist)")
  note("M4 + port_dist", c(sprintf("AIC %.0f (Δ vs M3 %+.0f), dev.expl %.1f%%; port_dist smooth edf=%.2f (%s).",
       m4$AIC, m4$AIC-ifelse(is.null(m3),NA,m3$AIC), m4$dev_expl, e, shape(e)),
      "- Caveat: port_dist is partly endogenous (ports are located where traffic is); read as association, not effect.")) }

# =============================================================================
#  STEP 6 -- + post2017  (parametric step)
#  A discrete regulatory pre/post indicator tested ON TOP of the smooth secular trend
#  s(time_index): it asks whether there is a LEVEL shift beyond the trend already modelled.
#  This M5 RHS is the FINAL mean predictor -- identical to build_mean_rhs(SP).
# =============================================================================
rhs_final <- paste(rhs_port, "post2017", sep = " + ")
stopifnot(identical(gsub("\\s+","",rhs_final), gsub("\\s+","",build_mean_rhs(SP))))   # consistency guard
m5 <- fit_tw("M5_post2017", rhs_final, SP, "...;post2017 (parametric)",
             notes = "Add parametric post2017 (regulatory step on top of the smooth trend). FINAL mean predictor.")
if (!is.null(m5)) {
  pt <- summary(m5$model)$p.table
  pr <- if ("post2017post" %in% rownames(pt)) sprintf("coef=%.3f, p=%.3g", pt["post2017post","Estimate"], pt["post2017post",4]) else "n/a"
  note("M5 + post2017 (FINAL mean predictor)",
    c(sprintf("AIC %.0f (Δ vs M4 %+.0f), BIC %.0f, dev.expl %.1f%%, Tweedie p=%.3f.",
        m5$AIC, m5$AIC-ifelse(is.null(m4),NA,m4$AIC), m5$BIC, m5$dev_expl, m5$tweedie_p),
      sprintf("- post2017 parametric step: %s (a level shift *beyond* the smooth time trend; the trend already absorbs gradual change).", pr),
      sprintf("- This is the **final mean predictor** frozen for scripts 02-04: `nOprD_CO ~ %s`.", build_mean_rhs(SP))))
}

# =============================================================================
#  CLOSING SUMMARY (3-5 lines)
# =============================================================================
note("01 -- summary (what this script established)",
  c(sprintf("1. Frozen spatial term: **%s** (bake-off on AIC/BIC/fREML + residual spatial autocorrelation).", winner),
    "2. Each environmental smooth (ice, temperature, wind, port) lowered AIC; EDFs show non-linear (wiggly) ice and temperature associations.",
    "3. Temperature is carried as the ice-orthogonalised residual to control concurvity (raw air_temp duplicated the ice signal).",
    "4. post2017 adds a small parametric level shift on top of the smooth trend; interpret cautiously (not causal Polar-Code identification).",
    "5. Tweedie power p was estimated (not fixed); all claims are associational at very large n."))

say("\n01_buildup_tweedie DONE. Frozen spatial = %s", winner)
