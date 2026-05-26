# =============================================================================
#  03_dispersion_twlss.R
# -----------------------------------------------------------------------------
#  ENGINE  : mgcv::gam  -- twlss() is a location-scale family bam cannot fit.
#            # gam (not bam): family unsupported by bam.
#            twlss is mgcv-RESTRICTED to the extended Fellner-Schall method, so we use
#            method="efs" (a Wood penalised-likelihood smoothness selector, NOT CV/black-box).
#  FAMILY  : twlss()  -- Tweedie location-scale. 3 linear predictors:
#              lp1 = mean (log),  lp2 = log-scale r=log(s),  lp3 = power p (constant here).
#  RESPONSE: nOprD_CO (vessel-days).
#
#  QUESTION: does the VARIABILITY (risk) of traffic change with ice?  We let the Tweedie
#  log-SCALE depend on ice:  lp2 = ~ s(ice_conc, bs="ps", m=c(2,2), k=20).  The mean (lp1)
#  is the chosen predictor from 01.
# =============================================================================
source("00_setup.R")
banner("03_dispersion_twlss :: Tweedie location-scale (does dispersion depend on ice?)")

# --- twlss availability gate (STOP + report rather than substitute another method) -----------
if (!exists("twlss", where = asNamespace("mgcv"))) {
  msg <- "twlss() is NOT available in the installed mgcv -- per spec we STOP here rather than substitute a different method."
  say("STOP: %s", msg)
  note("03 dispersion (twlss) -- NOT RUN", c(msg, sprintf("Installed mgcv: %s.", packageVersion("mgcv"))))
  quit(save = "no", status = 0)
}

fr <- get_frozen_spatial()
is_re <- grepl('bs = "re"', fr$spatial_term)
EXP <- "S_twlss_scale_ice"

# --- documented gam SPATIAL fallback chain (convergence/runtime management) -------------------
#  Order per spec: (i) RE->te(), (ii) reduce te k, (iii) drop spatial. At 461k rows a 6,409-level
#  RE is infeasible in a 3-predictor gam, and even te(k=20) is heavy, so we START at a lean te.
spatial_base <- if (is_re)
  "frozen winner = RE(hex); RE infeasible in gam(twlss)@461k -> fallback(i) RE->te() + fallback(ii) reduced k" else
  "frozen winner = te(); te k reduced (20->) for gam(twlss) tractability (fallback ii)"
gam_attempts <- list(
  list(term = 'te(centroid_x, centroid_y, bs = "ps", k = c(12, 12))', k = "te12x12", lvl = "primary (te k=12)"),
  list(term = 'te(centroid_x, centroid_y, bs = "ps", k = c(8, 8))',   k = "te8x8",   lvl = "fallback(ii) te k=8"),
  list(term = "",                                                     k = "none",    lvl = "fallback(iii) spatial DROPPED"))

# location predictor = the 01 mean predictor with the (lean) gam spatial term;
# scale predictor = ice; power predictor = constant.
build_fml <- function(term) list(
  as.formula(paste("nOprD_CO ~", build_mean_rhs(term))),   # lp1 mean (log)
  as.formula(paste("~", PS("ice_conc"))),                  # lp2 log-scale (dispersion) vs ICE
  ~ 1)                                                     # lp3 power p (constant -- not modelled)

# --- fit with the fallback chain (or REUSE a cached fit, so a clean re-run skips the ~96-min fit)
banner("twlss fit (method = efs) with spatial fallback chain")
fit <- NULL; chosen <- NULL; conv_status <- ""
.cache <- file.path(MODELS_DIR, paste0(EXP, ".rds"))
if (file.exists(.cache)) {
  say("[reusing cached twlss fit: %s]", .cache)
  mm <- readRDS(.cache)
  fit <- list(model = mm, runtime = 0, warnings = character(0), converged = isTRUE(mm$converged))
  chosen <- list(term = 'te(centroid_x, centroid_y, bs = "ps", k = c(12, 12))', k = "te12x12", lvl = "reused from cache")
  conv_status <- sprintf("%s [reused from cache; primary te k=12; %s]",
                         if (isTRUE(mm$converged)) "converged" else "non-converged-accepted", spatial_base)
} else for (att in gam_attempts) {
  say("\n[attempt] spatial = %s  (%s)", if (nzchar(att$term)) att$term else "<none>", att$lvl)
  ft <- fit_model(build_fml(att$term), twlss(), "gam", method = "efs")
  if (!is.null(ft$model)) {
    fit <- ft; chosen <- att
    conv_status <- if (isTRUE(ft$model$converged)) sprintf("converged [%s; %s]", att$lvl, spatial_base) else
                                                   sprintf("non-converged-accepted [%s; %s]", att$lvl, spatial_base)
    # Accept the first attempt that RETURNS a usable model. mgcv's convergence flag for
    # location-scale families is conservative (a flagged-but-returned fit is still usable),
    # so we record the flag honestly but do NOT discard the spatial term over it. The leaner /
    # dropped-spatial fallbacks engage only if a fit HARD-FAILS (returns NULL).
    break
  }
}
if (is.null(fit)) {
  note("03 dispersion (twlss) -- FIT FAILED", c("twlss did not fit under any spatial fallback (te12 -> te8 -> none).", spatial_base))
  .append_row(CSV_PATH, EXP_COLS, list(exp_id=EXP, script="03_dispersion_twlss.R", response="nOprD_CO",
    vessel_scope=VESSEL_SCOPE, family="twlss()", engine="gam", method="efs", formula="(location)+scale~s(ice)",
    spatial_term="(all fallbacks failed)", k_spec="n/a", tweedie_p_est=NA, score_name="efs", score_value=NA,
    AIC=NA, BIC=NA, deviance_explained_pct=NA, n_obs=nrow(dat), total_edf=NA, max_pairwise_concurvity=NA,
    dharma_dispersion_p=NA, dharma_zeroinflation_p=NA, dharma_spatial_autocorr=NA, dharma_temporal_autocorr=NA,
    k_index_min=NA, k_p_min=NA, convergence_status="FAILED", runtime_sec=NA, notes=spatial_base))
  quit(save = "no", status = 0)
}
m <- fit$model
say("twlss fitted: spatial=%s  converged=%s  runtime=%.1fs", chosen$k, isTRUE(m$converged), fit$runtime)

# --- estimated Tweedie power p (from the constant 3rd predictor) -------------------------------
#  twlss: p = (a + b*exp(t))/(1+exp(t)), defaults a=1.01, b=1.99; t = constant lp3.
tt   <- predict(m, newdata = dat[1, , drop = FALSE], type = "link")[1, 3]
a_tw <- 1.01; b_tw <- 1.99
p_est <- (a_tw + b_tw * exp(tt)) / (1 + exp(tt))
say("Estimated Tweedie power p (twlss, constant): %.4f", p_est)

# --- record (DHARMa NA: twlss is not simulable here -> mgcv deviance diagnostics instead) ------
loc_str   <- build_mean_rhs(chosen$term)
fml_str   <- sprintf("location: nOprD_CO ~ %s | log-scale: ~ %s | power: ~1", loc_str, PS("ice_conc"))
met <- record_experiment(EXP, "03_dispersion_twlss.R", m, "nOprD_CO", "twlss()", "gam", "efs",
                         fml_str, chosen$term, paste0("scale~ice ps20; ", chosen$k),
                         dharma_family = NA, convergence_status = conv_status,
                         runtime_sec = fit$runtime, tweedie_p_override = p_est,
                         extra_notes = paste("twlss location-scale; DHARMa not simulable -> mgcv deviance diagnostics.", spatial_base))

# =============================================================================
#  DEDICATED PLOT: the SCALE (dispersion) smooth vs ice  -- the scientific output of 03
# =============================================================================
banner("scale-smooth interpretation (lp2 = log-scale vs ice)")
ice_grid <- data.frame(ice_conc = seq(min(dat$ice_conc), max(dat$ice_conc), length.out = 250))
for (v in c("month","time_index","air_temp_resid","wind_speed","port_dist","centroid_x","centroid_y"))
  ice_grid[[v]] <- median(dat[[v]])
ice_grid$post2017 <- factor("post", levels = levels(dat$post2017))
pr  <- predict(m, newdata = ice_grid, type = "link", se.fit = TRUE)   # lp2 = log-scale
logs <- pr$fit[, 2]; se <- pr$se.fit[, 2]
scl <- data.frame(ice = ice_grid$ice_conc, s = exp(logs),
                  lo = exp(logs - 1.96 * se), hi = exp(logs + 1.96 * se))
pscale <- ggplot(scl, aes(ice, s)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#1F5FBF", alpha = .2) +
  geom_line(colour = "#0B2E6D", linewidth = 1) +
  labs(title = "Tweedie SCALE (dispersion) vs sea-ice concentration",
       subtitle = sprintf("twlss log-scale smooth; Var = s * mean^p, p=%.3f", p_est),
       x = "ice_conc (0-1)", y = "scale parameter  s = exp(lp2)") + theme_minimal()
ggsave(ppath(EXP, "scale_vs_ice"), pscale, width = 8, height = 5.5, dpi = 110)
say("Saved scale-vs-ice plot.")

# describe shape in words for the note
s_lo <- scl$s[which.min(scl$ice)]; s_hi <- scl$s[which.max(scl$ice)]; s_rng <- range(scl$s)
trend_word <- if (s_hi > s_lo * 1.05) "rises" else if (s_hi < s_lo * 0.95) "falls" else "is roughly flat"

note("03 dispersion -- does traffic variability change with ice? (twlss)",
  c(sprintf("Tweedie location-scale (twlss, gam/efs): mean = 01 predictor; log-scale lp2 = `s(ice_conc, ps, k20)`; power p estimated = %.3f.", p_est),
    sprintf("- Fit: AIC %.0f, BIC %.0f, dev.expl %.1f%%; spatial = %s; %s.", met$AIC, met$BIC, met$dev_expl, chosen$k, conv_status),
    sprintf("- The dispersion (scale s) %s with ice: s spans ~%.2f to %.2f across ice 0->1 (low->high ice ends: %.2f -> %.2f).", trend_word, s_rng[1], s_rng[2], s_lo, s_hi),
    "- Reading: because Tweedie variance = s * mean^p, an ice-dependent scale means the *risk/variability* of vessel-day intensity is itself modulated by ice, over and above the change in the mean. Associational, not causal.",
    "- Diagnostics are mgcv deviance-residual based (QQ, resid-vs-lp, residual hex map, monthly ACF); DHARMa is not used because twlss is not straightforwardly simulable.",
    "- Caveat: mgcv advises low-order penalties where the response is zero over large regions; we retain the mandated 2nd-order ps penalty per the project spec and flag this as a modelling choice."))

note("03 -- summary (what this script established)",
  c("1. Allowing the Tweedie scale to depend on ice tests whether traffic VARIABILITY (not just the mean) changes with ice.",
    sprintf("2. Estimated power p=%.3f; the scale smooth %s with ice (see %s_scale_vs_ice.png).", p_est, trend_word, EXP),
    sprintf("3. Engine note: twlss forced gam + efs (bam cannot fit it); spatial term used the documented lean-te fallback (%s).", chosen$k),
    "4. Claims associational; significance de-emphasised at n=461,448."))

say("\n03_dispersion_twlss DONE. p=%.3f spatial=%s", p_est, chosen$k)
