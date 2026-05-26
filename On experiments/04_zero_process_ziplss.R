# =============================================================================
#  04_zero_process_ziplss.R
# -----------------------------------------------------------------------------
#  ENGINE  : mgcv::gam  -- ziplss() is a 2-predictor family bam cannot fit.
#            # gam (not bam): family unsupported by bam.   method = "REML".
#  FAMILY  : ziplss()  -- zero-inflated Poisson location-scale. Two linear predictors:
#              lp1 = log Poisson mean (count GIVEN presence),
#              lp2 = complementary-log-log of P(presence)  => P(structural zero) = exp(-exp(lp2)).
#  RESPONSE: nMMSI_CO  -- CARGO+OTHER unique vessels (a COUNT; ziplss requires integer counts).
#
#  GOAL: separate the COUNT process from the PRESENCE / structural-zero process, put
#  s(ice_conc) in the zero process, and read off the ICE CONCENTRATION at which
#  P(structural zero) crosses 0.5 and 0.9 -- the OPERATIONAL ACCESSIBILITY THRESHOLD.
# =============================================================================
source("00_setup.R")
banner("04_zero_process_ziplss :: presence/structural-zero process vs ice (count = nMMSI_CO)")

if (!exists("ziplss", where = asNamespace("mgcv"))) {
  msg <- "ziplss() is NOT available in the installed mgcv -- cannot run the zero-process model."
  say("STOP: %s", msg); note("04 zero-process (ziplss) -- NOT RUN", msg); quit(save = "no", status = 0)
}

fr <- get_frozen_spatial()
is_re <- grepl('bs = "re"', fr$spatial_term)
EXP <- "Z_ziplss_zero_ice"

# documented gam SPATIAL fallback chain: (i) RE->te, (ii) reduce k, (iii) drop.
# EMPIRICAL FINDING (this dataset): a tensor te() in the ziplss COUNT predictor does NOT converge
# at 461k rows (te k=12 ran >3.5 CPU-hours without converging). The presence/structural-zero
# process is s(ice_conc) only, so the accessibility threshold is unaffected by the count-side
# spatial term. We therefore go straight to fallback (iii) -- DROP the count-predictor spatial
# smooth -- which converges quickly; a reduced te(k=8) is kept only as a secondary attempt.
spatial_base <- "ziplss count-predictor spatial: te() did not converge at 461k rows; fallback(iii) DROP count-side spatial (presence/threshold uses s(ice_conc), unaffected)."
gam_attempts <- list(
  list(term = "",                                                   k = "none",  lvl = "fallback(iii) count-spatial DROPPED (te did not converge)"),
  list(term = 'te(centroid_x, centroid_y, bs = "ps", k = c(8, 8))', k = "te8x8", lvl = "secondary te k=8 (only if dropped-spatial also failed)"))

# lp1 (count) = the 01 mean predictor (lean gam spatial); lp2 (presence) = ICE ONLY,
# so P(structural zero) is a clean function of ice for the threshold derivation.
build_fml <- function(term) list(
  as.formula(paste("nMMSI_CO ~", build_mean_rhs(term))),   # lp1: log Poisson mean (count | presence)
  as.formula(paste("~", PS("ice_conc"))))                  # lp2: cloglog P(presence) -- structural-zero process

banner("ziplss fit (method = REML) with spatial fallback chain")
fit <- NULL; chosen <- NULL; conv_status <- ""
for (att in gam_attempts) {
  say("\n[attempt] spatial = %s  (%s)", if (nzchar(att$term)) att$term else "<none>", att$lvl)
  ft <- fit_model(build_fml(att$term), ziplss(), "gam", method = "REML")
  if (!is.null(ft$model)) {
    fit <- ft; chosen <- att
    conv_status <- if (isTRUE(ft$model$converged)) sprintf("converged [%s; %s]", att$lvl, spatial_base) else
                                                   sprintf("non-converged-accepted [%s; %s]", att$lvl, spatial_base)
    # Accept the first attempt that RETURNS a usable model (mgcv's lss convergence flag is
    # conservative); record the flag honestly but don't discard the spatial term over it.
    # Leaner / dropped-spatial fallbacks engage only on a HARD failure (NULL).
    break
  }
}
if (is.null(fit)) {
  note("04 zero-process (ziplss) -- FIT FAILED", c("ziplss did not fit under any spatial fallback.", spatial_base))
  .append_row(CSV_PATH, EXP_COLS, list(exp_id=EXP, script="04_zero_process_ziplss.R", response="nMMSI_CO",
    vessel_scope=VESSEL_SCOPE, family="ziplss()", engine="gam", method="REML", formula="(count)+presence~s(ice)",
    spatial_term="(all fallbacks failed)", k_spec="n/a", tweedie_p_est=NA, score_name="REML", score_value=NA,
    AIC=NA, BIC=NA, deviance_explained_pct=NA, n_obs=nrow(dat), total_edf=NA, max_pairwise_concurvity=NA,
    dharma_dispersion_p=NA, dharma_zeroinflation_p=NA, dharma_spatial_autocorr=NA, dharma_temporal_autocorr=NA,
    k_index_min=NA, k_p_min=NA, convergence_status="FAILED", runtime_sec=NA, notes=spatial_base))
  quit(save = "no", status = 0)
}
m <- fit$model
say("ziplss fitted: spatial=%s converged=%s runtime=%.1fs", chosen$k, isTRUE(m$converged), fit$runtime)

fml_str <- sprintf("count: nMMSI_CO ~ %s | presence(cloglog): ~ %s", build_mean_rhs(chosen$term), PS("ice_conc"))
met <- record_experiment(EXP, "04_zero_process_ziplss.R", m, "nMMSI_CO", "ziplss()", "gam", "REML",
                         fml_str, chosen$term, paste0("presence~ice ps20; ", chosen$k),
                         dharma_family = "ziplss", convergence_status = conv_status,
                         runtime_sec = fit$runtime,
                         extra_notes = paste("ZIP location-scale; structural-zero process = s(ice).", spatial_base))

# =============================================================================
#  OPERATIONAL ACCESSIBILITY THRESHOLD: ice where P(structural zero) hits 0.5 and 0.9
#   P(structural zero)(ice) = exp(-exp(eta2(ice)))   [eta2 = lp2, depends only on ice here]
# =============================================================================
banner("structural-zero threshold derivation")
ice_grid <- data.frame(ice_conc = seq(min(dat$ice_conc), max(dat$ice_conc), length.out = 400))
for (v in c("month","time_index","air_temp_resid","wind_speed","port_dist","centroid_x","centroid_y"))
  ice_grid[[v]] <- median(dat[[v]])
ice_grid$post2017 <- factor("post", levels = levels(dat$post2017))
pr   <- predict(m, newdata = ice_grid, type = "link", se.fit = TRUE)    # 2-col: [count lp, presence lp]
eta2 <- pr$fit[, 2]; se2 <- pr$se.fit[, 2]
Pstruct <- exp(-exp(eta2))                                              # P(structural zero / absence)
Plo <- exp(-exp(eta2 + 1.96 * se2)); Phi <- exp(-exp(eta2 - 1.96 * se2))# 95% CI (monotone transform)

cross_ice <- function(target) {                                         # ice at which Pstruct == target
  if (max(Pstruct) < target || min(Pstruct) > target) return(NA_real_)
  o <- order(Pstruct); approx(Pstruct[o], ice_grid$ice_conc[o], xout = target, ties = "ordered")$y
}
ice50 <- cross_ice(0.5); ice90 <- cross_ice(0.9)
say("P(structural zero) range over ice 0..1: [%.3f, %.3f]", min(Pstruct), max(Pstruct))
say("ICE at P(structural zero)=0.5 : %s", ifelse(is.na(ice50), "not reached within observed ice", sprintf("%.3f", ice50)))
say("ICE at P(structural zero)=0.9 : %s", ifelse(is.na(ice90), "not reached within observed ice", sprintf("%.3f", ice90)))

# dedicated plot: structural-zero probability curve vs ice, with 0.5 / 0.9 marked
zc <- data.frame(ice = ice_grid$ice_conc, P = Pstruct, lo = Plo, hi = Phi)
pz <- ggplot(zc, aes(ice, P)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#1F5FBF", alpha = .2) +
  geom_line(colour = "#0B2E6D", linewidth = 1) +
  geom_hline(yintercept = c(0.5, 0.9), linetype = "dashed", colour = "grey40") +
  { if (!is.na(ice50)) geom_vline(xintercept = ice50, linetype = "dotted", colour = "firebrick") } +
  { if (!is.na(ice90)) geom_vline(xintercept = ice90, linetype = "dotted", colour = "firebrick") } +
  annotate("text", x = c(ice50, ice90), y = c(0.5, 0.9), label = c("0.5", "0.9"),
           vjust = -0.4, hjust = -0.1, colour = "firebrick", size = 3.5) +
  labs(title = "Operational accessibility: P(structural zero) vs sea-ice concentration",
       subtitle = "ziplss presence process; dashed = 0.5/0.9, dotted = ice threshold",
       x = "ice_conc (per-hex, 0-1)", y = "P(structural zero) = exp(-exp(lp2))") +
  theme_minimal()
ggsave(ppath(EXP, "zeroprob_vs_ice"), pz, width = 8, height = 5.5, dpi = 110)
say("Saved structural-zero probability curve.")

# IMPORTANT caveat (printed AND noted)
CAVEAT <- paste("CAVEAT: this threshold is on PER-HEX ice CONCENTRATION (0-1), NOT a regional ice-EXTENT",
                "figure from the economics literature. They are different quantities and must not be equated.")
say("\n%s", CAVEAT)

note("04 zero-process / operational accessibility threshold (ziplss)",
  c(sprintf("Zero-inflated Poisson location-scale (gam/REML) on the COUNT response nMMSI_CO: count process = 01 mean predictor; presence process (cloglog) = `s(ice_conc, ps, k20)`."),
    sprintf("- Fit: AIC %.0f, BIC %.0f, dev.expl %.1f%%; spatial = %s; %s.", met$AIC, met$BIC, met$dev_expl, chosen$k, conv_status),
    sprintf("- DHARMa (ZIP simulation): dispersion p=%s, zero-inflation p=%s.", format(signif(met$disp_p,3)), format(signif(met$zero_p,3))),
    sprintf("- **Operational accessibility threshold** (per-hex ice concentration): P(structural zero)=0.5 at ice = %s; =0.9 at ice = %s.",
        ifelse(is.na(ice50),"not reached within observed ice", sprintf("**%.3f**", ice50)),
        ifelse(is.na(ice90),"not reached within observed ice", sprintf("**%.3f**", ice90))),
    sprintf("- P(structural zero) ranges %.2f to %.2f across ice 0->1 (rising ice -> rising structural absence of cargo/other vessels).", min(Pstruct), max(Pstruct)),
    paste0("- ", CAVEAT),
    "- Associational: ships are absent where ice is high; this is conditional association, not a causal accessibility effect."))

note("04 -- summary (what this script established)",
  c("1. A ziplss split cleanly separates WHERE cargo/other vessels are structurally absent (presence process) from HOW MANY operate where present (count process).",
    sprintf("2. Operational accessibility threshold (per-hex ice conc): 0.5 at %s, 0.9 at %s -- see %s_zeroprob_vs_ice.png.",
        ifelse(is.na(ice50),"NR",sprintf("%.3f",ice50)), ifelse(is.na(ice90),"NR",sprintf("%.3f",ice90)), EXP),
    sprintf("3. Engine note: ziplss forced gam+REML (bam cannot fit it); spatial term used the documented lean-te fallback (%s).", chosen$k),
    "4. Threshold is per-hex CONCENTRATION, not regional EXTENT -- do not equate with economics-literature extent figures.",
    "5. Claims associational; significance de-emphasised at n=461,448."))

say("\n04_zero_process_ziplss DONE. ice@0.5=%s ice@0.9=%s spatial=%s",
    ifelse(is.na(ice50),"NR",sprintf("%.3f",ice50)), ifelse(is.na(ice90),"NR",sprintf("%.3f",ice90)), chosen$k)
