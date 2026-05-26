# =============================================================================
#  00_setup.R  --  shared setup + helpers for the Arctic-shipping GAM experiments
# -----------------------------------------------------------------------------
#  ROLE OF THIS FILE
#    Every modelling script (01..04) begins with  source("00_setup.R").  This file
#    therefore is the SINGLE place that:
#      * loads the READ-ONLY analytic dataset and derives the modelling responses,
#      * fixes the GLOBAL, NON-NEGOTIABLE smooth/engine conventions (one definition,
#        reused verbatim everywhere -- so the report can cite ONE specification),
#      * provides the logging helpers that write the comparison table
#        (experiments.csv), the per-term table (experiments_terms.csv), the saved
#        diagnostics (/plots) and the prose interpretation (experiments_interpretation.md).
#
#  WHAT THIS FILE NEVER DOES
#    * It NEVER writes, modifies or deletes any input DATA file.  readRDS() only.
#    * It NEVER installs packages (all required packages are assumed present).
# =============================================================================

# -----------------------------------------------------------------------------
# 0. SELF-LOCATE  --  resolve every path relative to THIS folder ("On experiments")
#    so the project is portable.  Data live in the sibling folder "On Data/Data".
# -----------------------------------------------------------------------------
.locate_dir <- function() {
  cand <- character(0)
  # 1. when source()'d (the real pipeline path): this file's own ofile
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { cand <- c(cand, dirname(normalizePath(of, mustWork = FALSE))); break }
  }
  # 2. Rscript --file, but ONLY trust it when it actually points at 00_setup.R
  a <- commandArgs(FALSE); fa <- grep("^--file=", a, value = TRUE)
  if (length(fa)) { f <- normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE)
                    if (basename(f) == "00_setup.R") cand <- c(cand, dirname(f)) }
  # 3. current working directory (scripts are run from inside "On experiments")
  cand <- c(cand, normalizePath(getwd(), mustWork = FALSE))
  # 4. last-resort hard fallback (known project location)
  cand <- c(cand, "/Users/pabloperezgonzalez/Downloads/Statistical Modeling Project/On experiments")
  for (d in cand) if (file.exists(file.path(d, "00_setup.R"))) return(d)
  cand[1]
}
OUT_DIR    <- .locate_dir()
PROJ_ROOT  <- normalizePath(file.path(OUT_DIR, ".."))              # project root
DATA_PATH  <- file.path(PROJ_ROOT, "On Data", "Data", "analytic_dataset_built.rds")
MODELS_DIR <- file.path(OUT_DIR, "models")
PLOTS_DIR  <- file.path(OUT_DIR, "plots")
dir.create(MODELS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR,  recursive = TRUE, showWarnings = FALSE)

CSV_PATH   <- file.path(OUT_DIR, "experiments.csv")
TERMS_PATH <- file.path(OUT_DIR, "experiments_terms.csv")
NOTE_PATH  <- file.path(OUT_DIR, "experiments_interpretation.md")
FROZEN_RDS <- file.path(OUT_DIR, "frozen_spatial.rds")            # spatial-term decision (frozen in 01)

# -----------------------------------------------------------------------------
# 1. PACKAGES  (library only -- never install)
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(mgcv)        # bam / gam, the whole modelling engine
  library(DHARMa)      # simulation-based residual diagnostics (dispersion / zero-inflation)
  library(ggplot2)     # plotting
})
HAS_GRATIA <- requireNamespace("gratia", quietly = TRUE)   # nicer smooth-effect plots (95% CI)
HAS_APE    <- requireNamespace("ape",    quietly = TRUE)   # Moran's I for residual spatial autocorrelation

options(width = 130)
set.seed(2025)                                   # determinism for every simulation (DHARMa etc.)
NTHREADS <- parallel::detectCores()              # ENGINE RULE: nthreads = <available>
if (is.na(NTHREADS) || NTHREADS < 1) NTHREADS <- 1L

# diagnostic sub-sampling knobs (n = 461k makes exhaustive simulation/Moran wasteful;
# all sub-samples are SEEDED so the run is reproducible -- documented in notes/CSV).
DHARMA_SUB  <- 40000L   # rows used for DHARMa residual tests
DHARMA_NSIM <- 100L     # simulations per DHARMa test
MORAN_SUB   <- 2500L    # hexes used for Moran's I

banner <- function(txt) cat("\n", strrep("=", 78), "\n## ", txt, "\n", strrep("=", 78), "\n", sep = "")
say    <- function(...) cat(sprintf(...), "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
# 2. DATA  --  load READ-ONLY, derive modelling responses, set types
# =============================================================================
banner("00_setup :: load analytic dataset (READ-ONLY) and derive responses")
if (!file.exists(DATA_PATH)) stop("Analytic dataset not found (read-only): ", DATA_PATH)
dat <- readRDS(DATA_PATH)                        # READ-ONLY: this object is never written back
say("Loaded %s", DATA_PATH)
say("raw dim: %d rows x %d cols", nrow(dat), ncol(dat))

# --- RESPONSE CONSTRUCTION (derived in-script from the data; source file untouched) ----------
#  Vessel scope = CARGO (C) + OTHER (O) ONLY.  We EXCLUDE fishing (F), tankers (T),
#  long-fast (L) and the all-combined (A) aggregate.  The C and O vessel classes are
#  DISJOINT partitions of the fleet, so summing their per-cell tallies is valid (no
#  double-counting): a cargo vessel cannot also be counted as an "other" vessel.
VESSEL_SCOPE <- "CARGO+OTHER"
dat$nOprD_CO <- dat$nOprD_C + dat$nOprD_O         # INTENSITY response: operating vessel-days (Tweedie / dispersion tracks)
dat$nMMSI_CO <- dat$nMMSI_C + dat$nMMSI_O         # COUNT response: unique vessels  (zero-inflation track)

# --- types --------------------------------------------------------------------
dat$hexID    <- factor(dat$hexID)                 # spatial random-effect grouping factor
dat$post2017 <- factor(ifelse(dat$post2017 == 1, "post", "pre"), levels = c("pre", "post"))  # parametric regulatory step
dat$month    <- as.numeric(dat$month)             # KEEP NUMERIC 1..12 for the cyclic seasonal term
# time_index: a numeric, monotone month counter 1..72 (months since study start 2015-01).
#  Distinct from `month`: `month` carries the wrap-around season, `time_index` the secular trend.
dat$time_index <- as.numeric((dat$year - min(dat$year)) * 12L + dat$month)

# --- observation summary (printed to the run log) -----------------------------
zp <- function(x) mean(x == 0)
say("\nObservation summary")
say("  rows (hex-months)        : %s", format(nrow(dat), big.mark = ","))
say("  unique hexes             : %s", format(nlevels(dat$hexID), big.mark = ","))
say("  unique months (dates)    : %d", length(unique(dat$date)))
say("  time span                : %s .. %s", as.character(min(dat$date)), as.character(max(dat$date)))
say("  time_index range         : %g .. %g", min(dat$time_index), max(dat$time_index))
say("  nOprD_CO  zero-proportion : %.4f   mean %.2f   var/mean %.1f   max %g",
    zp(dat$nOprD_CO), mean(dat$nOprD_CO), var(dat$nOprD_CO)/mean(dat$nOprD_CO), max(dat$nOprD_CO))
say("  nMMSI_CO  zero-proportion : %.4f   mean %.2f   var/mean %.1f   max %g",
    zp(dat$nMMSI_CO), mean(dat$nMMSI_CO), var(dat$nMMSI_CO)/mean(dat$nMMSI_CO), max(dat$nMMSI_CO))
say("  post2017 table           : pre=%d post=%d", sum(dat$post2017=="pre"), sum(dat$post2017=="post"))

# =============================================================================
# 3. GLOBAL SMOOTH SPECIFICATION  (do NOT deviate -- cited as a hard requirement)
# =============================================================================
#  Every CONTINUOUS covariate effect is a cubic B-spline (P-spline) basis with a
#  SECOND-ORDER difference penalty:  s(x, bs="ps", m=c(2,2), k=20).  One definition,
#  reused for ice/temperature/wind/port everywhere.
PS <- function(x, k = 20) sprintf('s(%s, bs = "ps", m = c(2, 2), k = %d)', x, k)

# Seasonality is the ONE labelled EXCEPTION to the ps rule:
# EXCEPTION to ps: cyclic spline for month (Dec->Jan continuity) -- justified in report.
SEASON <- 's(month, bs = "cc", k = 12)'
# cc boundary knots set so the period wraps with December(12) and January(1) ADJACENT
# (one step apart), not identical.  Passed as `knots=` to every fit carrying the season term.
MONTH_KNOTS <- list(month = c(0.5, 12.5))

# Secular trend: ps with second-order penalty (k=10 is ample for 72 monthly points).
TREND <- 's(time_index, bs = "ps", m = c(2, 2), k = 10)'

# Spatial term, two candidate forms (bake-off in 01, then FROZEN).
SPAT_RE <- 's(hexID, bs = "re")'                                           # (a) i.i.d. hex random effect
SPAT_TE <- 'te(centroid_x, centroid_y, bs = "ps", k = c(20, 20))'          # (b) tensor-product P-spline
#  NB: the "ps" basis family is PRESERVED for the 2D spatial term too (te marginals are ps),
#      so the whole model stays in one penalised-B-spline family bar the cyclic-month exception.

# FINAL mean predictor (RHS) given a chosen spatial term.  Built once here so that
# scripts 02/03/04 reconstruct the EXACT predictor selected in 01 (full consistency).
#  Temperature enters as air_temp_resid -- the ICE-ORTHOGONALISED residual chosen in 01
#  (raw air_temp is collinear with ice; the residual is carried forward, see 01 note).
#  Distance-to-port enters in km (port_dist); a log version (port_dist_log) exists in the
#  data but the ps spline already absorbs the right-skew, so km keeps the x-axis interpretable.
build_mean_rhs <- function(spatial_term = NULL) {
  # spatial_term="" / NULL omits the spatial smooth (used by the gam "drop spatial" fallback);
  # term ORDER is otherwise fixed so every script builds the identical predictor.
  parts <- c(SEASON, TREND,
             if (!is.null(spatial_term) && nzchar(spatial_term)) spatial_term,
             PS("ice_conc"), PS("air_temp_resid"), PS("wind_speed"), PS("port_dist"), "post2017")
  paste(parts, collapse = " + ")
}

# =============================================================================
# 4. ENGINE  --  the bam/gam switch (labelled; this is a reported design choice)
# =============================================================================
#  ENGINE RULE:
#    * bam(method="fREML", discrete=TRUE, nthreads=<available>) for every SINGLE-parameter
#      family: gaussian, Tweedie tw(), negative binomial nb().
#    * gam(method="REML") ONLY where bam cannot fit the family: ziplss(), twlss().
fit_model <- function(formula, family, engine = c("bam", "gam"), knots = NULL, method = NULL) {
  engine <- match.arg(engine)
  # auto-attach cyclic-month knots whenever the season term is present
  if (is.null(knots) && any(grepl("month", all.vars(if (is.list(formula)) formula[[1]] else formula))))
    knots <- MONTH_KNOTS
  warns <- character(0)
  t0 <- Sys.time()
  m <- withCallingHandlers(
    tryCatch({
      if (engine == "bam") {
        # bam: single-parameter families -- fREML + discretisation + threading
        bam(formula, family = family, data = dat,
            method = "fREML", discrete = TRUE, nthreads = NTHREADS, knots = knots)
      } else {
        # gam: location-scale families bam cannot fit. Default REML; twlss is mgcv-restricted
        # to the extended Fellner-Schall method, so 03 passes method="efs" (still a Wood
        # penalised-likelihood smoothness selector -- NOT cross-validation / black-box).
        gam(formula, family = family, data = dat, method = method %||% "REML", knots = knots)
      }
    }, error = function(e) { say("  [FIT ERROR] %s", conditionMessage(e)); NULL }),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  rt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(model = m, runtime = rt, warnings = warns,
       converged = !is.null(m) && isTRUE(m$converged))
}

# =============================================================================
# 5. METRIC HELPERS
# =============================================================================
# Estimated Tweedie power p (1<p<2).  We NEVER fix p; tw() estimates it -- here we read it back.
get_tweedie_p <- function(m) {
  f <- m$family$family
  if (!grepl("Tweedie", f, ignore.case = TRUE)) return(NA_real_)
  p <- tryCatch(as.numeric(m$family$getTheta(TRUE)), error = function(e) NA_real_)
  if (is.na(p) || p <= 1 || p >= 2) {                       # fallback: parse the family string
    g <- regmatches(f, regexpr("[0-9]+\\.[0-9]+", f)); if (length(g)) p <- as.numeric(g)
  }
  p
}

# Optimised smoothing-criterion score (fREML for bam, REML for gam).
get_score <- function(m, method) {
  v <- tryCatch(as.numeric(m$gcv.ubre), error = function(e) NA_real_)
  list(name = method, value = v)
}

# Maximum pairwise concurvity (worst-case "estimate" off-diagonal).  May be undefined
# for location-scale families -> returns NA, which is recorded honestly.
max_concurvity <- function(m) {
  # concurvity on a very high-dimensional term (e.g. the 6,409-level hex random effect) is both
  # scientifically meaningless (i.i.d. intercepts) and extremely slow, so skip it for such models.
  if (length(stats::coef(m)) > 2000) return(NA_real_)
  cc <- tryCatch(concurvity(m, full = FALSE), error = function(e) NULL)
  if (is.null(cc)) return(NA_real_)
  M <- cc$estimate; diag(M) <- NA
  suppressWarnings(max(M, na.rm = TRUE))
}

# k.check summary: smallest k-index and smallest basis-adequacy p across smooths.
kcheck_summary <- function(m) {
  kc <- tryCatch(mgcv::k.check(m), error = function(e) NULL)
  if (is.null(kc)) return(list(k_index_min = NA_real_, k_p_min = NA_real_, table = NULL))
  ki <- if ("k-index" %in% colnames(kc)) suppressWarnings(min(kc[, "k-index"], na.rm = TRUE)) else NA_real_
  kp <- if ("p-value" %in% colnames(kc)) suppressWarnings(min(kc[, "p-value"], na.rm = TRUE)) else NA_real_
  list(k_index_min = ki, k_p_min = kp, table = kc)
}

# -----------------------------------------------------------------------------
# DHARMa simulation-based dispersion + zero-inflation tests.
#  DHARMa cannot natively simulate from mgcv objects, so we simulate from the fitted
#  model BY HAND using the family's generative process, then createDHARMa().  Done on a
#  SEEDED row sub-sample for tractability at n=461k.  `fam` selects the simulator;
#  unsupported families (e.g. twlss) return NA and the caller falls back to mgcv
#  deviance-residual diagnostics.
# -----------------------------------------------------------------------------
dharma_tests <- function(m, fam) {
  if (is.na(fam) || !(fam %in% c("tweedie", "nb", "gaussian", "ziplss")))
    return(list(disp_p = NA_real_, zero_p = NA_real_))
  set.seed(2025)
  idx <- sort(sample.int(nrow(dat), min(DHARMA_SUB, nrow(dat))))
  res <- tryCatch({
    if (fam == "ziplss") {
      lp     <- predict(m, type = "link")[idx, , drop = FALSE]   # 2 columns: [count lp, zero lp]
      lambda <- exp(lp[, 1])                                     # Poisson mean (count process)
      pi_pre <- 1 - exp(-exp(lp[, 2]))                           # P(presence): cloglog 2nd lp (mgcv ziplss)
      mu     <- pi_pre * lambda                                  # expected observed count
      yobs   <- m$y[idx]
      sims   <- replicate(DHARMA_NSIM, rbinom(length(idx), 1, pi_pre) * rpois(length(idx), lambda))
      integerResp <- TRUE
    } else {
      mu   <- as.numeric(m$fitted.values)[idx]
      yobs <- m$y[idx]
      if (fam == "tweedie") {
        # NB: bam stores the estimated Tweedie dispersion in $sig2 (== summary()$dispersion);
        #     $scale is only a logical "scale estimated?" flag on bam objects.
        p <- get_tweedie_p(m); phi <- m$sig2
        sims <- replicate(DHARMA_NSIM, mgcv::rTweedie(mu, p = p, phi = phi)); integerResp <- FALSE
      } else if (fam == "nb") {
        theta <- m$family$getTheta(TRUE)
        sims <- replicate(DHARMA_NSIM, rnbinom(length(idx), mu = mu, size = theta)); integerResp <- TRUE
      } else { # gaussian
        sig <- sqrt(m$sig2)
        sims <- replicate(DHARMA_NSIM, rnorm(length(idx), mu, sig)); integerResp <- FALSE
      }
    }
    dh <- DHARMa::createDHARMa(simulatedResponse = sims, observedResponse = yobs,
                               fittedPredictedResponse = mu, integerResponse = integerResp)
    list(disp_p = DHARMa::testDispersion(dh, plot = FALSE)$p.value,
         zero_p = DHARMa::testZeroInflation(dh, plot = FALSE)$p.value)
  }, error = function(e) { say("  [DHARMa skipped: %s]", conditionMessage(e));
                           list(disp_p = NA_real_, zero_p = NA_real_) })
  res
}

# Residual SPATIAL autocorrelation: Moran's I on per-hex mean deviance residuals
#  (inverse-distance weights, seeded hex sub-sample).  Returns the I statistic + p-value.
resid_spatial_moran <- function(m) {
  if (!HAS_APE) return(list(I = NA_real_, p = NA_real_, agg = NULL))
  dr  <- residuals(m, type = "deviance")
  agg <- tapply(dr, dat$hexID, mean)                       # per-hex mean residual
  hx  <- !duplicated(dat$hexID)
  ids <- as.character(dat$hexID[hx]); cx <- dat$centroid_x[hx]; cy <- dat$centroid_y[hx]
  ord <- match(names(agg), ids); cx <- cx[ord]; cy <- cy[ord]; val <- as.numeric(agg)
  res <- tryCatch({
    set.seed(2025)
    if (length(val) > MORAN_SUB) { s <- sample.int(length(val), MORAN_SUB); val<-val[s]; cx<-cx[s]; cy<-cy[s] }
    D <- as.matrix(dist(cbind(cx, cy))); W <- 1/D; diag(W) <- 0; W[!is.finite(W)] <- 0
    mi <- ape::Moran.I(val, W)
    list(I = mi$observed, p = mi$p.value)
  }, error = function(e) { say("  [Moran skipped: %s]", conditionMessage(e)); list(I=NA_real_, p=NA_real_) })
  c(res, list(agg = agg))
}

# Residual TEMPORAL autocorrelation: ACF + Ljung-Box on monthly-aggregated residuals.
resid_temporal_acf <- function(m) {
  dr  <- residuals(m, type = "deviance")
  agg <- tapply(dr, dat$time_index, mean)
  agg <- agg[order(as.numeric(names(agg)))]
  bx  <- tryCatch(Box.test(agg, lag = 12, type = "Ljung-Box")$p.value, error = function(e) NA_real_)
  a1  <- tryCatch(as.numeric(acf(agg, plot = FALSE, lag.max = 1)$acf[2]), error = function(e) NA_real_)
  list(lag1 = a1, boxp = bx, series = as.numeric(agg))
}

# =============================================================================
# 6. OUTPUT WRITERS  (experiments.csv, experiments_terms.csv, interpretation md)
# =============================================================================
EXP_COLS <- c("exp_id","script","response","vessel_scope","family","engine","method",
              "formula","spatial_term","k_spec","tweedie_p_est","score_name","score_value",
              "AIC","BIC","deviance_explained_pct","n_obs","total_edf","max_pairwise_concurvity",
              "dharma_dispersion_p","dharma_zeroinflation_p","dharma_spatial_autocorr",
              "dharma_temporal_autocorr","k_index_min","k_p_min","convergence_status",
              "runtime_sec","notes")
TERM_COLS <- c("exp_id","term","edf","ref_df","statistic","p_value","k_prime","k_index","k_check_p")

# start_fresh(): called ONCE by 01 to (re)initialise the three output files with headers
#  ("Ignore any pre-existing experiments.csv; start fresh.").  02/03/04 only ever append.
start_fresh <- function() {
  writeLines(paste(EXP_COLS,  collapse = ","), CSV_PATH)
  writeLines(paste(TERM_COLS, collapse = ","), TERMS_PATH)
  writeLines(c("# Experiments -- interpretation log",
               sprintf("_Arctic shipping GAM sequence; generated %s; vessel scope %s._",
                       as.character(Sys.Date()), VESSEL_SCOPE),
               "", "Claims are ASSOCIATIONAL (conditional association, not causal).",
               "Significance is de-emphasised given the very large n (461,448).", ""),
             NOTE_PATH)
  say("Initialised fresh: experiments.csv, experiments_terms.csv, experiments_interpretation.md")
}

# robust one-row appender (auto-writes the header if a file is run out of order)
.append_row <- function(path, header, named_values) {
  if (!file.exists(path)) writeLines(paste(header, collapse = ","), path)
  df <- as.data.frame(as.list(named_values), stringsAsFactors = FALSE)[, header, drop = FALSE]
  write.table(df, path, sep = ",", append = TRUE, col.names = FALSE,
              row.names = FALSE, na = "NA", qmethod = "double")
}

# note(): append a labelled interpretation paragraph to the markdown log.
note <- function(title, lines) {
  cat(c(sprintf("\n### %s\n", title), paste0(lines, collapse = "\n"), "\n"),
      file = NOTE_PATH, append = TRUE, sep = "")
}

# append_terms(): per-term EDF + per-term k.check -> the long companion file.
append_terms <- function(exp_id, m) {
  s  <- summary(m)
  kc <- kcheck_summary(m)$table
  kc_df <- if (!is.null(kc)) data.frame(term = rownames(kc), kprime = kc[,1], kindex = kc[,"k-index"],
                                        kp = kc[,"p-value"], stringsAsFactors = FALSE) else
                             data.frame(term=character(0), kprime=numeric(0), kindex=numeric(0), kp=numeric(0))
  glk <- function(t) { i <- match(t, kc_df$term); if (is.na(i)) c(NA,NA,NA) else
                       c(kc_df$kprime[i], kc_df$kindex[i], kc_df$kp[i]) }
  # smooth terms
  if (!is.null(s$s.table) && nrow(s$s.table) > 0) {
    st <- s$s.table
    for (r in seq_len(nrow(st))) {
      k <- glk(rownames(st)[r])
      .append_row(TERMS_PATH, TERM_COLS, list(exp_id=exp_id, term=rownames(st)[r],
        edf=round(st[r,"edf"],3), ref_df=round(st[r,"Ref.df"],3),
        statistic=round(st[r,3],3), p_value=signif(st[r,4],4),
        k_prime=k[1], k_index=round(k[2],3), k_check_p=signif(k[3],4)))
    }
  }
  # parametric terms (e.g. post2017)
  if (!is.null(s$p.table) && nrow(s$p.table) > 0) {
    pt <- s$p.table
    for (r in seq_len(nrow(pt))) {
      .append_row(TERMS_PATH, TERM_COLS, list(exp_id=exp_id, term=rownames(pt)[r],
        edf=1, ref_df=1, statistic=round(pt[r,3],3), p_value=signif(pt[r,4],4),
        k_prime=NA, k_index=NA, k_check_p=NA))
    }
  }
}

# =============================================================================
# 7. DIAGNOSTIC PLOTS  (consistent <exp_id>_<type>.png naming, saved to /plots)
# =============================================================================
ppath <- function(exp_id, type) file.path(PLOTS_DIR, sprintf("%s_%s.png", exp_id, type))

save_diagnostics <- function(exp_id, m, moran_agg = NULL, acf_series = NULL) {
  # (a) smooth-effect plots with 95% CI -- one panel per smooth
  ok <- FALSE
  if (HAS_GRATIA) ok <- tryCatch({
    # scales="free": each smooth gets its OWN y-axis. With the default "fixed", the large
    # ice effect sets a common range and visually flattens the smaller smooths (wind/port/temp).
    p <- gratia::draw(m, ci_level = 0.95, scales = "free")
    ggplot2::ggsave(ppath(exp_id, "smooths"), p, width = 12, height = 9, dpi = 120); TRUE
  }, error = function(e) FALSE)
  if (!ok) tryCatch({                                   # fallback: base plot.gam
    grDevices::png(ppath(exp_id, "smooths"), width = 1100, height = 800, res = 110)
    op <- par(mfrow = c(2, 3)); plot(m, pages = 0, shade = TRUE, seWithMean = TRUE); par(op); dev.off()
  }, error = function(e) say("  [smooth plot skipped: %s]", conditionMessage(e)))

  # (b) mgcv residual diagnostics: QQ + residual-vs-linear-predictor (sampled for speed)
  tryCatch({
    set.seed(2025); dr <- residuals(m, type = "deviance")
    lp <- if (is.matrix(m$linear.predictors)) m$linear.predictors[,1] else m$linear.predictors
    s  <- sample.int(length(dr), min(20000L, length(dr)))
    grDevices::png(ppath(exp_id, "resid"), width = 1000, height = 500, res = 110)
    op <- par(mfrow = c(1, 2))
    qqnorm(dr[s], main = "QQ deviance residuals", pch = 16, cex = .3, col = "#0B2E6D"); qqline(dr[s], col = "firebrick")
    plot(lp[s], dr[s], pch = 16, cex = .3, col = "#0B2E6D33",
         xlab = "linear predictor", ylab = "deviance residual", main = "Resid vs linear predictor")
    abline(h = 0, col = "firebrick"); par(op); dev.off()
  }, error = function(e) say("  [resid plot skipped: %s]", conditionMessage(e)))

  # (c) residual SPATIAL map: per-hex mean deviance residual at the hex centroid
  if (!is.null(moran_agg)) tryCatch({
    hx <- !duplicated(dat$hexID); ids <- as.character(dat$hexID[hx])
    md <- data.frame(x = dat$centroid_x[hx], y = dat$centroid_y[hx],
                     r = as.numeric(moran_agg[match(ids, names(moran_agg))]))
    lim <- max(abs(md$r), na.rm = TRUE)
    p <- ggplot(md, aes(x, y, colour = r)) + geom_point(size = .5) + coord_fixed() +
      scale_colour_gradient2(low = "#B2182B", mid = "grey90", high = "#2166AC",
                             midpoint = 0, limits = c(-lim, lim), name = "mean dev. resid") +
      theme_minimal() + labs(title = sprintf("%s -- residual map (per hex)", exp_id),
                             x = "centroid_x (m)", y = "centroid_y (m)")
    ggsave(ppath(exp_id, "residmap"), p, width = 7, height = 7, dpi = 110)
  }, error = function(e) say("  [resid map skipped: %s]", conditionMessage(e)))

  # (d) residual TEMPORAL ACF on monthly-aggregated residuals
  if (!is.null(acf_series)) tryCatch({
    grDevices::png(ppath(exp_id, "acf"), width = 800, height = 500, res = 110)
    acf(acf_series, main = sprintf("%s -- ACF of monthly residuals", exp_id)); dev.off()
  }, error = function(e) say("  [acf plot skipped: %s]", conditionMessage(e)))
}

# =============================================================================
# 8. record_experiment()  --  fit-agnostic: takes a FITTED model + metadata, computes
#    every CSV metric, writes ONE comparison row + the per-term rows, saves the standard
#    diagnostic plot set, and RETURNS all metrics (so the caller can write a data-driven
#    interpretation note and make AIC/BIC/concurvity comparisons).
# =============================================================================
record_experiment <- function(exp_id, script, m, response, family_label, engine, method,
                               formula_str, spatial_term, k_spec, dharma_family = NA,
                               convergence_status = NULL, runtime_sec = NA, extra_notes = "",
                               do_spatial = TRUE, do_temporal = TRUE, save_plots = TRUE,
                               tweedie_p_override = NA) {
  s  <- summary(m)
  sc <- get_score(m, method)
  kc <- kcheck_summary(m)
  dh <- dharma_tests(m, dharma_family)
  sp <- if (do_spatial) resid_spatial_moran(m) else list(I = NA, p = NA, agg = NULL)
  tm <- if (do_temporal) resid_temporal_acf(m) else list(lag1 = NA, boxp = NA, series = NULL)
  if (is.null(convergence_status)) convergence_status <- if (isTRUE(m$converged)) "converged" else "check"

  row <- list(
    exp_id = exp_id, script = script, response = response, vessel_scope = VESSEL_SCOPE,
    family = family_label, engine = engine, method = method,
    formula = formula_str, spatial_term = spatial_term, k_spec = k_spec,
    tweedie_p_est = round(if (!is.na(tweedie_p_override)) tweedie_p_override else get_tweedie_p(m), 4),
    score_name = sc$name, score_value = round(sc$value, 2),
    AIC = round(AIC(m), 2), BIC = round(BIC(m), 2),
    deviance_explained_pct = round(100 * s$dev.expl, 2),
    n_obs = nrow(dat), total_edf = round(sum(m$edf), 2),
    max_pairwise_concurvity = round(max_concurvity(m), 4),
    dharma_dispersion_p = signif(dh$disp_p, 4), dharma_zeroinflation_p = signif(dh$zero_p, 4),
    dharma_spatial_autocorr = signif(sp$p, 4), dharma_temporal_autocorr = signif(tm$boxp, 4),
    k_index_min = round(kc$k_index_min, 3), k_p_min = signif(kc$k_p_min, 4),
    convergence_status = convergence_status, runtime_sec = round(runtime_sec, 1),
    notes = extra_notes)
  .append_row(CSV_PATH, EXP_COLS, row)
  append_terms(exp_id, m)
  if (save_plots) save_diagnostics(exp_id, m, moran_agg = sp$agg, acf_series = tm$series)
  saveRDS(m, file.path(MODELS_DIR, paste0(exp_id, ".rds")))

  say("  [%s] AIC=%.0f BIC=%.0f dev.expl=%.1f%% edf=%.1f score(%s)=%.0f conc=%.3f | disp.p=%s zero.p=%s moranI=%s(p=%s) box.p=%s",
      exp_id, AIC(m), BIC(m), 100*s$dev.expl, sum(m$edf), sc$name, sc$value,
      max_concurvity(m), format(dh$disp_p), format(dh$zero_p),
      format(round(sp$I,3)), format(signif(sp$p,3)), format(signif(tm$boxp,3)))

  invisible(list(exp_id=exp_id, AIC=AIC(m), BIC=BIC(m), score=sc$value, dev_expl=100*s$dev.expl,
                 edf=sum(m$edf), tweedie_p=get_tweedie_p(m), concurvity=max_concurvity(m),
                 disp_p=dh$disp_p, zero_p=dh$zero_p, moran_I=sp$I, moran_p=sp$p,
                 acf_lag1=tm$lag1, box_p=tm$boxp, k_index_min=kc$k_index_min, k_p_min=kc$k_p_min))
}

# =============================================================================
# 9. SPATIAL-TERM FREEZE  (01 writes the winner; 02/03/04 read it back)
# =============================================================================
freeze_spatial <- function(winner, spatial_term, k_spec) {
  saveRDS(list(winner = winner, spatial_term = spatial_term, k_spec = k_spec), FROZEN_RDS)
  say("FROZEN spatial term -> %s : %s", winner, spatial_term)
}
get_frozen_spatial <- function() {
  if (file.exists(FROZEN_RDS)) return(readRDS(FROZEN_RDS))
  # safe default if 01 has not run yet: the leaner tensor smooth
  list(winner = "te(default)", spatial_term = SPAT_TE, k_spec = "te:20x20")
}

# =============================================================================
# 10. SESSION INFO  (printed every run, per spec)
# =============================================================================
banner("00_setup :: session info")
si <- sessionInfo()
say("R           : %s", si$R.version$version.string)
say("platform    : %s", si$platform)
say("mgcv        : %s", as.character(packageVersion("mgcv")))
say("DHARMa      : %s", as.character(packageVersion("DHARMa")))
say("gratia      : %s", if (HAS_GRATIA) as.character(packageVersion("gratia")) else "not installed")
say("ape         : %s", if (HAS_APE)    as.character(packageVersion("ape"))    else "not installed")
say("nthreads    : %d", NTHREADS)
say("twlss avail : %s", exists("twlss", where = asNamespace("mgcv")))
say("ziplss avail: %s", exists("ziplss", where = asNamespace("mgcv")))
say("00_setup ready.\n")
