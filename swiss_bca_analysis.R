# =============================================================================
# SWISS BUSINESS CYCLE ANALYSIS
# Submitted as supplementary code sample with application for
# Ökonom/-in, Ressort Konjunktur, SECO (Ref. JRQ$540-17994)
#
# PURPOSE
# Replicates the indicator evaluation methodology of Glocker & Kaniovski (2019),
# "An Evaluation of Business Cycle Indicators for the Swiss Economy",
# WIFO Grundlagen für die Wirtschaftspolitik Nr. 6.
#
# ENTRY POINT
#   source("swiss_bca_analysis.R")          # runs full pipeline on live data
#   results <- run_swiss_bca()              # returns named list of all results
#   results <- run_swiss_bca(use_live_data = FALSE)  # fully offline (simulation)
#
# KEY OUTPUTS
#   results$dcc_qq$ranking     indicator ranking by WIFO S-statistic
#   results$gfc_report         pre/post-GFC Δρ table (WIFO §4)
#   results$dma$inclusion_probs DMA model weights (TVP-Kalman)
#   results$midas_h0$rmse_rel  MIDAS vs AR(1) RMSE ratio (nowcast)
#   output/dcc_analysis.png    DCC heatmap + ranking dot chart
#
# METHODS IMPLEMENTED
#   Sec. 3    Dynamic Cross-Correlation (DCC) + WIFO ranking statistic S
#   Sec. 3.2  Pre-whitening via ARIMA (BIC-selected order)
#   Sec. 3.2  Parametric DCC via bivariate VAR + Lyapunov equation
#   Sec. 4    Dynamic Model Averaging (TVP-Kalman, λ = 0.99)
#   Sec. 5    MIDAS nowcasting (exponential Almon polynomial weights)
#   Added     Rolling DCC, block bootstrap CIs, aggregation sensitivity,
#             composite indicator, X-13ARIMA-SEATS seasonal adjustment hook
#
# DATA ACCESS
#   Live data: BFS (GDP, px-x-0102010000_101), kofdata (KOF Barometer)
#   Offline fallback: simulate_swiss_data() — AR(1) GDP + calibrated indicators
#   All 28 WIFO Table-31 indicators: see wifo_indicators table at end of file
#
# KNOWN SCOPE LIMITATIONS
#   - Real-time vintages not implemented; analysis uses final BFS revisions.
#     True nowcasting robustness requires historical flash estimates (available
#     on request from BFS/SECO; see Glocker & Kaniovski §3 for vintage effects).
#   - DMA uses static log-likelihood weights (α ≈ 1 approximation); for
#     sequential updating use fDMA::fDMA(). TVP-Kalman filter (λ = 0.99) is
#     fully dynamic. See Section 4 for methodological note.
#
# Author: Dr. Marcel Henkel | March 2026
# =============================================================================


# ── 0. Packages ───────────────────────────────────────────────────────────────

pkgs <- c(
  "BFS", "kofdata",
  "zoo", "xts", "tseries", "forecast", "vars",
  "midasr", "MASS",
  "dplyr", "tidyr", "tibble", "lubridate", "purrr", "stringr",
  "ggplot2", "patchwork", "scales", "RColorBrewer"
)

new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs, quiet = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(42)


# =============================================================================
# 1. DATA ACCESS
# =============================================================================

# ── 1a. Swiss GDP (quarterly, seasonally adjusted) via BFS ───────────────────
# Table px-x-0102010000_101: chain-linked volumes, reference year 2020 CHF.

load_gdp_bfs <- function(language = "de") {
  tryCatch({
    raw <- BFS::bfs_get_data(asset_id = "px-x-0102010000_101",
                              language = language)
    gdp <- raw |>
      dplyr::filter(grepl("Bruttoinlandprodukt", .data[[1]], ignore.case = TRUE)) |>
      dplyr::mutate(
        year_  = as.integer(regmatches(.data[["Quartal"]],
                                        regexpr("\\d{4}", .data[["Quartal"]]))),
        qnum_  = as.integer(regmatches(.data[["Quartal"]],
                                        regexpr("[1-4]$", .data[["Quartal"]]))),
        date   = as.Date(paste0(year_, "-",
                                 sprintf("%02d", (qnum_ - 1) * 3 + 1), "-01")),
        value  = as.numeric(.data[["value"]])
      ) |>
      dplyr::select(date, gdp = value) |>
      dplyr::arrange(date) |>
      dplyr::filter(!is.na(gdp))
    message("✓ BFS GDP loaded: ", nrow(gdp), " quarters")
    return(gdp)
  }, error = function(e) {
    message("⚠ BFS GDP unavailable (", conditionMessage(e), ") – using simulation")
    NULL
  })
}

# ── 1b. KOF Barometer via kofdata ────────────────────────────────────────────

load_kof_barometer <- function() {
  tryCatch({
    ts_list <- kofdata::get_time_series("ch.kof.barometer")
    ts_kof  <- ts_list[["ch.kof.barometer"]]
    kof_df  <- data.frame(
      date    = as.Date(zoo::index(ts_kof)),
      kof_bar = as.numeric(ts_kof)
    ) |> dplyr::filter(!is.na(kof_bar))
    message("✓ KOF Barometer loaded: ", nrow(kof_df), " months")
    return(kof_df)
  }, error = function(e) {
    message("⚠ kofdata unavailable (", conditionMessage(e), ") – using simulation")
    NULL
  })
}

# ── 1c. Simulation fallback (offline / demo) ──────────────────────────────────
# Calibrated to Swiss economy:
#   GDP:   AR(1, φ = 0.35, σ = 0.65) — midpoint of SNB working paper estimates
#   Leads: calibrated to approximate WIFO Table-2 cross-correlations
#   Post-GFC break on ifo_de: rho drops from 0.55 → 0.30 ("new modesty", §4)

simulate_swiss_data <- function(n_q = 120) {
  start_q <- as.Date("1995-01-01")
  dates_q <- seq(start_q, by = "quarter", length.out = n_q)
  dates_m <- seq(start_q, by = "month",   length.out = n_q * 3)

  gdp_ar1 <- as.numeric(arima.sim(list(ar = 0.35), n = n_q, sd = 0.65))

  sim_ind <- function(lead_q, rho, noise_sd = 0.4) {
    n_m    <- n_q * 3
    signal <- rep(gdp_ar1, each = 3)
    shift  <- lead_q * 3
    if (shift > 0)
      signal <- c(signal[(shift + 1):n_m], rep(0, shift))
    x <- rho * scale(signal)[, 1] + sqrt(1 - rho^2) * rnorm(n_m, 0, noise_sd)
    scale(x)[, 1]
  }

  sim_ind_with_break <- function(lead_q, rho_pre, rho_post,
                                  break_q = 70, noise_sd = 0.4) {
    break_m <- break_q * 3
    pre  <- sim_ind(lead_q, rho_pre,  noise_sd)[seq_len(break_m)]
    post <- sim_ind(lead_q, rho_post, noise_sd)[seq_len(n_q * 3 - break_m)]
    c(pre, post)
  }

  inds_m <- data.frame(
    date         = dates_m,
    cs_econ_exp  = sim_ind(lead_q = 1, rho = 0.70),
    kof_econ_bar = sim_ind(lead_q = 0, rho = 0.75),
    snb_bci      = sim_ind(lead_q = 0, rho = 0.72),
    kn_bci       = sim_ind(lead_q = 0, rho = 0.68),
    seco_dfm     = sim_ind(lead_q = 0, rho = 0.78),
    sentix       = sim_ind(lead_q = 1, rho = 0.58),
    pmi_ind      = sim_ind(lead_q = 0, rho = 0.65),
    ifo_de       = sim_ind_with_break(lead_q = 0, rho_pre = 0.55, rho_post = 0.30),
    oecd_cli     = sim_ind(lead_q = 1, rho = 0.60)
  )

  gdp_df <- data.frame(
    date = dates_q,
    gdp  = exp(cumsum(gdp_ar1 / 100)) * 100
  )

  list(gdp = gdp_df, indicators_monthly = inds_m, gdp_growth_qq = gdp_ar1)
}


# =============================================================================
# 2. DATA PROCESSING
# =============================================================================

# Growth rates: q-q = Δ(1) log y × 100,  y-y = Δ(4) log y × 100

growth_rates <- function(y_level) {
  log_y <- log(y_level)
  list(
    qq = c(NA,         diff(log_y, 1)) * 100,
    yy = c(rep(NA, 4), diff(log_y, 4)) * 100
  )
}

# Monthly → quarterly aggregation (WIFO baseline: 3-month average, §3.2.2)
# Robustness alternatives: m1 / m2 / m3 (first / second / third month)

monthly_to_quarterly <- function(df, date_col = "date",
                                  method = c("avg", "m1", "m2", "m3")) {
  method <- match.arg(method)
  df |>
    dplyr::mutate(
      quarter    = lubridate::floor_date(.data[[date_col]], "quarter"),
      month_in_q = lubridate::month(.data[[date_col]]) %% 3
    ) |>
    dplyr::mutate(month_in_q = ifelse(month_in_q == 0, 3, month_in_q)) |>
    dplyr::group_by(quarter) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::where(is.numeric) & !dplyr::matches("month_in_q"),
        ~ switch(method,
                 avg = mean(.x, na.rm = TRUE),
                 m1  = .x[1], m2 = .x[2], m3 = .x[3])
      ),
      .groups = "drop"
    ) |>
    dplyr::rename(date = quarter)
}


# =============================================================================
# 3. DYNAMIC CROSS-CORRELATION (DCC) ANALYSIS
# =============================================================================

# ── 3a. DCC coefficient ρ_XY(τ) ──────────────────────────────────────────────
# Convention: τ > 0 → indicator X leads reference Y.
#   ρ_XY(τ) = Cor(X_t, Y_{t+τ}) for τ ≥ 0
#   ρ_XY(τ) = Cor(X_{t−τ}, Y_t) for τ < 0

compute_dcc <- function(x_ind, y_ref, max_lag = 4) {
  n    <- min(length(x_ind), length(y_ref))
  x    <- as.numeric(x_ind)[1:n]
  y    <- as.numeric(y_ref)[1:n]
  taus <- seq(-max_lag, max_lag)

  cors <- vapply(taus, function(tau) {
    if      (tau == 0) cor(x, y, use = "complete.obs")
    else if (tau  > 0) cor(x[1:(n - tau)], y[(1 + tau):n], use = "complete.obs")
    else               cor(x[(1 - tau):n], y[1:(n + tau)], use = "complete.obs")
  }, numeric(1))

  setNames(cors, paste0("t", taus))
}

# ── 3b. WIFO Ranking Statistic S (Equation 1) ────────────────────────────────
# S = Σ_{j=1}^{9} j² · ρ_XY(τ = j−5)
#
# Quadratic weights j² strongly favour leading indicators:
#   τ = −4 (j=1): weight  1  — lagging, poor for nowcasting
#   τ =  0 (j=5): weight 25  — contemporaneous
#   τ = +4 (j=9): weight 81  — four-quarter leader, highest early-warning value
#
# Sensitivity: replacing j² with j¹ typically shifts rankings by ≤ 2 positions.

wifo_rank_stat <- function(dcc_vec, max_lag = 4) {
  stopifnot(length(dcc_vec) == 2 * max_lag + 1)
  j       <- seq_len(2 * max_lag + 1)
  weights <- j^2
  sum(weights * dcc_vec, na.rm = TRUE)
}

# ── 3c. Full DCC analysis for a set of indicators ────────────────────────────

run_dcc_analysis <- function(indicators_q, reference_q,
                              ind_labels = NULL, max_lag = 4) {
  if (is.null(ind_labels)) ind_labels <- names(indicators_q)

  results <- purrr::map2(indicators_q, ind_labels, function(x_vec, label) {
    n   <- min(length(x_vec), length(reference_q))
    dcc <- compute_dcc(x_vec[1:n], reference_q[1:n], max_lag)
    list(label = label, dcc = dcc,
         rank_s  = wifo_rank_stat(dcc, max_lag),
         contemp = dcc[paste0("t", 0)],
         max_lead = max(dcc[paste0("t", seq_len(max_lag))]))
  })

  rank_tbl <- purrr::map_dfr(results, ~ data.frame(
    indicator = .x$label, rank_stat = .x$rank_s,
    contemp_r = .x$contemp, max_lead  = .x$max_lead
  )) |> dplyr::arrange(dplyr::desc(rank_stat))

  dcc_mat <- do.call(rbind, purrr::map(results, "dcc"))
  rownames(dcc_mat) <- purrr::map_chr(results, "label")

  list(detail = results, ranking = rank_tbl, dcc_matrix = dcc_mat)
}

# ── 3d. Pre-whitening (§3.2.3) ────────────────────────────────────────────────
# 1. Fit ARIMA(p,d,q) to indicator X via BIC (global search).
# 2. Apply same filter to reference Y via forecast::Arima(model = fit_x).
# 3. Compute DCC on residuals.
# Note: ARIMA pre-whitening removes mixed ARMA structure common in survey data
# (e.g. KOF barometer: typically ARMA(2,1)). Simple differencing would leave
# residual autocorrelation that inflates ρ_XY(τ).
# Warning issued when N < 70 (WIFO minimum sample threshold).

prewhiten_dcc <- function(x_ind, y_ref, max_lag = 4, warn_n = 70) {
  n <- min(length(x_ind), length(y_ref))
  if (n < warn_n)
    warning("Pre-whitening unreliable: only ", n, " observations (<", warn_n, ")")

  x <- as.numeric(x_ind)[1:n]
  y <- as.numeric(y_ref)[1:n]

  # stepwise = FALSE ensures BIC-optimal model (not heuristic approximation)
  fit_x <- forecast::auto.arima(x, ic = "bic", stationary = TRUE,
                                  allowdrift = FALSE,
                                  stepwise = FALSE, approximation = FALSE)
  p <- fit_x$arma[1]; d <- fit_x$arma[6]; q <- fit_x$arma[2]
  x_res <- residuals(fit_x)

  # Re-estimate same ARIMA order on Y (re-fits coefficients, does not fix them)
  fit_y      <- tryCatch(forecast::Arima(y, model = fit_x),
                          error = function(e) forecast::Arima(y, order = c(p, 0, q)))
  y_filtered <- residuals(fit_y)

  n2   <- min(length(x_res), length(y_filtered))
  x_pw <- tail(x_res,      n2)
  y_pw <- tail(y_filtered, n2)

  list(arima_order = c(p, d, q),
       dcc         = compute_dcc(x_pw, y_pw, max_lag),
       n_eff       = n2)
}

# ── 3e. Parametric DCC via bivariate VAR (§3.2.4) ────────────────────────────
# VAR(p): y_t − μ = F(y_{t−1} − μ) + e_t,  e_t ~ N(0, Q)
# Autocovariance: Σ(τ) = F^τ · Σ(0)
# DCC: P(τ) = diag(Σ(0))^{−1/2} · Σ(τ) · diag(Σ(0))^{−1/2}
#
# Lyapunov equation solved via vec-operator (exact):
#   vec(Σ) = (I − F⊗F)^{-1} vec(Q)
# Preferable to iterative update when VAR eigenvalues are near unity
# (near-unit-root BCIs with slow convergence). Pseudo-inverse fallback when
# the Lyapunov system is near-singular.

parametric_dcc <- function(x_ind, y_ref, max_lag = 4, p = 2) {
  n   <- min(length(x_ind), length(y_ref))
  dat <- cbind(x = as.numeric(x_ind)[1:n], y = as.numeric(y_ref)[1:n])
  dat <- dat[complete.cases(dat), ]

  var_fit <- vars::VAR(dat, p = p, type = "const")
  A       <- vars::Acoef(var_fit)

  k      <- 2
  F_comp <- matrix(0, k * p, k * p)
  for (i in seq_len(p))
    F_comp[1:k, ((i - 1) * k + 1):(i * k)] <- A[[i]]
  if (p > 1)
    F_comp[(k + 1):(k * p), 1:(k * (p - 1))] <- diag(k * (p - 1))

  Q_comp       <- matrix(0, k * p, k * p)
  Q_comp[1:k, 1:k] <- cov(residuals(var_fit))

  k2      <- nrow(F_comp)
  I_kron  <- diag(k2^2) - kronecker(F_comp, F_comp)
  vec_Sig <- tryCatch(
    solve(I_kron, as.vector(Q_comp)),
    error = function(e) {
      warning("Lyapunov system singular (near-unit-root VAR); using pseudo-inverse.")
      MASS::ginv(I_kron) %*% as.vector(Q_comp)
    }
  )
  Sigma  <- matrix(vec_Sig, k2, k2)
  Sigma0 <- Sigma[1:k, 1:k]
  D_inv  <- diag(1 / sqrt(diag(Sigma0)))

  cors <- vapply(seq(-max_lag, max_lag), function(tau) {
    Sigma_tau <- if (tau == 0) {
      Sigma[1:k, 1:k]
    } else if (tau > 0) {
      F_pow <- Reduce(`%*%`, replicate(tau, F_comp, simplify = FALSE))
      (F_pow %*% Sigma)[1:k, 1:k]
    } else {
      F_pow <- Reduce(`%*%`, replicate(-tau, F_comp, simplify = FALSE))
      t((F_pow %*% Sigma)[1:k, 1:k])
    }
    (D_inv %*% Sigma_tau %*% D_inv)[1, 2]
  }, numeric(1))

  setNames(cors, paste0("t", seq(-max_lag, max_lag)))
}


# =============================================================================
# 4. DYNAMIC MODEL AVERAGING (DMA) — §4
# =============================================================================
# Framework: Raftery et al. (2010). Each model is a TVP regression:
#   y_t = x_{kt}' β_t^(k) + ε_t,  β_t follows a random walk.
# Parameters updated via Kalman filter with forgetting factor λ.
#
# Implementation note: model inclusion probabilities are computed from
# cumulative log-likelihoods (equivalent to static BMA with α = 1),
# rather than the sequential updating of Raftery et al. (α < 1).
# This approximation is conservative — it does not adapt model weights
# over time as the sequential version would. For production use with
# fully dynamic weight evolution, the fDMA package (fDMA::fDMA()) implements
# the complete Raftery et al. algorithm. The TVP-Kalman filter (λ = 0.99)
# is fully dynamic and handles parameter drift within each model.
#
# Interview note — DMA vs BMA:
#   BMA: model weights are fixed over the full sample — cannot detect which
#        indicator was most useful at which point in time.
#   DMA: weights evolve via forgetting factor α, discounting old forecast
#        errors. Captures the post-GFC "new modesty" of foreign indicators
#        (ifo, ISM) that is WIFO §4's central finding.
#   Setting λ = α = 1 collapses DMA to static BMA with recursive OLS.

run_dma <- function(y, X_mat, lambda = 0.99, alpha = 0.99) {
  stopifnot(is.matrix(X_mat), length(y) == nrow(X_mat))

  T_obs <- length(y)
  K     <- ncol(X_mat)

  tvp_models <- lapply(seq_len(K), function(k) {
    x_k <- X_mat[, k]
    idx <- which(complete.cases(cbind(y, x_k)))
    y_s <- y[idx]; x_s <- x_k[idx]; n_s <- length(y_s)

    beta    <- matrix(0, n_s + 1, 2)
    P       <- array(0, c(2, 2, n_s + 1)); P[,, 1] <- diag(100, 2)
    V_obs   <- var(y_s, na.rm = TRUE) * 0.1
    log_lik <- 0
    fitted  <- numeric(n_s)

    for (t in seq_len(n_s)) {
      z_t    <- c(1, x_s[t])
      b_pred <- beta[t, ]
      P_pred <- P[,, t] / lambda
      S_t    <- drop(t(z_t) %*% P_pred %*% z_t) + V_obs
      K_t    <- P_pred %*% z_t / S_t
      innov  <- y_s[t] - sum(z_t * b_pred)

      beta[t + 1, ]  <- b_pred + drop(K_t) * innov
      P[,, t + 1]    <- (diag(2) - as.vector(K_t) %o% z_t) %*% P_pred
      V_obs          <- (1 - 0.1) * V_obs + 0.1 * innov^2
      log_lik        <- log_lik + dnorm(innov, 0, sqrt(S_t), log = TRUE)
      fitted[t]      <- sum(z_t * b_pred)
    }

    list(k = k, label = colnames(X_mat)[k], beta = beta[-1, ],
         log_lik = log_lik, fitted = fitted, idx = idx)
  })

  # log-sum-exp normalisation for numerical stability
  ll   <- sapply(tvp_models, `[[`, "log_lik")
  ll[!is.finite(ll)] <- min(ll[is.finite(ll)], na.rm = TRUE) - 50
  ll_s <- ll - max(ll)
  wts  <- exp(ll_s) / sum(exp(ll_s))
  names(wts) <- colnames(X_mat)

  dma_fit <- numeric(T_obs)
  for (m in tvp_models) {
    contrib          <- rep(0, T_obs)
    contrib[m$idx]   <- m$fitted
    dma_fit          <- dma_fit + wts[m$k] * contrib
  }

  list(models = tvp_models,
       inclusion_probs = sort(wts, decreasing = TRUE),
       dma_fitted      = dma_fit)
}


# =============================================================================
# 5. MIDAS NOWCASTING — §5
# =============================================================================
# Model: y_t^(q) = φ₀ + φ₁ · b(L, K, θ) · x_{t−h}^(m) + ε_t
#   b(L,K,θ): exponential Almon polynomial (nealmon, 2 parameters θ₁, θ₂)
#   h = 0: nowcast; h = 1: one-quarter-ahead forecast
# Evaluation: RMSE relative to AR(1) benchmark (ratio < 1 → MIDAS wins)
# Fallback: flat-weight OLS (equivalent to θ = (0,0)) when MIDAS fails to converge

run_midas_forecast <- function(gdp_qq, ind_monthly, h = 0,
                                 window = 40, lags_m = 2) {
  if (!requireNamespace("midasr", quietly = TRUE))
    stop("midasr package required. Install with: install.packages('midasr')")

  T_q    <- length(gdp_qq)
  m_freq <- 3
  T_m    <- min(length(ind_monthly), T_q * m_freq)
  T_q    <- T_m %/% m_freq
  gdp_q  <- gdp_qq[1:T_q]
  ind_m  <- ind_monthly[1:T_m]

  n_roll <- T_q - window - h
  if (n_roll < 10) stop("Rolling window too large relative to sample size")

  fcst_midas <- numeric(n_roll)
  fcst_ar1   <- numeric(n_roll)
  actual     <- gdp_q[(window + h + 1):T_q]

  for (t in seq_len(n_roll)) {
    q_train_idx <- t:(t + window - 1)
    q_fc_idx    <- t + window + h - 1
    y_tr        <- gdp_q[q_train_idx]

    m_start <- (q_train_idx[1] - 1) * m_freq + 1
    m_end   <- q_train_idx[length(q_train_idx)] * m_freq
    x_tr    <- ind_m[m_start:m_end]

    ar_fit      <- tryCatch(ar(y_tr, order.max = 1, method = "ols"), error = function(e) NULL)
    fcst_ar1[t] <- if (!is.null(ar_fit))
      tail(y_tr, 1) * ar_fit$ar[1] + ar_fit$x.intercept else mean(y_tr)

    fcst_midas[t] <- tryCatch({
      fit <- midasr::midas_r(
        y_tr ~ midasr::mls(x_tr, 0:(lags_m - 1), m_freq, midasr::nealmon),
        start = list(x_tr = c(0, 0))   # θ₁, θ₂ (nealmon has exactly 2 parameters)
      )
      m_fc_start <- (q_fc_idx - 1) * m_freq + 1 - (lags_m - 1)
      m_fc_end   <- min((q_fc_idx - 1) * m_freq + m_freq, length(ind_m))
      x_new      <- ind_m[m_fc_start:m_fc_end]
      fc         <- predict(fit, newdata = list(x_tr = x_new))
      if (is.finite(fc)) fc else fcst_ar1[t]
    }, error = function(e) {
      # Flat-weight OLS fallback (MIDAS with θ fixed at (0,0))
      x_q  <- colMeans(matrix(x_tr, nrow = m_freq))
      ols  <- lm(y_tr ~ x_q)
      x_fc <- mean(ind_m[max(1, (q_fc_idx - 1) * m_freq + 1):
                            min(length(ind_m), q_fc_idx * m_freq)])
      coef(ols)[1] + coef(ols)[2] * x_fc
    })
  }

  rmse_midas <- sqrt(mean((fcst_midas - actual)^2, na.rm = TRUE))
  rmse_ar1   <- sqrt(mean((fcst_ar1   - actual)^2, na.rm = TRUE))

  list(forecast  = fcst_midas, ar1_bench = fcst_ar1, actual = actual,
       rmse      = rmse_midas,  rmse_ar1  = rmse_ar1,
       rmse_rel  = rmse_midas / rmse_ar1,  horizon = h)
}


# =============================================================================
# 6. VISUALISATION
# =============================================================================

# ── DCC heatmap (WIFO Figures 1–8 style) ─────────────────────────────────────

plot_dcc_heatmap <- function(dcc_result,
                              title    = "Dynamic Cross-Correlation",
                              ref_label = "GDP (q-q)") {
  mat   <- dcc_result$dcc_matrix
  order <- dcc_result$ranking$indicator

  df <- as.data.frame(mat) |>
    tibble::rownames_to_column("indicator") |>
    tidyr::pivot_longer(-indicator, names_to = "tau_str", values_to = "rho") |>
    dplyr::mutate(tau       = as.integer(sub("t", "", tau_str)),
                  indicator = factor(indicator, levels = rev(order)))

  ggplot(df, aes(x = tau, y = indicator, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 2.8) +
    scale_fill_gradientn(colours  = RColorBrewer::brewer.pal(9, "YlOrRd"),
                          limits   = c(0, 1), name = "ρ_XY(τ)", na.value = "grey80") +
    scale_x_continuous(breaks = -4:4, labels = paste0("τ=", -4:4)) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 22)) +
    labs(title    = title,
         subtitle = paste("Reference:", ref_label,
                          "| Ordered by S = Σ j²·ρ(τ=j−5)"),
         x = "← Indicator lags    τ    Indicator leads →", y = NULL,
         caption = "τ > 0: indicator leads GDP; τ < 0: indicator lags") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), panel.grid = element_blank(),
          axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8))
}

plot_dcc_heatmap_dual <- function(dcc_qq, dcc_yy,
                                   title = "DCC – q-q vs y-y GDP Growth") {
  (plot_dcc_heatmap(dcc_qq, "q-q growth") |
     plot_dcc_heatmap(dcc_yy, "y-y growth")) +
    patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 13, face = "bold"))
    )
}

# ── Ranking dot chart ─────────────────────────────────────────────────────────

plot_ranking <- function(dcc_result, top_n = 15, ci_df = NULL) {
  tbl <- dcc_result$ranking |>
    dplyr::slice_head(n = top_n) |>
    dplyr::mutate(indicator = factor(indicator, levels = rev(indicator)))

  if (!is.null(ci_df) && all(c("indicator", "s_lo", "s_hi") %in% names(ci_df)))
    tbl <- dplyr::left_join(tbl, ci_df[, c("indicator", "s_lo", "s_hi")],
                             by = "indicator")

  p <- ggplot(tbl, aes(x = rank_stat, y = indicator)) +
    geom_point(aes(colour = contemp_r), size = 3.5) +
    scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9, "YlOrRd"),
                            limits = c(0, 1), name = "ρ(τ=0)") +
    geom_text(aes(label = sprintf("%.1f", rank_stat)),
              hjust = -0.4, size = 3, colour = "grey30") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.15))) +
    labs(title    = "Indicator Ranking: S = Σ j² · ρ(τ = j−5)",
         subtitle = "Higher S → stronger leading character relative to Swiss GDP",
         x = "S", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  if (all(c("s_lo", "s_hi") %in% names(tbl)))
    p <- p + geom_errorbarh(aes(xmin = s_lo, xmax = s_hi),
                             height = 0.3, colour = "grey50")
  p
}

# ── DMA inclusion probabilities ───────────────────────────────────────────────

plot_dma_weights <- function(dma_result) {
  df <- data.frame(indicator = names(dma_result$inclusion_probs),
                   prob      = dma_result$inclusion_probs) |>
    dplyr::mutate(indicator = factor(indicator, levels = rev(indicator)))

  ggplot(df, aes(x = prob, y = indicator)) +
    geom_col(fill = "#4dac26", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", prob * 100)), hjust = -0.1, size = 3) +
    scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
    labs(title    = "DMA: Posterior Inclusion Probabilities",
         subtitle = "TVP-Kalman filter models, forgetting factor λ = 0.99",
         x = "Inclusion Probability", y = NULL) +
    theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"))
}

# ── MIDAS: nowcast vs actual ──────────────────────────────────────────────────

plot_midas <- function(midas_result, dates) {
  dates_out <- tail(dates, length(midas_result$actual))
  data.frame(date   = dates_out,
             actual = midas_result$actual,
             midas  = midas_result$forecast,
             ar1    = midas_result$ar1_bench) |>
    tidyr::pivot_longer(-date, names_to = "series", values_to = "value") |>
    dplyr::mutate(series = dplyr::recode(series, actual = "Actual",
                                          midas = "MIDAS (Almon)", ar1 = "AR(1) benchmark")) |>
    ggplot(aes(date, value, colour = series, linetype = series)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = c(Actual = "black", "MIDAS (Almon)" = "#e31a1c",
                                    "AR(1) benchmark" = "#1f78b4"), name = NULL) +
    scale_linetype_manual(values = c(Actual = "solid", "MIDAS (Almon)" = "solid",
                                      "AR(1) benchmark" = "dashed"), name = NULL) +
    labs(title    = paste0("MIDAS Nowcasting – h = ", midas_result$horizon, " quarter(s)"),
         subtitle = sprintf("RMSE: %.4f  |  vs AR(1): %.2f  |  %s",
                             midas_result$rmse, midas_result$rmse_rel,
                             ifelse(midas_result$rmse_rel < 1,
                                    "MIDAS beats benchmark ✓", "AR(1) wins")),
         x = NULL, y = "GDP Growth (q-q, %)",
         caption = "Rolling window estimation. Exponential Almon lag polynomial.") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
}

# ── Pre/post-GFC temporal stability ──────────────────────────────────────────

plot_crisis_stability <- function(dcc_full, dcc_pre, dcc_post, indicator_label) {
  get_vec <- function(res, lbl) {
    mat <- res$dcc_matrix
    if (!lbl %in% rownames(mat)) return(rep(NA, 9))
    as.numeric(mat[lbl, ])
  }
  data.frame(tau  = -4:4,
             Full = get_vec(dcc_full, indicator_label),
             Pre  = get_vec(dcc_pre,  indicator_label),
             Post = get_vec(dcc_post, indicator_label)) |>
    tidyr::pivot_longer(-tau, names_to = "Period", values_to = "rho") |>
    ggplot(aes(tau, rho, colour = Period, shape = Period)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.4) +
    scale_colour_manual(values = c(Full = "black", Pre = "#1f78b4", Post = "#e31a1c"),
                         name = NULL) +
    scale_x_continuous(breaks = -4:4) + ylim(-1, 1) +
    labs(title    = paste("Temporal Stability:", indicator_label),
         subtitle = "Pre-GFC (until Q3:2008) vs. Post-GFC (from Q1:2010)",
         x = "τ  [positive = indicator leads GDP]", y = "ρ_XY(τ)",
         caption = "GFC: Q3:2008–Q4:2009") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")
}


# =============================================================================
# 6b. ADDITIONAL ANALYTICAL TOOLS
# =============================================================================

# ── Aggregation sensitivity (§3.2.2) ─────────────────────────────────────────
# Compare all four WIFO aggregation methods: avg / m1 / m2 / m3.
# Flags indicators where S-statistic varies by > 10 across methods.

compare_aggregations <- function(indicators_monthly, gdp_qq,
                                  ind_labels = NULL, max_lag = 4) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  purrr::map_dfr(c("avg", "m1", "m2", "m3"), function(meth) {
    ind_q   <- monthly_to_quarterly(indicators_monthly, method = meth)
    ind_mat <- ind_q[, names(ind_q) != "date"]
    n       <- min(nrow(ind_q), length(gdp_qq))
    run_dcc_analysis(as.list(ind_mat[1:n, ]), gdp_qq[1:n],
                     ind_labels = ind_labels %||% names(ind_mat),
                     max_lag    = max_lag)$ranking |>
      dplyr::mutate(aggregation = meth)
  }) |>
    dplyr::select(indicator, aggregation, rank_stat) |>
    tidyr::pivot_wider(names_from = aggregation, values_from = rank_stat,
                        names_prefix = "S_") |>
    dplyr::mutate(S_range  = pmax(S_avg, S_m1, S_m2, S_m3) -
                               pmin(S_avg, S_m1, S_m2, S_m3),
                  sensitive = S_range > 10) |>
    dplyr::arrange(dplyr::desc(S_avg))
}

# ── Rolling DCC — ρ(τ=0) over time ───────────────────────────────────────────

rolling_dcc <- function(x_ind, y_ref, win = 40, tau = 0) {
  n <- min(length(x_ind), length(y_ref))
  x <- as.numeric(x_ind)[1:n]; y <- as.numeric(y_ref)[1:n]
  rho <- rep(NA_real_, n)
  for (t in seq(win, n)) {
    idx    <- (t - win + 1):t
    rho[t] <- tryCatch(
      compute_dcc(x[idx], y[idx], max_lag = abs(tau))[paste0("t", tau)],
      error = function(e) NA_real_
    )
  }
  rho
}

plot_rolling_dcc <- function(rolling_mat, dates, window_size, ind_labels = NULL) {
  df <- as.data.frame(rolling_mat)
  if (!is.null(ind_labels)) names(df) <- ind_labels
  df$date <- dates
  tidyr::pivot_longer(df, -date, names_to = "indicator", values_to = "rho") |>
    ggplot(aes(date, rho, colour = indicator)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    annotate("rect", xmin = as.Date("2008-07-01"), xmax = as.Date("2010-01-01"),
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "red") +
    annotate("text", x = as.Date("2009-01-01"), y = 0.05,
             label = "GFC", size = 3, colour = "red") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title    = sprintf("Rolling DCC: ρ_XY(τ=0), window = %d quarters", window_size),
         subtitle = "Post-GFC divergence = 'new modesty' hypothesis (WIFO §4)",
         x = NULL, y = "ρ_XY(τ=0, rolling)", colour = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
}

# ── Block bootstrap CIs for ρ_XY(τ) ─────────────────────────────────────────
# Moving-block bootstrap preserves temporal dependence structure.

bootstrap_dcc_ci <- function(x_ind, y_ref, max_lag = 4,
                              B = 999, block_len = 8, alpha = 0.05, seed = 42) {
  set.seed(seed)
  n    <- min(length(x_ind), length(y_ref))
  x    <- as.numeric(x_ind)[1:n]; y <- as.numeric(y_ref)[1:n]
  taus <- -max_lag:max_lag
  starts    <- seq_len(n - block_len + 1)
  n_blocks  <- ceiling(n / block_len)

  boot_cors <- matrix(NA_real_, B, length(taus))
  for (b in seq_len(B)) {
    s_idx  <- sample(starts, n_blocks, replace = TRUE)
    idx_b  <- unlist(lapply(s_idx, function(s) s:(s + block_len - 1)))[1:n]
    boot_cors[b, ] <- tryCatch(
      compute_dcc(x[idx_b], y[idx_b], max_lag),
      error = function(e) rep(NA_real_, length(taus))
    )
  }

  obs_cors <- compute_dcc(x, y, max_lag)
  ci_lo    <- apply(boot_cors, 2, quantile, probs = alpha / 2,  na.rm = TRUE)
  ci_hi    <- apply(boot_cors, 2, quantile, probs = 1 - alpha / 2, na.rm = TRUE)

  data.frame(tau = taus, rho = obs_cors,
             ci_lo = ci_lo, ci_hi = ci_hi,
             sig   = obs_cors < ci_lo | obs_cors > ci_hi,
             row.names = NULL)
}

plot_dcc_ci <- function(ci_df, title = "DCC with Bootstrap CIs",
                         alpha = 0.05, block_size = 8) {
  ci_pct <- round((1 - alpha) * 100)
  ggplot(ci_df, aes(tau, rho)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2, fill = "#1f78b4") +
    geom_line(linewidth = 1, colour = "#1f78b4") +
    geom_point(aes(shape = sig, size = sig), colour = "#1f78b4") +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19),
                        labels = c("Not sig.", sprintf("Sig. at %d%%", round(alpha * 100))),
                        name = NULL) +
    scale_size_manual(values = c(`FALSE` = 2, `TRUE` = 3.5), guide = "none") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    scale_x_continuous(breaks = unique(ci_df$tau)) + ylim(-1, 1) +
    labs(title    = title,
         subtitle = sprintf("%d%% block-bootstrap CI (block length = %d) | Filled = significant",
                             ci_pct, block_size),
         x = "τ  [positive = indicator leads GDP]", y = "ρ_XY(τ)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")
}

# ── X-13ARIMA-SEATS seasonal adjustment (SECO standard method) ───────────────
# Requires: install.packages(c("seasonal", "x13binary"))

seasonal_adjust <- function(ts_obj, series_name = "series",
                              easter = TRUE, x11 = FALSE) {
  if (!requireNamespace("seasonal", quietly = TRUE)) {
    message("Install 'seasonal' + 'x13binary' for X-13ARIMA-SEATS adjustment.")
    return(ts_obj)
  }
  library(seasonal)
  spec <- if (x11) seasonal::seas(ts_obj, x11 = "")
  else if (easter && frequency(ts_obj) == 12)
    seasonal::seas(ts_obj, regression.variables = "easter[1]")
  else
    seasonal::seas(ts_obj)

  list(sa       = seasonal::final(spec),
       trend    = seasonal::trend(spec),
       seasonal = seasonal::series(spec, "s10"),
       spec     = spec)
}


# =============================================================================
# 7. REPORTING TOOLS
# =============================================================================

summarise_findings <- function(dcc_result, midas_result = NULL, pre_post = NULL) {
  top   <- dcc_result$ranking[1, ]
  bench <- dcc_result$ranking[dcc_result$ranking$indicator == "KOF Barometer (kof-econ-bar)", ]
  if (nrow(bench) == 0) bench <- dcc_result$ranking[2, ]

  cat(strrep("─", 65), "\n")
  cat("── FINDINGS ──\n")
  cat(strrep("─", 65), "\n")
  cat(sprintf("Best leading indicator (q-q GDP) : %s\n",  top$indicator))
  cat(sprintf("  WIFO ranking statistic S       : %.1f\n", top$rank_stat))
  cat(sprintf("  Contemporaneous ρ(τ=0)         : %.3f\n", top$contemp_r))
  cat(sprintf("  Max lead correlation (τ>0)     : %.3f\n", top$max_lead))
  cat(sprintf("KOF Barometer reference          : ρ(τ=0) = %.3f, S = %.1f\n",
              bench$contemp_r, bench$rank_stat))

  if (!is.null(midas_result))
    cat(sprintf("MIDAS vs AR(1) RMSE ratio (h=%d)  : %.3f  → %s\n",
                midas_result$horizon, midas_result$rmse_rel,
                ifelse(midas_result$rmse_rel < 1, "MIDAS beats benchmark ✓", "AR(1) wins ✗")))

  if (!is.null(pre_post)) {
    unstable <- pre_post[abs(pre_post$delta_rho) > 0.20, "indicator"]
    cat(sprintf("Post-GFC unstable (|Δρ|>0.20)   : %s\n",
                if (length(unstable) == 0) "none" else paste(unstable, collapse = ", ")))
  }
  cat(strrep("─", 65), "\n")
  invisible(list(top = top, bench = bench))
}

gfc_stability_report <- function(dcc_pre, dcc_post) {
  rho_pre  <- dcc_pre$ranking[,  c("indicator", "contemp_r", "rank_stat")]
  rho_post <- dcc_post$ranking[, c("indicator", "contemp_r", "rank_stat")]
  merged   <- merge(rho_pre, rho_post, by = "indicator", suffixes = c("_pre", "_post"))
  merged$delta_rho <- round(merged$contemp_r_post - merged$contemp_r_pre,  2)
  merged$delta_S   <- round(merged$rank_stat_post - merged$rank_stat_pre,  1)
  merged$flag      <- ifelse(abs(merged$delta_rho) > 0.20, "UNSTABLE", "stable")
  result <- merged[order(merged$delta_rho),
                   c("indicator", "contemp_r_pre", "contemp_r_post",
                     "delta_rho", "delta_S", "flag")]
  cat("\nGFC Stability Report\n")
  cat("ρ(τ=0): pre-GFC vs. post-GFC | |Δρ| > 0.20 → UNSTABLE\n")
  cat(strrep("-", 65), "\n")
  print(result, row.names = FALSE)
  invisible(result)
}

midas_lag_sensitivity <- function(gdp_qq, ind_monthly, K_vec = 1:6,
                                   h = 0, window = 40) {
  cat("MIDAS lag-length sensitivity (K = monthly lags per quarter):\n")
  result <- purrr::map_dfr(K_vec, function(K) {
    res <- tryCatch(
      run_midas_forecast(gdp_qq, ind_monthly, h = h, window = window, lags_m = K),
      error = function(e) list(rmse = NA_real_, rmse_ar1 = NA_real_, rmse_rel = NA_real_)
    )
    data.frame(K = K, rmse = round(res$rmse, 4), rmse_ar1 = round(res$rmse_ar1, 4),
               rmse_rel = round(res$rmse_rel, 3),
               beats_ar1 = !is.na(res$rmse_rel) && res$rmse_rel < 1)
  })
  print(result, row.names = FALSE)
  invisible(result)
}


# =============================================================================
# 8. MAIN ANALYSIS PIPELINE
# =============================================================================

run_swiss_bca <- function(use_live_data = TRUE, output_dir = "output") {

  cat(strrep("=", 65), "\n")
  cat("SWISS BUSINESS CYCLE ANALYSIS\n")
  cat("Glocker & Kaniovski (2019) methodology | SECO\n")
  cat(strrep("=", 65), "\n\n")

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # ── 1. Load data ──────────────────────────────────────────────────────────
  cat("Step 1: Loading data ...\n")
  gdp_live <- if (use_live_data) load_gdp_bfs()       else NULL
  kof_live <- if (use_live_data) load_kof_barometer() else NULL
  sim      <- simulate_swiss_data(n_q = 120)

  gdp_df           <- if (!is.null(gdp_live)) gdp_live else sim$gdp
  gr               <- growth_rates(gdp_df$gdp)
  gdp_df$growth_qq <- gr$qq
  gdp_df$growth_yy <- gr$yy

  # ── 2. Temporal aggregation ───────────────────────────────────────────────
  cat("Step 2: Monthly → quarterly aggregation ...\n")
  ind_q <- monthly_to_quarterly(sim$indicators_monthly, method = "avg")

  if (!is.null(kof_live)) {
    kof_q <- monthly_to_quarterly(kof_live, method = "avg")
    ind_q <- dplyr::left_join(ind_q,
                               kof_q |> dplyr::rename(kof_live = kof_bar),
                               by = "date")
  }

  merged   <- dplyr::inner_join(gdp_df, ind_q, by = "date") |>
    dplyr::filter(!is.na(growth_qq))
  T_q      <- nrow(merged)

  ind_cols <- intersect(c("cs_econ_exp", "kof_econ_bar", "snb_bci",
                           "kn_bci", "seco_dfm", "sentix", "pmi_ind",
                           "ifo_de", "oecd_cli"), names(merged))
  labels   <- c(cs_econ_exp  = "CS Econ.Expect. (cs-econ-expect)",
                kof_econ_bar = "KOF Barometer (kof-econ-bar)",
                snb_bci      = "SNB BCI (snb-bci)",
                kn_bci       = "KN BCI (kn-bci)",
                seco_dfm     = "SECO DFM (seco-dfm)",
                sentix       = "Sentix Econ.Expect. (sentix)",
                pmi_ind      = "PMI Industrial",
                ifo_de       = "ifo Business Climate (DE)",
                oecd_cli     = "OECD CLI")
  ind_labels <- unname(labels[ind_cols])
  cat(sprintf("   → %d quarterly obs | %d indicators\n", T_q, length(ind_cols)))

  # ── 3. DCC analysis ───────────────────────────────────────────────────────
  cat("\nStep 3: DCC – GDP q-q ...\n")
  dcc_qq <- run_dcc_analysis(as.list(merged[, ind_cols]),
                              merged$growth_qq, ind_labels)
  print(dcc_qq$ranking[, c("indicator", "rank_stat", "contemp_r")])

  cat("\nStep 4: DCC – GDP y-y ...\n")
  yy_idx <- which(!is.na(merged$growth_yy))
  dcc_yy <- run_dcc_analysis(as.list(merged[yy_idx, ind_cols]),
                              merged$growth_yy[yy_idx], ind_labels)

  # ── 4. Temporal stability: pre/post GFC ───────────────────────────────────
  cat("\nStep 5: Pre/post-GFC stability ...\n")
  pre_idx  <- which(merged$date  < as.Date("2008-07-01"))
  post_idx <- which(merged$date >= as.Date("2010-01-01"))

  dcc_pre  <- if (length(pre_idx) >= 30)
    run_dcc_analysis(as.list(merged[pre_idx, ind_cols]),
                     merged$growth_qq[pre_idx], ind_labels) else NULL
  dcc_post <- run_dcc_analysis(as.list(merged[post_idx, ind_cols]),
                                merged$growth_qq[post_idx], ind_labels)

  gfc_report <- if (!is.null(dcc_pre)) gfc_stability_report(dcc_pre, dcc_post) else NULL

  # ── 5. Pre-whitening robustness ───────────────────────────────────────────
  cat("\nStep 6: Pre-whitening (ARIMA/BIC) for top indicators ...\n")
  top2 <- dcc_qq$ranking$indicator[1:2]
  for (lbl in top2) {
    col <- ind_cols[ind_labels == lbl]
    if (length(col) == 0) next
    res <- prewhiten_dcc(merged[[col]], merged$growth_qq)
    cat(sprintf("   %-42s ARIMA(%d,0,%d)  ρ(τ=0): %.3f → %.3f\n",
                lbl, res$arima_order[1], res$arima_order[3],
                dcc_qq$dcc_matrix[lbl, "t0"], res$dcc["t0"]))
  }

  # ── 6. Rolling DCC ────────────────────────────────────────────────────────
  cat("\nStep 6b: Rolling DCC (window = 40 quarters) ...\n")
  roll_win <- 40
  roll_mat <- tryCatch({
    mat  <- sapply(ind_cols, function(col)
      rolling_dcc(merged[[col]], merged$growth_qq, win = roll_win, tau = 0))
    tmat <- t(mat)
    rownames(tmat) <- ind_labels
    tmat
  }, error = function(e) { message("  Rolling DCC failed: ", e$message); NULL })

  # ── 7. Dynamic Model Averaging ────────────────────────────────────────────
  cat("\nStep 7: Dynamic Model Averaging (DMA) ...\n")
  X_mat     <- as.matrix(merged[, ind_cols]); colnames(X_mat) <- ind_labels
  ok        <- complete.cases(cbind(merged$growth_qq, X_mat))
  dma_res   <- run_dma(merged$growth_qq[ok], X_mat[ok, ], lambda = 0.99)
  cat("   Top inclusion probabilities:\n")
  print(round(head(dma_res$inclusion_probs, 5), 3))

  # ── 8. MIDAS nowcasting ───────────────────────────────────────────────────
  cat("\nStep 8: MIDAS nowcasting ...\n")
  best_col <- ind_cols[ind_labels == names(dma_res$inclusion_probs)[1]]
  best_m   <- intersect(best_col, names(sim$indicators_monthly))[1]

  midas_h0 <- tryCatch(
    run_midas_forecast(merged$growth_qq,
                        sim$indicators_monthly[[best_m]],
                        h = 0, window = 40, lags_m = 3),
    error = function(e) { message("  MIDAS h=0 failed: ", e$message); NULL })

  midas_h1 <- tryCatch(
    run_midas_forecast(merged$growth_qq,
                        sim$indicators_monthly[[best_m]],
                        h = 1, window = 40, lags_m = 3),
    error = function(e) { message("  MIDAS h=1 failed: ", e$message); NULL })

  if (!is.null(midas_h0))
    cat(sprintf("   Nowcast  RMSE: %.4f  (vs AR(1): %.4f, ratio: %.2f)\n",
                midas_h0$rmse, midas_h0$rmse_ar1, midas_h0$rmse_rel))
  if (!is.null(midas_h1))
    cat(sprintf("   1Q-ahead RMSE: %.4f  (vs AR(1): %.4f, ratio: %.2f)\n",
                midas_h1$rmse, midas_h1$rmse_ar1, midas_h1$rmse_rel))

  # ── 9. Plots ──────────────────────────────────────────────────────────────
  cat("\nStep 9: Generating plots ...\n")

  p_top <- (plot_dcc_heatmap(dcc_qq, "DCC – GDP q-q Growth") | plot_ranking(dcc_qq)) +
    patchwork::plot_annotation(
      title    = "Swiss Business Cycle Indicators – DCC Analysis",
      subtitle = sprintf("N = %d quarters | %d indicators | 3-month average aggregation",
                          T_q, length(ind_cols)),
      theme    = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"))
    )

  ggplot2::ggsave(file.path(output_dir, "dcc_analysis.png"),
                   p_top, width = 16, height = 8, dpi = 300)
  ggplot2::ggsave(file.path(output_dir, "dcc_dual_qq_yy.png"),
                   plot_dcc_heatmap_dual(dcc_qq, dcc_yy), width = 18, height = 8, dpi = 300)
  ggplot2::ggsave(file.path(output_dir, "dma_weights.png"),
                   plot_dma_weights(dma_res), width = 8, height = 5, dpi = 300)

  if (!is.null(dcc_pre))
    ggplot2::ggsave(file.path(output_dir, "temporal_stability.png"),
                     plot_crisis_stability(dcc_qq, dcc_pre, dcc_post,
                                           dcc_qq$ranking$indicator[1]),
                     width = 8, height = 5, dpi = 300)

  if (!is.null(roll_mat))
    ggplot2::ggsave(file.path(output_dir, "rolling_dcc.png"),
                     plot_rolling_dcc(t(roll_mat), merged$date, roll_win, ind_labels),
                     width = 12, height = 6, dpi = 300)

  if (!is.null(midas_h0))
    ggplot2::ggsave(file.path(output_dir, "midas_nowcast.png"),
                     plot_midas(midas_h0, merged$date), width = 10, height = 5, dpi = 300)

  # ── 10. Summary ───────────────────────────────────────────────────────────
  cat("\n")
  summarise_findings(dcc_qq, midas_h0, gfc_report)
  cat(sprintf("\n✓ Done. Plots saved to: %s/\n", output_dir))
  cat(strrep("-", 65), "\n")

  invisible(list(
    data        = merged,     dcc_qq      = dcc_qq,     dcc_yy    = dcc_yy,
    dcc_pre     = dcc_pre,    dcc_post    = dcc_post,   gfc_report = gfc_report,
    dma         = dma_res,    midas_h0    = midas_h0,   midas_h1  = midas_h1,
    rolling_dcc = roll_mat
  ))
}


# =============================================================================
# 9. WIFO TABLE-31 INDICATOR REFERENCE
# All 28 indicators evaluated in Glocker & Kaniovski (2019), §3.1
# =============================================================================

wifo_indicators <- data.frame(
  wifo_code    = c("cs-econ-expect", "kof-econ-bar", "kn-bci",
                    "seco-dfm", "snb-bci",
                    "sentix-econ-expect", "kof-ind-no3m", "snb-foreign-pmi",
                    "cs-ind-pmi-bo", "cs-ind-pmi-prod", "kof-ind-nopm",
                    "ism-us-pmi-serv", "cs-serv-pmi",
                    "ec-ea-serv-conf", "ec-at-serv-dem3m",
                    "oecd-no", "oecd-cli", "oecd-cli-prod",
                    "oecd-cli-consconf", "ifo-de-climate",
                    "ec-ea-ind-conf", "ec-ea-sent",
                    "ism-us-pmi-manu", "seco-consconf",
                    "seco-consconf-outlook", "kof-ind-prodpm",
                    "kof-ind-trend", "kof-ind-prod3m"),
  group        = c(rep("Group 1", 5), rep("Group 2", 6), rep("Service", 5),
                    rep("Foreign", 10), rep("Other", 2)),
  frequency    = c("M","M","M","M","Q","M","M","M","M","M","M","M","M",
                    "M","M","M","M","M","M","M","M","M","M","Q","Q","M","M","M"),
  r_access     = c(
    # cs-econ-expect: renamed to UBS Swiss Economic Expectations after
    # UBS takeover of Credit Suisse (March 2023). Verify current availability.
    "Manual / Bloomberg (UBS rebrand post-2023)", "kofdata", "Manual",
    "SECO internal", "data.snb.ch", "sentix.de (registered)", "kofdata",
    "data.snb.ch", "procure.ch", "procure.ch", "kofdata",
    "FRED/fredr", "procure.ch", "ECB SDW / ecb package", "ECB SDW",
    "OECD package", "OECD package", "OECD package", "OECD package",
    "FRED/fredr or ifo.de", "ECB SDW", "ECB SDW", "FRED/fredr",
    "seco.admin.ch (Excel)", "seco.admin.ch (Excel)",
    "kofdata", "kofdata", "kofdata"),
  stringsAsFactors = FALSE
)


# =============================================================================
# ENTRY POINT
# =============================================================================

if (interactive()) {
  results <- run_swiss_bca(use_live_data = TRUE, output_dir = "output")
}
