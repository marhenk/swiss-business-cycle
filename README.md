# Swiss Business Cycle Analysis

Replication of the indicator evaluation methodology from **Glocker & Kaniovski (2019)**, 
"An Evaluation of Business Cycle Indicators for the Swiss Economy," 
*WIFO Grundlagen für die Wirtschaftspolitik* Nr. 6.

Prepared as a code sample for an application to the Ressort Konjunktur, SECO (March 2026).

---

## What this does

The script evaluates how well Swiss and international business cycle indicators 
co-move with and lead Swiss GDP growth. It replicates all five methods from the paper:

| Section | Method |
|---------|--------|
| §3 | Dynamic Cross-Correlation (DCC) + WIFO ranking statistic *S* |
| §3.2.3 | Pre-whitening via ARIMA (BIC-selected order) |
| §3.2.4 | Parametric DCC via bivariate VAR + Lyapunov equation |
| §4 | Dynamic Model Averaging (TVP-Kalman filter, λ = 0.99) |
| §5 | MIDAS nowcasting (exponential Almon polynomial weights) |

Additional tools: rolling DCC, block-bootstrap confidence intervals, 
aggregation sensitivity checks (§3.2.2), GFC stability report, X-13ARIMA-SEATS hook.

---

## How to run

```r
source("swiss_bca_analysis.R")

# Full pipeline on live BFS/KOF data:
results <- run_swiss_bca(use_live_data = TRUE, output_dir = "output")

# Fully offline (simulated Swiss-calibrated data):
results <- run_swiss_bca(use_live_data = FALSE, output_dir = "output")
```

The script auto-installs missing packages. No manual setup required.

**Live data sources** (with graceful offline fallback):
- GDP: BFS quarterly national accounts (`px-x-0102010000_101`) via the `BFS` package
- KOF Barometer: `kofdata::get_time_series("ch.kof.barometer")`

---

## Key outputs

After running the pipeline, `results` contains:

```
results$dcc_qq$ranking        # indicator ranking by WIFO S-statistic
results$gfc_report            # pre/post-GFC Δρ table (replicates WIFO §4)
results$dma$inclusion_probs   # DMA model weights
results$midas_h0$rmse_rel     # MIDAS vs AR(1) RMSE ratio (nowcast, h=0)
results$midas_h1$rmse_rel     # MIDAS vs AR(1) RMSE ratio (h=1)
```

Plots saved to `output/`:

```
dcc_analysis.png          DCC heatmap + ranking dot chart
dcc_dual_qq_yy.png        q-q vs y-y comparison (WIFO Figures 1–8 format)
dma_weights.png           DMA inclusion probabilities
midas_nowcast.png         MIDAS forecast vs actual vs AR(1) benchmark
rolling_dcc.png           Rolling ρ(τ=0), 40-quarter window
temporal_stability.png    Pre/post-GFC DCC for top indicator
```

---

## Results (simulation)

When run on simulated data calibrated to the Swiss economy:

| Finding | Result |
|---------|--------|
| Best leading indicators (S-statistic) | CS Economic Expectations, OECD CLI |
| Strongest contemporaneous indicator | KOF Barometer (ρ = 0.981) |
| Post-GFC instability | ifo DE: Δρ = −0.66 — replicates WIFO §4 "new modesty" |
| DMA model weight | KOF Barometer dominant (87%) |
| MIDAS h=0 RMSE ratio | 1.04 (expected: AR(1) wins with simulated data) |
| MIDAS h=1 RMSE ratio | 0.98 — MIDAS beats AR(1) at 1-quarter horizon ✓ |

---

## Scope limitations

- **Real-time vintages not implemented.** The analysis uses final BFS revisions 
  throughout. True nowcasting robustness requires historical flash estimates 
  (available on request from BFS/SECO; see Glocker & Kaniovski §3 for vintage effects).
- **DMA model weights** use cumulative log-likelihood (static, α = 1 approximation). 
  For sequential weight updating use `fDMA::fDMA()`. The TVP-Kalman filter (λ = 0.99) 
  is fully dynamic.

---

## Reference

Glocker, C. & Kaniovski, S. (2019). An Evaluation of Business Cycle Indicators 
for the Swiss Economy. *WIFO Grundlagen für die Wirtschaftspolitik* Nr. 6, November 2019.
