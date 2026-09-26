# Changelog

## v0.4.0 — 2026-09-25
- Froze the successful Phase-A historical replication and moved calendar-time interpretation into a separate Phase B1 path.
- Added exact-calendar autocorrelation: observations are paired only when their timestamps differ by the requested number of days.
- Added exact-calendar P10/P11 cross-correlation with an explicit dissertation-oriented lag convention (positive L pairs P10(t+L) with P11(t)).
- Calendar correlations report pair count, non-demeaned normalized dot product, Pearson r on raw flux, and Pearson r on log10 positive flux.
- Added single-spacecraft sample-lag→calendar-separation diagnostics, because the ACF lag axis is affected by missing dates too.
- Added Lomb–Scargle summaries using actual observation times and spectra of the observation-availability window.
- Added `run_phase_b_windows.bat`, `run_phase_b_unix.sh`, and `scripts/calendar_reanalysis.py`.
- Added explicit annual (330–400 d) versus ~1.3-year (430–520 d) comparisons and exact 340/341/365/469/470-day readouts.
- Test suite: **29 passing tests**.


## v0.3.3 — 2026-09-25
- Real SPDF `p10cpi.h24` inspection established the actual legacy daily layout: 15 numeric columns = `YY, DOY, HOUR` + the 12 CPI channels.
- The same real file contains exactly **7042 numeric records**, matching the dissertation's printed Pioneer 10 vector length before any fitting or filtering.
- Fixed the parser so the 15-column archive is no longer misclassified as decimal-year data.
- Added an explicit particle-only canonical product (`p10_cpi_daily.csv` / `p11_cpi_daily.csv`).
- Decoupled the locked ~470 d ACF and ~340 d unshifted CCF experiment from the later XYZ/solar-wind propagation correction.
- Full 17-column export is now emitted only when all particle, XYZ and solar-wind fields are genuinely available; missing geometry is not fabricated.
- Added a regression test using the exact observed SPDF record shape and values.
- Test suite: **22 passing tests**.

## v0.3.1
- Fixed NASA/PDS HAPI requests for the server's advertised HAPI 3.1 interface: native `dataset/start/stop` parameters are now used first.
- Added automatic fallback to legacy `id/time.min/time.max` syntax for compatibility.
- Corrected misleading diagnostics: HTTP 500 from the HAPI server is no longer reported as a DNS/Internet failure.
- Added response URL/body diagnostics for server-side HAPI failures.
- Added transient retry handling for HTTP 502/503/504.
- Added regression tests for HAPI 3 request syntax and legacy fallback.

## v0.3 — 2026-09-25
- Added PDS/PPI HAPI catalog configuration for P10/P11 CPI daily/hourly, PA hourly and PA trajectory products.
- Added metadata-driven HAPI client with year-chunked downloads.
- Added conservative channel, coordinate and plasma-speed field resolver; ambiguous mappings fail loudly.
- Added automatic daily CPI + PA speed reconstruction into the dissertation's 17-column schema.
- Added exact historical-style `.dat` export plus analysis-ready CSV.
- Added assembly reports with DOI/source IDs, mappings, missingness, row checks and SHA-256 provenance.
- Added `inspect_hapi.py`, `fetch_nasa.py`, and one-command `run_phase_a.py`.
- Added Windows and Unix launchers.
- Expanded test suite to 17 passing tests.
- Clarified that daily PA speed is a reconstruction assumption to be validated, not an established fact about the lost historical files.

## v0.2 — 2026-09-25
- Added P10/P11 common-date alignment and cross-correlation.
- Reconstructed propagation delay and backward P11 signal shift from dissertation code.
- Added locked target-band readouts, modern Lomb–Scargle and exploratory AR(1) testing.
- Added headless replication runner and provenance hashes.
## v0.3.2 — 2026-09-25

- Diagnosed the PDS/PPI HAPI daily-CPI failure from a real Windows run: the
  server returns HTTP 500 with `unable to parse time: 72Z`, i.e. the backend is
  choking on the legacy two-digit year representation rather than on local
  Internet/DNS access.
- Added a preferred NASA/SPDF legacy 24-hour CPI acquisition path.  This is also
  methodologically preferable for the 1998 replication because it is closer to
  the historical ASCII source than a modern HAPI transformation.
- Added a conservative legacy-ASCII parser supporting explicit YEAR/DOY,
  SCID/YEAR/DOY, and dissertation-style decimal-year layouts. Unsupported
  record widths fail loudly; no arbitrary column-window guessing is used.
- `scripts/fetch_nasa.py --cpi-source auto` now tries SPDF CPI first and HAPI
  only as a fallback. `--cpi-source spdf` and `--cpi-source hapi` can force a
  transport for diagnostics.
- HAPI `/info` metadata remains captured when available even when data bytes
  come directly from SPDF.
- If the daily CPI archive contains solar-wind speed, it is retained directly.
  Otherwise the pipeline still requires the official PA speed source and will
  not silently fabricate the field.
- Added parser regression tests for 17-column YEAR/DOY, 18-column
  SCID/YEAR/DOY, and decimal-year records with solar-wind speed.
- Test suite: 21 passing tests.

## v0.3.4
- Added exact dissertation benchmark diagnostics for the 11–20 MeV channel.
- Phase-A output now compares P10 ACF at 469 d and P11 ACF at 1697/467 d directly with the 1998 table.
- Cross-correlation output now reports both +315..365 and -365..-315 bands, exact thesis lags (1905, 341, 1524, -214), their opposite signs, and the strongest local peaks.
- Added common-date alignment/gap diagnostics to investigate the CCF discrepancy without changing the locked replication algorithm.
- Windows launcher now reuses existing processed P10/P11 CPI files automatically, avoiding unnecessary NASA re-downloads during diagnostics.

## v0.3.5
- Confirmed exact historical ACF/CCF benchmark reproduction for 11–20 MeV protons (to thesis rounding precision).
- Made cross-correlation sign convention explicit: thesis lag = -SciPy lag.
- Added calendar-separation diagnostics for gapped common-date series; historical xcorr lags remain sample-index lags.
- Added machine-readable benchmark replication checks and an all-four CCF replication verdict.
- Added tests for sign translation, calendar-lag diagnostics, and benchmark tolerance.

## v0.5.0 — Phase B2 multichannel uneven-sampling screen
- Added interpolation-free calendar DCF on actual daily lags, with log10 flux and configurable lag bins.
- Added generalized/floating-mean Lomb–Scargle screen on log10 flux at actual observation times.
- Added all-channel Phase B2 runner and compact CSV/JSON summaries for the fixed 330–400 d and 430–520 d bands.
- Sampling-window peaks are reported beside signal peaks; B2 explicitly makes no significance claim.
- Added regression tests for DCF lag recovery, zero-lag normalization, and irregular-time log-Lomb–Scargle recovery.

## v0.6.0 — Phase B3 surrogate significance
- Added fast FFT-based exact-calendar DCF engine for surrogate testing, regression-checked against the slower pairwise DCF implementation.
- Added daily AR(1) estimation using only genuinely adjacent calendar-day observations; gaps are never treated as one-day transitions.
- Added ACF band-max red-noise tests that simulate on the complete daily calendar and reapply the exact historical observation mask.
- Added P10/P11 circular-shift CCF null that preserves P11 values, autocorrelation structure, and missingness pattern while breaking absolute spacecraft alignment.
- Added raw-log10 and linearly detrended-log10 robustness variants.
- Added Holm family-wise and Benjamini-Hochberg FDR corrections across the 12 particle channels.
- Added `scripts/phase_b3_significance.py`, Windows/Unix launchers, JSON/CSV outputs, and explicit first-pass versus publication-resolution simulation guidance.
- Test suite: **36 passing tests**.

## v0.7.0 — publication-resolution B3 confirmation
- Added `run_phase_b3_publication_windows.bat` / Unix launcher for a full 12-channel rerun with **1999 surrogates per test**.
- Added `scripts/review_b3.py` to summarize cross-spacecraft candidates without creating new statistical tests.
- The review explicitly reports the Monte-Carlo p-value floor and the consequence for Holm correction. With 199 simulations, `p_min=0.005`; across 12 channels, even a zero-exceedance result has `Holm >= 0.06`, so the first-pass run cannot possibly cross a Holm 0.05 threshold.
- Candidate review distinguishes raw-log10 from linearly detrended robustness and flags P11 channels affected by the post-1980 CPI reliability mask.
- Publication-resolution outputs are written separately under `results/phase_b3_publication/` so the 199-simulation first-pass result remains frozen.

## v0.8.0
- Added Phase B4 full-lag morphology and sliding-window time localization.
- B4 reads publication-resolution B3 candidates and scans 200–700 calendar-day lags.
- Adds explicit band-edge and local-peak/prominence diagnostics so a band maximum is not automatically called a periodicity.
- Adds 5-year sliding windows (1-year step) for descriptive time localization.
- Writes full DCF profiles per candidate/variant plus `localization_summary.json` and `sliding_window_summary.csv`.

## v0.9.0
- Added Phase B5 to separate cross-spacecraft lag correlation from actual periodicity.
- Added individual P10/P11 calendar-ACF morphology and irregular-time Lomb–Scargle morphology.
- Added 5-year sliding Lomb–Scargle localization and descriptive P10/P11 period-concordance diagnostics.
- Added B5 Windows/Unix launchers and tests.

## v1.0.0
- Added Phase B6 all-channel individual-spacecraft periodicity atlas.
- Removed the B4 cross-spacecraft-survivor selection as a prerequisite for periodicity screening; all 12 CPI particle channels are now audited independently.
- Added channel-specific sampling-window spectra using the exact positive-flux mask used by log10 Lomb–Scargle.
- Added signal-vs-window peak separation in Rayleigh-frequency units; <=1 Rayleigh is explicitly flagged as unresolved from the observation-window spectrum, not automatically called an alias.
- Added exact-common-span P10/P11 Lomb–Scargle concordance for each channel and band.
- Attached frozen publication-resolution B3 calendar-ACF p/Holm/FDR results where available; B6 itself creates no new p-values.
- Added compact `candidate_flags.json` with explicit evidence flags and no opaque composite score.
- Added Windows/Unix launchers and regression tests.

## v1.1.0 — Phase B7 sampling-window co-tracking audit
- Added `scripts/phase_b7_sampling_tracking.py` and Windows/Unix launchers.
- B7 uses exact common valid calendar windows and compares each sliding-window Lomb–Scargle signal peak with the sampling-window peak in the same predeclared band.
- Reports Rayleigh-frequency separation, alias-warning fractions, signal/window frequency co-tracking, and cross-spacecraft concordant local-peak windows.
- B7 is explicitly descriptive: overlapping windows are not independent and no new significance p-values are created.

## v1.2.0 — Phase B8 Rayleigh-resolution stress test
- Added a multi-scale sampling-window audit at 5, 8, 10, 12 and 15-year common-calendar windows.
- B8 explicitly quantifies the first-order period width corresponding to one Rayleigh frequency resolution (`P^2/T`), addressing the fact that five-year windows are too coarse for a <=1-Rayleigh warning to be strongly discriminating in the 330–520 d range.
- Reports signal/window separation and alias-warning fractions as window length increases; longer windows narrow the spectral resolution at the cost of time localization.
- Channels truncated by the Pioneer 11 post-1980 reliability mask are retained but will naturally have no long-window results beyond their valid common span.
- B8 is descriptive and creates no new p-values.

## v1.3.0 — Phase B9 candidate-focused spectral significance
- Added a mask-preserving daily-AR(1) Monte-Carlo significance test for maximum classical Lomb–Scargle power inside a locked period band.
- B9 freezes two long-baseline ~1.3-year candidates after B8: He 20–24 MeV/n (primary) and He 11–20 MeV/n (secondary comparator).
- Daily AR(1) surrogates are generated on the complete calendar and then sampled through the exact positive-flux observation mask; no missing values are interpolated.
- The test statistic is the **maximum** spectral power anywhere inside 430–520 d, avoiding a single-frequency post-hoc test.
- Raw-log10 and linearly detrended-log10 variants are both tested for P10 and P11.
- Holm and Benjamini–Hochberg adjustments are applied across the locked four spacecraft/channel tests separately for each robustness variant.
- B9 is explicitly labelled post-selection evidence, not independent replication; surviving results still require time-frequency localization.
- Added Windows/Unix launchers and regression tests. Test suite: **50 passing tests**.

## v1.4.0 — Phase B10 intermittent joint spectral significance
- B9 found no full-record global 430–520 d Lomb–Scargle significance for either locked long-baseline helium candidate; B10 therefore tests the distinct remaining hypothesis of **intermittent** shared spectral structure rather than persistent periodicity.
- Added fixed 8-year calendar windows with a 2-year step, chosen before B10 results are examined.
- Added mask-preserving daily AR(1) Monte-Carlo calibration for the maximum over the entire predeclared time-frequency search, so the window/frequency look-elsewhere effect is included in the null statistic.
- Added a joint P10/P11 statistic: geometric-mean Lomb power at the **same period in the same window**, maximized across all locked windows and frequencies.
- P10 and P11 null processes are simulated independently on the exact common daily calendar and sampled through each spacecraft's exact positive-flux mask; no gaps are interpolated.
- The primary inferential family is the two locked joint-candidate tests, with Holm and BH/FDR adjustment separately for raw-log10 and linearly detrended-log10 variants.
- Added Windows/Unix launchers and regression tests.
