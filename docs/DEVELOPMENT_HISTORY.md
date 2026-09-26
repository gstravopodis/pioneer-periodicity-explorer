# Pioneer Periodicity Explorer

Local-first Python/Streamlit reconstruction and modern re-analysis of the 1998 University of Athens Physics dissertation on energetic-particle measurements from Pioneer 10 and 11.

## Status: v0.4.0 — Phase A frozen and reproduced; Phase B1 calendar-aware falsification added

### Implemented
- exact 17-column thesis schema from Table 3;
- thesis decimal-year parser (`72.1721` → day 63 of 1972, using DOY/365 even in leap years);
- documented Pioneer 11 CPI quality mask: columns 6–13 masked from day 239 of 1980 onward;
- reconstructed zero-padded FFT recipe (`7042 -> 8192`, DC removed, `2*abs(FFT)^2` power);
- MATLAB-like normalized `xcorr(...,'coeff')` without demeaning;
- common-date P10/P11 alignment;
- plain and propagation-shifted P10↔P11 cross-correlation;
- reconstructed propagation correction (`dt = round((dr/v)*1731)`) and backward P11 index shift;
- pre-declared target bands around ~470 d and ~340 d;
- gap-preserving Lomb–Scargle and an exploratory AR(1) red-noise test;
- metadata-driven PDS/PPI HAPI client;
- automatic fetch of the official legacy SPDF 24-hour CPI particle products for P10/P11;
- conservative metadata/name mapping of NASA fields to the thesis schema (ambiguities stop the run);
- exact parser for the observed SPDF `YY DOY HOUR + 12 CPI channels` daily layout;
- particle-only Phase-A files are emitted immediately; full 17-column reconstruction is explicitly deferred until separate official XYZ/speed sources are validated;
- machine-readable assembly reports with source IDs, DOI, mappings, missingness and SHA-256 hashes;
- one-command NASA acquisition → assembly → Phase-A replication runner;
- explicit SciPy↔1998 cross-correlation lag-sign translation (thesis lag = -SciPy lag);
- calendar-separation diagnostics showing what an xcorr sample lag means when common-date series contain gaps;
- machine-readable benchmark replication verdicts;
- 29 automated tests.

### Important scientific status
The official legacy Pioneer 10 CPI file has now been reached in a real Windows run and contains **exactly 7042 numeric records**, matching the row count printed in the 1998 dissertation. Its observed record layout is `YY DOY HOUR + 12 CPI particle channels`. This establishes the particle backbone needed for the locked ~470 d ACF and ~340 d unshifted P10/P11 CCF experiment. A real Phase-A run has now reproduced the printed 1998 11–20 MeV proton ACF/CCF benchmarks to rounding precision. The historical computation is therefore considered reproduced. Physical interpretation is now separated into Phase B because the historical lag axis is a sample-index axis, not always literal calendar time.

## Official NASA data path

The project uses the PDS/PPI HAPI service (`https://pds-ppi.igpp.ucla.edu/hapi`) and keeps particle, plasma and trajectory provenance separate.

### Pioneer 10
- CPI daily: `urn:nasa:pds:p10-cpi-jup-cal:data-1day` — DOI `10.17189/yfb9-ve02`
- CPI hourly: `urn:nasa:pds:p10-cpi-jup-cal:data-1hr` — DOI `10.17189/1mq1-dk04`
- PA hourly: `urn:nasa:pds:p10-pa:data-avg-1hr` — DOI `10.17189/75z7-gw37`
- PA trajectory: `urn:nasa:pds:p10-pa:data-traj` — DOI `10.17189/99wc-jw74`

### Pioneer 11
- CPI daily: `urn:nasa:pds:p11-cpi:data-1day` — DOI `10.17189/jqm4-yq88`
- CPI hourly: `urn:nasa:pds:p11-cpi:data-1hr` — DOI `10.17189/k63r-df32`
- PA hourly: `urn:nasa:pds:p11-pa:data-avg-1hr` — DOI `10.17189/rxvv-eb83`
- PA trajectory: `urn:nasa:pds:p11-pa:data-traj` — DOI `10.17189/cjw6-7b69`

The legacy SPDF 24-hour CPI files supply the 12 particle channels. The real Pioneer 10 archive record observed in v0.3.2 does **not** contain XYZ or solar-wind speed; therefore those fields are no longer assumed to be part of that file. They are needed only for the later propagation-corrected extension and must be recovered from separate official trajectory/plasma products.

## Fastest route on Windows

Double-click:

```text
run_phase_a_windows.bat
```

It creates a virtual environment, installs dependencies, downloads the official legacy CPI daily files, writes the particle-only Phase-A CSVs, and runs the locked ACF/CCF analysis for the 11–20 MeV proton channel. A full 17-column reconstruction is not required for these first two historical targets.

## Manual route

```bash
python -m venv .venv
# Windows
.venv\Scripts\activate
# macOS/Linux
# source .venv/bin/activate
pip install -r requirements.txt

# Inspect metadata before downloading, if desired
python scripts/inspect_hapi.py --spacecraft p10 --dataset cpi_daily

# Fetch the exact daily CPI particle backbone for both spacecraft
python scripts/fetch_nasa.py --spacecraft both

# Run the first locked Phase-A target channel
python scripts/run_phase_a.py --channel p_11_20_mev --skip-download

# Open the local UI
streamlit run app.py
```


## Phase B1 — calendar-aware falsification

After a successful Phase A, run:

```text
run_phase_b_windows.bat
```

This does **not** modify the frozen historical replication. It tests whether the historical sample-index features survive when elapsed time is measured in actual calendar days. It writes:

```text
results/phase_b_calendar/calendar_reanalysis_summary.json
results/phase_b_calendar/*_calendar_acf.csv
results/phase_b_calendar/*_calendar_ccf.csv
results/phase_b_calendar/*_lomb_scargle.csv
results/phase_b_calendar/*_sampling_window.csv
```

For each exact calendar lag it reports pair count, a non-demeaned normalized dot product, raw Pearson correlation, and log10-flux Pearson correlation. It also computes Lomb–Scargle periodograms on the actual observation times and spectra of the observation-availability window. The purpose is to distinguish a physical ~1-year/~1.3-year feature from an index-lag or sampling-window effect.

## Outputs

After acquisition:

```text
data/nasa/metadata/    # exact HAPI info JSON used for field mapping
data/nasa/raw/         # unmodified HAPI CSV downloads
data/nasa/processed/   # CPI daily Phase-A CSVs + provenance reports; thesis17 only when complete
results/phase_a/       # frozen FFT/ACF/CCF outputs and summary JSON
```

For Pioneer 10, `p10_cpi_report.json` explicitly checks the dissertation's stated **7042-row** reference. The real archive file has already shown 7042 numeric records; the successful parser/run will now verify the post-parse row count as well.

## Test

```bash
pytest -q
```

Expected for v0.5.0: `32 passed`.

## Research principle

First reproduce or falsify the 1998 outputs. Only then run modern inference. Modern methods are not allowed to redefine the historical targets after seeing the data.

### HAPI compatibility note

The PDS/PPI endpoint currently advertises HAPI 3.1. Version 3 renamed the request parameters to `dataset`, `start`, and `stop`. The downloader uses those native HAPI-3 names first and automatically falls back to the older `id`, `time.min`, and `time.max` form if necessary. Server-side HTTP errors are reported with the attempted URL and response body instead of being mislabeled as Internet/DNS failures.

### NASA/PDS HAPI daily-CPI issue and v0.3.2 acquisition path

A real Phase-A run on 2026-09-25 established that the PDS/PPI HAPI server
advertises HAPI 3.1 and serves metadata correctly, but its `/data` conversion
for the Pioneer daily CPI collection returns HTTP 500 with the server message
`unable to parse time: 72Z`. This is an upstream conversion problem involving
legacy two-digit years, not a local network problem.

For that reason the default CPI path is now:

1. download the official legacy NASA/SPDF 24-hour ASCII CPI product;
2. parse it conservatively into the canonical dissertation columns;
3. retain HAPI metadata/provenance when available;
4. use HAPI CPI data only as a fallback/diagnostic path.

This path is also preferable for the replication objective because the legacy
SPDF daily file is closer to the data product used in 1998. To force a path:

```bash
python scripts/fetch_nasa.py --spacecraft both --cpi-source spdf
python scripts/fetch_nasa.py --spacecraft both --cpi-source hapi
```

The one-command launchers continue to use `auto`, which prefers SPDF.

### v0.3.3: particle backbone discovery

The v0.3.2 Windows run reached the real SPDF `p10cpi.h24` file. It contained
exactly 7042 numeric records, and every record had 15 numeric fields. Inspection
of the actual rows shows the layout to be:

```text
YY  DOY  HOUR  RID2P RID2HE RID3P RID3HE RID4P RID4HE RID5P RID5HE RID5E1 RID5E2 RID7 RID7ZG5
```

That is the 12-channel energetic-particle sequence documented in the thesis and
in OMNIWeb, preceded by the archival time triplet. The previous parser wrongly
interpreted the first `72` as a decimal-year field and therefore expected XYZ
columns that are not in this archive file. v0.3.3 fixes this without inventing
coordinates.

Crucially, the two first locked historical tests do not need XYZ or solar-wind
speed:

1. P10/P11 autocorrelation target near 469–470 days;
2. unshifted P10↔P11 common-date cross-correlation target near 339–341 days.

The propagation-shift experiment remains a later extension and will run only
when genuine XYZ and solar-wind fields have been recovered.

### v0.3.4 cross-correlation diagnostic
The 11–20 MeV Phase-A run now prints exact 1998 benchmark values and both lag signs. This is diagnostic only: the locked xcorr implementation and data pipeline are unchanged.

### v0.3.5: Phase-A historical benchmark reproduced

A real run against the recovered NASA/SPDF daily CPI files reproduces the printed 1998 benchmarks for the 11–20 MeV proton channel to rounding precision:

- P10 ACF: 469 samples, 0.1899745 (thesis: 469, 0.19)
- P11 ACF: 1697 samples, 0.2062353 and 467 samples, 0.1066832 (thesis: 1697/0.21 and 467/0.11)
- P10↔P11 CCF amplitudes and absolute lags: 1905/0.8742, 341/0.23870, 1524/0.12930, -214/0.08092, matching the thesis 1905/0.88, 341/0.24, 1524/0.13, -214/0.08 after translating the lag-sign convention (`thesis_lag = -scipy_lag`).

This establishes historical computational reproducibility. It does **not** yet establish that the quoted lag numbers are literal calendar-day periodicities. The common-date P10/P11 series contains missing days, so MATLAB/SciPy xcorr lags are sample-index lags. v0.3.5 therefore reports, for each benchmark lag, the actual distribution of calendar-day separations between rows that are that many common samples apart. This distinction becomes a central target of the modern re-analysis.

## Phase B2 — uneven-sampling multichannel screening (v0.5.0)

Phase B1 showed that the historical sample-index lags cannot automatically be interpreted as elapsed calendar days. Phase B2 therefore adds two interpolation-free screening tools before any expensive significance simulations:

1. a calendar-lag **Discrete Correlation Function (DCF)** on log10 positive flux, using the actual observation dates and 10-day lag bins;
2. a generalized/floating-mean **Lomb–Scargle** screen on log10 flux at the actual observation times.

Run all 12 particle channels with:

```text
run_phase_b2_windows.bat
```

Outputs:

```text
results/phase_b2_screening/multichannel_screening_summary.json
results/phase_b2_screening/multichannel_screening_summary.csv
```

The same two predeclared calendar bands are retained: 330–400 d and 430–520 d. The summary also reports the strongest sampling-window period in each band so a spectral peak close to the observation-window peak is explicitly flagged for later testing. B2 is a **screen**, not a significance claim; red-noise/surrogate inference is deferred to Phase B3 and will be run only for candidates that survive B2.

## Phase B3 — red-noise and alignment-surrogate significance (v0.6.0)

B2 is deliberately only a screen. Phase B3 adds explicit null-model testing while keeping the historical Phase A result frozen.

Run the first-pass all-channel test with:

```text
run_phase_b3_windows.bat
```

The default launcher uses 199 surrogates per channel and robustness variant so it can be completed on a normal local workstation in minutes. It is intended to identify which signals deserve a higher-resolution run. For any publication-facing candidate rerun, for example:

```text
.venv\Scripts\python.exe scripts\phase_b3_significance.py --simulations 1999
```

Phase B3 uses two complementary nulls:

1. **Single-spacecraft ACF:** a one-calendar-day AR(1) coefficient is estimated only from genuinely adjacent observed days. Red-noise surrogates are simulated on the full daily calendar and then sampled through the exact historical missing-data mask. The test statistic is the **maximum 10-day-binned DCF anywhere inside the predeclared band**, so the p-value already accounts for choosing the strongest bin within that band.
2. **P10/P11 CCF:** the complete P11 value+missingness pattern is randomly circularly shifted in calendar time by at least 700 days. This preserves P11's internal red-noise structure and sampling pattern while breaking its absolute alignment to P10.

Both raw log10 flux and a linearly detrended log10 robustness variant are tested. Across the 12 particle channels, both strict Holm family-wise adjusted p-values and Benjamini-Hochberg FDR q-values are reported separately for each band/statistic/variant.

Outputs:

```text
results/phase_b3_significance/significance_summary.json
results/phase_b3_significance/significance_summary.csv
```

A small surrogate p-value is not, by itself, a physical-causality claim. In particular, a cross-spacecraft result still requires time localization and a physically meaningful treatment of propagation before it can be interpreted heliophysically.

### Publication-resolution significance rerun
After the 199-simulation B3 screen, run:

```bat
run_phase_b3_publication_windows.bat
```

This reruns the complete 12-channel B3 family with 1999 surrogates per test and creates `results/phase_b3_publication/candidate_review.json`. The larger run is required because the 199-simulation p-value floor (0.005) makes Holm <0.05 mathematically impossible across 12 channels.

### Phase B4: morphology and localization
After `run_phase_b3_publication_windows.bat`, run:

```bat
run_phase_b4_windows.bat
```

This phase does **not** add a new significance claim. It checks whether B3 survivors are distinct lag peaks and whether they persist through time. Outputs are under `results/phase_b4_localization/`.

## Phase B5 — individual-spacecraft periodicity morphology (v0.9.0)

B4 demonstrated why a statistically unusual P10↔P11 lag must not be called a
period. B5 therefore returns to the individual spacecraft records. For each B4
candidate it computes, separately for P10 and P11 and for raw/detrended log10
flux:

- calendar-time ACF morphology over 200–700 d;
- Lomb–Scargle morphology over 200–700 d using the actual irregular times;
- 5-year sliding Lomb–Scargle localization;
- descriptive P10/P11 period concordance in common 5-year windows.

No missing values are interpolated, and B5 creates no new significance claim.
Its purpose is to decide which B3/B4 lag features also look like periodic
structure in the individual records.

Run:

```bat
run_phase_b5_windows.bat
```

Primary output:

```text
results/phase_b5_periodicity/periodicity_summary.json
```

## Phase B6 — all-channel periodicity atlas + sampling-window audit (v1.0.0)

B5 followed only the cross-spacecraft B4 survivors. That is appropriate for
following up unusual P10/P11 lags, but it is not a valid selection rule for a
survey of periodicity inside the individual spacecraft records. A periodic
feature can exist in P10 or P11 without producing a strong cross-spacecraft CCF.

B6 therefore removes that selection and audits **all 12 CPI particle channels**
in both predeclared bands. For P10 and P11, raw-log10 and linearly detrended
log10, it reports calendar-ACF morphology, irregular-time Lomb–Scargle
morphology, the exact positive-flux sampling-window spectrum, and the frequency
separation between signal and sampling-window peaks in Rayleigh-resolution
units. It also attaches the already-frozen publication-resolution B3 ACF
p/Holm/FDR values and compares P10/P11 spectra over the exact same valid calendar
span for each channel.

B6 adds **no new p-values**. Its purpose is to avoid selection-by-CCF and to
identify a small, defensible set for later dedicated Lomb–Scargle/time-frequency
significance testing.

Run:

```bat
run_phase_b6_windows.bat
```

Primary outputs:

```text
results/phase_b6_periodicity_atlas/periodicity_atlas_summary.json
results/phase_b6_periodicity_atlas/periodicity_atlas_summary.csv
results/phase_b6_periodicity_atlas/candidate_flags.json
```

### Phase B7: sampling-window co-tracking
After B6, run `run_phase_b7_windows.bat` (Windows) or `run_phase_b7_unix.sh` (Unix). The main outputs are `results/phase_b7_sampling_tracking/sampling_tracking_summary.json` and `sampling_tracking_flags.json`. B7 compares signal-period tracks with the observation-window spectrum in the same 5-year windows; it is a descriptive alias audit and does not create new p-values.

## Phase B8 — Rayleigh-resolution stress test (v1.2.0)

B7's five-year windows are useful for time localization but have coarse
frequency resolution. Around a period `P`, one Rayleigh frequency interval
corresponds approximately to `P^2/T` days in period space, where `T` is the
window span. At ~365–500 d, a five-year window can therefore make a
`<=1 Rayleigh` warning cover a large fraction of the entire predefined band.

B8 repeats the signal-vs-sampling-window audit at 5, 8, 10, 12 and 15-year
windows. The purpose is not to create another significance test, but to ask
whether sampling-window proximity survives as spectral resolution improves.

Run:

```bat
run_phase_b8_windows.bat
```

Primary outputs:

```text
results/phase_b8_resolution_stress/resolution_stress_summary.json
results/phase_b8_resolution_stress/resolution_stress_flags.json
```

## Phase B9 — candidate-focused red-noise spectral significance (v1.3.0)

B8 showed that the universal five-year `<=1 Rayleigh` warnings were partly a
resolution effect: some peaks separate from the sampling-window spectrum in
12–15 year windows. B9 therefore freezes a small long-baseline candidate set
before adding any new p-values:

- He 20–24 MeV/n, 430–520 d (primary after B8)
- He 11–20 MeV/n, 430–520 d (secondary long-baseline comparator)

For each candidate, P10/P11 and raw-log10/linear-detrended-log10 are tested
against a daily AR(1) null. Surrogates are simulated on the complete daily
calendar and then sampled through the exact historical positive-flux mask. The
statistic is the maximum classical Lomb–Scargle power anywhere inside the
locked 430–520 d band, so the within-band frequency search is included in the
Monte-Carlo null. B9 applies Holm and Benjamini–Hochberg adjustments across the
four locked spacecraft/channel tests separately for each robustness variant.

B9 is intentionally labelled **post-selection evidence** because the candidate
set was chosen after earlier phases examined the same Pioneer data. It is not an
independent replication and does not by itself establish a physical origin.

Run:

```bat
run_phase_b9_windows.bat
```

Primary outputs:

```text
results/phase_b9_spectral_significance/spectral_significance_summary.json
results/phase_b9_spectral_significance/spectral_significance_summary.csv
```

### Phase B10 — intermittent same-window/same-period significance
After B9, use B10 only to test the distinct intermittent-signal hypothesis, not to rescue a non-significant persistent-periodicity claim. The launcher scans fixed 8-year common-calendar windows (2-year step) inside the locked 430–520 d band and calibrates the maximum over all windows/frequencies with 1999 mask-preserving daily AR(1) surrogates. The joint statistic requires P10 and P11 power at the same period in the same window.

Windows:
```bat
run_phase_b10_windows.bat
```
Output: `results/phase_b10_intermittent_significance/intermittent_significance_summary.json`.
