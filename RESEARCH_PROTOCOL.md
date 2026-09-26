# Research protocol — Pioneer Periodicity Explorer

**Protocol version: 0.5 — 2026-09-26**

## Purpose
Reconstruct the 1998 Pioneer 10/11 analysis before any modern re-analysis.  The software is an instrument for testing the old findings, not for assuming they are correct.

## Phase A — replication lock
1. Use the dissertation's 17-column schema exactly.
2. Preserve the thesis-style daily sampling and original missing/invalid-data rules where documented.
3. Reconstruct its FFT recipe: zero-padding to the next power of two, DC removal, phase, and `2*abs(FFT)^2` power.
4. Reconstruct MATLAB `xcorr(...,'coeff')` autocorrelation/cross-correlation behavior.
5. Match P10 and P11 on common observation dates.
6. Reconstruct radial solar-wind propagation correction from spacecraft heliocentric coordinates and solar-wind speed.
7. Primary replication targets, declared in advance:
   - recurring autocorrelation feature around 469–470 days;
   - P10/P11 cross-correlation feature around 339–341 days.

No modern detrending or significance machinery is allowed to alter the Phase A target analysis.

## Phase B — robustness / falsification
Only after Phase A is frozen. Phase B1 is now pre-specified as the calendar-time test:
- map the historical sample-index ACF/CCF lags to their actual calendar-day separation distributions;
- recompute ACF/CCF at exact integer calendar-day offsets with missing dates left missing;
- report raw normalized overlap, Pearson correlation, and log10-flux Pearson correlation;
- compare the annual band (330–400 d) with the ~1.3-year band (430–520 d);
- compute Lomb–Scargle spectra at actual observation times;
- compute the observation-window spectrum to test annual/gap artefacts.

Subsequent Phase B work:
- Lomb–Scargle periodograms;
- linear/constant detrending sensitivity;
- windowing and gap sensitivity;
- red-noise null models;
- phase-randomised / block-surrogate tests;
- wavelet time-frequency localisation;
- solar-cycle sub-window analysis;
- channel-by-channel and radial-distance dependence;
- leave-one-epoch-out sensitivity;
- tests of annual/sampling/geometry artefacts around ~340 days.

## Paper-oriented hypotheses
H1: the ~470-day feature is a statistically defensible, intermittent ~1.3-year heliospheric quasi-periodicity rather than a stationary clock.

H2: the ~340-day P10/P11 feature either (a) survives sampling, propagation, geometry and red-noise controls, or (b) can be explained as a methodological/observational artefact. Both outcomes are scientifically useful if established rigorously.

## Guardrails
- Keep replication and modern analysis outputs visibly separate.
- Never tune preprocessing merely to maximize the 340- or 470-day peaks.
- Record every transformation and exclusion.
- Prefer NASA/PDS/SPDF archived products and preserve provenance/DOIs.
- Freeze analysis code used for paper figures with a version/tag.

## Phase B2 — pre-specified uneven-sampling screen
After B1, run every available CPI particle channel through:
- a log10-flux Discrete Correlation Function on true calendar lags, no interpolation, 10-day bins;
- generalized/floating-mean Lomb–Scargle on actual observation times;
- the same annual (330–400 d) and ~1.3-year (430–520 d) bands;
- direct comparison of Lomb–Scargle band peaks with the sampling-window spectrum.

B2 is only a screening stage. It cannot establish significance. Only features that survive calendar timing, appear in physically coherent channels/spacecraft, and are not trivially coincident with the sampling window advance to B3 red-noise/surrogate testing.

## Phase B3 protocol lock — 2026-09-25

B3 is defined before examining its surrogate p-values.

- Fixed bands remain **330–400 calendar days** and **430–520 calendar days**.
- The primary DCF statistic is the **maximum positive 10-day-binned DCF within each fixed band**. A single peak bin is never tested as though it had been pre-specified.
- P10 and P11 ACF null: daily AR(1) red noise, estimated from actual one-day observation pairs, simulated on a complete daily calendar, then sampled with the historical observation mask.
- P10/P11 CCF null: random circular calendar shifts of the full P11 values plus missingness mask, with an absolute shift of at least 700 days.
- Primary transform: log10 of strictly positive particle flux.
- Robustness transform: the same log10 series after linear detrending in calendar time.
- Multiplicity: Holm family-wise and Benjamini-Hochberg FDR adjustments across all available particle channels, separately by band, statistic, and robustness variant.
- The 199-surrogate launcher is **screening resolution only**. Any result considered for the manuscript must be rerun with at least 1999 surrogates and then subjected to a time-localization stage.
- No candidate will be called a physical periodicity solely because B3 rejects a null model.

## Phase B3 publication-resolution confirmation
The 199-surrogate B3 run is a screening pass only. Its Monte-Carlo p-value floor is 1/(199+1)=0.005. Because multiple-testing correction is performed over 12 channels, the minimum possible single-step Holm-adjusted value is 12*0.005=0.06; therefore no first-pass result can satisfy Holm < 0.05, regardless of effect size. Surviving findings must be rerun with at least 1999 surrogates across the full 12-channel family before any publication-level significance claim. The first-pass candidate selection is not treated as an independent replication.

Interpretation must also distinguish full-span P11 channels from columns 6-13, which are masked after day 239 of 1980 because the dissertation records CPI reliability loss. Cross-spacecraft findings in those affected channels therefore refer to the pre-failure overlap only.

## Phase B4 — morphology and time localization
Publication-resolution B3 candidates must be characterized before any physical-period claim. B4 scans a wider 200–700 d lag range, checks whether the predeclared-band maximum is an actual local peak or merely a band-edge value on a broad correlation tail, and localizes the feature in 5-year sliding windows. Windowed outputs are descriptive and do not add new significance claims.

## Phase B5 — individual-spacecraft periodicity morphology
A cross-spacecraft CCF lag is not itself an oscillation period. B5 therefore
requires candidate timescales to be examined independently within P10 and P11.
It uses calendar-time ACF morphology and Lomb–Scargle spectra on actual
observation times over 200–700 d, with raw-log10 and linearly detrended variants.
Five-year sliding spectra are descriptive localization only; adjacent windows
overlap and must not be counted as independent replications. A publication-level
periodicity claim requires a distinct within-spacecraft spectral/ACF feature and
subsequent red-noise/time-frequency confirmation, not only a B3 CCF p-value.


## Phase B6 selection rule

Periodicity analysis must not be conditioned on survival of the P10↔P11 CCF.
Cross-spacecraft lag correlation and within-spacecraft periodic structure are
distinct hypotheses. Phase B6 therefore audits all 12 particle channels in the
predeclared 330–400 d and 430–520 d bands. Any later publication-facing period
claim should require, at minimum: (i) calendar-aware single-spacecraft evidence,
(ii) a non-edge local spectral feature, (iii) robustness to linear detrending,
(iv) an explicit sampling-window/alias audit, and (v) dedicated red-noise/time-
frequency significance testing on a candidate set defined before that final test.

## Phase B7 — time-localized sampling-window co-tracking
B6 showed that no candidate cleared every desired publication-facing criterion because every candidate retained at least one one-Rayleigh sampling-window warning. B7 therefore does not add another significance claim. Instead, on exact common valid calendar spans, it tracks Lomb–Scargle signal peaks and sampling-window peaks through overlapping 5-year windows. A feature that repeatedly moves with the sampling-window frequency is treated as sampling-sensitive; this is a warning rather than proof of aliasing. Only candidates that remain morphologically coherent without persistent window co-tracking should proceed to expensive dedicated periodogram/time-frequency significance work.

## Phase B8: multi-scale Rayleigh-resolution stress test
B7 uses five-year sliding windows, which localize changes in time but have broad
Rayleigh frequency resolution. A <=1-Rayleigh signal/window separation in such
a short window is therefore a warning with limited discriminating power. B8
repeats the audit at progressively longer common-calendar windows (default 5,
8, 10, 12, 15 years) and reports the approximate period resolution `P^2/T`.
No p-values are created. A feature that remains unresolved from the sampling
window as T increases deserves greater sampling concern; a feature that
separates at longer T has passed only this specific alias check and still
requires dedicated red-noise spectral/time-frequency significance testing.

## Phase B9: locked candidate spectral significance
B8 is used only to freeze the next candidate set before B9. The locked set is:

- **Primary:** He 20–24 MeV/n, 430–520 calendar days.
- **Secondary comparator:** He 11–20 MeV/n, 430–520 calendar days.

B9 uses the same candidate definitions for P10 and P11 and for raw-log10 and
linearly detrended-log10 variants. The test statistic is the maximum classical
Lomb–Scargle power anywhere in the locked 430–520 d band. Daily AR(1) red-noise
surrogates are generated on the complete calendar and sampled through the exact
positive-flux observation mask. Thus sampling gaps are preserved and no missing
values are interpolated. Publication launcher resolution is 1999 surrogates.
Holm and Benjamini–Hochberg adjustments are applied over the four locked
spacecraft/channel tests separately for each robustness variant.

Because the candidate set was selected after B2–B8 inspection of these same data,
B9 p-values are **post-selection evidence** and are not described as independent
confirmation. A surviving B9 candidate must still undergo time-frequency
localization and sensitivity analysis before a physical quasi-periodicity claim.

### Phase B10 — intermittent shared-signal falsification
B9 is the full-record global spectral test. A non-significant B9 result rejects a persistent/coherent global 430–520 d spectral claim under the adopted AR(1) null, but does not by itself exclude an intermittent quasi-periodicity. B10 therefore freezes a separate alternative before inspection: He 20–24 MeV/n (primary) and He 11–20 MeV/n (comparator), 430–520 d, fixed 8-year windows stepped by 2 calendar years. The statistic is maximized over all windows and frequencies and is calibrated by the same maximization on mask-preserving daily AR(1) surrogates. The primary cross-spacecraft statistic requires P10 and P11 power at the same tested period in the same window. A non-significant joint B10 result is evidence against the remaining shared intermittent-signal hypothesis under this null model; a significant result would still be post-selection evidence and would require independent/external validation before a physical claim.
