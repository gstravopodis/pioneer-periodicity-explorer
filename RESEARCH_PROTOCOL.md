# Research protocol — Pioneer Periodicity Explorer

**Protocol version: 0.5 — 2026-09-26**

## Purpose
Reconstruct the 1998 Pioneer 10/11 analysis before any modern re-analysis. The software is an instrument for testing the old findings, not for assuming they are correct.

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
Only after Phase A is frozen. Phase B1 is pre-specified as the calendar-time test: map historical sample-index ACF/CCF lags to actual calendar-day separation distributions; recompute ACF/CCF at exact calendar-day offsets; compare annual and ~1.3-year bands; compute Lomb–Scargle spectra and the observation-window spectrum.

Subsequent phases add multichannel screening, mask-preserving AR(1) and alignment-surrogate tests, morphology/localization audits, sampling-window stress tests, candidate-focused spectral significance, and one final intermittent same-window/same-period test.

## Guardrails
- Keep replication and modern analysis outputs visibly separate.
- Never tune preprocessing merely to maximize the 340- or 470-day peaks.
- Record every transformation and exclusion.
- Prefer NASA/PDS/SPDF archived products and preserve provenance/DOIs.
- Freeze analysis code used for public figures/results.
- Treat B9/B10 as post-selection evidence, not independent confirmation.
- Stop after B10 rather than change bands/windows/nulls to seek a lower p-value.

## Final locked status
The historical computation is reproduced. The historical lag coordinate is a retained-sample index in a gapped sequence and cannot be read directly as elapsed calendar days. The modern AI-assisted robustness pipeline did not provide robust confirmation of a persistent or shared intermittent ~1.3-year signal under the implemented null models. A detrended He 20–24 MeV/n intermittent feature was nominally significant but did not survive the locked family-wise correction and was not reproduced in the raw-log variant.
