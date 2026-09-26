# Analysis freeze - 26 September 2026

This file records the stopping state associated with the public technical report v1.0 and software v1.4.0.

## Frozen conclusions

- The principal historical 11-20 MeV proton ACF/CCF benchmarks reproduce to rounding precision.
- The historical cross-correlation lag axis is a retained-sample index on a gapped common-date series, not a uniformly elapsed calendar-time axis.
- The historical 341-sample feature therefore cannot be interpreted directly as a fixed 341-day physical lag.
- Calendar-aware and spectral screening identify descriptive structure in several channels, but morphology alone is not treated as periodicity or significance.
- Candidate-focused full-record ~1.3-year Lomb-Scargle tests (B9) are non-significant under the implemented mask-preserving daily AR(1) nulls.
- The final intermittent same-window/same-period test (B10) produces a detrended He 20-24 MeV/n feature near 459 d in 1973-1981 with nominal p=0.027 and Holm-adjusted p=0.054; the corresponding raw result is non-significant.
- The stopping rule was invoked after B10: no further tuning of bands, windows, detrending, candidate family, or null model is used to search the same data for a more favourable p-value.

## Interpretation boundary

Part I (historical reconstruction and time-axis audit) is directly auditable against the dissertation and archival data. Part II is an AI-assisted computational robustness assessment. It is reproducible but not independently specialist-validated signal analysis.
