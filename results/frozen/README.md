# Frozen result summaries

These JSON files are the machine-readable summaries retained for the public report. They are included to make the computational audit trail inspectable without requiring a complete rerun.

- `phase_a_replication_summary.json`: historical 11-20 MeV proton ACF/CCF replication and sample-lag/calendar-separation diagnostics.
- `phase_b2_*`: uneven-sampling multichannel screening.
- `phase_b3_*`: publication-resolution surrogate significance and candidate review.
- `phase_b4_*` to `phase_b8_*`: morphology, individual-spacecraft periodicity, sampling-window, and resolution stress tests.
- `phase_b9_*`: locked full-record ~1.3-year spectral significance.
- `phase_b10_*`: final intermittent same-window/same-period significance test.

The absence of raw NASA particle files is intentional. The acquisition code and official source identifiers are preserved separately so inputs can be re-fetched from NASA archives.
