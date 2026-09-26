# Pioneer Periodicity Explorer

**Public analysis freeze: v1.4.0 - 26 September 2026**

This repository reconstructs the principal numerical analyses reported in a 1998 University of Athens BSc Physics dissertation on Pioneer 10 and Pioneer 11 energetic-particle measurements, then records an AI-assisted modern robustness reassessment.

## Interpretation boundary

The project intentionally separates two layers:

1. **Historical reconstruction.** The principal 11-20 MeV proton ACF/CCF benchmarks printed in the 1998 dissertation reproduce to rounding precision from the recovered archival CPI particle products.
2. **Modern reassessment.** Calendar-aware correlation, Lomb-Scargle screening, sampling-window diagnostics, AR(1) surrogates, multiple-testing corrections, and intermittent same-window/same-period tests were developed with substantial AI assistance. These methods are reproducible from the code and frozen outputs but have not been independently validated by a statistical signal-analysis specialist.

The central robust finding is methodological: the historical `xcorr` lag coordinate is a retained-sample index in a gapped sequence and therefore is not automatically an elapsed calendar-day period. The modern stress tests did not provide robust confirmation of the candidate periods under the implemented null models.

## Accompanying report

See `report/Pioneer_10_11_Reproducibility_Report_v1.0.pdf` in the public-release package. The report is a non-peer-reviewed technical research note, not a journal article.

## Data provenance

The code downloads or reconstructs inputs from official NASA archival products; raw NASA data are **not redistributed in this repository**.

Primary daily CPI collections:

- Pioneer 10 CPI one-day: DOI `10.17189/yfb9-ve02`
- Pioneer 11 CPI one-day: DOI `10.17189/jqm4-yq88`

See `docs/DATA_PROVENANCE.md` for source identifiers and acquisition details.

## Quick start

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
pytest -q
```

The frozen public codebase has **52 passing tests** in the release environment.

Historical acquisition/replication:

```bash
python scripts/fetch_nasa.py --spacecraft both
python scripts/run_phase_a.py --channel p_11_20_mev --skip-download
```

Windows launchers are also provided for the sequential Phase A/B analyses.

## Frozen audit outputs

`results/frozen/` contains the machine-readable summaries used for the public technical report, from the historical Phase A replication through the final B10 intermittent test. These are retained as an audit trail; they are not independent expert validation of the Part II methods.

## Research protocol and history

- `RESEARCH_PROTOCOL.md` - analysis locks, guardrails, and stopping rule.
- `CHANGELOG.md` - implementation history.
- `docs/DEVELOPMENT_HISTORY.md` - the longer development-era README retained for traceability.
- `docs/ANALYSIS_FREEZE.md` - final interpretation and stopping status for this release.

## AI assistance

Advanced time-series methodology, software implementation, iterative robustness testing, figure generation, and substantial drafting assistance were developed with ChatGPT (OpenAI). The project does not present the Part II signal-analysis methodology as an independently specialist-validated contribution.

## License

Software: MIT License. The accompanying technical report is distributed separately and is recommended for release under CC BY 4.0.
