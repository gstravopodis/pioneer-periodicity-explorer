# Data provenance lock

## Historical source described by the dissertation

The dissertation used two 17-column daily files named `p10cpi.dat` and `p11cpi.dat`. Table 3 defines:

1. decimal sampling date
2. protons 11–20 MeV
3. alpha particles 11–20 MeV
4. protons 20–24 MeV
5. alpha particles 20–24 MeV
6. protons 24–29 MeV
7. alpha particles 24–29 MeV
8. protons 29–67 MeV
9. alpha particles 29–67 MeV
10. electrons 7–17 MeV
11. 2 × minimum ionizing
12. ions >67 MeV/nucleon
13. ions Z>5 and >67 MeV/nucleon
14. heliocentric x (AU)
15. heliocentric y (AU)
16. heliocentric z (AU)
17. solar-wind plasma speed (km/s)

The dissertation states that P11 columns 6–13 are unreliable from day 239 of 1980 onward because of CPI failure. It also states that the loaded Pioneer 10 matrix had 7042 samples; this is treated as a validation target, not forced by filtering.

## Current NASA source strategy

### Particle/coordinate backbone — CPI daily
- P10: `urn:nasa:pds:p10-cpi-jup-cal:data-1day`, DOI `10.17189/yfb9-ve02`
- P11: `urn:nasa:pds:p11-cpi:data-1day`, DOI `10.17189/jqm4-yq88`

### Solar-wind speed — PA hourly
- P10: `urn:nasa:pds:p10-pa:data-avg-1hr`, DOI `10.17189/75z7-gw37`
- P11: `urn:nasa:pds:p11-pa:data-avg-1hr`, DOI `10.17189/rxvv-eb83`

The PA products contain hourly bulk velocity. For reconstruction of the thesis's one-row-per-day file, finite positive hourly values are reduced by arithmetic daily mean. This is explicitly a reconstruction decision and is written into the generated provenance JSON.

### Trajectory fallback / validation
- P10: `urn:nasa:pds:p10-pa:data-traj`, DOI `10.17189/99wc-jw74`
- P11: `urn:nasa:pds:p11-pa:data-traj`, DOI `10.17189/cjw6-7b69`

The current assembler first uses daily heliocentric coordinates available with CPI. PA trajectory data are retained as an independent source for a later cross-check rather than silently mixed into the baseline.

## Mapping rule

Field names and HAPI metadata are saved before transformation. Known SPDF aliases (`RID2P`, `RID2HE`, … `RID7ZG5`) are preferred. Description/units matching is used only as a fallback. Any top-score tie raises an error; ambiguous fields are never guessed.

## Reproducibility rule

No transformed dataset overwrites raw NASA data. Each assembled file has a JSON report with:
- source collection ID and DOI;
- HAPI parameter mapping;
- daily aggregation rule;
- row count and date coverage;
- missing/zero/finite counts per thesis field;
- SHA-256 hashes of raw and processed files;
- explicit 7042-row Pioneer-10 historical check.
