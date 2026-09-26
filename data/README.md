# Data directory

Do not commit large raw NASA data products by default.

Expected local inputs can be any of:

1. Historical 17-column whitespace files matching the dissertation (`p10cpi.dat`, `p11cpi.dat`).
2. CSV files with the dissertation field names.
3. OMNI/SPDF-style CSV exports containing particle fields such as `RID2P`, `RID2HE`, `RID3P`, `RID3HE`, `RID4P`, `RID4HE`, `RID5P`, `RID5HE`, `RID5E1`, `RID5E2`, `RID7`, `RID7ZG5`, plus date information.

Official daily PDS4 collection inventories can be downloaded into `data/pds_manifests/` with:

```bash
python scripts/fetch_pds_manifests.py
```

PDS provenance:
- Pioneer 10 daily CPI collection DOI: 10.17189/yfb9-ve02
- Pioneer 11 daily CPI collection DOI: 10.17189/jqm4-yq88
