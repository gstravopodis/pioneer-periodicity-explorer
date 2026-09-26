from __future__ import annotations

HAPI_BASE_URL = "https://pds-ppi.igpp.ucla.edu/hapi"

# PDS/PPI HAPI collection IDs.  These are intentionally kept separate by
# measurement family so the historical 17-column file can be rebuilt from
# independently archived particle, plasma and trajectory products.
HAPI_DATASETS = {
    "p10": {
        "cpi_daily": {
            "id": "urn:nasa:pds:p10-cpi-jup-cal:data-1day",
            "doi": "10.17189/yfb9-ve02",
            "coverage": ("1972-03-03", "1992-09-01"),
        },
        "cpi_hourly": {
            "id": "urn:nasa:pds:p10-cpi-jup-cal:data-1hr",
            "doi": "10.17189/1mq1-dk04",
            "coverage": ("1972-03-03", "1992-09-01"),
        },
        "pa_hourly": {
            "id": "urn:nasa:pds:p10-pa:data-avg-1hr",
            "doi": "10.17189/75z7-gw37",
            "coverage": ("1972-04-18", "1995-09-07"),
        },
        "pa_trajectory": {
            "id": "urn:nasa:pds:p10-pa:data-traj",
            "doi": "10.17189/99wc-jw74",
            "coverage": ("1972-03-24", "1995-12-31"),
        },
    },
    "p11": {
        "cpi_daily": {
            "id": "urn:nasa:pds:p11-cpi:data-1day",
            "doi": "10.17189/jqm4-yq88",
            "coverage": ("1973-04-06", "1992-11-28"),
        },
        "cpi_hourly": {
            "id": "urn:nasa:pds:p11-cpi:data-1hr",
            "doi": "10.17189/k63r-df32",
            "coverage": ("1973-04-06", "1992-11-28"),
        },
        "pa_hourly": {
            "id": "urn:nasa:pds:p11-pa:data-avg-1hr",
            "doi": "10.17189/rxvv-eb83",
            "coverage": ("1973-04-21", "1992-05-31"),
        },
        "pa_trajectory": {
            "id": "urn:nasa:pds:p11-pa:data-traj",
            "doi": "10.17189/cjw6-7b69",
            "coverage": ("1974-01-01", "1999-01-01"),
        },
    },
}

PDS_SOURCES = {
    "Pioneer 10 CPI daily": {
        "collection": "urn:nasa:pds:p10-cpi-jup-cal:data-1day::1.0",
        "doi": "10.17189/yfb9-ve02",
        "coverage": "1972-03-03 to 1992-08-31",
        "release": "p10-cpi-jup-cal-20250814",
        "manifest": "https://pds.nasa.gov/data/pds4/releases/ppi/p10-cpi-jup-cal-20250814/collection-data-1day-1.0.csv",
    },
    "Pioneer 11 CPI daily": {
        "collection": "urn:nasa:pds:p11-cpi:data-1day::1.0",
        "doi": "10.17189/jqm4-yq88",
        "coverage": "1973-04-06 to 1992-11-27",
        "release": "p11-cpi-20260422",
        "manifest": "https://pds.nasa.gov/data/pds4/releases/ppi/p11-cpi-20260422/collection-data-1day-1.0.csv",
    },
}

OMNI_SOURCES = {
    "Pioneer 10 CPI hourly": "https://omniweb.gsfc.nasa.gov/ftpbrowser/p10_part_cpi_1h.html",
    "Pioneer 11 CPI hourly": "https://omniweb.gsfc.nasa.gov/ftpbrowser/p11_part_cpi_1h.html",
}


# Historical SPDF daily CPI files.  These are the legacy 24-hour ASCII
# products closest to the files used in the 1998 dissertation.  Pioneer 10's
# p10cpi.h24 is explicitly cited as the daily CPI source in Vogt et al. (2018,
# A&A 613, A28).  Pioneer 11 uses the symmetric SPDF archive layout exposed by
# the official OMNIWeb CPI page.
SPDF_CPI_DAILY = {
    "p10": {
        "url": "https://spdf.gsfc.nasa.gov/pub/data/pioneer/pioneer10/particle/cpi/ip_1day_ascii/p10cpi.h24",
        "archive_dir": "https://spdf.gsfc.nasa.gov/pub/data/pioneer/pioneer10/particle/cpi/ip_1day_ascii/",
        "note": "Legacy SPDF 24-hour CPI ASCII; preferred for 1998 replication fidelity.",
    },
    "p11": {
        "url": "https://spdf.gsfc.nasa.gov/pub/data/pioneer/pioneer11/particle/cpi/ip_1day_ascii/p11cpi.h24",
        "archive_dir": "https://spdf.gsfc.nasa.gov/pub/data/pioneer/pioneer11/particle/cpi/ip_1day_ascii/",
        "note": "Legacy SPDF 24-hour CPI ASCII in the official Pioneer 11 CPI archive.",
    },
}
