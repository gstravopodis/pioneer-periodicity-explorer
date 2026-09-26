from __future__ import annotations

# Exact 17-column layout documented in Table 3 of the 1998 dissertation.
THESIS_COLUMNS = [
    "time_decimal_year",
    "p_11_20_mev",
    "he_11_20_mev_n",
    "p_20_24_mev",
    "he_20_24_mev_n",
    "p_24_29_mev",
    "he_24_29_mev_n",
    "p_29_67_mev",
    "he_29_67_mev_n",
    "e_7_17_mev",
    "minimum_ionizing_x2",
    "ions_gt67_mev_n",
    "ions_zgt5_gt67_mev_n",
    "x_au",
    "y_au",
    "z_au",
    "solar_wind_speed_km_s",
]

CHANNEL_LABELS = {
    "p_11_20_mev": "Protons 11–20 MeV",
    "he_11_20_mev_n": "Helium 11–20 MeV/n",
    "p_20_24_mev": "Protons 20–24 MeV",
    "he_20_24_mev_n": "Helium 20–24 MeV/n",
    "p_24_29_mev": "Protons 24–29 MeV",
    "he_24_29_mev_n": "Helium 24–29 MeV/n",
    "p_29_67_mev": "Protons 29–67 MeV",
    "he_29_67_mev_n": "Helium 29–67 MeV/n",
    "e_7_17_mev": "Electrons 7–17 MeV",
    "minimum_ionizing_x2": "2 × minimum-ionizing",
    "ions_gt67_mev_n": "Ions >67 MeV/n",
    "ions_zgt5_gt67_mev_n": "Z>5 ions >67 MeV/n",
}

PARTICLE_COLUMNS = list(CHANNEL_LABELS)

# Dissertation footnote: P11 columns 6–13 are unreliable from day 239 of 1980 onward.
P11_AFFECTED_COLUMNS = THESIS_COLUMNS[5:13]

# Common field names used by SPDF/OMNIWeb exports.
NASA_FIELD_ALIASES = {
    "RID2P": "p_11_20_mev",
    "RID2HE": "he_11_20_mev_n",
    "RID3P": "p_20_24_mev",
    "RID3HE": "he_20_24_mev_n",
    "RID4P": "p_24_29_mev",
    "RID4HE": "he_24_29_mev_n",
    "RID5P": "p_29_67_mev",
    "RID5HE": "he_29_67_mev_n",
    "RID5E1": "e_7_17_mev",
    "RID5E2": "minimum_ionizing_x2",
    "RID7": "ions_gt67_mev_n",
    "RID7ZG5": "ions_zgt5_gt67_mev_n",
    "X": "x_au",
    "Y": "y_au",
    "Z": "z_au",
    "SW_SPEED": "solar_wind_speed_km_s",
    "SOLAR_WIND_SPEED": "solar_wind_speed_km_s",
    "VSW": "solar_wind_speed_km_s",
}
