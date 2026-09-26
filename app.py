from __future__ import annotations

from pathlib import Path
import pandas as pd
import plotly.express as px
import streamlit as st

from src.io import load_auto, apply_p11_cpi_quality_mask, data_quality_summary
from src.schema import PARTICLE_COLUMNS, CHANNEL_LABELS
from src.replication import (
    thesis_fft,
    autocorrelation,
    strongest_positive_lags,
    crosscorrelation_common_dates,
    propagation_shifted_crosscorrelation,
    strongest_in_band,
)
from src.modern import lomb_scargle_periodogram, band_peak, red_noise_band_test
from src.provenance import PDS_SOURCES

st.set_page_config(page_title="Pioneer Periodicity Explorer", layout="wide")
st.title("Pioneer Periodicity Explorer")
st.caption("Local-first reconstruction and modern falsification of the 1998 Pioneer 10/11 periodicity analysis.")

with st.sidebar:
    mode = st.radio("Analysis mode", ["1998 Replication", "Modern Re-analysis"])
    st.divider()
    st.markdown("**Inputs**")
    p10_file = st.file_uploader("Pioneer 10 — 17-column .dat or CSV", type=["dat", "txt", "csv"], key="p10")
    p11_file = st.file_uploader("Pioneer 11 — 17-column .dat or CSV", type=["dat", "txt", "csv"], key="p11")
    apply_quality = st.checkbox("Apply documented P11 CPI quality mask", value=True)

st.markdown("### Research status")
st.markdown(
    "Replication targets are fixed **before** modern analysis: recurring ~469–470 d autocorrelation "
    "and ~339–341 d P10↔P11 cross-correlation. Modern methods are used to test, not tune, those findings."
)

if p10_file is None and p11_file is None:
    st.info("Load at least one spacecraft file. For P10↔P11 correlation and propagation correction, load both.")
    with st.expander("Official NASA/PDS provenance"):
        st.dataframe(pd.DataFrame(PDS_SOURCES).T, use_container_width=True)
    st.stop()

frames: dict[str, pd.DataFrame] = {}
try:
    if p10_file is not None:
        frames["Pioneer 10"] = load_auto(p10_file, p10_file.name)
    if p11_file is not None:
        p11 = load_auto(p11_file, p11_file.name)
        frames["Pioneer 11"] = apply_p11_cpi_quality_mask(p11) if apply_quality else p11
except Exception as e:
    st.error(str(e))
    st.stop()

spacecraft = st.selectbox("Primary spacecraft", list(frames))
df = frames[spacecraft]
available = [c for c in PARTICLE_COLUMNS if c in df.columns]
if not available:
    st.error("No recognized CPI particle channels were found.")
    st.stop()
channel = st.selectbox("Particle / energy channel", available, format_func=lambda c: CHANNEL_LABELS[c])

valid_dates = df.date.dropna()
if valid_dates.empty:
    st.error("No valid dates found.")
    st.stop()
start, stop = st.slider(
    "Date window",
    min_value=valid_dates.min().date(),
    max_value=valid_dates.max().date(),
    value=(valid_dates.min().date(), valid_dates.max().date()),
)

def crop(x: pd.DataFrame) -> pd.DataFrame:
    return x[(x.date.dt.date >= start) & (x.date.dt.date <= stop)].copy()

q = crop(df)
series = q[["date", channel]].dropna()

c1, c2, c3, c4 = st.columns(4)
c1.metric("Rows", f"{len(q):,}")
c2.metric("Finite channel values", f"{len(series):,}")
c3.metric("Start", str(q.date.min().date()))
c4.metric("End", str(q.date.max().date()))

fig = px.line(series, x="date", y=channel, title=f"{spacecraft} — {CHANNEL_LABELS[channel]}")
st.plotly_chart(fig, use_container_width=True)

with st.expander("Data quality summary"):
    st.dataframe(data_quality_summary(q), use_container_width=True)

if mode == "1998 Replication":
    st.subheader("FFT — reconstructed 1998 recipe")
    fft = thesis_fft(series[channel].to_numpy())
    if not fft.empty:
        f = px.line(fft[(fft.period_days >= 2) & (fft.period_days <= 4000)], x="period_days", y="power", log_x=True)
        f.add_vline(x=340, line_dash="dash")
        f.add_vline(x=470, line_dash="dash")
        st.plotly_chart(f, use_container_width=True)

    st.subheader("Autocorrelation")
    acf = autocorrelation(series[channel].to_numpy())
    positive = acf[(acf.lag_days >= 0) & (acf.lag_days <= min(7000, max(0, len(series)-1)))]
    fa = px.scatter(positive, x="lag_days", y="coefficient")
    fa.add_vline(x=340, line_dash="dash")
    fa.add_vline(x=470, line_dash="dash")
    st.plotly_chart(fa, use_container_width=True)
    target470 = strongest_in_band(acf, 430, 520)
    st.markdown(f"**Pre-declared ~470 d band:** strongest lag = `{target470['lag_days']}` d, coefficient = `{target470['coefficient']:.4f}`")
    st.markdown("**Local positive maxima with coefficient ≥ 0.1**")
    st.dataframe(strongest_positive_lags(acf, threshold=0.1), use_container_width=True)

    if "Pioneer 10" in frames and "Pioneer 11" in frames:
        st.divider()
        st.subheader("Pioneer 10 ↔ Pioneer 11 cross-correlation")
        p10 = crop(frames["Pioneer 10"])
        p11 = crop(frames["Pioneer 11"])
        if channel not in p10.columns or channel not in p11.columns:
            st.warning("The selected channel is not available in both inputs.")
        else:
            cc = crosscorrelation_common_dates(p10, p11, channel)
            fc = px.scatter(cc, x="lag_days", y="coefficient", title="Common-date cross-correlation")
            fc.add_vline(x=340, line_dash="dash")
            st.plotly_chart(fc, use_container_width=True)
            t340 = strongest_in_band(cc, 315, 365)
            st.markdown(f"**Pre-declared ~340 d band:** strongest lag = `{t340['lag_days']}` d, coefficient = `{t340['coefficient']:.4f}`")

            required_geom = {"x_au", "y_au", "z_au", "solar_wind_speed_km_s"}
            if required_geom.issubset(p10.columns) and required_geom.issubset(p11.columns):
                try:
                    cc_shift, shifted = propagation_shifted_crosscorrelation(p10, p11, channel, thesis_exact=True)
                    fs = px.scatter(cc_shift, x="lag_days", y="coefficient", title="After dissertation propagation-time shift")
                    fs.add_vline(x=340, line_dash="dash")
                    st.plotly_chart(fs, use_container_width=True)
                    st.caption(
                        f"Common aligned samples: {len(shifted):,}; median propagation correction: "
                        f"{shifted.delay_days.median():.0f} d. Uses the thesis factor 1731 and index-shift algorithm."
                    )
                except Exception as e:
                    st.warning(f"Propagation correction could not be computed: {e}")
            else:
                st.info("Propagation correction requires x/y/z heliocentric coordinates and solar-wind speed in both files.")
else:
    st.subheader("Lomb–Scargle periodogram — gaps preserved")
    detrend = st.selectbox("Detrending", ["linear", "constant"])
    pg = lomb_scargle_periodogram(q, channel, detrend=detrend)
    if not pg.empty:
        fp = px.line(pg, x="period_days", y="power", log_x=True)
        fp.add_vrect(x0=430, x1=520, opacity=0.1, line_width=0)
        fp.add_vrect(x0=315, x1=365, opacity=0.1, line_width=0)
        st.plotly_chart(fp, use_container_width=True)
        p470 = band_peak(pg, 430, 520)
        p340 = band_peak(pg, 315, 365)
        st.dataframe(pd.DataFrame([
            {"target": "~1.3 y / 470 d", **p470},
            {"target": "~340 d", **p340},
        ]), use_container_width=True)

        with st.expander("Exploratory AR(1) red-noise falsification test"):
            sims = st.slider("Monte Carlo simulations", 100, 2000, 500, 100)
            if st.button("Run red-noise tests"):
                with st.spinner("Running AR(1) surrogates..."):
                    r470 = red_noise_band_test(q, channel, 430, 520, simulations=sims, detrend=detrend)
                    r340 = red_noise_band_test(q, channel, 315, 365, simulations=sims, detrend=detrend)
                st.dataframe(pd.DataFrame([
                    {"target": "430–520 d", **r470},
                    {"target": "315–365 d", **r340},
                ]), use_container_width=True)
                st.caption("Exploratory only: this is a simple AR(1) null model, not yet the final publication-grade significance framework.")
