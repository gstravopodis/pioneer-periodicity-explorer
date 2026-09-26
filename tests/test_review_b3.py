from scripts.review_b3 import build_review


def _node(p, q, h, lag=334.5, dcf=.2, ccf=False):
    key = "p_value_band_max_circular_shift" if ccf else "p_value_band_max_ar1"
    return {"status":"ok", "bands": {
        "annual_330_400": {"observed":{"lag_calendar_days":lag,"dcf":dcf}, key:p, "p_holm_across_12_channels":h, "q_bh_fdr_across_12_channels":q},
        "one_point_three_year_430_520": {"observed":{"lag_calendar_days":lag,"dcf":dcf}, key:p, "p_holm_across_12_channels":h, "q_bh_fdr_across_12_channels":q},
    }}


def test_review_exposes_mc_floor_and_robust_ccf():
    obj={"simulations_per_test":199,"channels":[{"channel":"he_20_24_mev_n","label":"Helium 20–24 MeV/n","variants":{}}]}
    for v in ("raw_log10","linear_detrended_log10"):
        obj["channels"][0]["variants"][v]={
            "p10_acf":_node(.005,.01,.06),
            "p11_acf":_node(.005,.01,.06),
            "p10_p11_ccf":_node(.005,.02,.06,ccf=True),
        }
    review, df=build_review(obj)
    assert review["monte_carlo_p_floor"] == 0.005
    assert abs(review["holm_floor_if_12_equal_min_p"] - 0.06) < 1e-12
    assert len(review["robust_cross_spacecraft_fdr05_both_variants"]) == 2  # one per fixed band
    assert df["ccf_robust_both_variants_fdr05"].all()
