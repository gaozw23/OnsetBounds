from onsetbounds.inference import aggregate_contributions, patient_contribution


def test_ppv_denominator_is_all_alerts():
    rows = [
        patient_contribution(True, True, 10, 12, 4),
        patient_contribution(False, False, None, None, 2),
    ]
    result = aggregate_contributions(rows)
    assert result["alerts"] == 2
    assert result["ppv_lower"] == result["ppv_upper"] == 0.5


def test_missing_denominator_status():
    result = aggregate_contributions([patient_contribution(False, False, None, None, float("inf"))])
    assert result["status"] == "MISSING_BOTH_DENOMINATORS"
    assert result["ppv_lower"] is None
