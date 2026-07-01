"""ACMG code-string -> points parsing (plan Step 0.2 / Gate G1)."""
from rules.acmg_codes import code_points, contributions_from_codes
from rules.point_engine import classify


def test_defaults_and_modifiers():
    assert code_points("PVS1") == (8.0, False)
    assert code_points("PS3") == (4.0, False)
    assert code_points("PM2") == (2.0, False)
    assert code_points("PP4") == (1.0, False)
    assert code_points("PVS1_Strong") == (4.0, False)
    assert code_points("PM3_Very Strong") == (8.0, False)   # eRepo uses a space
    assert code_points("PP4_Moderate") == (2.0, False)
    assert code_points("PM2_Supporting") == (1.0, False)
    assert code_points("BS1") == (-4.0, False)
    assert code_points("BS3_Supporting") == (-1.0, False)
    assert code_points("BP4") == (-1.0, False)
    assert code_points("") == (0.0, False)
    assert code_points("FOO") == (0.0, False)


def test_ba1_flagged_not_summed():
    pts, is_ba1 = code_points("BA1")
    assert is_ba1 and pts == 0.0


def test_end_to_end_pathogenic_call():
    # PAH c.1A>G example from eRepo: PS3, PM3, PP4_Moderate, PM2 -> 4+2+2+2 = 10 -> Pathogenic
    contribs, ba1 = contributions_from_codes(["PS3", "PM3", "PP4_Moderate", "PM2"])
    from core.schemas import PointsLedger, Variant
    led = PointsLedger(Variant("12", 1, "A", "G", gene="PAH"), contributions=contribs)
    assert led.points == 10.0
    assert not ba1
    assert classify(led).name == "P"
