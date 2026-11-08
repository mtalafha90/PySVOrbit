"""Regression tests for real-world INP failure modes (no student data committed)."""
from pathlib import Path
import numpy as np
import pytest

import matplotlib
matplotlib.use("Agg")

import backend
import app


def _input(tmp_path, name="target.inp", *, pos_count=10, rv_rows=(), extra="", p=10.0, t=2000.0):
    lines = [
        "Object: Test Binary", "R.A.: 12 30 00", "Dec: -24 30 00",
        extra, f"P {p}", f"T {t}", "e 0.3", "a 0.2",
        "W 30", "w 40", "i 50",
    ]
    for k in range(pos_count):
        lines.append(f"{2000 + k * 0.5:.5f} 100.0 0.20 0.01 I1 TEST")
    lines.extend(rv_rows)
    file = tmp_path / name
    file.write_text("\n".join(lines) + "\n")
    return file


def test_reader_resets_metadata_and_understands_whitespace_coordinates(tmp_path):
    first = _input(tmp_path, "one.inp", extra="Parallax: 25.0")
    second = _input(tmp_path, "two.inp")
    backend.readinp(str(first))
    assert backend.orb.obj["parallax"] == pytest.approx(25.0)
    assert backend.orb.obj["radeg"] == pytest.approx(187.5)
    assert backend.orb.obj["dedeg"] == pytest.approx(-24.5)
    backend.readinp(str(second))
    assert backend.orb.obj["parallax"] == 0.0
    assert backend.orb.obj["name"] == "Test Binary"


def test_visual_only_file_automatically_holds_velocity_elements(tmp_path):
    file = _input(tmp_path, pos_count=12, extra="Period_Unit: years", p=421.98, t=1969.97)
    backend.readinp(str(file))
    assert backend.orb.obj["period_unit"] == "years"
    for nm in ("K1", "K2", "V0"):
        assert backend.orb.fixel[backend.orb.elname.index(nm)] == 0
    assert backend._period_in_years(421.98, T=1969.97) == 421.98
    assert backend.calculate_total_mass(421.98, 0.5, 50.0, T=1969.97) == pytest.approx(
        10.0 ** 3 / 421.98 ** 2
    )


def test_file_marked_fixed_parameter_is_not_accidentally_free(tmp_path):
    file = _input(tmp_path)
    file.write_text(file.read_text().replace("P 10.0", "*P 10.0"))
    backend.readinp(str(file))
    assert backend.orb.fixel[backend.orb.elname.index("P")] == 0


def test_sb1_has_no_unconstrained_k2(tmp_path):
    file = _input(tmp_path, rv_rows=["2001.0 4.2 0.2 Va COR"])
    backend.readinp(str(file))
    assert backend.orb.fixel[backend.orb.elname.index("K1")] == 1
    assert backend.orb.fixel[backend.orb.elname.index("K2")] == 0
    assert backend.orb.rv1_source == ["COR"]


def test_inp_preserves_rv_sources_and_exposes_offsets(tmp_path):
    file = _input(tmp_path, rv_rows=[
        "2000.0 4.2 0.2 Va COR",
        "2001.0 5.2 0.2 Va SBX",
        "2001.5 6.2 0.2 Vb COR",
    ])
    backend.readinp(str(file))
    assert backend.orb.rv1_source == ["COR", "SBX"]
    assert backend.orb.rv2_source == ["COR"]
    assert "SBX" in backend.orb.rv_offset_index


def test_duplicate_rv_rows_are_rejected_before_fit(tmp_path):
    file = _input(tmp_path, rv_rows=[
        "2001.0 4.2 0.2 Va SBX",
        "2001.0 4.2 0.2 Va SBX",
    ])
    backend.readinp(str(file))
    with pytest.raises(ValueError, match="exactly duplicated"):
        backend.fitorb()


def test_incomplete_marker_refuses_fit(tmp_path):
    file = _input(tmp_path, name="HIP_123_INCOMPLETE.inp")
    backend.readinp(str(file))
    with pytest.raises(ValueError, match="INCOMPLETE_DO_NOT_FIT"):
        backend.fitorb()


def test_invalid_position_uncertainty_is_refused(tmp_path):
    file = _input(tmp_path)
    text = file.read_text().replace("0.20 0.01 I1", "0.20 0 I1")
    file.write_text(text)
    backend.readinp(str(file))
    with pytest.raises(ValueError, match="positive finite uncertainty"):
        backend.fitorb()


def test_julian_date_is_not_mistaken_for_mjd():
    arr = np.array([[2451544.5, 1.0, 1.0]], dtype=float)
    backend.correct(arr, 2000.0)
    assert abs(arr[0, 0] - 2000.0) < 0.03
    arr2 = np.array([[2000.0, 1.0, 1.0]], dtype=float)
    backend.correct(arr2, 2451544.5)
    assert abs(arr2[0, 0] - 2451544.5) < 10


def test_derived_mass_uses_long_period_in_years():
    elements = [
        {"name": "P", "value": 421.98},
        {"name": "T", "value": 1969.97},
        {"name": "a", "value": 0.5},
        {"name": "e", "value": 0.1},
        {"name": "i", "value": 60},
    ]
    derived = app.compute_derived(elements, {"parallax": 50.0, "period_unit": "years"})
    assert derived["P_yr"] == pytest.approx(421.98)
    assert derived["M_total_Msun"] == pytest.approx(1000 / 421.98**2)


def test_explicit_day_period_is_converted():
    assert backend._period_in_years(365.25, T=2450000.0, unit="days") == pytest.approx(1.0)
