"""Tests for the PySVOrbit fitting engine.

Several of these lock in behaviour that was found to be wrong during
review: the '*' fixed-element convention in .inp files, the weighting of
the Jacobian used for the covariance, and the finite-difference step for
elements sitting at zero. They are regression tests as much as unit tests.
"""
import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")

import backend  # noqa: E402

EXAMPLES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "static", "examples")
TEACHING = os.path.join(EXAMPLES, "example_visual_binary.inp")
GL765 = os.path.join(EXAMPLES, "GL765_Test1.inp")


@pytest.fixture(autouse=True)
def _quiet(monkeypatch, capsys):
    """fitorb() draws diagnostic figures; skip them so tests stay fast."""
    monkeypatch.setattr(backend, "orbplot", lambda *a, **k: None)


# --------------------------------------------------------------------------
# Kepler solver and ephemeris
# --------------------------------------------------------------------------

@pytest.mark.parametrize("e", [0.0, 0.1, 0.35, 0.6, 0.9])
@pytest.mark.parametrize("frac", [0.0, 0.13, 0.5, 0.77, 0.99])
def test_kepler_solver_satisfies_its_own_equation(e, frac):
    """E returned by _solve_kepler must satisfy M = E - e sin E."""
    M = 2 * math.pi * frac
    E = backend._solve_kepler(M, e)
    assert abs((E - e * math.sin(E)) - M) < 1e-4


def test_kepler_solver_terminates_on_unphysical_eccentricity():
    """A trial e >= 1 must not hang the solver (it is capped by max_iter)."""
    E = backend._solve_kepler(1.0, 1.2)
    assert np.isfinite(E)


@pytest.mark.parametrize("diff,expected", [
    (0.0, 0.0), (10.0, 10.0), (-10.0, -10.0),
    (350.0, -10.0), (-350.0, 10.0), (359.0, -1.0), (181.0, -179.0),
])
def test_wrap_deg_diff(diff, expected):
    assert backend._wrap_deg_diff(diff) == pytest.approx(expected, abs=1e-9)


def test_circular_face_on_orbit_has_constant_separation():
    """e=0, i=0 => the projected separation is the semi-major axis."""
    el = np.array([10.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    t = np.linspace(0, 10, 25)
    res = backend.eph(el, t, rho=True)
    assert np.allclose(res[:, 1], 0.5, atol=1e-6)


def test_radial_velocity_amplitude_of_a_circular_orbit():
    """For e=0 the primary RV should span V0 +/- K1."""
    K1, V0 = 8.0, -3.0
    el = np.array([10.0, 0.0, 0.0, 0.5, 0.0, 0.0, 90.0, K1, 5.0, V0])
    rv = backend.eph(el, np.linspace(0, 10, 400), rv=True)[:, 0]
    assert rv.max() == pytest.approx(V0 + K1, abs=0.05)
    assert rv.min() == pytest.approx(V0 - K1, abs=0.05)


# --------------------------------------------------------------------------
# Input parsing
# --------------------------------------------------------------------------

@pytest.mark.parametrize("text,expected", [
    ("12.5762", 12.5762),
    ("12 34 33.1", 12 + 34 / 60 + 33.1 / 3600),
    ("-45:30:00", -(45 + 30 / 60)),
    ("+45 30", 45.5),
])
def test_getcoord_accepts_decimal_and_sexagesimal(text, expected):
    assert backend.getcoord(text) == pytest.approx(expected, abs=1e-9)


def test_star_prefix_holds_an_element_fixed(tmp_path):
    """Regression: '*K1 0.0' must set the value AND clear the fit flag.

    Previously the '*' was left on the name token, so the element failed the
    name lookup, kept its default of zero and stayed free -- an
    astrometry-only file then fitted K1/K2/V0 against no velocities at all.
    """
    p = tmp_path / "fixed.inp"
    p.write_text(
        "Object: Test\nParallax: 20.0\n"
        "P 10.0\nT 2000.0\ne 0.3\na 0.2\nW 100.0\nw 50.0\ni 80.0\n"
        "*K1 1.5\n* K2 2.5\n*V0 -7.0\n"
        "2000.0 100.0 0.2 0.01 I1\n2001.0 110.0 0.2 0.01 I1\n"
    )
    backend.readinp(str(p))
    names = list(backend.orb.elname)
    for nm, val in (("K1", 1.5), ("K2", 2.5), ("V0", -7.0)):
        k = names.index(nm)
        assert backend.orb.el[k] == pytest.approx(val), f"{nm} value not read"
        assert backend.orb.fixel[k] == 0, f"{nm} should be held fixed"
    assert backend.orb.fixel[names.index("P")] == 1


def test_metadata_key_variants_are_recognised(tmp_path):
    p = tmp_path / "meta.inp"
    p.write_text("Object: Some Star\nR.A.:  12 00 00\nDec: -30 00 00\nPlx: 25.5\nP 5.0\n")
    backend.readinp(str(p))
    assert backend.orb.obj["name"] == "Some Star"
    assert backend.orb.obj["parallax"] == pytest.approx(25.5)
    assert backend.orb.obj["radeg"] == pytest.approx(180.0)
    assert backend.orb.obj["dedeg"] == pytest.approx(-30.0)


def test_parallax_is_not_swallowed_by_the_ra_key(tmp_path):
    """Regression: 'Parallax' contains 'ra', which once matched the RA branch."""
    p = tmp_path / "plx.csv"
    p.write_text("Object,Star\nRA,12.0\nParallax,31.0\nP,10.0\n")
    backend.readcsv_custom(str(p))
    assert backend.orb.obj["parallax"] == pytest.approx(31.0)
    assert backend.orb.obj["radeg"] == pytest.approx(180.0)


# --------------------------------------------------------------------------
# Derived physical quantities
# --------------------------------------------------------------------------

def test_spectroscopic_mass_matches_the_classical_coefficient():
    """M(1+2)sin^3 i must agree with 1.0385e-7 (1-e^2)^{3/2}(K1+K2)^3 P_days."""
    P, e, i, K1, K2 = 11.7295, 0.2489, 81.909, 7.9420, 7.6956
    m12, _, _ = backend.calculate_spectroscopic_masses(P, e, i, K1, K2)
    classical = 1.0385e-7 * (1 - e ** 2) ** 1.5 * (K1 + K2) ** 3 * (P * 365.25)
    assert m12 == pytest.approx(classical, rel=0.01)


def test_spectroscopic_masses_sum_and_ratio():
    P, e, i, K1, K2 = 11.7295, 0.2489, 81.909, 7.9420, 7.6956
    m12, M1, M2 = backend.calculate_spectroscopic_masses(P, e, i, K1, K2)
    assert (M1 + M2) == pytest.approx(m12 / math.sin(math.radians(i)) ** 3, rel=1e-9)
    assert (M2 / M1) == pytest.approx(K1 / K2, rel=1e-9)


def test_total_mass_follows_keplers_third_law():
    P_yr, a_arcsec, plx = 10.0, 0.2, 25.0          # d = 40 pc -> a = 8 AU
    M = backend.calculate_total_mass(P_yr, a_arcsec, plx)
    assert M == pytest.approx(8.0 ** 3 / 10.0 ** 2, rel=1e-9)


def test_total_mass_is_zero_without_a_parallax():
    assert backend.calculate_total_mass(10.0, 0.2, 0.0) == 0.0


@pytest.mark.parametrize("P,expected", [(11.7, 11.7), (4284.0, 4284.0 / 365.25)])
def test_period_unit_heuristic(P, expected):
    assert backend._period_in_years(P) == pytest.approx(expected)


# --------------------------------------------------------------------------
# Fitting: end-to-end behaviour and the covariance
# --------------------------------------------------------------------------

def test_recovers_the_elements_of_the_synthetic_example():
    """The bundled teaching file is generated from known elements.

    The data carry simulated noise, so recovery is judged against the fit's
    own uncertainties rather than against arbitrary absolute tolerances:
    every element must land within 3 sigma of the value used to build it.
    """
    backend.readinp(TEACHING)
    truth = backend.orb.el.copy()
    backend.orb.el[0] *= 1.02          # start slightly off
    backend.orb.el[2] += 0.02
    backend.fitorb()
    names = list(backend.orb.elname)
    for nm in ("P", "T", "e", "a", "W", "w", "i", "K1", "K2", "V0"):
        k = names.index(nm)
        sigma = backend.orb.elerr[k]
        assert sigma > 0, f"{nm} has no uncertainty"
        pull = abs(backend.orb.el[k] - truth[k]) / sigma
        assert pull < 3.0, f"{nm} recovered {pull:.1f} sigma from the input value"


def test_zero_valued_free_element_does_not_poison_the_jacobian():
    """Regression: a step of exactly zero used to divide by zero -> NaN."""
    backend.readinp(TEACHING)
    backend.orb.el[7] = 0.0            # K1 = 0 with K1 still free
    n = 2 * backend.orb.obj["npos"] + backend.orb.obj["nrv1"] + backend.orb.obj["nrv2"]
    par = backend.orb.el[np.where(backend.orb.fixel > 0)[0]]
    rows = np.array([backend.alleph(par, i) for i in range(min(n, 20))])
    assert np.isfinite(rows).all()


@pytest.mark.skipif(not os.path.exists(GL765), reason="GL765 example not present")
def test_covariance_uses_the_error_weighted_jacobian():
    """Regression: an unweighted J^T J inflated sigma(a) by a factor of ~29.

    With heterogeneous data (arcsec astrometry alongside km/s velocities) the
    unweighted normal matrix gave sigma(a) = 0.377 on a value of 0.215. The
    weighted form gives ~0.013, consistent with a delta-chi^2 = 1 scan.
    """
    backend.readinp(GL765)
    backend.fitorb()
    names = list(backend.orb.elname)
    a, sa = backend.orb.el[names.index("a")], backend.orb.elerr[names.index("a")]
    assert sa / a < 0.15, f"sigma(a)/a = {sa/a:.2f}; covariance is probably unweighted"
    assert sa == pytest.approx(0.0131, abs=0.004)


@pytest.mark.skipif(not os.path.exists(GL765), reason="GL765 example not present")
def test_gl765_reproduces_the_published_solution():
    """End-to-end: the combined fit should land on the paper's solution."""
    backend.readinp(GL765)
    backend.fitorb()
    names = list(backend.orb.elname)
    expected = {"P": (11.7295, 0.01), "e": (0.2489, 0.005),
                "i": (81.909, 0.5), "K1": (7.942, 0.05), "K2": (7.696, 0.05)}
    for nm, (val, tol) in expected.items():
        assert backend.orb.el[names.index(nm)] == pytest.approx(val, abs=tol), nm
    dof = 2 * backend.orb.obj["npos"] + backend.orb.obj["nrv1"] + backend.orb.obj["nrv2"] - 10
    assert backend.orb.obj["chi2"] / dof == pytest.approx(1.05, abs=0.1)


@pytest.mark.skipif(not os.path.exists(GL765), reason="GL765 example not present")
def test_restrict_points_drops_observations():
    backend.readinp(GL765)
    n0 = backend.orb.obj["npos"]
    backend.restrict_points(pos_keep=list(range(n0 - 2)))
    assert backend.orb.obj["npos"] == n0 - 2
    assert backend.orb.pos.shape[0] == n0 - 2


@pytest.mark.skipif(not os.path.exists(GL765), reason="GL765 example not present")
def test_visual_only_mode_holds_spectroscopic_elements_fixed():
    backend.readinp(GL765)
    backend.restrict_points(rv1_keep=[], rv2_keep=[])
    names = list(backend.orb.elname)
    for nm in ("K1", "K2", "V0"):
        backend.orb.fixel[names.index(nm)] = 0
    before = backend.orb.el.copy()
    backend.fitorb()
    for nm in ("K1", "K2", "V0"):
        k = names.index(nm)
        assert backend.orb.el[k] == pytest.approx(before[k]), f"{nm} moved despite being fixed"
