"""Regression checks for invalid-e sqrt warnings and high-eccentricity orbits."""
import warnings

import numpy as np
import pytest

import backend


@pytest.mark.parametrize("eccentricity", [-1.01, -0.01, 1.0, 1.01, 2.0, np.nan, np.inf])
def test_ephemeris_rejects_invalid_eccentricity_without_numpy_warnings(eccentricity):
    elements = np.array([10., 2000., eccentricity, 0.2, 40., 70.,
                         80., 7., 8., -3.])
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        with pytest.raises(ValueError, match="eccentricity|finite"):
            backend.eph(elements, [2000., 2002.], rho=True)
        with pytest.raises(ValueError, match="eccentricity|finite"):
            backend.eph(elements, [2000., 2002.], rv=True)


@pytest.mark.parametrize("eccentricity", [0.0, 0.23, 0.9754, 0.999999])
def test_ephemeris_remains_finite_for_valid_elliptical_orbits(eccentricity):
    elements = np.array([10., 2000., eccentricity, 0.2, 40., 70.,
                         80., 7., 8., -3.])
    times = 2000. + 10 * np.array([0., 0.00001, 0.001, 0.02, 0.25, 0.5, 0.9])
    with np.errstate(invalid="raise", divide="raise", over="raise"):
        positions = backend.eph(elements, times, rho=True)
        velocities = backend.eph(elements, times, rv=True)
    assert np.isfinite(positions).all()
    assert np.isfinite(velocities).all()
    assert (positions[:, 1] >= 0).all()


@pytest.mark.parametrize("eccentricity", [0.0, 0.2, 0.65, 0.97])
def test_ephemeris_is_algebraically_equivalent_to_true_anomaly_form(eccentricity):
    """Check the numerical implementation against the older true-anomaly equations."""
    p, t0, a = 11.7, 2000., 0.23
    W, w, inc = np.radians([40., 70., 80.])
    k1, k2, systemic = 7., 8., -3.
    el = np.array([p, t0, eccentricity, a, np.degrees(W),
                   np.degrees(w), np.degrees(inc), k1, k2, systemic])
    times = t0 + p * np.array([0.01, 0.17, 0.42, 0.71, 0.99])
    AA = a * (np.cos(w)*np.cos(W) - np.sin(w)*np.sin(W)*np.cos(inc))
    BB = a * (np.cos(w)*np.sin(W) + np.sin(w)*np.cos(W)*np.cos(inc))
    FF = a * (-np.sin(w)*np.cos(W) - np.cos(w)*np.sin(W)*np.cos(inc))
    GG = a * (-np.sin(w)*np.sin(W) + np.cos(w)*np.cos(W)*np.cos(inc))
    positions = []
    velocities = []
    for epoch in times:
        mean = 2*np.pi * ((epoch-t0)/p % 1)
        E = backend._solve_kepler(mean, eccentricity)
        nu = 2*np.arctan(
            np.sqrt((1+eccentricity)/(1-eccentricity))*np.tan(E/2)
        )
        radius = (1-eccentricity**2)/(1+eccentricity*np.cos(nu))
        x = radius*np.cos(nu)
        y = radius*np.sin(nu)
        positions.append([AA*x + FF*y, BB*x + GG*y])
        factor = eccentricity*np.cos(w) + np.cos(nu+w)
        velocities.append([systemic+k1*factor, systemic-k2*factor])
    assert backend.eph(el, times) == pytest.approx(positions, abs=1e-10)
    assert backend.eph(el, times, rv=True) == pytest.approx(velocities, abs=1e-10)


def test_optimizer_invalid_trial_is_penalized_without_ephemeris_evaluation(monkeypatch):
    """Even an unconstrained LM trial may not evaluate sqrt(1-e**2) at e>1."""
    from pathlib import Path
    import matplotlib
    matplotlib.use("Agg")
    source = Path(__file__).resolve().parents[1] / "static/examples/example_visual_binary.inp"
    backend.readinp(str(source))
    original_fit = backend.least_squares
    original_plot = backend.orbplot
    attempted = []

    def probe_optimizer(fun, par, **kwargs):
        k = list(np.where(backend.orb.fixel > 0)[0]).index(2)
        invalid = par.copy()
        invalid[k] = 1.2
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            result = fun(invalid)
        assert np.isfinite(result).all()
        attempted.append(True)
        return original_fit(fun, par, **kwargs)

    try:
        monkeypatch.setattr(backend, "least_squares", probe_optimizer)
        monkeypatch.setattr(backend, "orbplot", lambda: None)
        backend.fitorb()
    finally:
        monkeypatch.setattr(backend, "least_squares", original_fit)
        monkeypatch.setattr(backend, "orbplot", original_plot)
    assert attempted
