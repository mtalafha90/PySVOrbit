"""Regression tests for the statistical definition of fitted uncertainties."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")

import backend  # noqa: E402

EXAMPLE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "static", "examples", "GL765_Test1.inp")


def test_covariance_uses_absolute_measurement_errors_without_chi2_rescaling(monkeypatch):
    """Formal errors must be sqrt(diag((J_r^T J_r)^-1)), with no chi2_nu factor.

    ``least_squares`` receives residuals already divided by their quoted 1-sigma
    measurement errors, so its Jacobian is the Jacobian J_r of the standardized
    residual vector.  Reduced chi-square is a goodness-of-fit diagnostic only.
    """
    monkeypatch.setattr(backend, "orbplot", lambda *args, **kwargs: None)
    real_least_squares = backend.least_squares
    captured = {}

    def capture(*args, **kwargs):
        result = real_least_squares(*args, **kwargs)
        captured["result"] = result
        return result

    monkeypatch.setattr(backend, "least_squares", capture)
    backend.readinp(EXAMPLE)
    backend.fitorb()

    result = captured["result"]
    expected = np.sqrt(np.diag(np.linalg.inv(result.jac.T @ result.jac)))
    free = np.where(backend.orb.fixel > 0)[0]

    assert np.allclose(backend.orb.elerr[free], expected, rtol=1e-10, atol=1e-12)

    n = 2 * backend.orb.obj["npos"] + backend.orb.obj["nrv1"] + backend.orb.obj["nrv2"]
    reduced_chi2 = backend.orb.obj["chi2"] / (n - len(free))
    rescaled = expected * np.sqrt(reduced_chi2)
    assert not np.allclose(backend.orb.elerr[free], rescaled, rtol=1e-4, atol=1e-8)
