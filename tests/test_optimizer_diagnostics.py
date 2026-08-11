import os

import pytest

from scripts.optimizer_diagnostics import EXAMPLE, run_diagnostics


@pytest.mark.skipif(not os.path.exists(EXAMPLE), reason="GL765 example not present")
def test_optimizer_effort_is_measured_not_inferred():
    report = run_diagnostics()

    assert report["scipy_nfev"] > 0
    # For LM with a numerically estimated Jacobian, current SciPy reports
    # njev=None; the important regression is that we separately count the
    # actual residual-vector calls instead of deriving them from nfev.
    assert report["scipy_njev"] is None
    assert report["residual_vector_calls"] >= report["scipy_nfev"]
    assert report["ephemeris_calls"] >= report["residual_vector_calls"]
