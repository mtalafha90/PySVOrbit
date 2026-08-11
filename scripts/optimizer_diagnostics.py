"""Reproduce optimizer-effort diagnostics for the bundled GJ 765.2 fit.

This script distinguishes SciPy's ``OptimizeResult.nfev`` from the actual
number of calls made to the residual-vector callable when a numerical
Jacobian is used.  It also instruments the PySVOrbit ephemeris function so
that any ephemeris-call count reported in the manuscript is measured rather
than inferred from solver metadata.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import backend  # noqa: E402


EXAMPLE = ROOT / "static" / "examples" / "GL765_Test1.inp"


def run_diagnostics(example: os.PathLike[str] | str = EXAMPLE) -> dict[str, int | None | str]:
    counts = {"residual_vector_calls": 0, "ephemeris_calls": 0}
    captured = {}

    real_least_squares = backend.least_squares
    real_eph = backend.eph
    real_orbplot = backend.orbplot

    def counted_eph(*args, **kwargs):
        counts["ephemeris_calls"] += 1
        return real_eph(*args, **kwargs)

    def counted_least_squares(fun, *args, **kwargs):
        def counted_fun(*fargs, **fkwargs):
            counts["residual_vector_calls"] += 1
            return fun(*fargs, **fkwargs)

        result = real_least_squares(counted_fun, *args, **kwargs)
        captured["result"] = result
        return result

    try:
        backend.eph = counted_eph
        backend.least_squares = counted_least_squares
        # Plotting performs additional ephemeris evaluations that are not part
        # of the fit itself, so suppress it for a clean optimizer benchmark.
        backend.orbplot = lambda *args, **kwargs: None
        backend.readinp(str(example))
        backend.fitorb()
    finally:
        backend.eph = real_eph
        backend.least_squares = real_least_squares
        backend.orbplot = real_orbplot

    result = captured["result"]
    return {
        "scipy_nfev": int(result.nfev),
        "scipy_njev": None if result.njev is None else int(result.njev),
        "residual_vector_calls": counts["residual_vector_calls"],
        "ephemeris_calls": counts["ephemeris_calls"],
        "termination": str(result.message),
    }


def main() -> None:
    report = run_diagnostics()
    print("=== OPTIMIZER EFFORT DIAGNOSTICS ===")
    print(f"SciPy nfev: {report['scipy_nfev']}")
    print(f"SciPy njev: {report['scipy_njev']}")
    print(f"Instrumented residual-vector calls: {report['residual_vector_calls']}")
    print(f"Instrumented ephemeris calls during fit: {report['ephemeris_calls']}")
    print(f"Termination: {report['termination']}")


if __name__ == "__main__":
    main()
