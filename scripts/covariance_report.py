import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import backend

backend.orbplot = lambda *args, **kwargs: None
EXAMPLE = "static/examples/GL765_Test1.inp"

_real_least_squares = backend.least_squares
_last_result = None


def _capture_least_squares(*args, **kwargs):
    global _last_result
    _last_result = _real_least_squares(*args, **kwargs)
    return _last_result


backend.least_squares = _capture_least_squares


def report(label):
    print(f"\n=== {label} ===")
    names = list(backend.orb.elname)
    for name, value, error, fixed in zip(names, backend.orb.el, backend.orb.elerr, backend.orb.fixel):
        print(f"{name:>4s}  {value:.8f}  {error:.8f}  fixed={int(fixed == 0)}")
    n = 2 * backend.orb.obj["npos"] + backend.orb.obj["nrv1"] + backend.orb.obj["nrv2"]
    nfree = int(sum(backend.orb.fixel > 0))
    dof = n - nfree
    print(f"chi2={backend.orb.obj['chi2']:.8f}")
    print(f"dof={dof}")
    print(f"reduced_chi2={backend.orb.obj['chi2']/dof:.8f}")


def report_combined_derived():
    cov = np.linalg.inv(_last_result.jac.T @ _last_result.jac)
    rng = np.random.default_rng(20260811)
    draws = rng.multivariate_normal(backend.orb.el[:10], cov, size=200000)
    # The local Gaussian approximation is very tight for this fit; retain only
    # physically valid samples in case a tail crosses a Keplerian boundary.
    draws = draws[(draws[:, 0] > 0) & (draws[:, 2] >= 0) & (draws[:, 2] < 1)]

    P = draws[:, 0]
    e = draws[:, 2]
    a = draws[:, 3]
    inc = np.deg2rad(draws[:, 6])
    K1 = draws[:, 7]
    K2 = draws[:, 8]

    P_sec = P * backend.SEC_PER_YEAR
    ksum = (K1 + K2) * 1000.0
    msini_kg = (ksum**3 * P_sec * (1 - e**2)**1.5) / (2 * np.pi * backend.G_NEWTON)
    msini = msini_kg / backend.M_SUN_KG
    mtotal = msini / np.sin(inc)**3
    q = K1 / K2
    m1 = mtotal / (1 + q)
    m2 = q * m1

    distance = (mtotal * P**2 / a**3)**(1 / 3)
    parallax = 1000.0 / distance

    print("\n=== COMBINED_DERIVED_MC ===")
    for name, values in (
        ("Msin3i", msini),
        ("Mtotal", mtotal),
        ("M1", m1),
        ("M2", m2),
        ("parallax_mas", parallax),
        ("distance_pc", distance),
    ):
        print(f"{name:>14s} mean={np.mean(values):.8f} std={np.std(values, ddof=1):.8f}")
    print(f"mc_samples={len(draws)}")


backend.readinp(EXAMPLE)
backend.fitorb()
report("COMBINED")
report_combined_derived()

backend.readinp(EXAMPLE)
backend.restrict_points(rv1_keep=[], rv2_keep=[])
for name in ("K1", "K2", "V0"):
    backend.orb.fixel[list(backend.orb.elname).index(name)] = 0
backend.fitorb()
report("VISUAL_ONLY")

backend.readinp(EXAMPLE)
backend.restrict_points(pos_keep=[])
for name in ("a", "W", "i"):
    backend.orb.fixel[list(backend.orb.elname).index(name)] = 0
backend.fitorb()
report("SPECTROSCOPIC_ONLY")
