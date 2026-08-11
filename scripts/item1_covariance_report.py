import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import backend

EXAMPLE = os.path.join(ROOT, "static", "examples", "GL765_Test1.inp")

plt.show = lambda *args, **kwargs: None
backend.plt.show = lambda *args, **kwargs: None

captured = {}
real_least_squares = backend.least_squares

def capture_least_squares(*args, **kwargs):
    result = real_least_squares(*args, **kwargs)
    captured["result"] = result
    return result

backend.least_squares = capture_least_squares
backend.readinp(EXAMPLE)
backend.fitorb()

result = captured["result"]
Jr = result.jac
cov = np.linalg.inv(Jr.T @ Jr)
mu = backend.orb.el.copy()

print("=== COMBINED COVARIANCE ===")
print(f"P={mu[0]:.10f} sigmaP={np.sqrt(cov[0,0]):.10f}")
print(f"a={mu[3]:.10f} sigmaa={np.sqrt(cov[3,3]):.10f}")
print(f"cov_P_a={cov[0,3]:.12e}")
print(f"corr_P_a={cov[0,3]/np.sqrt(cov[0,0]*cov[3,3]):.8f}")

rng = np.random.default_rng(20260811)
samples = rng.multivariate_normal(mu, cov, size=300000)
P = samples[:, 0]
e = samples[:, 2]
a = samples[:, 3]
i = samples[:, 6]
K1 = samples[:, 7]
K2 = samples[:, 8]
plx = rng.normal(31.0, 0.5, size=len(samples))

# Astrometric total mass using the input/published parallax.
d_pc = 1000.0 / plx
a_au = a * d_pc
M_ast = a_au**3 / P**2
print("=== ASTROMETRIC ROUTE MC ===")
print(f"M_ast_mean={np.mean(M_ast):.8f} M_ast_std={np.std(M_ast, ddof=1):.8f}")

# Spectroscopic masses, q, orbital parallax and distance.
P_sec = P * backend.SEC_PER_YEAR
Ksum_ms = (K1 + K2) * 1000.0
Msin3 = (Ksum_ms**3 * P_sec * (1-e**2)**1.5) / (2*np.pi*backend.G_NEWTON) / backend.M_SUN_KG
sin3i = np.sin(np.deg2rad(i))**3
Mtot = Msin3 / sin3i
q = K1 / K2
M1 = Mtot / (1+q)
M2 = q*M1
orb_plx = 1000.0 * a / np.cbrt(Mtot * P**2)
distance = 1000.0 / orb_plx
print("=== SPECTROSCOPIC DERIVED MC ===")
for name, x in [("Msin3",Msin3),("Mtot",Mtot),("M1",M1),("M2",M2),("q",q),("orb_plx",orb_plx),("distance",distance)]:
    print(f"{name}_mean={np.mean(x):.8f} {name}_std={np.std(x, ddof=1):.8f}")

# Sigma separation between the two mass routes using the MC standard deviations.
delta = abs(np.mean(M_ast) - np.mean(Mtot))
sigma_delta = np.sqrt(np.var(M_ast, ddof=1) + np.var(Mtot, ddof=1))
print(f"mass_routes_delta_sigma={delta/sigma_delta:.8f}")
