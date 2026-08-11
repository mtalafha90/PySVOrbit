import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import backend

backend.orbplot = lambda *args, **kwargs: None
EXAMPLE = "static/examples/GL765_Test1.inp"


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


backend.readinp(EXAMPLE)
backend.fitorb()
report("COMBINED")

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
