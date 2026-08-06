# PySVOrbit

**A Python package for combined spectroscopic and visual orbit fitting for
binary stars.**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21830266.svg)](https://doi.org/10.5281/zenodo.21830266)
[![Tests](https://github.com/mtalafha90/PySVOrbit/actions/workflows/tests.yml/badge.svg)](https://github.com/mtalafha90/PySVOrbit/actions/workflows/tests.yml)
[![Licence: MIT](https://img.shields.io/badge/licence-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)

PySVOrbit determines the Keplerian orbits of binary stars from relative
astrometry and radial-velocity measurements. The two kinds of observation
may be fitted jointly or on their own, by weighted non-linear least
squares, and the dynamical masses that follow from the solution are derived
along with their uncertainties.

An interactive version runs in the browser at
**<https://pysvorbit.onrender.com/>** — no installation required.

## What it does

- **Three fitting modes**: visual-only, spectroscopic-only, or a combined
  solution using both data types at once.
- **Weighted least squares**: each observation contributes according to its
  own measurement error, and the parameter uncertainties come from the
  resulting covariance matrix.
- **Derived masses**: total and individual dynamical masses, with parallax,
  distance and physical semi-major axis all linked and solvable from any
  one known quantity.
- **Multiple instruments**: radial velocities from different spectrographs
  are combined by fitting a zero-point offset per instrument alongside the
  orbital elements.
- **Interactive refitting**: individual observations can be excluded and
  the fit re-run without editing the input file.
- **Publication outputs**: orbit, radial-velocity and residual plots, a
  results table, and a PDF report.

## Installation

Python 3.10 or newer is required.

```bash
git clone https://github.com/mtalafha90/PySVOrbit.git
cd PySVOrbit
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

### Web interface

```bash
python app.py
```

Then open <http://127.0.0.1:8000> and upload an input file. Two examples
are bundled in `static/examples/`:

| File | Description |
| --- | --- |
| `example_visual_binary.inp` | Synthetic visual binary with known elements |
| `GL765_Test1.inp` | Combined visual and spectroscopic data for GL 765.2 |

### Scripted use

The fitting engine can be driven directly, which is the route to take for
batch work or for reproducing a published solution:

```python
import backend

backend.readinp("static/examples/GL765_Test1.inp")
backend.fitorb()

for name, value, error in zip(backend.orb.elname,
                              backend.orb.el,
                              backend.orb.elerr):
    print(f"{name:>6} = {value:12.6f} +/- {error:.6f}")

backend.orbsave()   # writes <name>_output.csv to the working directory
```

Note that `backend` holds the current system in a single module-level
object, `backend.orb`. One fit is therefore in progress at any one time:
call `readinp()` again to start a new system, and do not run fits
concurrently in the same process.

## Input formats

Two formats are accepted, both documented by the bundled examples.

**`.inp`** — the classical fixed-format file. Orbital elements may be held
fixed during the fit by prefixing the element name with an asterisk, for
example `*e 0.0` to force a circular orbit.

**CSV** — a more flexible variant in which each radial-velocity row may
carry a source label in its final column, identifying the instrument it
came from.

In both formats, `RA` and `Dec` accept either decimal values (hours and
degrees) or sexagesimal, separated by spaces or colons — `12 34 33.1` and
`-45:30:00` are both understood.

## Notes on orbital fitting

- **Excluding individual data points**: the results page lists every
  position/RV1/RV2 observation with a checkbox (rows more than 3× the
  fit's typical residual are highlighted). Uncheck any points and click
  "Refit with selected points" to re-run the fit on just the remaining
  ones — no need to edit and re-upload the file. Excluded points stay
  visible (struck through) so they can be brought back later. If
  excluding points removes every observation from one instrument in a
  multi-instrument CSV, that instrument's `dV0` zero-point offset
  parameter is automatically dropped too.
- **Solving derived quantities from any known value**: parallax,
  distance, physical semi-major axis (AU), total mass, M1, and M2 are all
  linked by the same physics (parallax ↔ distance ↔ a(AU) ↔ M_total ↔
  M1/M2). By default the parallax from the uploaded file is the known
  starting point and everything else is derived from it, but the results
  page has a "Know one of these instead?" panel — enter any *one* of
  these quantities (e.g. a mass from spectral typing) and the rest are
  solved from it, without re-running the orbital fit. This only affects
  the derived quantities; it never changes the fitted orbital elements
  themselves.
- **Period units**: the orbital period `P` may be supplied in days or in
  years; values above 200 are treated as days when computing physical
  masses. Everywhere else (the fit itself, phase-folding) `P` is used in
  whatever unit the epoch of periastron `T` is given in, so keep `P` and
  `T` in a consistent day/year convention within one input file.
- **Combining RVs from multiple instruments**: for CSV uploads, each RV
  point carries a source/instrument label (the last column). If more than
  one distinct label is present, the fit automatically adds a `dV0_<source>`
  zero-point offset parameter per extra instrument (relative to one
  reference instrument), fitted simultaneously with the orbital elements.
  RV plots and the results table show these offsets; `.inp` uploads have no
  source labels, so this doesn't apply to them.

## Tests

```bash
pip install -r requirements.txt pytest
pytest tests/ -v
```

The suite covers the Kepler solver, checked against its own defining
equation; the ephemeris, checked against cases with analytic answers; input
parsing; the derived masses, checked against Kepler's third law; the
covariance estimate; and end-to-end fits of both bundled examples.

Several are regression tests for faults found during review — the `*`
prefix that holds an element fixed, the error-weighted Jacobian used for
the covariance, and the finite-difference step for elements whose value is
zero — so those faults cannot silently return. The suite runs on every push
and pull request through GitHub Actions against Python 3.10, 3.11 and 3.12.

## Citation

If you use PySVOrbit in work that leads to a publication, please cite
**both** the archived software release and the accompanying paper.

The archived release is the citation that pins the exact version of the
code you ran; the repository continues to develop, so it should not be
cited in its place.

> Al-Wardat, M. A., Talafha, M. H., & Alshamsi, S. N. (2026).
> *PySVOrbit: A Python Package for Combined Spectroscopic and Visual Orbit
> Fitting for Binary Stars* (version 1.0.0). Zenodo.
> <https://doi.org/10.5281/zenodo.21830266>

```bibtex
@software{alwardat_pysvorbit_software_2026,
  author    = {Al-Wardat, Mashhoor A. and
               Talafha, Mohammed H. and
               Alshamsi, Shaikha N.},
  title     = {{PySVOrbit}: A Python Package for Combined Spectroscopic
               and Visual Orbit Fitting for Binary Stars},
  version   = {1.0.0},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.21830266},
  url       = {https://doi.org/10.5281/zenodo.21830266}
}

@article{alwardat_pysvorbit_paper_2026,
  author  = {Al-Wardat, Mashhoor A. and
             Talafha, Mohammed H. and
             Alshamsi, Shaikha N.},
  title   = {{PySVOrbit}: A Python Package for Combined Spectroscopic
             and Visual Orbit Fitting for Binary Stars},
  journal = {Astronomy \& Computing},
  year    = {2026},
  note    = {Submitted}
}
```

The file [`CITATION.cff`](CITATION.cff) carries the same information in
machine-readable form. GitHub reads it to populate the **Cite this
repository** panel in the sidebar, which will export the citation in
BibTeX or APA on request.

## Deployment

`render.yaml` configures the repository as a Render Blueprint: it builds
with `pip install -r requirements.txt && kaleido_get_chrome`, serves with
`gunicorn -w 2 -b 0.0.0.0:$PORT app:app`, and generates a
`FLASK_SECRET_KEY` automatically. In the Render dashboard choose
**New > Blueprint** and point it at this repository.

To serve it behind your own reverse proxy instead:

```bash
gunicorn -w 2 -b 127.0.0.1:8000 app:app
```

Uploaded orbits, plots, CSVs and PDF reports are written to a local `runs/`
directory. On Render's free and standard instances this disk is ephemeral —
wiped on redeploy or restart, and not shared between instances — so it
suits interactive single-instance use, but is not a permanent results
store.

### PDF report plots

The PDF report renders the same Plotly figures shown on the results page,
using [Kaleido](https://github.com/plotly/Kaleido), so the two match
exactly. This needs a headless Chrome, which Kaleido manages itself but
does not fetch on demand; the Render build step runs `kaleido_get_chrome`
once at deploy time. Run that command yourself once if you are working
locally and want matching plots. If it is skipped, or Chrome is
unavailable, the report falls back to Matplotlib plots rather than failing.
Set `KALEIDO_CHROME_PATH` to point at an existing Chrome or Chromium
executable instead of Kaleido's managed copy.

## Licence

Released under the MIT Licence — see [`LICENSE`](LICENSE) for the full
text. You may use, modify and redistribute the software, including in
commercial work, provided the copyright notice is retained.
