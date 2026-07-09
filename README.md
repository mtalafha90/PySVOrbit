# Orbital Fitting Normal Web App (Flask)

This is a normal server-rendered web app (Flask + HTML) built on your original fitting code.

## Run locally
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python app.py
```
Open: http://127.0.0.1:8000

## Production (university server)
```bash
source .venv/bin/activate
gunicorn -w 2 -b 127.0.0.1:8000 app:app
```
Reverse proxy with Nginx/Apache to 127.0.0.1:8000.

## Notes on orbital fitting

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
- **Coordinates**: `RA`/`Dec` accept either plain decimal (hours/degrees)
  or sexagesimal, space- or colon-separated (e.g. `12 34 33.1` or
  `-45:30:00`).
