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

## Deploy on Render

`render.yaml` at the repo root configures this as a Render Blueprint:
build with `pip install -r requirements.txt && kaleido_get_chrome` (the
second step downloads the headless Chrome build the PDF report uses to
render plots — see below), run with `gunicorn -w 2 -b 0.0.0.0:$PORT
app:app`, and auto-generate a `FLASK_SECRET_KEY`. In the Render dashboard:
**New > Blueprint**, point it at this repo, and deploy — no manual
configuration needed. To configure a Web Service by hand instead, use the
same build/start commands above.

Note: uploaded orbits, generated plots, CSVs and PDF reports are written
to a local `runs/` directory. Render's free/standard instance disk is
ephemeral (wiped on redeploy or restart) and not shared across multiple
instances, so download links from past runs won't survive either — fine
for a single-instance deployment used interactively, but not a permanent
results store.

## PDF report plots

The PDF report renders the same Plotly figures shown on the results page
(via [Kaleido](https://github.com/plotly/Kaleido)), so they match exactly
rather than falling back to a separately-styled Matplotlib version. This
needs a headless Chrome, which Kaleido downloads and manages itself, but
does not fetch on demand — `kaleido_get_chrome` in the Render build step
does this once at deploy time. Running locally, run `kaleido_get_chrome`
yourself once beforehand to get matching plots; if it's skipped, or
Chrome isn't available for any reason, the report falls back to the
Matplotlib plots instead of failing. To point at an existing
Chrome/Chromium install rather than Kaleido's managed one, set
`KALEIDO_CHROME_PATH` to its executable path.

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
