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
