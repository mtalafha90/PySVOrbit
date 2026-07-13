import os
import uuid
import shutil
import traceback
import io
from pathlib import Path
import math
import plotly.graph_objects as go
import json
from plotly.utils import PlotlyJSONEncoder
import numpy as np
from flask import Flask, render_template, request, redirect, url_for, send_from_directory, flash
from flask import send_file
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib.utils import ImageReader
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, Image as RLImage, PageBreak
from werkzeug.utils import secure_filename
from plotly.subplots import make_subplots

import matplotlib
matplotlib.use("Agg")

import backend  # your original code renamed to backend.py

APP_DIR = Path(__file__).resolve().parent
RUNS_DIR = APP_DIR / "runs"
RUNS_DIR.mkdir(exist_ok=True)

ALLOWED_EXTENSIONS = {".inp", ".csv", ".txt"}


def allowed_file(filename: str) -> bool:
    return Path(filename).suffix.lower() in ALLOWED_EXTENSIONS


def save_fig(fig, path: Path):
    fig.savefig(path, dpi=160, bbox_inches="tight")
    try:
        import matplotlib.pyplot as plt
        plt.close(fig)
    except Exception:
        pass


def default_stats():
    return {
        "name": "",
        "npos": 0,
        "nrv1": 0,
        "nrv2": 0,
        "parallax": 0.0,
        "chi2": 0.0,
        "chi2n_theta": 0.0,
        "chi2n_rho": 0.0,
        "chi2n_rv1": 0.0,
        "chi2n_rv2": 0.0,
        "rms_theta": 0.0,
        "rms_rho": 0.0,
        "rms_rv1": 0.0,
        "rms_rv2": 0.0,
    }


def _get_el(elements, name):
    """Return element value by name from elements list (dicts)."""
    for r in elements:
        if r["name"] == name:
            return float(r["value"])
    return None

# The physical quantities describing the system size/mass form a single
# chain, each link solvable in either direction:
#   parallax_mas <-> distance_pc <-> a_AU <-> M_total_Msun <-> (M1_Msun, M2_Msun)
# DERIVED_QUANTITY_LABELS names every quantity a user can supply as the
# known value, with the rest solved from it.
DERIVED_QUANTITY_LABELS = {
    "parallax_mas": "Parallax (mas)",
    "distance_pc": "Distance (pc)",
    "a_AU": "Physical semi-major axis a (AU)",
    "M_total_Msun": "Total mass M_total (Msun)",
    "M1_Msun": "Primary mass M1 (Msun)",
    "M2_Msun": "Secondary mass M2 (Msun)",
}


def compute_derived(elements, stats, override_key=None, override_value=None):
    """
    Compute derived quantities (parallax, distance, physical semi-major
    axis, total mass, M1, M2) by solving the chain

        parallax_mas <-> distance_pc <-> a_AU <-> M_total_Msun <-> (M1, M2)

    from whichever single quantity is known. By default that's the
    parallax read from the uploaded file (the original one-directional
    behaviour: parallax -> ... -> masses). Passing override_key (one of
    DERIVED_QUANTITY_LABELS) and override_value seeds the chain from that
    quantity instead - e.g. a known mass from spectral typing, solved
    backwards to the parallax/distance it implies.

    Returns a dict of the resolved values (only those the chain could
    reach) plus a "warnings" list explaining anything that couldn't be
    computed.
    """
    out = {}
    warn = []

    P = _get_el(elements, "P")
    a_arcsec = _get_el(elements, "a")
    i_deg = _get_el(elements, "i")
    e = _get_el(elements, "e")
    K1 = _get_el(elements, "K1")
    K2 = _get_el(elements, "K2")

    # Detect P unit (years vs days) - best-effort heuristic
    P_yr = None
    if P is not None:
        if P > 200:  # likely days
            P_yr = P / 365.25
            warn.append("Detected P looks like days; converted to years for mass calculations.")
        else:
            P_yr = P
    else:
        warn.append("P is missing; cannot compute masses.")

    # Mass ratio q = M2/M1 = K1/K2 (RV amplitude is inversely proportional
    # to mass), needed to split M_total into M1/M2 in either direction.
    q = None
    if K1 is not None and K2 is not None and K2 != 0:
        q = K1 / K2
    else:
        warn.append("K1/K2 not available; cannot relate individual masses to the total mass.")

    # Seed the chain: an explicit override, or (default) the file's parallax.
    v = {k: None for k in DERIVED_QUANTITY_LABELS}
    if override_key in v and override_value is not None:
        v[override_key] = float(override_value)
    else:
        plx_from_file = float(stats.get("parallax", 0.0) or 0.0) if isinstance(stats, dict) else 0.0
        if plx_from_file > 0:
            v["parallax_mas"] = plx_from_file
        else:
            warn.append("Parallax is missing/zero; cannot compute distance or physical masses "
                         "(enter a known value below to work out the rest instead).")

    # Relax the chain until nothing new is resolved. 4 passes is enough to
    # cover the longest possible path (parallax <-> ... <-> M1/M2, 4 hops)
    # regardless of which single quantity was seeded.
    for _ in range(4):
        if v["distance_pc"] is None and v["parallax_mas"]:
            v["distance_pc"] = 1000.0 / v["parallax_mas"]
        if v["parallax_mas"] is None and v["distance_pc"]:
            v["parallax_mas"] = 1000.0 / v["distance_pc"]

        if v["a_AU"] is None and v["distance_pc"] is not None and a_arcsec:
            v["a_AU"] = a_arcsec * v["distance_pc"]
        if v["distance_pc"] is None and v["a_AU"] is not None and a_arcsec:
            v["distance_pc"] = v["a_AU"] / a_arcsec

        if v["M_total_Msun"] is None and v["a_AU"] is not None and P_yr:
            v["M_total_Msun"] = v["a_AU"] ** 3 / P_yr ** 2
        if v["a_AU"] is None and v["M_total_Msun"] is not None and P_yr:
            v["a_AU"] = (v["M_total_Msun"] * P_yr ** 2) ** (1 / 3)

        if q is not None and q > 0:
            if v["M_total_Msun"] is None and v["M1_Msun"] is not None:
                v["M_total_Msun"] = v["M1_Msun"] * (1.0 + q)
            if v["M_total_Msun"] is None and v["M2_Msun"] is not None:
                v["M_total_Msun"] = v["M2_Msun"] * (1.0 + q) / q
            if v["M_total_Msun"] is not None:
                if v["M1_Msun"] is None:
                    v["M1_Msun"] = v["M_total_Msun"] / (1.0 + q)
                if v["M2_Msun"] is None:
                    v["M2_Msun"] = v["M_total_Msun"] * q / (1.0 + q)

    if not a_arcsec:
        warn.append("a is missing; cannot relate distance to the physical semi-major axis.")

    for key in DERIVED_QUANTITY_LABELS:
        if v[key] is not None:
            out[key] = v[key]
    if q is not None:
        out["q_M2_over_M1"] = q
    if P_yr is not None:
        out["P_yr"] = P_yr
    if a_arcsec is not None:
        out["a_arcsec"] = a_arcsec
    if i_deg is not None:
        out["inclination_deg"] = i_deg
    if e is not None:
        out["eccentricity"] = e
    if override_key in v and override_value is not None:
        out["override_key"] = override_key

    out["warnings"] = warn
    return out
def residuals_plotly(name: str):
    """
    Plotly version of backend.residual_plots():
    Δρ and Δθ vs epoch, plus side boxplots (same structure as PDF).
    """
    if int(backend.orb.obj.get("npos", 0)) <= 0:
        return None

    # Observations
    t_obs = backend.orb.pos[:, 0]
    theta_obs = backend.orb.pos[:, 1]  # deg
    rho_obs = backend.orb.pos[:, 2]    # arcsec

    # Fitted values at observed epochs: columns [theta_fit, rho_fit]
    res_obs = backend.eph(backend.orb.el, t_obs, rho=True)
    theta_fit, rho_fit = res_obs[:, 0], res_obs[:, 1]

    # Residuals (match backend.residual_plots)
    dtheta = theta_obs - theta_fit
    dtheta = (dtheta + 180) % 360 - 180   # wrap to [-180,180)
    drho = rho_obs - rho_fit

    fig = make_subplots(
        rows=2, cols=2,
        column_widths=[0.78, 0.22],
        row_heights=[0.5, 0.5],
        horizontal_spacing=0.06,
        vertical_spacing=0.10,
        shared_xaxes=True,
        specs=[
            [{"type": "scatter"}, {"type": "box"}],
            [{"type": "scatter"}, {"type": "box"}],
        ],
        subplot_titles=(None, "ρ residuals", None, "θ residuals")
    )

    # Δρ vs epoch (top-left)
    fig.add_trace(
        go.Scatter(
            x=t_obs, y=drho,
            mode="markers",
            name="Δρ",
            marker=dict(size=8, symbol="triangle-up")
        ),
        row=1, col=1
    )
    fig.add_hline(y=0, row=1, col=1)

    # Δθ vs epoch (bottom-left)
    fig.add_trace(
        go.Scatter(
            x=t_obs, y=dtheta,
            mode="markers",
            name="Δθ",
            marker=dict(size=8, symbol="star")
        ),
        row=2, col=1
    )
    fig.add_hline(y=0, row=2, col=1)

    # Boxplots (right column)
    fig.add_trace(go.Box(y=drho, name="Δρ", boxpoints=False), row=1, col=2)
    fig.add_trace(go.Box(y=dtheta, name="Δθ", boxpoints=False), row=2, col=2)

    # Axis labels (match PDF intention)
    fig.update_yaxes(title_text="Δρ (arcsec)", row=1, col=1)
    fig.update_yaxes(title_text="Δθ (deg)", row=2, col=1)
    fig.update_xaxes(title_text="Epoch (year)", row=2, col=1)

    # Keep right boxplots aligned vertically with left panels
    # by matching y-ranges
    y1 = [float(np.nanmin(drho)), float(np.nanmax(drho))]
    y2 = [float(np.nanmin(dtheta)), float(np.nanmax(dtheta))]
    if np.isfinite(y1).all() and y1[0] != y1[1]:
        pad = 0.08 * (y1[1] - y1[0])
        fig.update_yaxes(range=[y1[0]-pad, y1[1]+pad], row=1, col=1)
        fig.update_yaxes(range=[y1[0]-pad, y1[1]+pad], row=1, col=2)
    if np.isfinite(y2).all() and y2[0] != y2[1]:
        pad = 0.08 * (y2[1] - y2[0])
        fig.update_yaxes(range=[y2[0]-pad, y2[1]+pad], row=2, col=1)
        fig.update_yaxes(range=[y2[0]-pad, y2[1]+pad], row=2, col=2)

    fig.update_layout(
        title=f"Residuals (Δρ, Δθ) vs Epoch: {name}",
        height=620,
        legend=dict(orientation="h"),
        margin=dict(l=40, r=20, t=70, b=40),
    )

    return fig

def orbplot_plotly():
    """
    Plotly replacement for orbplot_streamlit().
    Returns a list of dicts: [{title, fig_json}, ...]
    """
    figs = []
    gr = 180 / np.pi
    name = backend.orb.obj.get("name", backend.orb.obj.get("fname", "Orbit"))

    # ---------- Visual Orbit ----------
    if int(backend.orb.obj.get("npos", 0)) > 0:
        # model orbit
        time = np.linspace(0, backend.orb.el[0], 200) + backend.orb.el[1]
        xye = backend.eph(backend.orb.el, time)
        # observed
        xobs = -backend.orb.pos[:, 2] * np.sin(backend.orb.pos[:, 1] / gr)
        yobs =  backend.orb.pos[:, 2] * np.cos(backend.orb.pos[:, 1] / gr)

        # initial orbit overlay (your code keeps backend.orb.initial_el)
        fig = go.Figure()
        # initial orbit overlay (only if available)
        init_el = getattr(backend.orb, "initial_el", None)
        if init_el is not None:
            xye_init = backend.eph(init_el, time)
            fig.add_trace(go.Scatter(
                x=-xye_init[:, 1], y=xye_init[:, 0],
                mode="lines", name="Initial Orbit",
                line=dict(dash="dot")
            ))
        fig.add_trace(go.Scatter(
            x=-xye[:, 1], y=xye[:, 0],
            mode="lines", name="Fitted Orbit"
        ))
        t_obs = backend.orb.pos[:, 0]
        labels = [f"{t:.2f}" for t in t_obs]
        fig.add_trace(go.Scatter(
            x=xobs, y=yobs,
            mode="markers+text",
            name="Observed",
            marker=dict(size=8, symbol="circle"),
            text=labels,                         # time labels
            textposition="top center",
            textfont=dict(size=11),
            hovertemplate="t=%{text}<br>x=%{x:.4f}<br>y=%{y:.4f}<extra></extra>"
        ))
        fig.update_layout(
            title=f"Visual Orbit: {name}",
            xaxis_title="ΔRA (arcsec)",
            yaxis_title="ΔDec (arcsec)",
            width=None, height=520,
            legend=dict(orientation="h")
        )
        fig.update_yaxes(scaleanchor="x", scaleratio=1)  # equal aspect
        figs.append({"title": "Visual Orbit", "fig_json": json.dumps(fig, cls=PlotlyJSONEncoder)})

    # ---------- RV vs Time ----------
    if int(backend.orb.obj.get("nrv1", 0)) > 0 or int(backend.orb.obj.get("nrv2", 0)) > 0:
        # Build model times over the actual observed RV time span
        tmin, tmax = None, None
        if int(backend.orb.obj.get("nrv1", 0)) > 0 and getattr(backend.orb, "rv1", None) is not None and len(backend.orb.rv1) > 0:
            tmin = backend.orb.rv1[:, 0].min()
            tmax = backend.orb.rv1[:, 0].max()
        if int(backend.orb.obj.get("nrv2", 0)) > 0 and getattr(backend.orb, "rv2", None) is not None and len(backend.orb.rv2) > 0:
            tmin2 = backend.orb.rv2[:, 0].min()
            tmax2 = backend.orb.rv2[:, 0].max()
            tmin = tmin2 if tmin is None else min(tmin, tmin2)
            tmax = tmax2 if tmax is None else max(tmax, tmax2)
        # Fallback: if something is odd, use old method
        if tmin is None or tmax is None or not np.isfinite(tmin) or not np.isfinite(tmax) or tmax <= tmin:
            t = np.linspace(0, backend.orb.el[0], 600) + backend.orb.el[1]
        else:
            pad = 0.03 * (tmax - tmin)
            t = np.linspace(tmin - pad, tmax + pad, 800)
        rv = backend.eph(backend.orb.el, t, rv=True)
        fig2 = go.Figure()
        if int(backend.orb.obj.get("nrv1", 0)) > 0 and getattr(backend.orb, "rv1", None) is not None and len(backend.orb.rv1) > 0:
            rv1_corr = backend.orb.rv1[:, 1] - backend.rv_offset_for_points(backend.orb.rv1_source, len(backend.orb.rv1))
            fig2.add_trace(go.Scatter(
                x=backend.orb.rv1[:, 0], y=rv1_corr,
                mode="markers", name="Primary RV",
                error_y=dict(type="data", array=backend.orb.rv1[:, 2], visible=True),
                marker=dict(size=8, symbol="circle")
            ))
            fig2.add_trace(go.Scatter(
                x=t, y=rv[:, 0],
                mode="lines", name="Primary Fit"
            ))
            if int(backend.orb.obj.get("nrv2", 0)) > 0 and getattr(backend.orb, "rv2", None) is not None and len(backend.orb.rv2) > 0:
                rv2_corr = backend.orb.rv2[:, 1] - backend.rv_offset_for_points(backend.orb.rv2_source, len(backend.orb.rv2))
                fig2.add_trace(go.Scatter(
                    x=backend.orb.rv2[:, 0], y=rv2_corr,
                    mode="markers", name="Secondary RV",
                    error_y=dict(type="data", array=backend.orb.rv2[:, 2], visible=True),
                    marker=dict(size=8, symbol="diamond")
                ))
                fig2.add_trace(go.Scatter(
                    x=t, y=rv[:, 1],
                    mode="lines", name="Secondary Fit",
                    line=dict(dash="dash")
                ))
            fig2.update_layout(
                title=f"RV vs Time: {name}",
                xaxis_title="Time",
                yaxis_title="Radial Velocity (km/s)",
                height=520,
                legend=dict(orientation="h")
            )
            figs.append({"title": "RV vs Time", "fig_json": json.dumps(fig2, cls=PlotlyJSONEncoder)})

    # ---------- RV vs Phase ----------
    if int(backend.orb.obj.get("nrv1", 0)) > 0 or int(backend.orb.obj.get("nrv2", 0)) > 0:
        phases = np.linspace(0, 1, 200)
        t_phase = phases * backend.orb.el[0] + backend.orb.el[1]
        rv_phase = backend.eph(backend.orb.el, t_phase, rv=True)

        fig3 = go.Figure()
        if int(backend.orb.obj.get("nrv1", 0)) > 0:
            phase1 = ((backend.orb.rv1[:, 0] - backend.orb.el[1]) / backend.orb.el[0]) % 1
            rv1_corr = backend.orb.rv1[:, 1] - backend.rv_offset_for_points(backend.orb.rv1_source, len(backend.orb.rv1))
            fig3.add_trace(go.Scatter(
                x=phase1, y=rv1_corr,
                mode="markers", name="Primary RV",
                error_y=dict(type="data", array=backend.orb.rv1[:, 2], visible=True),
                marker=dict(size=8, symbol="circle")
            ))
            fig3.add_trace(go.Scatter(
                x=phases, y=rv_phase[:, 0],
                mode="lines", name="Primary Fit"
            ))

        if int(backend.orb.obj.get("nrv2", 0)) > 0:
            phase2 = ((backend.orb.rv2[:, 0] - backend.orb.el[1]) / backend.orb.el[0]) % 1
            rv2_corr = backend.orb.rv2[:, 1] - backend.rv_offset_for_points(backend.orb.rv2_source, len(backend.orb.rv2))
            fig3.add_trace(go.Scatter(
                x=phase2, y=rv2_corr,
                mode="markers", name="Secondary RV",
                error_y=dict(type="data", array=backend.orb.rv2[:, 2], visible=True),
                marker=dict(size=8, symbol="diamond")
            ))
            fig3.add_trace(go.Scatter(
                x=phases, y=rv_phase[:, 1],
                mode="lines", name="Secondary Fit",
                line=dict(dash="dash")
            ))

        fig3.update_layout(
            title=f"RV vs Phase: {name}",
            xaxis_title="Phase",
            yaxis_title="Radial Velocity (km/s)",
            height=520,
            legend=dict(orientation="h")
        )
        figs.append({"title": "RV vs Phase", "fig_json": json.dumps(fig3, cls=PlotlyJSONEncoder)})

    # ---------- Residuals (Interactive) ----------
    fig_res = residuals_plotly(name)
    if fig_res is not None:
        figs.append({"title": "Residuals (Δρ & Δθ)", "fig_json": json.dumps(fig_res, cls=PlotlyJSONEncoder)})
    return figs


def _snapshot_points():
    """Capture the full position/RV1/RV2 arrays and source labels as read
    from the file, before restrict_points() (if any) shrinks them. Row i
    here is always "original file row i" - the stable index used by the
    point-selection checkboxes and by point_filter."""
    def _copy(name, ncols):
        n = int(backend.orb.obj.get(f"n{name}", 0))
        if n == 0:
            return np.zeros((0, ncols)), []
        return getattr(backend.orb, name).copy(), list(getattr(backend.orb, f"{name}_source"))
    pos, pos_src = _copy("pos", 6)
    rv1, rv1_src = _copy("rv1", 3)
    rv2, rv2_src = _copy("rv2", 3)
    return {"pos": pos, "pos_source": pos_src, "rv1": rv1, "rv1_source": rv1_src, "rv2": rv2, "rv2_source": rv2_src}


def _build_point_tables(full, point_filter):
    """Per-observation data + residuals (observed minus fitted, using the
    just-fitted backend.orb.el) for EVERY point read from the file - not
    just the ones used in the fit - so a previously-excluded point can be
    re-included later. `point_filter` (as passed to run_fit_to_dir) marks
    which are currently included; a missing/None key means "all included".
    Call after fitorb() so the residuals reflect the fit just performed.
    """
    tables = {"pos": [], "rv1": [], "rv2": []}
    keep = point_filter or {}

    pos = full["pos"]
    if pos.shape[0] > 0:
        t_obs, theta_obs, rho_obs, err = pos[:, 0], pos[:, 1], pos[:, 2], pos[:, 3]
        res = backend.eph(backend.orb.el, t_obs, rho=True)
        dtheta = (theta_obs - res[:, 0] + 180) % 360 - 180
        drho = rho_obs - res[:, 1]
        included = set(keep["pos"]) if keep.get("pos") is not None else None
        src = full["pos_source"]
        for i in range(pos.shape[0]):
            tables["pos"].append({
                "idx": i, "epoch": float(t_obs[i]), "theta": float(theta_obs[i]),
                "rho": float(rho_obs[i]), "err": float(err[i]), "source": src[i] if i < len(src) else "",
                "dtheta": float(dtheta[i]), "drho": float(drho[i]),
                "included": included is None or i in included,
            })

    for key, col in (("rv1", 0), ("rv2", 1)):
        arr = full[key]
        if arr.shape[0] > 0:
            t, rv_obs, err = arr[:, 0], arr[:, 1], arr[:, 2]
            src = full[f"{key}_source"]
            offset = backend.rv_offset_for_points(src, arr.shape[0])
            rv_fit = backend.eph(backend.orb.el, t, rv=True)[:, col] + offset
            drv = rv_obs - rv_fit
            included = set(keep[key]) if keep.get(key) is not None else None
            for i in range(arr.shape[0]):
                tables[key].append({
                    "idx": i, "epoch": float(t[i]), "rv": float(rv_obs[i]),
                    "err": float(err[i]), "source": src[i] if i < len(src) else "", "drv": float(drv[i]),
                    "included": included is None or i in included,
                })

    return tables


def run_fit_to_dir(run_dir: Path, uploaded_path: Path, fixel_map: dict[str, int],
                    point_filter: dict | None = None) -> dict:
    """
    readinp/readcsv_custom -> [restrict_points] -> apply fixel -> fitorb ->
    orbsave -> plots. All outputs written into run_dir.
    Returns dict for templating (ALWAYS includes stats).

    point_filter, if given, is a dict with any of "pos"/"rv1"/"rv2" mapped
    to a list of row indices (into the freshly-read file, before any
    filtering) to keep; omitted keys keep all rows of that type.
    """
    cwd = os.getcwd()
    os.chdir(run_dir)
    try:
        ext = uploaded_path.suffix.lower()
        if ext == ".inp":
            backend.readinp(str(uploaded_path))
        else:
            if hasattr(backend, "readcsv_custom"):
                backend.readcsv_custom(str(uploaded_path))
            else:
                raise RuntimeError("CSV upload not supported: backend.readcsv_custom not found in backend.py")

        full_points = _snapshot_points()

        if point_filter:
            backend.restrict_points(
                pos_keep=point_filter.get("pos"),
                rv1_keep=point_filter.get("rv1"),
                rv2_keep=point_filter.get("rv2"),
            )

        # Apply fit/fix flags (1=fit, 0=fixed)
        for i, nm in enumerate(list(backend.orb.elname)):
            backend.orb.fixel[i] = int(fixel_map.get(nm, 1))

        backend.fitorb()
        backend.orbsave()

        # Output file naming
        stem = Path(backend.orb.obj.get("fname", uploaded_path.name)).stem
        output_csv = f"{stem}_output.csv"
        output_csv_path = run_dir / output_csv

        # Plots
        # ---------- PLOTS: (A) PNG for PDF + (B) Plotly for interactive web ----------
        # (A) Save Matplotlib PNG plots for PDF / archiving
        plot_files = []
        plot_func = None
        for cand in ("orbplot_streamlit", "orbplot", "make_plots"):
            if hasattr(backend, cand):
                plot_func = getattr(backend, cand)
                break
        if plot_func is not None:
            figs = plot_func()
            for idx, fig in enumerate(figs, start=1):
                p = run_dir / f"plot_{idx}.png"
                save_fig(fig, p)
                plot_files.append(p.name)
        else:
            print("No plot function found in backend (expected orbplot_streamlit/orbplot/make_plots).")
        # (B) Build Plotly interactive plots for the web UI
        plotly_figs = []
        try:
            plotly_figs = orbplot_plotly()   # your Plotly builder function
        except Exception:
            print("Plotly plot build failed:\n", traceback.format_exc())


        residual_file = None
        if int(backend.orb.obj.get("npos", 0)) > 0 and hasattr(backend, "residual_plots"):
            fig = backend.residual_plots()
            p = run_dir / "residuals.png"
            save_fig(fig, p)
            residual_file = p.name

        # Elements table
        elements = []
        elerr = getattr(backend.orb, "elerr", [0.0] * len(list(backend.orb.elname)))
        for i, nm in enumerate(list(backend.orb.elname)):
            elements.append(
                {
                    "name": nm,
                    "value": float(backend.orb.el[i]),
                    "err": float(elerr[i]) if i < len(elerr) else 0.0,
                    "fit": int(backend.orb.fixel[i]),
                }
            )

        # Stats (ALWAYS dict)
        obj = getattr(backend.orb, "obj", {}) or {}
        stats = default_stats()
        stats.update(
            {
                "name": obj.get("name", "") or "",
                "npos": int(obj.get("npos", 0) or 0),
                "nrv1": int(obj.get("nrv1", 0) or 0),
                "nrv2": int(obj.get("nrv2", 0) or 0),
                "parallax": float(obj.get("parallax", 0.0) or 0.0),
                "chi2": float(obj.get("chi2", 0.0) or 0.0),
            }
        )
        try:
            chi2n = obj.get("chi2n", [0, 0, 0, 0])
            rms = obj.get("rms", [0, 0, 0, 0])
            stats.update(
                {
                    "chi2n_theta": float(chi2n[0]),
                    "chi2n_rho": float(chi2n[1]),
                    "chi2n_rv1": float(chi2n[2]),
                    "chi2n_rv2": float(chi2n[3]),
                    "rms_theta": float(rms[0]),
                    "rms_rho": float(rms[1]),
                    "rms_rv1": float(rms[2]),
                    "rms_rv2": float(rms[3]),
                }
            )
        except Exception:
            pass
        derived = compute_derived(elements, stats)
        point_tables = _build_point_tables(full_points, point_filter)

        # Persist the numeric results, plus enough to re-run a fit later
        # (fixel choices, the uploaded filename, current point selection),
        # so /report/<run_id> and /refit/<run_id> (separate later requests)
        # work without depending on the shared, possibly-since-overwritten
        # backend.orb global.
        try:
            (run_dir / "summary.json").write_text(json.dumps({
                "run_id": run_dir.name, "elements": elements, "stats": stats, "derived": derived,
                "plotly_figs": plotly_figs,  # so the PDF report can render the same plots
                "point_tables": point_tables,
                "fixel_map": fixel_map,
                "uploaded_filename": uploaded_path.name,
                "point_filter": point_filter,
            }))
        except Exception:
            pass

        return {
            "run_id": run_dir.name,
            "elements": elements,
            "stats": stats,
            "derived": derived,
            "plotly_figs": plotly_figs,     # web interactive
            "plot_files": plot_files,       # optional (kept for debugging)
            "residual_file": residual_file, # still useful
            "output_csv": output_csv if output_csv_path.exists() else None,
            "point_tables": point_tables,
        }

    finally:
        os.chdir(cwd)


app = Flask(__name__)
app.secret_key = os.environ.get("FLASK_SECRET_KEY", "change-me")

# Plain-language reference for the orbital elements, shown as tooltips,
# under-checkbox captions, and the expandable glossary on the upload page.
ELEMENT_FULL_NAMES = {
    "P": "Orbital period",
    "T": "Time of periastron",
    "e": "Eccentricity",
    "a": "Semi-major axis",
    "W": "Longitude of ascending node (Ω)",
    "w": "Argument of periastron (ω)",
    "i": "Inclination",
    "K1": "RV semi-amplitude (primary)",
    "K2": "RV semi-amplitude (secondary)",
    "V0": "Systemic velocity",
}
ELEMENT_GLOSSARY = {
    "P": "How long one full orbit takes (years or days, matching your file).",
    "T": "The epoch when the two stars are closest together (same time system as P).",
    "e": "How elongated the orbit is: 0 = perfect circle, close to 1 = very stretched ellipse.",
    "a": "The apparent size of the orbit on the sky, in arcseconds (for astrometric/visual data).",
    "W": "Which way the orbit is tilted on the sky, in degrees — orientation of the line of nodes.",
    "w": "Where periastron (closest approach) sits within the orbital plane, in degrees.",
    "i": "How tilted the orbital plane is to our line of sight: 0° = face-on, 90° = edge-on.",
    "K1": "How fast the primary star's radial velocity swings up and down, in km/s.",
    "K2": "How fast the secondary star's radial velocity swings up and down, in km/s.",
    "V0": "The whole system's velocity through space, in km/s (the RV curves' midline).",
}
ELEMENT_SHORT_LABELS = {
    "P": "period", "T": "periastron time", "e": "eccentricity", "a": "semi-major axis",
    "W": "node angle", "w": "periastron angle", "i": "inclination",
    "K1": "RV amp. 1", "K2": "RV amp. 2", "V0": "systemic v",
}


@app.route("/", methods=["GET"])
def index():
    elnames = list(getattr(backend, "orb").elname)
    defaults = {nm: 1 for nm in elnames}  # checked=fit
    return render_template(
        "index.html", elnames=elnames, defaults=defaults,
        full_names=ELEMENT_FULL_NAMES, glossary=ELEMENT_GLOSSARY, short_labels=ELEMENT_SHORT_LABELS,
    )


@app.route("/run", methods=["POST"])
def run():
    run_id = uuid.uuid4().hex[:12]  # create early so error page can show it
    try:
        if "datafile" not in request.files:
            flash("No file selected.", "error")
            return redirect(url_for("index"))

        f = request.files["datafile"]
        if not f.filename:
            flash("No file selected.", "error")
            return redirect(url_for("index"))

        if not allowed_file(f.filename):
            flash("Unsupported file type. Upload .inp, .csv, or .txt", "error")
            return redirect(url_for("index"))

        run_dir = RUNS_DIR / run_id
        run_dir.mkdir(parents=True, exist_ok=True)

        filename = secure_filename(f.filename)
        uploaded_path = run_dir / filename
        f.save(uploaded_path)

        elnames = list(getattr(backend, "orb").elname)
        fixel_map = {nm: (1 if request.form.get(f"fit_{nm}") == "on" else 0) for nm in elnames}

        summary = run_fit_to_dir(run_dir, uploaded_path, fixel_map)
        return render_template("result.html", glossary=ELEMENT_GLOSSARY, derived_labels=DERIVED_QUANTITY_LABELS, **summary)

    except Exception as e:
        tb = traceback.format_exc()
        # error.html displays `error` prominently itself, so no flash() here
        # to avoid showing the same message twice.
        # IMPORTANT: render error page with safe defaults (stats always exists)
        return render_template(
            "error.html",
            error=str(e),
            traceback=tb,
            run_id=run_id,
            stats=default_stats(),
            elements=[],
            plot_files=[],
            residual_file=None,
            output_csv=None,
        )


@app.route("/refit/<run_id>", methods=["POST"])
def refit(run_id):
    """Re-run the fit for an existing run using only the selected data
    points, keeping the same fit/fix parameter choices as the original run
    (change those from the upload page instead)."""
    run_dir = RUNS_DIR / run_id
    summary_path = run_dir / "summary.json"
    if not run_dir.exists() or not summary_path.exists():
        flash("Run not found or its saved results have expired.", "error")
        return redirect(url_for("index"))

    try:
        data = json.loads(summary_path.read_text())
    except Exception:
        flash("Could not read this run's saved results.", "error")
        return redirect(url_for("index"))

    uploaded_filename = data.get("uploaded_filename")
    uploaded_path = (run_dir / uploaded_filename) if uploaded_filename else None
    if not uploaded_path or not uploaded_path.exists():
        flash("The originally uploaded file is no longer available for this run.", "error")
        return redirect(url_for("index"))

    fixel_map = data.get("fixel_map") or {}
    point_tables = data.get("point_tables") or {}
    point_filter = {}
    for key in ("pos", "rv1", "rv2"):
        rows = point_tables.get(key) or []
        if not rows:
            continue
        point_filter[key] = [row["idx"] for row in rows if request.form.get(f"{key}_keep_{row['idx']}") == "on"]

    try:
        summary = run_fit_to_dir(run_dir, uploaded_path, fixel_map, point_filter=point_filter)
    except Exception as e:
        tb = traceback.format_exc()
        return render_template(
            "error.html", error=str(e), traceback=tb, run_id=run_id,
            stats=default_stats(), elements=[], plot_files=[], residual_file=None, output_csv=None,
        )

    return render_template("result.html", glossary=ELEMENT_GLOSSARY, derived_labels=DERIVED_QUANTITY_LABELS, **summary)


@app.route("/runs/<run_id>/<path:filename>")
def run_file(run_id, filename):
    run_dir = RUNS_DIR / run_id
    return send_from_directory(run_dir, filename, as_attachment=False)


@app.route("/download/<run_id>/<path:filename>")
def download(run_id, filename):
    run_dir = RUNS_DIR / run_id
    return send_from_directory(run_dir, filename, as_attachment=True)


@app.route("/cleanup/<run_id>", methods=["POST"])
def cleanup(run_id):
    run_dir = RUNS_DIR / run_id
    if run_dir.exists():
        shutil.rmtree(run_dir, ignore_errors=True)
    flash("Run directory cleaned up.", "info")
    return redirect(url_for("index"))


@app.route("/recompute/<run_id>", methods=["POST"])
def recompute(run_id):
    """Re-solve the derived quantities (parallax/distance/a_AU/masses) from
    a user-supplied value instead of the file's parallax, without re-running
    the orbital fit itself."""
    run_dir = RUNS_DIR / run_id
    summary_path = run_dir / "summary.json"
    if not run_dir.exists() or not summary_path.exists():
        flash("Run not found or its saved results have expired.", "error")
        return redirect(url_for("index"))

    try:
        data = json.loads(summary_path.read_text())
    except Exception:
        flash("Could not read this run's saved results.", "error")
        return redirect(url_for("index"))

    override_key = request.form.get("override_key") or None
    override_raw = request.form.get("override_value", "").strip()
    override_value = None
    if override_key and override_key not in DERIVED_QUANTITY_LABELS:
        override_key = None
    if override_key and override_raw:
        try:
            override_value = float(override_raw)
            if override_value <= 0:
                flash("Please enter a positive number.", "error")
                override_key = None
        except ValueError:
            flash("That doesn't look like a number.", "error")
            override_key = None
    elif override_key and not override_raw:
        flash("Enter a value to use for that quantity.", "error")
        override_key = None

    elements = data.get("elements", [])
    stats = data.get("stats", {})
    derived = compute_derived(elements, stats, override_key, override_value)

    # Persist so a subsequent PDF report reflects the same override.
    data["derived"] = derived
    try:
        summary_path.write_text(json.dumps(data))
    except Exception:
        pass

    output_csv = None
    for f in run_dir.iterdir():
        if f.name.endswith("_output.csv"):
            output_csv = f.name
            break

    return render_template(
        "result.html", glossary=ELEMENT_GLOSSARY, derived_labels=DERIVED_QUANTITY_LABELS,
        run_id=run_id, elements=elements, stats=stats, derived=derived,
        plotly_figs=data.get("plotly_figs", []), output_csv=output_csv,
        point_tables=data.get("point_tables", {}),
    )


_PDF_STYLES = getSampleStyleSheet()
_PDF_TITLE = ParagraphStyle("PdfTitle", parent=_PDF_STYLES["Title"], fontSize=18, spaceAfter=2)
_PDF_META = ParagraphStyle("PdfMeta", parent=_PDF_STYLES["Normal"], fontSize=9, textColor=colors.HexColor("#5a6578"))
_PDF_H2 = ParagraphStyle("PdfH2", parent=_PDF_STYLES["Heading2"], fontSize=13, spaceBefore=10, spaceAfter=4,
                         textColor=colors.HexColor("#1a1a2e"))
_PDF_H3 = ParagraphStyle("PdfH3", parent=_PDF_STYLES["Heading3"], fontSize=10.5, spaceBefore=6, spaceAfter=2)
_PDF_BODY = ParagraphStyle("PdfBody", parent=_PDF_STYLES["Normal"], fontSize=9.5, leading=13)
_PDF_CAPTION = ParagraphStyle("PdfCaption", parent=_PDF_STYLES["Normal"], fontSize=8.5,
                               textColor=colors.HexColor("#5a6578"), spaceBefore=6)
_PDF_WARN = ParagraphStyle("PdfWarn", parent=_PDF_STYLES["Normal"], fontSize=9, textColor=colors.HexColor("#8a5a00"))


def _pdf_table(rows, col_widths=None):
    t = Table(rows, colWidths=col_widths, repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#7c5cff")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#dcdce6")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f6f6fb")]),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
    ]))
    return t


# Optional override for sandboxed/offline environments where Kaleido can't
# auto-detect or download its own Chrome; normally left unset so Kaleido
# manages its own browser (see the Render build step in render.yaml).
_KALEIDO_CHROME_PATH = os.environ.get("KALEIDO_CHROME_PATH")


def _plotly_json_to_png(fig_json: str, width=1000, height=620, scale=2) -> bytes:
    """Render a Plotly figure (as saved by orbplot_plotly) to PNG bytes via
    Kaleido, so the PDF report shows the exact same plot as the web UI."""
    import kaleido
    fig_dict = json.loads(fig_json)
    kopts = {"path": _KALEIDO_CHROME_PATH} if _KALEIDO_CHROME_PATH else {}
    return kaleido.calc_fig_sync(
        fig_dict, opts={"width": width, "height": height, "scale": scale, "format": "png"}, kopts=kopts,
    )


@app.route("/report/<run_id>")
def report(run_id):
    run_dir = RUNS_DIR / run_id
    if not run_dir.exists():
        flash("Run not found.", "error")
        return redirect(url_for("index"))

    data = None
    summary_path = run_dir / "summary.json"
    if summary_path.exists():
        try:
            data = json.loads(summary_path.read_text())
        except Exception:
            data = None

    pdf_path = run_dir / f"{run_id}_report.pdf"
    doc = SimpleDocTemplate(str(pdf_path), pagesize=A4,
                             topMargin=1.8*cm, bottomMargin=1.8*cm, leftMargin=2*cm, rightMargin=2*cm)
    story = []

    story.append(Paragraph("Binary-Star Orbit Fitting Report", _PDF_TITLE))
    story.append(Paragraph(f"Run ID: {run_id}", _PDF_META))
    obj_name = (data or {}).get("stats", {}).get("name")
    if obj_name:
        story.append(Paragraph(f"Object: {obj_name}", _PDF_META))
    story.append(Paragraph(
        "Python Implementation of Tokovinin’s original OrbitX code. This version, enhanced by "
        "Prof. Mashhoor Al-Wardat and Dr. Mohammed Hussin Talafha, analyzes orbits of visual "
        "binaries (VB) using data from the 4th Catalogue of Interferometric Measurements or SOAR "
        "telescope, spectroscopic binaries (SB) using SB9 catalogue data, or any binary system "
        "with both interferometric and spectroscopic measurements.", _PDF_META))
    story.append(Spacer(1, 0.5*cm))

    if data:
        story.append(Paragraph(
            "This report summarises a least-squares fit of a binary-star orbit to the "
            "uploaded observations. The fit adjusts the checked orbital elements until the "
            "predicted orbit matches the data as closely as possible; the quoted errors show "
            "how tightly each element is constrained by the data.", _PDF_BODY))

        elements = data.get("elements", [])
        if elements:
            story.append(Paragraph("Orbital Elements", _PDF_H2))
            rows = [["Parameter", "Value", "Error", "Fit?"]]
            for row in elements:
                err = f'{row["err"]:.3g}' if row.get("fit") == 1 else "—"
                rows.append([row["name"], f'{row["value"]:.10g}', err, "Yes" if row.get("fit") == 1 else "No"])
            story.append(_pdf_table(rows, col_widths=[3.5*cm, 5*cm, 4*cm, 2.5*cm]))

            story.append(Paragraph("What these mean", _PDF_H3))
            gloss_rows = [["Symbol", "Meaning"]]
            for row in elements:
                nm = row["name"]
                if nm in ELEMENT_GLOSSARY:
                    gloss_rows.append([nm, ELEMENT_GLOSSARY[nm]])
            story.append(_pdf_table(gloss_rows, col_widths=[2*cm, 13*cm]))

        stats = data.get("stats", {})
        if stats.get("npos", 0) > 0 or stats.get("nrv1", 0) > 0 or stats.get("nrv2", 0) > 0:
            story.append(Paragraph("Fit Quality", _PDF_H2))
            story.append(Paragraph(
                "Reduced chi-squared (χ²) near 1 means the fit matches the data about as well "
                "as the quoted measurement errors predict; much larger suggests a poor fit or "
                "underestimated errors, much smaller suggests overestimated errors. RMS is the "
                "typical size of the fit's miss, in real units.", _PDF_BODY))
            rows = [["Data type", "Reduced χ²", "RMS residual"]]
            if stats.get("npos", 0) > 0:
                rows.append(["Position angle (θ)", f'{stats["chi2n_theta"]:.3g}', f'{stats["rms_theta"]:.4g}°'])
                rows.append(["Separation (ρ)", f'{stats["chi2n_rho"]:.3g}', f'{stats["rms_rho"]:.4g}"'])
            if stats.get("nrv1", 0) > 0:
                rows.append(["RV1 (primary)", f'{stats["chi2n_rv1"]:.3g}', f'{stats["rms_rv1"]:.4g} km/s'])
            if stats.get("nrv2", 0) > 0:
                rows.append(["RV2 (secondary)", f'{stats["chi2n_rv2"]:.3g}', f'{stats["rms_rv2"]:.4g} km/s'])
            story.append(_pdf_table(rows, col_widths=[6*cm, 4.5*cm, 4.5*cm]))

        derived = data.get("derived", {})
        if derived:
            story.append(Paragraph("Derived Quantities (Masses)", _PDF_H2))
            for w in derived.get("warnings", []):
                story.append(Paragraph(f"Note: {w}", _PDF_WARN))
            labels = [
                ("parallax_mas", "Parallax (mas)"), ("distance_pc", "Distance (pc)"),
                ("a_AU", "a (AU)"), ("P_yr", "P (yr used)"),
                ("M_total_Msun", "M_total (Msun)"), ("q_M2_over_M1", "q = M2/M1"),
                ("M1_Msun", "M1 (Msun)"), ("M2_Msun", "M2 (Msun)"),
            ]
            rows = [["Quantity", "Value"]]
            for key, label in labels:
                if key in derived:
                    rows.append([label, f'{derived[key]:.6g}'])
            if len(rows) > 1:
                story.append(_pdf_table(rows, col_widths=[8*cm, 7*cm]))
    else:
        story.append(Paragraph(
            "Numeric results were not available when this report was generated; showing plots only.",
            _PDF_BODY))

    # Plots. Prefer rendering the same Plotly figures shown on the results
    # page (via Kaleido) so the PDF visually matches the web UI; fall back
    # to the separate Matplotlib PNGs if that's unavailable (e.g. Kaleido's
    # browser isn't installed) or for older runs saved before this existed.
    max_w = A4[0] - 4*cm
    max_h = A4[1] - 6*cm

    def _add_image(target, caption, img_reader_source, is_bytesio=False):
        img = ImageReader(img_reader_source)
        iw, ih = img.getSize()
        draw_w, draw_h = max_w, ih * (max_w / iw)
        if draw_h > max_h:
            draw_h, draw_w = max_h, iw * (max_h / ih)
        target.append(Paragraph(caption, _PDF_CAPTION))
        if is_bytesio:
            img_reader_source.seek(0)
        target.append(RLImage(img_reader_source, width=draw_w, height=draw_h))

    plotly_figs = (data or {}).get("plotly_figs") or []
    plots_added = False
    if plotly_figs:
        plotly_story = [
            PageBreak(),
            Paragraph("Plots", _PDF_H2),
            Paragraph("These match the interactive plots in the web UI.", _PDF_BODY),
        ]
        try:
            for item in plotly_figs:
                png_bytes = _plotly_json_to_png(item["fig_json"])
                _add_image(plotly_story, item.get("title", ""), io.BytesIO(png_bytes), is_bytesio=True)
            story.extend(plotly_story)
            plots_added = True
        except Exception:
            traceback.print_exc()  # fall through to the file-based plots below

    if not plots_added:
        plot_files = sorted(run_dir.glob("plot_*.png"))
        if (run_dir / "residuals.png").exists():
            plot_files.append(run_dir / "residuals.png")
        if plot_files:
            story.append(PageBreak())
            story.append(Paragraph("Plots", _PDF_H2))
            story.append(Paragraph("Click-to-enlarge and marker styling are available in the web UI; this PDF embeds the same images.", _PDF_BODY))
            for p in plot_files:
                try:
                    _add_image(story, p.name, str(p))
                except Exception:
                    story.append(Paragraph(f"(Could not embed {p.name})", _PDF_BODY))

    doc.build(story)
    return send_file(str(pdf_path), as_attachment=True, download_name=pdf_path.name)


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 8000))
    debug = os.environ.get("FLASK_DEBUG", "0") == "1"
    app.run(host="0.0.0.0", port=port, debug=debug)
