import streamlit as st
from pathlib import Path
import tempfile
import shutil
import json
import io

import plotly.graph_objects as go
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib.utils import ImageReader
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Image as RLImage,
    PageBreak,
)

from app import (
    run_fit_to_dir,
    backend,
    compute_derived,
    DERIVED_QUANTITY_LABELS,
    ELEMENT_GLOSSARY,
    _pdf_table,
    _plotly_json_to_png,
    _PDF_TITLE,
    _PDF_META,
    _PDF_H2,
    _PDF_H3,
    _PDF_CAPTION,
    _PDF_WARN,
    ELEMENT_SHORT_LABELS,
)


def build_streamlit_pdf(summary, run_dir):
    """Build a PDF from the current Streamlit summary without refitting."""
    pdf_buffer = io.BytesIO()
    doc = SimpleDocTemplate(
        pdf_buffer,
        pagesize=A4,
        topMargin=1.8 * cm,
        bottomMargin=1.8 * cm,
        leftMargin=2 * cm,
        rightMargin=2 * cm,
    )
    story = []

    stats = summary.get("stats", {})
    elements = summary.get("elements", [])
    derived = summary.get("derived", {})

    story.append(Paragraph("Binary-Star Orbit Fitting Report", _PDF_TITLE))
    story.append(Paragraph(f"Run ID: {summary.get('run_id', 'Unknown')}", _PDF_META))
    story.append(Paragraph(f"Object: {stats.get('name') or 'Unknown'}", _PDF_META))
    story.append(Spacer(1, 0.5 * cm))

    story.append(Paragraph("Results Summary", _PDF_H2))
    story.append(
        _pdf_table(
            [
                ["Quantity", "Value"],
                ["Positions", str(stats.get("npos", 0))],
                ["RV1", str(stats.get("nrv1", 0))],
                ["RV2", str(stats.get("nrv2", 0))],
                ["Parallax (mas)", str(stats.get("parallax", 0.0))],
                ["Chi²", f"{stats.get('chi2', 0.0):.4g}"],
            ],
            col_widths=[8 * cm, 7 * cm],
        )
    )

    if elements:
        story.append(Paragraph("Orbital Elements", _PDF_H2))
        element_rows = [["Parameter", "Value", "Unit", "Error", "Fit?"]]
        for row in elements:
            element_rows.append(
                [
                    row["name"],
                    f'{row["value"]:.10g}',
                    row.get("unit", ""),
                    f'{row.get("err", 0.0):.3g}' if row.get("fit") == 1 else "-",
                    "Yes" if row.get("fit") == 1 else "No",
                ]
            )
        story.append(
            _pdf_table(
                element_rows,
                col_widths=[3 * cm, 4.3 * cm, 2 * cm, 3.2 * cm, 2 * cm],
            )
        )
        story.append(Paragraph("What these mean", _PDF_H3))
        glossary_rows = [["Symbol", "Meaning"]]
        for row in elements:
            if row["name"] in ELEMENT_GLOSSARY:
                glossary_rows.append([row["name"], ELEMENT_GLOSSARY[row["name"]]])
        story.append(_pdf_table(glossary_rows, col_widths=[2 * cm, 13 * cm]))

    story.append(Paragraph("Fit Quality", _PDF_H2))
    statistics_rows = [["Data type", "Reduced Chi²", "RMS residual"]]
    if stats.get("npos", 0) > 0:
        statistics_rows.extend(
            [
                [
                    "Position angle (theta)",
                    f'{stats.get("chi2n_theta", 0.0):.3g}',
                    f'{stats.get("rms_theta", 0.0):.4g} deg',
                ],
                [
                    "Separation (rho)",
                    f'{stats.get("chi2n_rho", 0.0):.3g}',
                    f'{stats.get("rms_rho", 0.0):.4g} arcsec',
                ],
            ]
        )
    if stats.get("nrv1", 0) > 0:
        statistics_rows.append(
            [
                "RV1 (primary)",
                f'{stats.get("chi2n_rv1", 0.0):.3g}',
                f'{stats.get("rms_rv1", 0.0):.4g} km/s',
            ]
        )
    if stats.get("nrv2", 0) > 0:
        statistics_rows.append(
            [
                "RV2 (secondary)",
                f'{stats.get("chi2n_rv2", 0.0):.3g}',
                f'{stats.get("rms_rv2", 0.0):.4g} km/s',
            ]
        )
    story.append(_pdf_table(statistics_rows, col_widths=[6 * cm, 4.5 * cm, 4.5 * cm]))

    if derived:
        story.append(Paragraph("Derived Quantities", _PDF_H2))
        for warning in derived.get("warnings", []):
            story.append(Paragraph(f"Note: {warning}", _PDF_WARN))
        derived_rows = [["Quantity", "Value"]]
        for key, value in derived.items():
            if key != "warnings":
                display_value = f"{value:.6g}" if isinstance(value, (int, float)) else str(value)
                derived_rows.append([key, display_value])
        if len(derived_rows) > 1:
            story.append(_pdf_table(derived_rows, col_widths=[8 * cm, 7 * cm]))

    point_table_labels = {
        "pos": ("Positions", ["idx", "epoch", "theta", "rho", "err", "source", "dtheta", "drho", "included"]),
        "rv1": ("RV1 (primary)", ["idx", "epoch", "rv", "err", "source", "drv", "included"]),
        "rv2": ("RV2 (secondary)", ["idx", "epoch", "rv", "err", "source", "drv", "included"]),
    }
    for table_name, (title, columns) in point_table_labels.items():
        rows = summary.get("point_tables", {}).get(table_name, [])
        if not rows:
            continue
        story.append(Paragraph(title, _PDF_H2))
        table_rows = [columns]
        table_rows.extend([[str(row.get(column, "")) for column in columns] for row in rows])
        story.append(_pdf_table(table_rows))

    plotly_figs = summary.get("plotly_figs", [])
    plot_images = []
    plot_render_failed = False
    for item in plotly_figs:
        try:
            plot_images.append(
                (item.get("title", ""), io.BytesIO(_plotly_json_to_png(item["fig_json"])))
            )
        except Exception:
            plot_render_failed = True

    if plot_images:
        story.append(PageBreak())
        story.append(Paragraph("Plots", _PDF_H2))
        for caption, image_source in plot_images:
            image = ImageReader(image_source)
            image_width, image_height = image.getSize()
            max_width = A4[0] - 4 * cm
            max_height = A4[1] - 6 * cm
            draw_width = max_width
            draw_height = image_height * (max_width / image_width)
            if draw_height > max_height:
                draw_height = max_height
                draw_width = image_width * (max_height / image_height)
            story.append(Paragraph(caption, _PDF_CAPTION))
            story.append(RLImage(image_source, width=draw_width, height=draw_height))

    if not plot_images or plot_render_failed:
        fallback_names = []
        fallback_names.extend(summary.get("plot_files", []))
        fallback_names.extend(summary.get("residual_files", []))
        if summary.get("residual_file"):
            fallback_names.append(summary["residual_file"])
        fallback_names = list(dict.fromkeys(fallback_names))
        fallback_paths = [run_dir / name for name in fallback_names if (run_dir / name).exists()]
        if fallback_paths:
            story.append(PageBreak())
            story.append(Paragraph("Saved Plot Files", _PDF_H2))
            for path in fallback_paths:
                story.append(Paragraph(path.name, _PDF_CAPTION))
                story.append(RLImage(str(path), width=15 * cm, height=9 * cm))

    doc.build(story)
    return pdf_buffer.getvalue()

# ---------------------------------------------------------
# Page configuration
# ---------------------------------------------------------

st.set_page_config(
    page_title="PySVOrbit",
    layout="wide",
)


# ---------------------------------------------------------
# Render website styling
# ---------------------------------------------------------

st.markdown(
    """
    <style>

    :root {
    --bg: #f4f6fa;
    --panel: #ffffff;
    --panel2: #f8fafc;
    --text: #1f2937;
    --muted: #5f6b7a;
    --border: #d9dee7;
    --accent: #6d4aff;
    --accent2: #14b8a6;
    --danger: #dc2626;
}

/* ---------- Main application ---------- */

.stApp {
    background:
        radial-gradient(
            1000px 600px at 10% -10%,
            rgba(109, 74, 255, 0.08),
            transparent 60%
        ),
        radial-gradient(
            900px 600px at 110% 10%,
            rgba(20, 184, 166, 0.07),
            transparent 55%
        ),
        var(--bg);
    color: var(--text);
}

/* ---------- Header ---------- */

header[data-testid="stHeader"] {
    background: transparent;
}

/* ---------- Main content ---------- */

.block-container {
    max-width: 1400px;
    padding-top: 1.5rem;
    padding-right: 2rem;
    padding-bottom: 2rem;
    padding-left: 2rem;
}

@media (max-width: 900px) {
    .block-container {
        padding-right: 1rem;
        padding-left: 1rem;
    }
}

/* ---------- Text ---------- */

h1,
h2,
h3,
h4 {
    color: var(--text) !important;
}

p {
    color: var(--text);
}

label {
    color: var(--text) !important;
}

.pysvorbit-subtitle {
    color: var(--muted) !important;
    max-width: 750px;
    line-height: 1.5;
}

/* ---------- Navigation ---------- */

.pysvorbit-nav {
    width: 100%;
    border-bottom: 1px solid var(--border);
    padding: 12px 0;
    margin-bottom: 10px;
}

.pysvorbit-nav-inner {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 12px;
}

.pysvorbit-brand {
    color: var(--text);
    font-weight: 800;
    letter-spacing: .2px;
    font-size: 16px;
}

.pysvorbit-badge {
    display: inline-block;
    margin-left: 8px;
    padding: 4px 10px;
    border-radius: 999px;
    border: 1px solid var(--border);
    color: var(--muted) !important;
    background: #ffffff;
    font-size: 12px;
    font-weight: 600;
}

/* ---------- Hero ---------- */

.pysvorbit-hero {
    padding: 18px 0 20px;
}

.pysvorbit-hero h1 {
    margin: 0 0 8px 0;
    font-size: 32px;
    color: var(--text) !important;
}

/* ---------- Cards ---------- */

div[data-testid="stVerticalBlockBorderWrapper"] {
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 16px;
    box-shadow: 0 6px 20px rgba(15, 23, 42, 0.06);
}

/* ---------- Step numbers ---------- */

.step-num {
    display: inline-flex;
    width: 26px;
    height: 26px;
    align-items: center;
    justify-content: center;
    border-radius: 50%;
    background: linear-gradient(
        135deg,
        var(--accent),
        var(--accent2)
    );
    color: white !important;
    font-size: 13px;
    font-weight: 800;
    margin-right: 8px;
    vertical-align: middle;
}

/* ---------- Callout ---------- */

.callout {
    border: 1px solid var(--border);
    border-left: 4px solid var(--accent);
    border-radius: 10px;
    padding: 10px 14px;
    background: rgba(109, 74, 255, 0.06);
    color: var(--text);
    font-size: 14px;
    margin: 12px 0;
}

/* ---------- Radio buttons ---------- */

div[data-testid="stRadio"] label {
    color: var(--text) !important;
}

div[data-testid="stRadio"] p {
    color: var(--text) !important;
}

/* ---------- Select boxes ---------- */

div[data-baseweb="select"] > div {
    background: #ffffff !important;
    border: 1px solid var(--border) !important;
    border-radius: 12px;
}

div[data-baseweb="select"] span {
    color: var(--text) !important;
}

div[data-baseweb="select"] input {
    color: var(--text) !important;
}

div[data-baseweb="select"] svg {
    fill: var(--text) !important;
}

/* ---------- Selectbox dropdown menu ---------- */

div[role="listbox"] {
    background: #ffffff !important;
    border: 1px solid var(--border) !important;
}

div[role="option"] {
    background: #ffffff !important;
    color: var(--text) !important;
}

div[role="option"]:hover {
    background: #f1efff !important;
    color: var(--text) !important;
}

/* ---------- File uploader ---------- */

section[data-testid="stFileUploaderDropzone"] {
    background: #ffffff !important;
    border: 1px dashed #b8c0cc !important;
    border-radius: 12px;
}

section[data-testid="stFileUploaderDropzone"] * {
    color: var(--text) !important;
}

section[data-testid="stFileUploaderDropzone"] small {
    color: var(--muted) !important;
}

/* Upload button */

section[data-testid="stFileUploaderDropzone"] button {
    background: #ffffff !important;
    color: var(--text) !important;
    border: 1px solid var(--border) !important;
}

/* ---------- Success messages ---------- */

div[data-testid="stAlert"] {
    border-radius: 10px;
}

/* ---------- Buttons ---------- */

div.stButton > button {
    border-radius: 12px;
    border: 1px solid var(--border);
    background: #ffffff;
    color: var(--text);
    font-weight: 700;
    padding: 10px 14px;
}

div.stButton > button:hover {
    border-color: var(--accent);
    color: var(--accent);
}

div[data-testid="stDownloadButton"] > button {
    border-radius: 12px;
    border: 1px solid var(--border);
    background: #ffffff;
    color: var(--text);
    font-weight: 700;
    padding: 10px 14px;
}

div[data-testid="stDownloadButton"] > button:hover {
    border-color: var(--accent);
    color: var(--accent);
}

div[data-testid="stDownloadButton"] > button[kind="primary"] {
    background: linear-gradient(
        135deg,
        var(--accent),
        var(--accent2)
    );
    border-color: transparent;
    color: white !important;
}

/* ---------- Primary button ---------- */

div.stButton > button[kind="primary"] {
    background: linear-gradient(
        135deg,
        var(--accent),
        var(--accent2)
    );
    border-color: transparent;
    color: white !important;
}

div.stButton > button[kind="primary"]:hover {
    color: white !important;
    opacity: 0.92;
}

/* ---------- Divider ---------- */

hr {
    border-color: var(--border);
}

/* ---------- Footer ---------- */

.pysvorbit-footer {
    color: var(--muted);
    font-size: 12px;
    padding: 24px 0 10px;
    line-height: 1.5;
}

    </style>
    """,
    unsafe_allow_html=True,
)


# ---------------------------------------------------------
# Application paths
# ---------------------------------------------------------

APP_DIR = Path(__file__).resolve().parent
EXAMPLES_DIR = APP_DIR / "static" / "examples"


# ---------------------------------------------------------
# Sample files
# ---------------------------------------------------------

SAMPLE_FILES = {
    "Synthetic Visual Binary": "example_visual_binary.inp",
    "GL765": "GL765_Test1.inp",
    "HD New": "HD_New.inp",
    "HIP12780": "HIP12780 _FIN379_Parallax1.inp",
    "HIP51360": "HIP51360.inp",
    "HIP53206": "HIP53206.inp",
    "HIP72217 - Parallax 1": "HIP72217_parallax1.inp",
    "HIP72217 - Parallax 2": "HIP72217_parallax2.inp",
}


# ---------------------------------------------------------
# Navigation
# ---------------------------------------------------------

st.markdown(
    """
    <div class="pysvorbit-nav">
        <div class="pysvorbit-nav-inner">
            <div class="pysvorbit-brand">
                🛰️ Orbit Analysis for Visual & Spectroscopic Binaries
                <span class="pysvorbit-badge">RV + Visual</span>
            </div>
        </div>
    </div>
    """,
    unsafe_allow_html=True,
)


def start_new_fit():
    """Clear the current Streamlit run and return to the input form."""
    for state_key in (
        "summary",
        "run_dir",
        "uploaded_path",
        "fixel_map",
        "pdf_bytes",
        "plot_style_signature",
    ):
        st.session_state.pop(state_key, None)

    for state_key in list(st.session_state):
        if state_key.startswith(("point_", "fit_", "plot_")):
            del st.session_state[state_key]


def render_plot_card(item, plot_index):
    """Render one Flask-style plot card and return its customized figure data."""
    marker_symbols = ["circle", "square", "diamond", "triangle-up", "star"]
    fig = go.Figure(json.loads(item["fig_json"]))
    editable_trace_indexes = [
        trace_index
        for trace_index, trace in enumerate(fig.data)
        if "markers" in (getattr(trace, "mode", None) or "")
        and (
            getattr(trace, "meta", None) is None
            or not isinstance(getattr(trace, "meta", None), dict)
            or trace.meta.get("marker_editable", True)
        )
    ]

    style_values = (plot_index, None, None, None)
    if editable_trace_indexes:
        control_col1, control_col2, control_col3 = st.columns([1, 1.35, 1.2])
        color_key = f"plot_color_{plot_index}"
        size_key = f"plot_size_{plot_index}"
        symbol_key = f"plot_symbol_{plot_index}"

        with control_col1:
            color = st.color_picker(
                "Color",
                value=st.session_state.get(color_key, "#7c5cff"),
                key=color_key,
            )
        with control_col2:
            size = st.slider(
                "Size",
                min_value=3,
                max_value=18,
                value=st.session_state.get(size_key, 8),
                key=size_key,
            )
        with control_col3:
            selected_symbol = st.session_state.get(symbol_key, "circle")
            symbol = st.selectbox(
                "Symbol",
                marker_symbols,
                index=(
                    marker_symbols.index(selected_symbol)
                    if selected_symbol in marker_symbols
                    else 0
                ),
                key=symbol_key,
            )

        style_values = (plot_index, color, size, symbol)
        for trace_index in editable_trace_indexes:
            fig.data[trace_index].marker.color = color
            fig.data[trace_index].marker.size = size
            fig.data[trace_index].marker.symbol = symbol

    return {
        **item,
        "fig_json": fig.to_json(),
        "style_values": style_values,
    }


if "summary" not in st.session_state:
    if st.button("New Fit", key="new_fit"):
        start_new_fit()
        st.rerun()


# ---------------------------------------------------------
# Hero
# ---------------------------------------------------------

st.markdown(
    """
    <div class="pysvorbit-hero">
        <h1>Binary-Star Orbit Fitting</h1>
        <div class="pysvorbit-subtitle">
            This tool works out how two stars orbit each other from
            your observations — and, if you provide a parallax,
            how much each star weighs. New to this? Start with the
            example file below and the glossary in step 2.
        </div>
    </div>
    """,
    unsafe_allow_html=True,
)

# ---------------------------------------------------------
# Input selection
# ---------------------------------------------------------

# ---------------------------------------------------------
# Input selection
# ---------------------------------------------------------

left_col, right_col = st.columns(2, gap="medium")

with left_col:
    with st.container(border=True):

        st.markdown(
            '<h2><span class="step-num">1</span>Upload your data</h2>',
            unsafe_allow_html=True,
        )

        st.markdown(
            '<p class="pysvorbit-subtitle">'
            "Choose one of the provided sample files or upload your own data file."
            "</p>",
            unsafe_allow_html=True,
        )

        st.markdown(
            '<div class="callout">'
            "<strong>Tip:</strong> Start with a sample file if you are new to PySVOrbit."
            "</div>",
            unsafe_allow_html=True,
        )

        input_mode = st.radio(
            "Choose input method:",
            ["Sample file", "Upload file"],
            horizontal=True,
        )

        # Default values
        selected_path = None
        selected_filename = None
        sample_name = None
        uploaded_file = None

        # -------------------------------------------------
        # Sample file
        # -------------------------------------------------

        if input_mode == "Sample file":

            sample_name = st.selectbox(
                "Select a sample file:",
                list(SAMPLE_FILES.keys()),
            )

            selected_filename = SAMPLE_FILES[sample_name]
            selected_path = EXAMPLES_DIR / selected_filename

            if selected_path.exists():
                st.success(f"Selected: {selected_filename}")
            else:
                st.error(f"Sample file not found: {selected_filename}")

        # -------------------------------------------------
        # User-uploaded file
        # -------------------------------------------------

        else:

            uploaded_file = st.file_uploader(
                "Upload your input file",
                type=["inp", "csv"],
            )

            if uploaded_file is not None:
                st.success(f"Uploaded: {uploaded_file.name}")


with right_col:
    with st.container(border=True):

        st.markdown(
            "<h2>What you'll get</h2>",
            unsafe_allow_html=True,
        )

        st.markdown(
            """
            <ul>
                <li>Orbital elements and fitted parameters</li>
                <li>Radial velocity and orbit plots</li>
                <li>Residuals and diagnostic information</li>
                <li>Derived physical quantities</li>
                <li>Downloadable results and PDF report</li>
            </ul>
            """,
            unsafe_allow_html=True,
        )


# ---------------------------------------------------------
# Sample-file selection
# ---------------------------------------------------------


# ---------------------------------------------------------
# File upload
# ---------------------------------------------------------


# ---------------------------------------------------------
# Temporary test button
# ---------------------------------------------------------
# ---------------------------------------------------------
# Choose what to fit
# ---------------------------------------------------------

st.divider()

with st.container(border=True):

    st.markdown(
        '<h2><span class="step-num">2</span>Choose what to fit</h2>',
        unsafe_allow_html=True,
    )

    st.markdown(
        '<p class="pysvorbit-subtitle">'
        'Each orbital element below can be <b>fit</b> '
        '(checked — the tool solves for its best value) or '
        '<b>fixed</b> (unchecked — kept at the value from your file). '
        'Tip: if your file has <b>no radial velocities</b>, fix '
        '<code>K1</code>, <code>K2</code> and <code>V0</code>, '
        'since there is nothing to constrain them.'
        '</p>',
        unsafe_allow_html=True,
    )

    # Get the orbital elements used by the existing backend
    elnames = list(getattr(backend, "orb").elname)

    # Store the user's fit/fix choices
    fit_choices = {}

    # Display the orbital-element checkboxes in columns
    fit_columns = st.columns(4, gap="medium")

    for index, name in enumerate(elnames):

        with fit_columns[index % 4]:

            fit_choices[name] = st.checkbox(
                name,
                value=True,
                help=ELEMENT_GLOSSARY.get(name, ""),
                key=f"fit_{name}",
            )

            short_label = ELEMENT_SHORT_LABELS.get(name, "")

            if short_label:
                st.caption(short_label)
st.divider()


if st.button("Run", type="primary"):

    # -----------------------------------------------------
    # Check that an input file has been selected
    # -----------------------------------------------------

    if input_mode == "Sample file":

        selected_filename = SAMPLE_FILES[sample_name]
        source_path = EXAMPLES_DIR / selected_filename

        if not source_path.exists():
            st.error(f"Sample file not found: {selected_filename}")
            st.stop()

    else:

        if uploaded_file is None:
            st.warning("Please select or upload an input file.")
            st.stop()

    # -----------------------------------------------------
    # Create a temporary run directory
    # -----------------------------------------------------

    run_dir = Path(tempfile.mkdtemp(prefix="pysvorbit_"))

    try:

        # -------------------------------------------------
        # Prepare the input file
        # -------------------------------------------------

        if input_mode == "Sample file":

            uploaded_path = run_dir / selected_filename

            shutil.copy2(source_path, uploaded_path)

        else:

            uploaded_path = run_dir / uploaded_file.name

            with open(uploaded_path, "wb") as f:
                f.write(uploaded_file.getbuffer())

        # -------------------------------------------------
        # Initial fit/fix settings
        # -------------------------------------------------

        fixel_map = {
            name: 1 if fit_choices[name] else 0
            for name in elnames
        }

        # -------------------------------------------------
        # Run the existing PySVOrbit fitting code
        # -------------------------------------------------

        with st.spinner("Running orbit fitting..."):

            summary = run_fit_to_dir(
                run_dir=run_dir,
                uploaded_path=uploaded_path,
                fixel_map=fixel_map,
                point_filter=None,
            )

        st.success("Orbit fitting completed successfully.")

        # Keep the current run available across Streamlit reruns.
        st.session_state["summary"] = summary
        st.session_state["run_dir"] = str(run_dir)
        st.session_state["uploaded_path"] = str(uploaded_path)
        st.session_state["fixel_map"] = fixel_map
        st.session_state.pop("pdf_bytes", None)
        for state_key in list(st.session_state):
            if state_key.startswith("point_"):
                del st.session_state[state_key]

        st.write(f"Results generated in: `{run_dir}`")

    except Exception as exc:

        st.error("Orbit fitting failed.")
        st.exception(exc)

# ---------------------------------------------------------
# Display results
# ---------------------------------------------------------

if "summary" in st.session_state:

    summary = st.session_state["summary"]
    run_dir = Path(st.session_state["run_dir"])

    try:

        st.divider()
        st.header("Fit Results")
        st.caption(f"Run ID: {summary.get('run_id', 'Unknown')}")

        # -------------------------------------------------
        # Results summary
        # -------------------------------------------------

        stats = summary.get("stats", {})
        elements = summary.get("elements", [])
        derived = summary.get("derived", {})
        point_tables = summary.get("point_tables", {})

        with st.container(border=True):
            st.subheader("Summary")

            st.markdown(
                f"**Object:** {stats.get('name') or 'Unknown'} &nbsp;•&nbsp; "
                f"**Positions:** {stats.get('npos', 0)} &nbsp;•&nbsp; "
                f"**RV1:** {stats.get('nrv1', 0)} &nbsp;•&nbsp; "
                f"**RV2:** {stats.get('nrv2', 0)} &nbsp;•&nbsp; "
                f"**Parallax:** {stats.get('parallax', 0.0)} mas &nbsp;•&nbsp; "
                f"**Chi²:** {stats.get('chi2', 0.0):.4g}"
            )

            action_col1, action_col2, action_col3 = st.columns([1, 1.35, 0.8])

            with action_col1:
                output_csv = summary.get("output_csv")
                output_csv_path = run_dir / output_csv if output_csv else None
                if output_csv_path and output_csv_path.exists():
                    st.download_button(
                        "Download CSV",
                        data=output_csv_path.read_bytes(),
                        file_name=output_csv_path.name,
                        key="summary_download_csv",
                        type="primary",
                    )

            with action_col2:
                if "pdf_bytes" not in st.session_state:
                    if st.button(
                        "Generate PDF Report",
                        key="summary_generate_pdf",
                        type="primary",
                    ):
                        try:
                            st.session_state["pdf_bytes"] = build_streamlit_pdf(
                                summary,
                                run_dir,
                            )
                            st.rerun()
                        except Exception as exc:
                            st.error("Could not generate the PDF report.")
                            st.exception(exc)
                else:
                    st.download_button(
                        "Download PDF Report",
                        data=st.session_state["pdf_bytes"],
                        file_name=f"{summary.get('run_id', 'pysvorbit')}_report.pdf",
                        mime="application/pdf",
                        key="summary_download_pdf",
                        type="primary",
                    )

            with action_col3:
                if st.button("New Fit", key="new_fit_results"):
                    start_new_fit()
                    st.rerun()

            with st.expander("❔ How do I tell if this is a good fit?"):
                st.markdown(
                    "Look at **reduced chi² per data type** below: values close to "
                    "**1** mean the fit matches your data about as well as your "
                    "quoted measurement errors would predict. Also check **RMS** "
                    "(root-mean-square residual), which shows how far the fit "
                    "typically misses each data point in real units."
                )

        st.markdown(
            """
            <div class="pysvorbit-subtitle">
                The following results were calculated from the selected input data
                using the PySVOrbit orbit-fitting backend.
            </div>
            """,
            unsafe_allow_html=True,
        )

        # -------------------------------------------------
        # Orbital elements
        # -------------------------------------------------

        if elements:
            st.subheader("Orbital Elements")
            st.dataframe(elements, width="stretch", hide_index=True)

        # -------------------------------------------------
        # Fit statistics
        # -------------------------------------------------

        with st.container(border=True):

            st.markdown(
                "<h2>Fit Statistics</h2>",
                unsafe_allow_html=True,
            )

            stat_col1, stat_col2, stat_col3, stat_col4 = st.columns(4)

            with stat_col1:
                st.metric(
                    "Data points",
                    stats.get("npos", 0),
                )

            with stat_col2:
                st.metric(
                    "RV1 points",
                    stats.get("nrv1", 0),
                )

            with stat_col3:
                st.metric(
                    "RV2 points",
                    stats.get("nrv2", 0),
                )

            with stat_col4:
                st.metric(
                    "Parallax (mas)",
                    f"{stats.get('parallax', 0.0):.6g}",
                )

            st.markdown("### Goodness of fit")

            statistics_rows = []

            if stats.get("npos", 0) > 0:
                statistics_rows.append(
                    {
                        "Data type": "Position angle (θ)",
                        "Reduced χ²": f"{stats.get('chi2n_theta', 0.0):.3g}",
                        "RMS residual": f"{stats.get('rms_theta', 0.0):.4g}°",
                    }
                )

                statistics_rows.append(
                    {
                        "Data type": "Separation (ρ)",
                        "Reduced χ²": f"{stats.get('chi2n_rho', 0.0):.3g}",
                        "RMS residual": f"{stats.get('rms_rho', 0.0):.4g}″",
                    }
                )

            if stats.get("nrv1", 0) > 0:
                statistics_rows.append(
                    {
                        "Data type": "RV1 (primary)",
                        "Reduced χ²": f"{stats.get('chi2n_rv1', 0.0):.3g}",
                        "RMS residual": f"{stats.get('rms_rv1', 0.0):.4g} km/s",
                    }
                )

            if stats.get("nrv2", 0) > 0:
                statistics_rows.append(
                    {
                        "Data type": "RV2 (secondary)",
                        "Reduced χ²": f"{stats.get('chi2n_rv2', 0.0):.3g}",
                        "RMS residual": f"{stats.get('rms_rv2', 0.0):.4g} km/s",
                    }
                )

            if statistics_rows:
                st.table(statistics_rows)
            else:
                st.info("No fit statistics are available.")

        # -------------------------------------------------
        # Derived quantities
        # -------------------------------------------------

        if derived:
            st.subheader("Derived Quantities")
            derived_rows = [
                {
                    "Quantity": key,
                    "Value": f"{value:.6g}"
                    if isinstance(value, (int, float))
                    else str(value),
                }
                for key, value in derived.items()
                if key != "warnings"
            ]
            if derived_rows:
                st.dataframe(derived_rows, width="stretch", hide_index=True)
            for warning in derived.get("warnings", []):
                st.warning(warning)

        # -------------------------------------------------
        # Recalculate derived quantities
        # -------------------------------------------------

        with st.container(border=True):
            st.subheader("Recalculate Derived Quantities")

            recompute_options = [""] + list(DERIVED_QUANTITY_LABELS)
            override_key = st.selectbox(
                "Quantity to use as the known value:",
                recompute_options,
                format_func=lambda key: (
                    "Use parallax from file"
                    if not key
                    else DERIVED_QUANTITY_LABELS[key]
                ),
                key="recompute_key",
            )
            override_raw = st.text_input(
                "Value:",
                placeholder="Enter a positive number",
                key="recompute_value",
            )

            if st.button("Recalculate", key="recalculate_derived"):
                override_value = None

                if override_key:
                    if not override_raw.strip():
                        st.warning("Enter a value to use for that quantity.")
                    else:
                        try:
                            override_value = float(override_raw.strip())
                        except ValueError:
                            st.warning("That doesn't look like a number.")
                        else:
                            if override_value <= 0:
                                st.warning("Please enter a positive number.")
                                override_value = None

                if not override_key or override_value is not None:
                    current_summary = st.session_state["summary"]
                    recalculated_derived = compute_derived(
                        current_summary.get("elements", []),
                        current_summary.get("stats", {}),
                        override_key or None,
                        override_value,
                    )
                    updated_summary = dict(current_summary)
                    updated_summary["derived"] = recalculated_derived
                    st.session_state["summary"] = updated_summary
                    st.rerun()

        # -------------------------------------------------
        # Point tables
        # -------------------------------------------------

        point_table_specs = {
            "pos": (
                "Positions",
                ["Epoch", "PA θ (°)", "Sep ρ (arcsec)", "Err", "Source", "Δθ", "Δρ"],
            ),
            "rv1": (
                "RV1 (primary)",
                ["Epoch", "RV (km/s)", "Err", "Source", "ΔRV"],
            ),
            "rv2": (
                "RV2 (secondary)",
                ["Epoch", "RV (km/s)", "Err", "Source", "ΔRV"],
            ),
        }

        for table_name, rows in point_tables.items():
            if rows:
                for row in rows:
                    state_key = f"point_{table_name}_{row['idx']}"
                    if state_key not in st.session_state:
                        st.session_state[state_key] = row.get("included", True)

        for table_name, (table_title, headers) in point_table_specs.items():
            rows = point_tables.get(table_name, [])
            if not rows:
                continue

            st.subheader(f"{table_title} ({len(rows)})")
            select_col, none_col = st.columns(2)
            with select_col:
                if st.button("Select All", key=f"select_all_{table_name}"):
                    for row in rows:
                        st.session_state[f"point_{table_name}_{row['idx']}"] = True
                    st.rerun()
            with none_col:
                if st.button("Select None", key=f"select_none_{table_name}"):
                    for row in rows:
                        st.session_state[f"point_{table_name}_{row['idx']}"] = False
                    st.rerun()

            header_columns = st.columns(len(headers) + 1)
            header_columns[0].caption("Keep")
            for column, header in zip(header_columns[1:], headers):
                column.caption(header)

            for row in rows:
                row_columns = st.columns(len(headers) + 1)
                row_columns[0].checkbox(
                    "Keep point",
                    key=f"point_{table_name}_{row['idx']}",
                    label_visibility="collapsed",
                )
                if table_name == "pos":
                    values = (
                        f"{row['epoch']:.4f}",
                        f"{row['theta']:.4f}",
                        f"{row['rho']:.4f}",
                        f"{row['err']:.4g}",
                        row["source"],
                        f"{row['dtheta']:.4f}",
                        f"{row['drho']:.4f}",
                    )
                else:
                    values = (
                        f"{row['epoch']:.4f}",
                        f"{row['rv']:.4f}",
                        f"{row['err']:.4g}",
                        row["source"],
                        f"{row['drv']:.4f}",
                    )
                for column, value in zip(row_columns[1:], values):
                    column.write(value)

        # -------------------------------------------------
        # Refit with selected points
        # -------------------------------------------------

        if st.button("Refit", type="primary", key="refit"):
            point_filter = {
                table_name: [
                    row["idx"]
                    for row in point_tables.get(table_name, [])
                    if st.session_state.get(
                        f"point_{table_name}_{row['idx']}",
                        row.get("included", True),
                    )
                ]
                for table_name in ("pos", "rv1", "rv2")
            }

            if not any(point_filter.values()):
                st.warning("Select at least one observation before refitting.")
            else:
                try:
                    with st.spinner("Refitting with selected points..."):
                        refit_summary = run_fit_to_dir(
                            run_dir=Path(st.session_state["run_dir"]),
                            uploaded_path=Path(st.session_state["uploaded_path"]),
                            fixel_map=st.session_state["fixel_map"],
                            point_filter=point_filter,
                        )
                    st.session_state["summary"] = refit_summary
                    st.session_state["run_dir"] = str(run_dir)
                    st.session_state["uploaded_path"] = str(
                        st.session_state["uploaded_path"]
                    )
                    st.session_state["fixel_map"] = st.session_state["fixel_map"]
                    st.session_state.pop("pdf_bytes", None)
                    st.rerun()
                except Exception as exc:
                    st.error("Refit failed.")
                    st.exception(exc)

        # -------------------------------------------------
        # Display Plotly figures generated by PySVOrbit
        # -------------------------------------------------

        plotly_figs = summary.get("plotly_figs", [])

        if plotly_figs:
            st.subheader("Interactive Plots")
            st.caption("Zoom/pan and adjust marker style. Click a plot to enlarge.")
            customized_plotly_figs = []

            for row_start in range(0, len(plotly_figs), 2):
                plot_columns = st.columns(2, gap="medium")
                for column, (plot_index, item) in zip(
                    plot_columns,
                    enumerate(plotly_figs[row_start : row_start + 2], start=row_start),
                ):
                    with column:
                        with st.container(border=True):
                            title = item.get("title", "Plot")
                            st.markdown(f"### {title}")
                            customized_item = render_plot_card(item, plot_index)
                            customized_plotly_figs.append(customized_item)
                            fig = go.Figure(
                                json.loads(customized_item["fig_json"])
                            )
                            st.plotly_chart(fig, width="stretch")

            style_signature = tuple(
                item["style_values"] for item in customized_plotly_figs
            )
            if st.session_state.get("plot_style_signature") != style_signature:
                st.session_state["plot_style_signature"] = style_signature
                st.session_state.pop("pdf_bytes", None)

            summary["plotly_figs"] = customized_plotly_figs
            st.session_state["summary"] = summary
        else:
            st.info("No Plotly figures were found.")

    except Exception as exc:
        st.error("Could not display the generated results.")
        st.exception(exc)

else:
    st.info("Run a fit to display results.")
