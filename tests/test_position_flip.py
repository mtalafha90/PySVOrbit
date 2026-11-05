import json
import shutil
from pathlib import Path

import pytest

import app
import backend


def test_flip_is_applied_before_filtering_and_can_be_undone(tmp_path):
    source = Path(__file__).resolve().parents[1] / "static/examples/example_visual_binary.inp"
    uploaded = tmp_path / source.name
    shutil.copy(source, uploaded)
    backend.readinp(str(uploaded))
    original = backend.orb.pos[:, 1].copy()
    keep = list(range(1, len(original)))

    result = app.run_fit_to_dir(tmp_path, uploaded, {}, {"pos": keep}, [0, 1])
    rows = result["point_tables"]["pos"]
    assert rows[0]["theta"] == pytest.approx((original[0] + 180) % 360)
    assert rows[0]["flipped"] and not rows[0]["included"]
    assert rows[1]["theta"] == pytest.approx((original[1] + 180) % 360)
    assert backend.orb.pos[0, 1] == pytest.approx((original[1] + 180) % 360)
    assert json.loads((tmp_path / "summary.json").read_text())["flipped_positions"] == [0, 1]

    result = app.run_fit_to_dir(tmp_path, uploaded, {}, {"pos": keep}, [0])
    assert result["point_tables"]["pos"][1]["theta"] == pytest.approx(original[1])
    assert not result["point_tables"]["pos"][1]["flipped"]


def test_flask_flip_and_refit_preserve_selection_and_flip_state(tmp_path, monkeypatch):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "input.inp").write_text("input")
    data = {"uploaded_filename": "input.inp", "fixel_map": {},
            "flipped_positions": [1], "point_tables": {"pos": [
                {"idx": 0}, {"idx": 1}], "rv1": [], "rv2": []}}
    (run_dir / "summary.json").write_text(json.dumps(data))
    monkeypatch.setattr(app, "RUNS_DIR", tmp_path)
    captured = {}

    def fake_fit(*args, **kwargs):
        captured.update(kwargs)
        return {}

    monkeypatch.setattr(app, "run_fit_to_dir", fake_fit)
    monkeypatch.setattr(app, "render_template", lambda *args, **kwargs: "ok")
    response = app.app.test_client().post("/refit/run", data={"pos_keep_1": "on", "flip_position": "0"})
    assert response.status_code == 200
    assert captured["point_filter"]["pos"] == [1]
    assert captured["flipped_positions"] == [0, 1]
