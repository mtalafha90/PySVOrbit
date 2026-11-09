"""Ensure the course builder is truly self-contained in this repository."""
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1] / "scripts" / "course24"

@pytest.mark.parametrize("script", ["build_24_course_orbits.py", "rebuild_revised_from_48.py"])
def test_entrypoint_imports_from_repo(script):
    result = subprocess.run(
        [sys.executable, str(ROOT / script), "--help"],
        cwd=ROOT.parents[1],
        text=True, capture_output=True, timeout=20,
    )
    assert result.returncode == 0, result.stderr
    assert "--output" in result.stdout
