from pathlib import Path

p = Path('README.md')
t = p.read_text(encoding='utf-8')

old_install = '''## Installation

Python 3.10 or newer is required.

```bash
git clone https://github.com/mtalafha90/PySVOrbit.git
cd PySVOrbit
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```
'''
new_install = '''## Installation

Python 3.10 or newer is required. PySVOrbit now uses standard Python package
metadata in `pyproject.toml`, so the numerical library can be installed directly
from a checkout:

```bash
git clone https://github.com/mtalafha90/PySVOrbit.git
cd PySVOrbit
python -m venv .venv
source .venv/bin/activate
pip install .
```

For the browser interface and PDF/report dependencies, install the optional web
extra instead:

```bash
pip install ".[web]"
```

The historical `requirements.txt` is retained for deployment compatibility,
but `pyproject.toml` is the package metadata used by `pip install .`.
'''
if t.count(old_install) != 1:
    raise RuntimeError(f'installation block count={t.count(old_install)}')
t = t.replace(old_install, new_install, 1)

old_script = '''```python
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
'''
new_script = '''```python
import pysvorbit

pysvorbit.readinp("static/examples/GL765_Test1.inp")
pysvorbit.fitorb()

for name, value, error in zip(pysvorbit.orb.elname,
                              pysvorbit.orb.el,
                              pysvorbit.orb.elerr):
    print(f"{name:>6} = {value:12.6f} +/- {error:.6f}")

pysvorbit.orbsave()   # writes <name>_output.csv to the working directory
```

The legacy `import backend` interface remains available for backward
compatibility. The public `pysvorbit` module is a thin facade over that same
engine, so both names operate on the same module-level orbit object. One fit is
therefore in progress at any one time: call `readinp()` again to start a new
system, and do not run fits concurrently in the same process.
'''
if t.count(old_script) != 1:
    raise RuntimeError(f'scripted-use block count={t.count(old_script)}')
t = t.replace(old_script, new_script, 1)

old_tests = '''```bash
pip install -r requirements.txt pytest
pytest tests/ -v
```
'''
new_tests = '''```bash
pip install ".[test]"
pytest tests/ -v
```
'''
if t.count(old_tests) != 1:
    raise RuntimeError(f'test-install block count={t.count(old_tests)}')
t = t.replace(old_tests, new_tests, 1)

p.write_text(t, encoding='utf-8')
print('Updated README for installable package metadata.')
