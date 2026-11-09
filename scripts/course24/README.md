# Reproducible 24-system binary-star course builder

This directory contains the **source-driven, review-only** INT4 + ORB6 + SBX + SOAR dataset generator. All required Python modules are included here; no separate checkout of BinaryStarCatalogFusion is needed. Raw catalogues are **not** bundled or redistributed.

## Regenerate from local catalogues

From the PySVOrbit repository root:

```bash
BASE="/media/talafha/Disk_1/BinaryStarCatalogFusion"
python3 scripts/course24/build_24_course_orbits.py \
  --int4 "$BASE/catalog_run/data/int4_all_downloaded.dat" \
  --orb6 "$BASE/catalog_run/data/orb6orbits.txt" \
  --sb9-dir "$BASE/catalog_run/data/sbx_catalog_official" \
  --soar-dir "$BASE/catalog_run/data/soar/2015_2025" \
  --count 48 \
  --output "$BASE/course_48_from_main"

python3 scripts/course24/rebuild_revised_from_48.py \
  --pool "$BASE/course_48_from_main" \
  --output "$BASE/course_24_revised_from_main"
```

The 24 resulting `.inp` candidates are in `course_24_revised_from_main/fit_inputs`. Choose **new, nonexistent** output directories: these tools intentionally refuse to overwrite existing results.

The revised roster replaces Projects 1, 3, 7 and 20 using candidate-pool indices 25, 28, 29 and 38 respectively. Project 12 has the published periastron epoch fixed (`*T`) for a conditional nearly-circular fit. These choices require the same catalogue versions as the previously verified 48-system pool.

## Scientific safeguards and limits

- All observations and available catalogue uncertainties retain their provenance; missing uncertainties are never manufactured.
- Uncertain source rows remain in audit files but are excluded from the weighted `.inp` fit.
- Parallax and spectroscopic subsystem identities require verification; numerical acceptance is not scientific certification.
- SOAR internal errors are not automatically treated as isotropic separation uncertainties.
- These are **review candidates**, not publication-certified student assignments. The PySVOrbit backend and UI are unchanged by this addition.

Regression import check: `python -m pytest tests/test_course24_builder_imports.py`.
