#!/usr/bin/env python3
"""Split the USNO/GSU Fourth Interferometric Catalog (INT4) into star CSVs.

Only fields contained in INT4 are exported; the .inp file, if supplied,
provides search identifiers ONLY and is never used as a measurement source.

Format reference: https://www.astro.gsu.edu/wds/int4/format.html

Examples
--------
# All catalogued systems, one CSV per system inside a ZIP file:
python int4_to_star_csv.py

# Only the example binary MLR 224 / GJ 765.2:
python int4_to_star_csv.py --inp 'GL765_Test1(2).inp' --output gl765_int4.zip

# Export as ordinary CSV files, filtering by several names:
python int4_to_star_csv.py --star 'MLR 224' --star 'HIP 96656' --output-dir stars

# Keep only systems with at least four catalogue observations:
python int4_to_star_csv.py --min-points 4 --output-dir stars_4plus

# Keep only systems with at least four resolved astrometric measurements:
python int4_to_star_csv.py --resolved-only --min-points 4 --output-dir resolved_4plus

# Use a catalog file already downloaded, keep all INT4 columns:
python int4_to_star_csv.py --catalog int4_all.txt.gz --full --output all_int4.zip
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import re
import sys
import tempfile
import urllib.error
import urllib.request
import zipfile
from collections import Counter
from pathlib import Path
from typing import Iterator

# A compressed file is preferable to the ~44 MB uncompressed distribution.
CATALOG_URLS = (
    'https://www.astro.gsu.edu/wds/int4/int4_all.txt.gz',
    'https://crf.usno.navy.mil/data_products/WDS/int4/int4_all.txt.gz',
    'https://www.astro.gsu.edu/wds/int4/int4_all.txt',
    'https://crf.usno.navy.mil/data_products/WDS/int4/int4_all.txt',
)

# INT4 fixed-column positions are 1-based in the official documentation;
# Python indices here are 0-based and END-EXCLUSIVE.
ID_COLUMNS = {
    'coordinates_j2000': (0, 18),
    'name_1': (20, 46),
    'name_2': (46, 72),
    'hd_dm': (72, 85),
    'other_catalog': (85, 104),
    'wds': (104, 114),
    'general_flag': (115, 116),
    'orbit_flag': (117, 118),
}

OBS_COLUMNS = {
    'epoch_flag': (1, 2),
    'epoch_besselian': (2, 11),
    'theta_flag': (13, 14),
    'theta_deg': (14, 21),
    'theta_error_flag': (22, 23),
    'theta_error_deg': (23, 28),
    'rho_flag': (28, 29),
    'rho_original': (29, 39),
    'rho_error_flag': (40, 41),
    'rho_error_original': (41, 49),
    'mag_primary_flag': (50, 51),
    'mag_primary': (51, 57),
    'mag_primary_error_flag': (58, 59),
    'mag_primary_error': (59, 64),
    'mag_secondary_flag': (65, 66),
    'mag_secondary_or_delta': (66, 72),
    'mag_secondary_error_flag': (73, 74),
    'mag_secondary_error': (74, 79),
    'filter_wavelength': (82, 86),
    'filter_fwhm': (86, 90),
    'filter_flag': (90, 91),
    'telescope_aperture': (92, 96),
    'telescope_aperture_code': (96, 97),
    'number_of_nights': (98, 100),
    'reference': (102, 110),
    'technique': (111, 114),
}

# The first nine columns resemble the visual/interferometric section in .inp.
COMPACT_COLUMNS = [
    'epoch_besselian', 'theta_deg', 'rho_arcsec',
    'theta_error_deg', 'rho_error_arcsec',
    'theta_flag', 'rho_flag', 'technique', 'reference',
]
FULL_COLUMNS = (
    list(ID_COLUMNS)
    + list(OBS_COLUMNS)
    + ['rho_arcsec', 'rho_error_arcsec']
)

COORD_RE = re.compile(r'^\d{6}\.\d{2}[+-]\d{6}\.\d')
WDS_RE = re.compile(r'^\d{5}[+-]\d{4}$')
MISSING = {'', '.'}


def clean(value: str) -> str:
    value = value.strip()
    return '' if value in MISSING else value


def fixed_fields(line: str, positions: dict[str, tuple[int, int]]) -> dict[str, str]:
    return {key: clean(line[start:end]) for key, (start, end) in positions.items()}


def normalize_name(value: str) -> str:
    """Case/space/punctuation insensitive; GL -> GJ for Gliese catalogue IDs."""
    value = re.sub(r'[^A-Z0-9]', '', value.upper())
    # Gliese, GL and GJ are often alternative prefixes for the same object.
    if value.startswith('GL') and re.match(r'^GL\d', value):
        value = 'GJ' + value[2:]
    return value


def identifiers_from_inp(path: Path) -> list[str]:
    """Extract names from an input's Object: line; no numerical data imported."""
    found = []
    with path.open(encoding='utf-8', errors='replace') as stream:
        for line in stream:
            if line.lower().startswith('object:'):
                found.extend(part.strip() for part in line.split(':', 1)[1].split('='))
    if not found:
        raise ValueError(f'No "Object:" line found in {path}')
    return [item for item in found if item]


def match_system(ids: dict[str, str], wanted: set[str]) -> bool:
    if not wanted:
        return True
    names = [ids[key] for key in ('name_1', 'name_2', 'hd_dm', 'other_catalog', 'wds')]
    # Other_catalog may be 'HIP     96656', WDS may be '19391+7625'.
    keys = {normalize_name(name) for name in names if name}
    keys |= {'WDS' + normalize_name(ids['wds'])} if ids['wds'] else set()
    return bool(keys & wanted)


def is_identification(line: str) -> bool:
    # Every normal ID record begins with the precise epoch-2000 coordinates.
    # The fallback handles catalog entries with unspecified precise coordinates.
    return bool(COORD_RE.match(line) or
                (len(line) >= 114 and WDS_RE.fullmatch(line[104:114].strip())))


def read_catalog(path: Path) -> Iterator[tuple[dict[str, str], list[dict[str, str]]]]:
    """Stream system blocks from the fixed-width ASCII file or .gz file."""
    with path.open('rb') as f:
        compressed = f.read(2) == b'\x1f\x8b'
    opener = gzip.open if compressed else open
    with opener(path, mode='rt', encoding='latin-1', errors='replace') as f:
        system: dict[str, str] | None = None
        observations: list[dict[str, str]] = []
        for raw in f:
            line = raw.rstrip('\r\n')
            if not line.strip():
                continue
            # Allow occasional leading whitespace on identifier records.
            candidate = line.lstrip() if COORD_RE.match(line.lstrip()) else line
            if is_identification(candidate):
                if system is not None:
                    yield system, observations
                system = fixed_fields(candidate, ID_COLUMNS)
                observations = []
                continue
            if system is None:
                continue  # Any title/comment before the first ID record.
            row = fixed_fields(line, OBS_COLUMNS)
            # A real observation line has a numeric Besselian-year epoch.
            try:
                float(row['epoch_besselian'])
            except ValueError:
                continue
            observations.append(row)
        if system is not None:
            yield system, observations


def arcseconds(raw: str, flag: str) -> str:
    """Normalize unusual INT4 separation units; don't fabricate absent data."""
    if not raw:
        return ''
    factors = {'m': 1e-3, 'M': 60.0, 'D': 3600.0}
    if flag not in factors:
        return raw  # In all other cases the documented unit is arcseconds.
    try:
        return format(float(raw) * factors[flag], '.12g')
    except ValueError:
        return ''


def make_row(ids: dict[str, str], obs: dict[str, str], full: bool) -> dict[str, str]:
    row = dict(ids) if full else {}
    row.update(obs)
    row['rho_arcsec'] = arcseconds(obs['rho_original'], obs['rho_flag'])
    row['rho_error_arcsec'] = arcseconds(obs['rho_error_original'], obs['rho_flag'])
    return row


def filename_for(ids: dict[str, str], used: Counter[str]) -> str:
    # WDS + discoverer/pair is usually stable and more descriptive than RA.
    prefix = ids['wds'] or ids['coordinates_j2000'] or 'unknown'
    suffix = ids['name_2'] or ids['name_1'] or 'star'
    name = f'{prefix}_{suffix}'
    name = re.sub(r'[^A-Za-z0-9_.+-]+', '_', name).strip('._')[:110]
    used[name] += 1
    if used[name] > 1:
        name += f'_{used[name]}'
    return f'{name}.csv'


def to_csv(records: list[dict[str, str]], header: list[str]) -> str:
    out = io.StringIO(newline='')
    writer = csv.DictWriter(out, fieldnames=header, extrasaction='ignore', lineterminator='\n')
    writer.writeheader()
    writer.writerows(records)
    return out.getvalue()


def download_catalog(destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    errors = []
    for url in CATALOG_URLS:
        tmp = None
        print(f'Downloading {url}', file=sys.stderr)
        try:
            request = urllib.request.Request(url, headers={
                'User-Agent': 'Mozilla/5.0 (compatible; INT4-scientific-data-export/1.0)'
            })
            with urllib.request.urlopen(request, timeout=90) as source, \
                    tempfile.NamedTemporaryFile(mode='w+b', dir=destination.parent,
                                                delete=False) as temporary:
                tmp = Path(temporary.name)
                first_bytes = source.read(512)
                if (b'<!doctype html' in first_bytes.lower()
                        or b'<html' in first_bytes.lower()):
                    raise ValueError('Received an HTML page, not the INT4 data file')
                temporary.write(first_bytes)
                while True:
                    chunk = source.read(1024 * 1024)
                    if not chunk:
                        break
                    temporary.write(chunk)
            if tmp.stat().st_size < 1000:
                tmp.unlink(missing_ok=True)
                raise ValueError('Downloaded content is unexpectedly small')
            tmp.replace(destination)
            print(f'Saved catalogue ({destination.stat().st_size:,} bytes): {destination}',
                  file=sys.stderr)
            return
        except (OSError, ValueError, urllib.error.URLError) as exc:
            errors.append(f'{url}: {exc}')
            if tmp is not None:
                tmp.unlink(missing_ok=True)
    raise RuntimeError('Could not download INT4. Download int4_all.txt manually '
                       'and use --catalog PATH. Details:\n' + '\n'.join(errors))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--catalog', type=Path, help='Local INT4 .txt or .txt.gz; no download')
    parser.add_argument('--output', type=Path, default=Path('int4_stars.zip'),
                        help='ZIP archive holding a separate CSV per system')
    parser.add_argument('--output-dir', type=Path,
                        help='Write individual CSVs to this directory instead of ZIP')
    parser.add_argument('--star', action='append', default=[],
                        help='Star identifier to export; repeat for multiple stars')
    parser.add_argument('--inp', type=Path,
                        help='Read only star identifiers from an .inp Object: line')
    parser.add_argument('--stars-file', type=Path,
                        help='Text file containing one star identifier per line')
    parser.add_argument('--full', action='store_true',
                        help='Include ALL published INT4 measurement and ID columns')
    parser.add_argument('--resolved-only', action='store_true',
                        help='Keep only rows with numeric theta and separation (exclude limits)')
    parser.add_argument('--min-points', type=int, default=1, metavar='N',
                        help='Export a star only if at least N observation rows remain after filtering (default: 1)')
    args = parser.parse_args()
    if args.min_points < 1:
        parser.error('--min-points must be a positive integer')

    selected = list(args.star)
    if args.inp:
        selected.extend(identifiers_from_inp(args.inp))
    if args.stars_file:
        selected.extend(line.strip() for line in args.stars_file.read_text().splitlines()
                        if line.strip() and not line.lstrip().startswith('#'))
    wanted = {normalize_name(item) for item in selected if normalize_name(item)}
    catalog = args.catalog or Path('int4_all_downloaded.dat')
    if not catalog.exists():
        if args.catalog:
            parser.error(f'Catalog not found: {catalog}')
        download_catalog(catalog)
    headings = FULL_COLUMNS if args.full else COMPACT_COLUMNS
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)

    systems_read = files_written = observations_written = 0
    used: Counter[str] = Counter()
    zip_stream = (zipfile.ZipFile(args.output, 'w', compression=zipfile.ZIP_DEFLATED,
                                  compresslevel=6) if not args.output_dir else None)
    try:
        for ids, measurements in read_catalog(catalog):
            systems_read += 1
            if not match_system(ids, wanted):
                continue
            if args.resolved_only:
                measurements = [m for m in measurements
                                if m['theta_deg'] and m['rho_original']
                                and m['rho_flag'] not in ('<', '>', 'U', 'R', 'G', 'O', 'S', 'V', 'X')]
            # Require N observations per system after the requested row filters.
            if len(measurements) < args.min_points:
                continue
            rows = [make_row(ids, obs, args.full) for obs in measurements]
            name = filename_for(ids, used)
            content = to_csv(rows, headings)
            if zip_stream:
                zip_stream.writestr(name, content)
            else:
                (args.output_dir / name).write_text(content, encoding='utf-8', newline='')
            files_written += 1
            observations_written += len(rows)
    finally:
        if zip_stream:
            zip_stream.close()
    print(f'Systems scanned: {systems_read:,}', file=sys.stderr)
    print(f'CSV files created: {files_written:,}', file=sys.stderr)
    print(f'INT4 measurement rows exported: {observations_written:,}', file=sys.stderr)
    print(f'Output: {args.output_dir or args.output}', file=sys.stderr)
    if wanted and not files_written:
        print('WARNING: No identifiers matched. Try --star "MLR 224" or --star "19391+7625".',
              file=sys.stderr)


if __name__ == '__main__':
    main()
