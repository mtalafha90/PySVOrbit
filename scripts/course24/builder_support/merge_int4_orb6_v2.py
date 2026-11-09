#!/usr/bin/env python3
"""Match per-system INT4 CSVs with ORB6 elements; make merged CSV and .inp.

Uses the CSVs produced by int4_to_star_csv.py (compact or --full).  No
artificial measurements, parallax, RV amplitudes, or uncertainties are added.
Run with --help for options.  Python 3.10+, standard library only.

Compared with v1 this version provides an audit of equal-grade orbit ties,
optional exclusion of ambiguous orbit selections, and an explicit opt-in
assumed rho uncertainty for INP export. None of these are physical observations.

ORB6 format: https://crf.usno.navy.mil/data_products/WDS/orb6/orb6.html
ORB6 fixed columns: https://crf.usno.navy.mil/showtxthtml.php?fileHead=&source_path=/data_products/WDS/orb6/orb6format.txt
"""
from __future__ import annotations

import argparse
import csv
import math
import io
import zipfile
import re
import sys
import tempfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

ORB6_URLS = (
    'https://www.astro.gsu.edu/wds/orb6/orb6orbits.txt',
    'https://crf.usno.navy.mil/data_products/WDS/orb6/orb6orbits.txt',
)
WDS_RE = re.compile(r'^\d{5}[+-]\d{4}$')
# ORB6 Fortran T positions (column numbers) converted to Python [start:end].
# Official field widths are important: splitting on whitespace DOES NOT work
# for flags and tightly adjacent fields.
ORB6_FIELDS = {
    'coord': (0, 18), 'wds': (19, 29), 'discoverer': (30, 44),
    'ads': (45, 50), 'hd': (51, 57), 'hip': (58, 64),
    'mag1': (66, 72), 'mag2': (73, 79),
    'period': (81, 92), 'period_unit': (92, 93),
    'period_error': (94, 104),
    'semimajor': (105, 114), 'semimajor_unit': (114, 115),
    'semimajor_error': (116, 124),
    'inclination': (125, 133), 'inclination_error': (134, 142),
    'node': (143, 151), 'node_flag': (151, 152),
    'node_error': (153, 161),
    'periastron': (162, 174), 'periastron_unit': (174, 175),
    'periastron_error': (176, 186),
    'eccentricity': (187, 195), 'eccentricity_error': (196, 204),
    'omega': (205, 213), 'omega_flag': (213, 214),
    'omega_error': (214, 222),
    'equinox': (223, 227), 'last_used': (228, 232),
    'grade': (233, 234), 'notes_flag': (235, 236),
    'reference': (237, 245), 'figure': (246, 264),
}
BAD_RHO_FLAGS = {'<', '>', 'U', 'R', 'G', 'O', 'S', 'V', 'X', 'e'}
MERGED_FIELDS = [
    'wds', 'discoverer', 'hip', 'hd', 'orbit_reference', 'orbit_grade',
    'P_years', 'T_besselian', 'e', 'a_arcsec', 'W_deg', 'w_deg', 'i_deg',
    'node_flag', 'omega_flag', 'epoch_besselian', 'theta_deg', 'rho_arcsec',
    'theta_error_deg', 'rho_error_arcsec', 'theta_flag', 'rho_flag',
    'technique', 'reference', 'parallax_mas', 'K1', 'K2', 'V0'
]


def clean(s: str) -> str:
    s = str(s).strip()
    return '' if s in ('', '.') else s


def num(s: str) -> float | None:
    try:
        x = float(clean(s))
        return x if math.isfinite(x) else None
    except (ValueError, TypeError):
        return None


def fmt(x: float | None, places: int = 8) -> str:
    return '' if x is None else f'{x:.{places}f}'.rstrip('0').rstrip('.')


def designation(s: str) -> str:
    return re.sub(r'[^A-Z0-9]', '', s.upper())


def epoch_from_jd(jd: float) -> float:
    # IAU definition of a Besselian year; do not use 365.25 here.
    return 1900.0 + (jd - 2415020.31352) / 365.242198781


def convert_period(period: float | None, unit: str) -> float | None:
    # Empty flags on some older records denote the default years unit.
    units = {'': 1.0, 'y': 1.0, 'd': 1 / 365.242198781,
             'c': 100.0, 'h': 1 / (24 * 365.242198781),
             'm': 1 / (1440 * 365.242198781)}
    return None if period is None or unit not in units else period * units[unit]


def convert_a(a: float | None, unit: str) -> float | None:
    units = {'': 1.0, 'a': 1.0, 'm': 0.001, 'M': 60.0, 'u': 0.000001}
    return None if a is None or unit not in units else a * units[unit]


def convert_T(t0: float | None, unit: str) -> float | None:
    if t0 is None:
        return None
    if unit in ('y', ''):
        return t0
    if unit == 'c':
        return t0 * 100.0
    if unit == 'd':
        return epoch_from_jd(2400000.0 + t0)
    if unit == 'm':
        return epoch_from_jd(2400000.5 + t0)
    return None


def parse_orbit(line: str) -> dict[str, str] | None:
    row = {k: clean(line[s:e]) for k, (s, e) in ORB6_FIELDS.items()}
    if not WDS_RE.fullmatch(row['wds']):
        return None
    row['P_years'] = fmt(convert_period(num(row['period']), row['period_unit']))
    row['a_arcsec'] = fmt(convert_a(num(row['semimajor']), row['semimajor_unit']))
    row['T_besselian'] = fmt(convert_T(num(row['periastron']), row['periastron_unit']))
    return row


def parse_orb6(path: Path) -> dict[str, list[dict[str, str]]]:
    orbits: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(encoding='latin-1', errors='replace') as f:
        for line in f:
            record = parse_orbit(line.rstrip('\r\n'))
            if record:
                orbits[record['wds']].append(record)
    if not orbits:
        raise ValueError('ORB6 file is empty or does not match the official fixed-width layout.')
    return orbits


def download_orb6(dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    failures = []
    for url in ORB6_URLS:
        tmp: Path | None = None
        try:
            print(f'Downloading {url}', file=sys.stderr)
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=90) as response, \
                 tempfile.NamedTemporaryFile(dir=dest.parent, delete=False) as out:
                tmp = Path(out.name)
                head = response.read(512)
                if b'<html' in head.lower() or b'<!doctype' in head.lower():
                    raise ValueError('Server returned HTML, not the orbital elements text')
                out.write(head)
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            if tmp.stat().st_size < 10000:
                raise ValueError('Catalogue download too small; check source')
            tmp.replace(dest)
            return
        except Exception as exc:
            failures.append(f'{url}: {exc}')
            if tmp:
                tmp.unlink(missing_ok=True)
    raise RuntimeError('Failed to download ORB6; save orb6orbits.txt manually and '
                       'use --orb6 PATH.\n' + '\n'.join(failures))


def read_optional_metadata(path: Path | None) -> dict[str, dict[str, str]]:
    """Optional independently sourced columns: wds,parallax_mas,K1,K2,V0."""
    if not path:
        return {}
    by_id = {}
    with path.open(newline='', encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            wds = clean(row.get('wds', ''))
            if not WDS_RE.fullmatch(wds):
                raise ValueError(f'Invalid WDS identifier in metadata CSV: {wds!r}')
            if wds in by_id:
                raise ValueError(f'More than one metadata row for WDS {wds}; '
                                 'do not attach pair-specific RV parameters by WDS alone')
            by_id[wds] = {k: clean(row.get(k, '')) for k in ('parallax_mas', 'K1', 'K2', 'V0')}
    return by_id


def csv_identifiers(path: Path, rows: list[dict[str, str]]) -> tuple[str, set[str]]:
    # Default filenames created by int4_to_star_csv.py: WDS_DISCARDNAME.csv
    m = re.match(r'^(\d{5}[+-]\d{4})(?:_(.*))?$', path.stem)
    wds = clean(rows[0].get('wds', '')) if rows else ''
    if m and wds and m.group(1) != wds:
        raise ValueError(f'WDS mismatch: filename {path.name}, header {wds}')
    wds = wds or (m.group(1) if m else '')
    if not WDS_RE.fullmatch(wds):
        return '', set()
    names = {designation(clean(rows[0].get(k, ''))) for k in ('name_1', 'name_2')}
    if m and m.group(2):
        names.add(designation(m.group(2)))
    return wds, {s for s in names if s}


def best_matching_orbit(choices: list[dict[str, str]], names: set[str],
                        max_grade: int) -> tuple[dict[str, str] | None, str]:
    eligible = []
    for o in choices:
        g = num(o['grade'])
        if g is None or g > max_grade or g < 1:
            continue
        if not all(num(o[k]) is not None for k in
                   ('P_years', 'T_besselian', 'eccentricity', 'a_arcsec',
                    'node', 'omega', 'inclination')):
            continue
        eligible.append(o)
    if not eligible:
        return None, 'no eligible complete ORB6 orbit'
    # Safe association: match WDS *and* discoverer/component designation.
    # Never silently take an AB orbit for an Aa,Ab measurement series.
    matches = [o for o in eligible if designation(o['discoverer']) in names]
    if not matches:
        return None, 'no matching ORB6 component designation'
    matches.sort(key=lambda o: (num(o['grade']) or 99, o['reference']))
    if len(matches) > 1 and matches[0]['grade'] == matches[1]['grade']:
        return matches[0], 'multiple matching solutions with equal grade: first by reference used'
    return matches[0], ''


def valid_astrometry(row: dict[str, str]) -> bool:
    theta = num(row.get('theta_deg', ''))
    rho = num(row.get('rho_arcsec', ''))
    epoch = num(row.get('epoch_besselian', ''))
    return (epoch is not None and theta is not None and 0 <= theta <= 360
            and rho is not None and rho > 0
            and clean(row.get('rho_flag', '')) not in BAD_RHO_FLAGS
            and clean(row.get('theta_flag', '')) != 'V')


def radec(coord: str) -> tuple[str, str]:
    if not re.fullmatch(r'\d{6}\.\d{2}[+-]\d{6}\.\d', coord):
        return '', ''
    h, mi, sec = int(coord[:2]), int(coord[2:4]), float(coord[4:9])
    sign, deg, amin, asec = coord[9], int(coord[10:12]), int(coord[12:14]), float(coord[14:])
    right = h + mi / 60 + sec / 3600  # decimal HOURS
    decl = (deg + amin / 60 + asec / 3600) * (-1 if sign == '-' else 1)
    return fmt(right, 9), fmt(decl, 9)


def make_inp(o: dict[str, str], measurements: list[dict[str, str]],
             aux: dict[str, str], assumed_sigma: float | None = None
             ) -> tuple[str | None, int, int]:
    """Return input text, # rejected missing sigmas, # adopted assumptions.

    A missing uncertainty must never be silently treated as a measurement.
    The caller must explicitly opt into assumed_sigma, which will be recorded
    in the output file and in the match summary for a sensitivity analysis.
    """
    usable = []
    no_sigma = 0
    adopted = 0
    for r in measurements:
        sigma = num(r.get('rho_error_arcsec', ''))
        sigma_flag = clean(r.get('rho_error_flag', ''))
        published = (sigma is not None and sigma > 0 and
                     sigma_flag not in ('<', '>', 'V', 'P', 'F', 'G'))
        if published:
            usable.append((r, sigma))
        elif assumed_sigma is not None:
            adopted += 1
            usable.append((r, assumed_sigma))
        else:
            no_sigma += 1
    if not usable:
        return None, no_sigma, adopted
    ra, dec = radec(o['coord'])
    lines = [f"Object: {o['discoverer']} = WDS {o['wds']}"]
    if ra:
        lines.append(f'R.A.:   {ra}')
    if dec:
        lines.append(f'Dec:     {dec}')
    if aux.get('parallax_mas'):
        lines.append(f"Parallax: {aux['parallax_mas']}    # mas; supplied metadata")
    for label, field in [('P', 'P_years'), ('T', 'T_besselian'),
                         ('e', 'eccentricity'), ('a', 'a_arcsec'),
                         ('W', 'node'), ('w', 'omega'), ('i', 'inclination')]:
        lines.append(f'{label:<7} {o[field]}')
    for label in ('K1', 'K2', 'V0'):
        if aux.get(label):
            lines.append(f'{label:<7} {aux[label]}')
    lines.extend([
        f'C ORB6 orbital solution: {o["reference"]}; grade {o["grade"]}',
        'C Orbital P and T are in Besselian years; a in arcsec; angles in degrees',
        'C Original ORB6 flags: node=' + (o['node_flag'] or 'none')
        + ', omega=' + (o['omega_flag'] or 'none'),
        'C Radial velocities not supplied by ORB6 or INT4',
    ])
    if adopted:
        lines.append(f'C WARNING: {adopted} rho uncertainty values below are ASSUMED, '
                     f'not published; adopted sigma_rho={fmt(assumed_sigma)} arcsec.')
        lines.append('C Use these adopted weights only after a justified error-model '
                     'and sensitivity analysis.')
    lines.append('C Visual and interferometric data')
    for r, sigma in usable:
        ref = clean(r.get('reference', '')) or 'INT4'
        lines.append(f"{r['epoch_besselian']} {r['theta_deg']} {r['rho_arcsec']} "
                     f'{fmt(sigma)} I1 {ref}')
    return '\n'.join(lines) + '\n', no_sigma, adopted


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--int4-dir', type=Path, required=True,
                    help='Directory OR ZIP of INT4 star CSV files from int4_to_star_csv.py')
    ap.add_argument('--orb6', type=Path, default=Path('orb6orbits.txt'),
                    help='Local ORB6 fixed-width ASCII file (downloaded if missing)')
    ap.add_argument('--output-dir', type=Path, default=Path('merged_int4_orb6'))
    ap.add_argument('--min-points', type=int, default=4,
                    help='Minimum resolved INT4 theta+rho records per star (default 4)')
    ap.add_argument('--max-grade', type=int, default=5,
                    help='Maximum permitted ORB6 grade (1-5, default 5)')
    ap.add_argument('--tie-policy', choices=('reference', 'skip'), default='reference',
                    help='For equal-grade orbits, choose alphabetically first reference '
                         '(original behavior) or skip ambiguous system (safer).')
    ap.add_argument('--assumed-rho-error', type=float, default=None, metavar='ARCSEC',
                    help='Optional, explicit assumed rho sigma (arcsec) when INT4 error '
                         'is absent/unusable; labels assumption in each .inp. '
                         'Use only with scientific justification.')
    ap.add_argument('--metadata', type=Path,
                    help='Optional independently sourced CSV: wds,parallax_mas,K1,K2,V0')
    args = ap.parse_args()
    if not (args.int4_dir.is_dir() or (args.int4_dir.is_file() and zipfile.is_zipfile(args.int4_dir))):
        ap.error(f'INT4 folder or ZIP not found: {args.int4_dir}')
    if not 1 <= args.max_grade <= 5:
        ap.error('--max-grade must be between 1 and 5')
    if args.min_points < 1:
        ap.error('--min-points must be positive')
    if args.assumed_rho_error is not None and not (math.isfinite(args.assumed_rho_error) and args.assumed_rho_error > 0):
        ap.error('--assumed-rho-error must be a finite positive arcsecond value')
    if not args.orb6.exists():
        download_orb6(args.orb6)
    orbits = parse_orb6(args.orb6)
    metadata = read_optional_metadata(args.metadata)
    if args.int4_dir.is_dir():
        paths = sorted(args.int4_dir.glob('*.csv'))
    else:
        with zipfile.ZipFile(args.int4_dir) as archive:
            paths = [Path(name) for name in sorted(archive.namelist())
                     if name.lower().endswith('.csv') and not name.endswith('/')]
    if not paths:
        ap.error(f'No *.csv files in {args.int4_dir}')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    summary = []
    ties = []
    is_zip = args.int4_dir.is_file()
    archive = zipfile.ZipFile(args.int4_dir) if is_zip else None
    try:
        for path in paths:
            counts['input_csvs'] += 1
            if archive:
                with archive.open(path.as_posix()) as f:
                    rows = list(csv.DictReader(io.TextIOWrapper(f, encoding='utf-8-sig', newline='')))
            else:
                with path.open(encoding='utf-8-sig', newline='') as f:
                    rows = list(csv.DictReader(f))
            if not rows:
                counts['empty_input'] += 1
                continue
            wds, names = csv_identifiers(path, rows)
            if not wds:
                counts['bad_wds'] += 1
                continue
            resolved = [r for r in rows if valid_astrometry(r)]
            if len(resolved) < args.min_points:
                counts['too_few_observations'] += 1
                continue
            if wds not in orbits:
                counts['not_in_orb6'] += 1
                continue
            o, warning = best_matching_orbit(orbits[wds], names, args.max_grade)
            if o is None:
                counts['no_matching_orbit'] += 1
                continue
            if warning:
                counts['equal_grade_ambiguous'] += 1
                # Preserve *all* orbit candidates of the best available grade.
                best_grade = num(o['grade'])
                for alt in orbits[wds]:
                    if (num(alt['grade']) == best_grade and
                            designation(alt['discoverer']) in names and
                            all(num(alt[k]) is not None for k in
                                ('P_years', 'T_besselian', 'eccentricity',
                                 'a_arcsec', 'node', 'omega', 'inclination'))):
                        ties.append({'int4_file': path.name, 'wds': wds,
                                     'discoverer': alt['discoverer'],
                                     'grade': alt['grade'],
                                     'orbit_reference': alt['reference'],
                                     'P_years': alt['P_years'],
                                     'T_besselian': alt['T_besselian'],
                                     'e': alt['eccentricity'],
                                     'a_arcsec': alt['a_arcsec'],
                                     'W_deg': alt['node'],
                                     'w_deg': alt['omega'],
                                     'i_deg': alt['inclination']})
                if args.tie_policy == 'skip':
                    counts['ambiguous_orbits_skipped'] += 1
                    print(f'WARNING {path.name}: equal-grade orbit tie; '
                          'skipping (--tie-policy skip)', file=sys.stderr)
                    continue
                print(f'WARNING {path.name}: {warning}; '
                      'see equal_grade_orbit_candidates.csv', file=sys.stderr)
            aux = metadata.get(wds, {})
            combined = []
            for r in resolved:
                out = dict.fromkeys(MERGED_FIELDS, '')
                out.update({
                    'wds': wds, 'discoverer': o['discoverer'], 'hip': o['hip'],
                    'hd': o['hd'], 'orbit_reference': o['reference'],
                    'orbit_grade': o['grade'], 'P_years': o['P_years'],
                    'T_besselian': o['T_besselian'], 'e': o['eccentricity'],
                    'a_arcsec': o['a_arcsec'], 'W_deg': o['node'],
                    'w_deg': o['omega'], 'i_deg': o['inclination'],
                    'node_flag': o['node_flag'], 'omega_flag': o['omega_flag'],
                })
                for k in ('epoch_besselian', 'theta_deg', 'rho_arcsec',
                          'theta_error_deg', 'rho_error_arcsec', 'theta_flag',
                          'rho_flag', 'technique', 'reference'):
                    out[k] = clean(r.get(k, ''))
                out.update(aux)
                combined.append(out)
            dest_csv = args.output_dir / f'{path.stem}_merged.csv'
            with dest_csv.open('w', encoding='utf-8', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=MERGED_FIELDS)
                writer.writeheader()
                writer.writerows(combined)
            inp, skipped, adopted = make_inp(o, resolved, aux, args.assumed_rho_error)
            inp_created = bool(inp and (len(resolved) - skipped) >= args.min_points)
            if inp_created:
                (args.output_dir / f'{path.stem}.inp').write_text(inp, encoding='utf-8')
                counts['inp_written'] += 1
            else:
                counts['insufficient_rho_errors_for_inp'] += 1
            counts['csv_written'] += 1
            counts['observations_exported'] += len(combined)
            counts['inp_observations_omitted_missing_uncertainty'] += skipped
            counts['inp_observations_using_assumed_uncertainty'] += adopted if inp_created else 0
            summary.append({'wds': wds, 'discoverer': o['discoverer'],
                            'int4_file': path.name, 'orbit_reference': o['reference'],
                            'grade': o['grade'], 'astrometric_rows': len(resolved),
                            'inp_rows': len(resolved) - skipped if inp_created else 0,
                            'inp_excluded_no_rho_sigma': skipped,
                            'inp_assumed_rho_sigma_count': adopted if inp_created else 0,
                            'inp_created': 'yes' if inp_created else 'no'})
    finally:
        if archive:
            archive.close()
    if ties:
        with (args.output_dir / 'equal_grade_orbit_candidates.csv').open(
                'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=list(ties[0]))
            writer.writeheader()
            writer.writerows(ties)
    if summary:
        with (args.output_dir / 'match_summary.csv').open('w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=list(summary[0]))
            writer.writeheader()
            writer.writerows(summary)
    print('Summary:', file=sys.stderr)
    for k, v in counts.items():
        print(f'  {k}: {v:,}', file=sys.stderr)
    print(f'  systems_in_orb6: {sum(map(len, orbits.values())):,} orbits', file=sys.stderr)
    print(f'Output: {args.output_dir}', file=sys.stderr)
    if not summary:
        print('No matched systems. Confirm both IDs and component designations;',
              'inspect --min-points and --max-grade.', file=sys.stderr)


if __name__ == '__main__':
    main()
