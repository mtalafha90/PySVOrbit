#!/usr/bin/env python3
"""Add observed SBX radial velocities to INT4+ORB6, with RA/Dec crossmatching.

Read the ``*_merged.csv`` files made by merge_int4_orb6_v2.py, and the
four fixed-width files in CDS B/sbx (alias.dat, systems.dat, orbits.dat,
velocities.dat).  Produce a combined, provenance-labelled CSV per matched
system; optionally also .inp files using published RV/astrometry uncertainties.

Important: Matching a common HIP/HD name does not prove the same binary orbit
in hierarchical multiples.  A consistent ORB6/SBX period and eccentricity
is REQUIRED; ambiguous orbit matches are skipped by default.  Uncertainties
are never invented. Dates that cannot be interpreted safely are kept in CSV
but excluded from .inp. All source RV measurements are retained in CSV.

Python >=3.10; standard library only.
CDS format: https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/B/sbx
"""
from __future__ import annotations

import argparse
import csv
import io
import math
import re
import sys
import shutil
import tempfile
import urllib.request
import urllib.parse
from collections import Counter, defaultdict
from pathlib import Path

CDS_MIRRORS = [
    'https://cdsarc.cds.unistra.fr/ftp/cats/B/sbx',
    'https://cdsarc.cds.unistra.fr/ftp/B/sbx',
    'https://vizier.cds.unistra.fr/ftp/cats/B/sbx',
]
TAP_MIRRORS = [
    'https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync',
    'https://tapvizier1.cds.unistra.fr/TAPVizieR/tap/sync',
]

# VizieR TAP exposes the same official B/sbx tables through SQL-like ADQL.
# Select only columns that the original fixed-width parser actually needs.
# Names are VizieR labels (not underlying SQL/database property names).
TAP_FIELDS = {
    'alias': {'seq':'Seq','catalog':'Catalog','version':'Version','ident':'Ident'},
    'systems': {'seq':'Seq','ra_deg':'RA_ICRS','dec_deg':'DE_ICRS',
                'component':'Comp','parallax':'Plx','parallax_error':'e_Plx',
                'parallax_ref':'r_Plx'},
    'orbits': {'seq':'Seq','orbit':'O','period_d':'Per','T0':'T0',
               'ecc':'ecc','omega':'Periastron','K1':'K1','K1_error':'e_K1',
               'K2':'K2','K2_error':'e_K2','V0':'velocity',
               'V0_error':'e_velocity','n_rv1':'o_rms1','n_rv2':'o_rms2',
               'grade':'Grade','bibcode':'BibCode','time_ref':'epoch_reference'},
    'velocities': {'seq':'Seq','orbit':'O','epoch':'Epoch','rv':'RV',
                   'rv_error':'e_RV','weight':'Weight','component':'Comp',
                   'vel_ref':'VelRef','legacy_offset':'offLegacy',
                   'comment':'comment'},
}
TAP_MIN_ROWS = {'alias': 50000, 'systems': 3000, 'orbits': 4000, 'velocities': 100000}

CAT_FILES = ('alias.dat', 'systems.dat', 'orbits.dat', 'velocities.dat')
WDS_PATTERN = re.compile(r'^\d{5}[+-]\d{4}$')

# Official CDS B/sbx byte offsets (1-based inclusive converted to 0-based slices)
FIELDS = {
    'alias': {'seq': (1, 4), 'catalog': (6, 15), 'version': (17, 19),
              'ident': (21, 49)},
    'systems': {'seq': (1, 4), 'ra_deg': (6, 26), 'dec_deg': (43, 64),
                'component': (109, 115), 'parallax': (259, 273),
                'parallax_error': (275, 287), 'parallax_ref': (289, 307)},
    'orbits': {'seq': (1, 4), 'orbit': (6, 6), 'period_d': (8, 25),
               'T0': (45, 61), 'ecc': (84, 94), 'omega': (105, 114),
               'K1': (127, 137), 'K1_error': (139, 146),
               'K2': (148, 158), 'K2_error': (160, 167),
               'V0': (169, 180), 'V0_error': (182, 189),
               'n_rv1': (212, 218), 'n_rv2': (241, 247),
               'grade': (249, 252), 'bibcode': (259, 309),
               'time_ref': (326, 328)},
    'velocities': {'seq': (1, 4), 'orbit': (6, 6), 'epoch': (8, 25),
                   'rv': (27, 37), 'rv_error': (39, 58),
                   'weight': (60, 75), 'component': (77, 82),
                   'vel_ref': (84, 102), 'legacy_offset': (104, 117),
                   'comment': (119, 140)},
}

META_FIELDS = [
    'wds','discoverer','hip','hd','orbit_reference','orbit_grade',
    'P_years','T_besselian','e','a_arcsec','W_deg','w_deg','i_deg',
    'node_flag','omega_flag','parallax_mas','K1','K2','V0',
    'sbx_system','sbx_orbit','sbx_orbit_reference','match_identifier',
    'period_fractional_difference','eccentricity_difference',
    'sbx_parallax_mas','sbx_parallax_error_mas','sbx_parallax_reference',
    'sbx_ra_deg','sbx_dec_deg',
    'match_method','match_separation_arcsec','match_coordinate_source','match_component_check','match_confidence','match_evidence',
    'sbx_K1_kms','sbx_K2_kms','sbx_V0_kms',
]
ASTRO_FIELDS = ['epoch_besselian','theta_deg','rho_arcsec','theta_error_deg',
                'rho_error_arcsec','theta_flag','rho_flag','technique','reference']
RV_FIELDS = ['rv_epoch_raw','rv_epoch_reference','rv_epoch_jd_minus_2400000',
             'rv_epoch_conversion','rv_kms','rv_error_kms','rv_weight',
             'rv_component','rv_component_raw','rv_source','rv_legacy_offset','rv_comment']
OUT_FIELDS = META_FIELDS + ['measurement_type'] + ASTRO_FIELDS + RV_FIELDS
SUMMARY_FIELDS = ['file','wds','discoverer','match_status','detail',
                  'sbx_system','sbx_orbit','match_identifier','orb6_period_yr',
                  'sbx_period_yr','period_fractional_difference',
                  'eccentricity_difference','astrometry_rows','rv_rows',
                  'rv_a_rows','rv_b_rows','rv_inp_rows','astro_inp_rows',
                  'inp_written','match_method','match_separation_arcsec',
                  'match_coordinate_source','match_component_check','match_confidence','match_evidence']


def clean(v: object) -> str:
    s = str(v or '').strip()
    return '' if s in {'', '.', '--', 'nan', 'None'} else s


def number(s: object) -> float | None:
    try:
        x = float(clean(s))
        return x if math.isfinite(x) else None
    except (ValueError, TypeError):
        return None


def nfmt(v: float | None, digits: int = 9) -> str:
    if v is None:
        return ''
    return f'{v:.{digits}f}'.rstrip('0').rstrip('.')


def field_row(line: str, kind: str) -> dict[str, str]:
    return {k: clean(line[a-1:b]) for k, (a,b) in FIELDS[kind].items()}


def norm_catalog(name: str) -> str:
    return re.sub(r'[^A-Z0-9]+','',name.upper())


def norm_id(name: str) -> str:
    return re.sub(r'[^A-Z0-9+-]+','',name.upper())


def digits(s: str) -> str:
    return str(int(s)) if s.isdigit() else ''


def identifier_key(catalog: str, ident: str) -> tuple[str,str] | None:
    cat, val = norm_catalog(catalog), clean(ident)
    # Most SBX aliases separate catalog name and ID; tolerate repeated prefixes.
    if cat in {'HIP', 'HIPPARCOS', 'HD', 'HENRYDRAPER'}:
        prefix = 'HIP' if cat.startswith('HIP') else 'HD'
        val = re.sub(r'^(HIP|HD)\s*[-:]?\s*', '', val, flags=re.I)
        matched = re.search(r'\d+', val)
        if matched and not re.search(r'[A-Za-z]', val[matched.end():]):
            return prefix, str(int(matched.group()))
    if cat in {'WDS', 'WASHINGTONDOUBLESTAR'}:
        m = re.search(r'\d{5}[+-]\d{4}', val)
        if m:
            return 'WDS', m.group()
    return None


def _is_html(head: bytes) -> bool:
    sample = head.lstrip().lower()
    return (sample.startswith((b'<!doctype', b'<html', b'<?xml', b'<error'))
            or b'<html' in sample[:400])


def download_archive_one(folder: Path, filename: str) -> None:
    errors = []
    for urlbase in CDS_MIRRORS:
        url = f'{urlbase}/{filename}'
        tmp: Path | None = None
        try:
            print(f'Downloading {url}', file=sys.stderr)
            request = urllib.request.Request(url, headers={
                'User-Agent': 'Mozilla/5.0 (compatible; SBX research export; Python)'})
            with urllib.request.urlopen(request, timeout=90) as resp, \
                 tempfile.NamedTemporaryFile(dir=folder, delete=False) as dst:
                tmp = Path(dst.name)
                head = resp.read(512)
                if _is_html(head):
                    raise ValueError('archive server returned HTML rather than SBX data')
                dst.write(head)
                shutil.copyfileobj(resp, dst, length=1024 * 1024)
            if tmp.stat().st_size < 1024:
                raise ValueError('download too small')
            tmp.replace(folder / filename)
            return
        except Exception as exc:
            if tmp:
                tmp.unlink(missing_ok=True)
            errors.append(f'{url}: {exc}')
    raise RuntimeError('\n'.join(errors))


def _tap_select(kind: str, first_seq: int, last_seq: int) -> str:
    columns = ', '.join('"'+col+'"' for col in TAP_FIELDS[kind].values())
    return (f'SELECT {columns} FROM "B/sbx/{kind}" '
            f'WHERE "Seq" BETWEEN {first_seq} AND {last_seq}')


def _fetch_tap_csv(kind: str, low: int, high: int) -> list[dict[str, str]]:
    query = _tap_select(kind, low, high)
    post_data = urllib.parse.urlencode({
        'REQUEST':'doQuery', 'LANG':'ADQL','FORMAT':'csv', 'QUERY':query,
        'MAXREC':'1000000'}).encode('ascii')
    errors=[]
    for endpoint in TAP_MIRRORS:
      for transport in ('GET','POST'):
        try:
            if transport == 'GET':
                url = endpoint + '?' + post_data.decode('ascii')
                req = urllib.request.Request(url, headers={
                    'User-Agent':'Mozilla/5.0 (SBX research catalogue ingestion)',
                    'Accept':'text/csv'})
            else:
                req = urllib.request.Request(endpoint, data=post_data, headers={
                    'User-Agent':'Mozilla/5.0 (SBX research catalogue ingestion)',
                    'Content-Type':'application/x-www-form-urlencoded',
                    'Accept':'text/csv'})
            with urllib.request.urlopen(req, timeout=240) as response:
                blob = response.read()
            if _is_html(blob[:512]):
                raise ValueError('VizieR TAP returned HTML/XML, not CSV: '+repr(blob[:140]))
            decoded = blob.decode('utf-8-sig', errors='replace')
            reader = csv.DictReader(io.StringIO(decoded, newline=''))
            actual = set((reader.fieldnames or []))
            required = set(TAP_FIELDS[kind].values())
            if not required.issubset(actual):
                raise ValueError('VizieR TAP output columns differ from expected '+
                                 f'{kind}: missing {sorted(required-actual)}; received '+
                                 f'{sorted(actual)}; preview {decoded[:200]!r}')
            result = list(reader)
            if any(r.get(None) for r in result):
                raise ValueError('Malformed CSV record in TAP response')
            if result and not all(str(r.get('Seq','')).strip().isdigit() for r in result):
                raise ValueError('Non-numeric Seq in TAP response')
            return result
        except Exception as exc:
            errors.append(f'{endpoint} {transport}: {type(exc).__name__}: {exc}')
    raise RuntimeError('TAP retrieval failed for '+f'{kind} Seq {low}-{high}:\n'+'\n'.join(errors))


def _fixed_field(value: object, width: int, numeric: bool) -> str:
    s = clean(value)
    if len(s) > width and numeric:
        try:
            x = float(s)
            for decimals in range(min(12, width-3), -1, -1):
                candidate = f'{x:.{decimals}f}'
                if len(candidate) <= width:
                    s = candidate
                    break
            if len(s) > width:
                s = f'{x:.{max(1,width-7)}g}'
        except (ValueError, OverflowError):
            pass
    if len(s) > width:
        raise ValueError(f'TAP field does not fit official SBX fixed-width allocation: {s!r} > {width}')
    return s.rjust(width) if numeric else s.ljust(width)


def tap_row_to_fixed(kind: str, record: dict[str,str]) -> str:
    spec = FIELDS[kind]
    output = [' '] * max(end for start,end in spec.values())
    for field, label in TAP_FIELDS[kind].items():
        start, end = spec[field]
        # Character fields in the original CDS ReadMe; others are numeric.
        char = field in {'catalog','version','ident','component','parallax_ref',
                         'bibcode','time_ref','vel_ref','comment'}
        value = _fixed_field(record.get(label,''), end-start+1, not char)
        output[start-1:end] = value
    return ''.join(output) + '\n'


def download_tap_one(folder: Path, kind: str, chunk_size: int = 400) -> None:
    out = folder / f'{kind}.dat'
    tmp: Path | None = None
    count = 0
    try:
        with tempfile.NamedTemporaryFile(dir=folder, mode='w', encoding='latin-1',
                                         newline='\n', delete=False) as dst:
            tmp = Path(dst.name)
            # Seq indexes may include gaps. This upper limit covers published 4080 systems.
            for low in range(1, 4201, chunk_size):
                high = min(low+chunk_size-1, 4200)
                print(f'VizieR TAP {kind}: SBX system numbers {low}-{high}',file=sys.stderr)
                batch = _fetch_tap_csv(kind,low,high)
                for rec in batch:
                    dst.write(tap_row_to_fixed(kind,rec))
                    count += 1
        if count < TAP_MIN_ROWS[kind]:
            raise RuntimeError(f'TAP returned only {count} records for {kind}; '
                               f'expected at least {TAP_MIN_ROWS[kind]}. '
                               'This may be a server-side row limit or incomplete data.')
        tmp.replace(out)
        print(f'Saved {count:,} validated VizieR TAP rows -> {out}',file=sys.stderr)
    except Exception:
        if tmp:
            tmp.unlink(missing_ok=True)
        raise


def download_files(folder: Path, method: str = 'auto') -> None:
    folder.mkdir(parents=True,exist_ok=True)
    for filename in CAT_FILES:
        out = folder / filename
        if out.is_file() and out.stat().st_size > 0:
            with out.open('rb') as f:
                head = f.read(512)
            if not _is_html(head):
                print(f'Using cached {out}',file=sys.stderr)
                continue
            print(f'Ignoring invalid cached HTML file: {out}',file=sys.stderr)
            out.unlink()
        kind = filename.removesuffix('.dat')
        if method in ('archive','auto'):
            try:
                download_archive_one(folder,filename)
                continue
            except RuntimeError as exc:
                if method == 'archive':
                    raise
                print(f'CDS archive inaccessible for {filename}; '
                      'trying VizieR TAP fallback.\n'+str(exc),file=sys.stderr)
        download_tap_one(folder,kind)


def index_sbx(folder: Path):
    alias_index: dict[tuple[str,str],set[str]] = defaultdict(set)
    systems = {}
    orbits: dict[str,list[dict[str,str]]] = defaultdict(list)
    rv_counts = Counter()
    for kind in ('alias','systems','orbits'):
        with (folder/f'{kind}.dat').open(encoding='latin-1',errors='replace') as inp:
            for line in inp:
                r = field_row(line,kind)
                seq = r.get('seq','')
                if not seq.isdigit(): continue
                if kind == 'alias':
                    key = identifier_key(r['catalog'],r['ident'])
                    if key: alias_index[key].add(seq)
                elif kind == 'systems':
                    systems[seq] = r
                else:
                    if number(r['period_d']) is not None:
                        orbits[seq].append(r)
    # Process velocities just once; store only (system,orbit) -> list of small dicts.
    velocities = defaultdict(list)
    with (folder/'velocities.dat').open(encoding='latin-1',errors='replace') as inp:
        for line in inp:
            r = field_row(line,'velocities')
            if (not r['seq'].isdigit() or not r['orbit'].isdigit()
                    or number(r['rv']) is None or number(r['epoch']) is None):
                continue
            velocities[(r['seq'],r['orbit'])].append(r)
            rv_counts['velocity_rows'] += 1
    if not systems or not orbits or not velocities:
        raise ValueError('SBX catalog files appear empty or not in documented fixed-width format')
    return alias_index, systems, orbits, velocities, rv_counts


def target_keys(row: dict[str,str]) -> list[tuple[str,str]]:
    keys = []
    for cat, column in [('HIP','hip'),('HD','hd')]:
        val = re.search(r'\d+',clean(row.get(column,'')))
        if val: keys.append((cat,str(int(val.group()))))
    wds = clean(row.get('wds',''))
    if WDS_PATTERN.fullmatch(wds): keys.append(('WDS',wds))
    return keys


def make_candidates(first: dict[str,str], index, systems, orbits, velocities, args):
    available_keys = target_keys(first)
    # Weak WDS association alone is excluded unless explicitly opted in.
    priorities = [key for key in available_keys if key[0] in ('HIP','HD')]
    if args.allow_wds_only:
        priorities.extend(key for key in available_keys if key[0]=='WDS')
    target_period = number(first.get('P_years'))
    target_e = number(first.get('e'))
    if not target_period or target_period <= 0 or target_e is None:
        return [], 'Missing ORB6 P or e; cannot verify same binary orbit'
    ids_to_keys = defaultdict(list)
    for key in priorities:
        for seq in index.get(key,set()):
            ids_to_keys[seq].append(key)
    if not ids_to_keys:
        return [], 'No shared HIP/HD alias in SBX' + (' (or WDS)' if args.allow_wds_only else '')
    matches = []
    for seq, found_keys in ids_to_keys.items():
        for orbit in orbits.get(seq,[]):
            period_days = number(orbit['period_d'])
            ecc = number(orbit['ecc'])
            if not period_days or period_days <= 0 or ecc is None: continue
            period_yr = period_days / 365.242198781
            delta_p = abs(period_yr-target_period)/target_period
            delta_e = abs(ecc-target_e)
            if delta_p > args.max_period_fraction or delta_e > args.max_ecc_difference:
                continue
            raw_rvs = velocities.get((seq,orbit['orbit']),[])
            if not raw_rvs: continue
            # Prefer corroborating HIP and HD over one alias, and those over WDS only.
            has_hip = any(k[0]=='HIP' for k in found_keys)
            has_hd = any(k[0]=='HD' for k in found_keys)
            if has_hip and has_hd: quality = 3
            elif has_hip: quality = 2
            elif has_hd: quality = 1
            else: quality = 0
            matches.append({'seq':seq,'orbit':orbit,'rvs':raw_rvs,'keys':found_keys,
                            'quality':quality,'delta_p':delta_p,'delta_e':delta_e,
                            'period_yr':period_yr,
                            'component':systems.get(seq,{}).get('component','')})
    if not matches:
        return [], 'Identifier found but no RV orbit passes period/eccentricity consistency checks'
    # Keep all possible identifier orbits; hybrid reconciliation decides using
    # *independent* sky and pair evidence instead of deleting weaker matches.
    matches.sort(key=lambda c:(c['delta_p'],c['delta_e'],-len(c['rvs']),c['seq'],c['orbit']['orbit']))
    return matches,''


# ORB6 has precise J2000 coordinates for each published orbit. The WDS name
# encodes a far coarser position (0.1 min in RA, 1 arcmin in Dec). SBX positions
# are ICRS at the SBX system's tabulated position epoch (usually 1991-2016).
ORB6_URLS = (
    'https://www.astro.gsu.edu/wds/orb6/orb6orbits.txt',
    'https://crf.usno.navy.mil/data_products/WDS/orb6/orb6orbits.txt',
)
COORD_PATTERN = re.compile(r'^(\d{2})(\d{2})(\d{2}\.\d{2})([+-])(\d{2})(\d{2})(\d{2}\.\d)$')


def precise_orb6_position(coord: str) -> tuple[float,float] | None:
    match=COORD_PATTERN.fullmatch(coord.strip())
    if not match: return None
    hh,mm,ss,sign,dd,dm,ds=match.groups()
    h,m,sec=int(hh),int(mm),float(ss)
    d,mi,dsec=int(dd),int(dm),float(ds)
    if h>=24 or m>=60 or sec>=60 or d>90 or mi>=60 or dsec>=60: return None
    ra=15*(h+m/60+sec/3600)
    dec=(d+mi/60+dsec/3600)*(1 if sign=='+' else -1)
    if abs(dec)>90: return None
    return ra,dec


def wds_approx_position(wds: str) -> tuple[float,float] | None:
    match = re.fullmatch(r'(\d{2})(\d{3})([+-])(\d{2})(\d{2})',wds)
    if not match: return None
    hh,tenth_min,sign,dd,dm=match.groups()
    hh,tenth_min,dd,dm=int(hh),int(tenth_min),int(dd),int(dm)
    if hh>23 or tenth_min>=600 or dd>90 or dm>=60:return None
    ra=15*(hh+tenth_min/600)
    dec=(dd+dm/60)*(1 if sign=='+' else -1)
    return ra,dec


def angular_separation_arcsec(ra1:float,dec1:float,ra2:float,dec2:float)->float:
    # Haversine is numerically stable for subarcsecond separations and RA wrap.
    d_ra=math.radians((ra2-ra1+180)%360-180)
    d_dec=math.radians(dec2-dec1)
    a=math.sin(d_dec/2)**2 + math.cos(math.radians(dec1))*math.cos(math.radians(dec2))*math.sin(d_ra/2)**2
    return math.degrees(2*math.asin(math.sqrt(min(1,max(0,a)))))*3600


def _download_orb6_positions(path: Path) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    errors=[]
    for url in ORB6_URLS:
        tmp=None
        try:
            print('Downloading ORB6 coordinates '+url, file=sys.stderr)
            req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'})
            with urllib.request.urlopen(req,timeout=120) as response, tempfile.NamedTemporaryFile(dir=path.parent,delete=False) as dest:
                tmp=Path(dest.name)
                first=response.read(512)
                if _is_html(first): raise ValueError('HTML returned instead of ORB6')
                dest.write(first)
                shutil.copyfileobj(response,dest)
            if tmp.stat().st_size<10000:raise ValueError('ORB6 file too small')
            tmp.replace(path)
            return
        except Exception as exc:
            errors.append(f'{url}: {exc}')
            if tmp:tmp.unlink(missing_ok=True)
    raise RuntimeError('ORB6 file required for precise position matching. Download '+
                       'orb6orbits.txt manually or use --orb6 PATH.\n'+'\n'.join(errors))


def load_orb6_positions(path:Path) -> dict[tuple[str,str],list[dict]]:
    if not path.exists(): _download_orb6_positions(path)
    index=defaultdict(list)
    with path.open(encoding='latin-1',errors='replace') as f:
        for line in f:
            if len(line)<246: continue
            wds=line[19:29].strip()
            if not WDS_PATTERN.fullmatch(wds): continue
            pos=precise_orb6_position(line[0:18])
            if pos is None: continue
            discoverer=line[30:44].strip()
            ref=line[237:245].strip()
            index[(wds,norm_catalog(discoverer))].append({'pos':pos,'reference':ref})
    if not index:raise ValueError('ORB6 file has no valid precise J2000 coordinates: '+str(path))
    return index


def lookup_target_position(first:dict, orb6_index:dict, enable_wds_fallback:bool):
    wds=clean(first.get('wds',''))
    name=norm_catalog(clean(first.get('discoverer','')))
    rows=orb6_index.get((wds,name),[])
    ref=clean(first.get('orbit_reference',''))
    exact=[r for r in rows if r['reference']==ref]
    if exact: rows=exact
    # If more than one ORB6 row remains, require agreement within one arcsecond;
    # otherwise the selected orbit does not have an unambiguous precise position.
    if rows:
        coords=[r['pos'] for r in rows]
        if all(angular_separation_arcsec(*coords[0],*c)<=1 for c in coords[1:]):
            return coords[0], 'ORB6_J2000_precise'
    if enable_wds_fallback:
        pos=wds_approx_position(wds)
        if pos:return pos,'WDS_J2000_approximate'
    return None,'coordinate_unavailable'


def make_sky_index(systems:dict)->dict[int,list[tuple[str,float,float]]]:
    bins=defaultdict(list)
    for seq,s in systems.items():
        ra=number(s.get('ra_deg'))
        dec=number(s.get('dec_deg'))
        if ra is not None and dec is not None and 0<=ra<360 and -90<=dec<=90:
            bins[math.floor(dec)].append((seq,ra,dec))
    return bins


def neighbor_systems(point:tuple[float,float], bins:dict, radius_arcsec:float):
    ra,dec=point
    pad=radius_arcsec/3600
    lower=math.floor(max(-90,dec-pad))
    upper=math.floor(min(90,dec+pad))
    results=[]
    for key in range(lower,upper+1):
        for seq,rs,ds in bins.get(key,[]):
            sep=angular_separation_arcsec(ra,dec,rs,ds)
            if sep<=radius_arcsec: results.append((seq,sep))
    return sorted(results,key=lambda x:(x[1],x[0]))


def orb6_pair_label(discoverer:str)->str:
    # ORB6 format: absence of a component suffix implies the AB pair.
    m=re.match(r'^\s*[A-Za-z]{1,5}\s*\d+\s*(.*)$',discoverer)
    suffix=norm_catalog(m.group(1)) if m else ''
    return suffix or 'AB'


def component_comparison(discoverer: str, sbx_comp: str) -> str:
    """Conservative component comparison: conflicts require two explicit pairs.

    A single SBX component label A/B is not itself evidence that a visual
    AB orbit represents different bodies. Hierarchical labels are not reduced
    to the same pair; uncertain cases are labelled unknown, never agree.
    """
    target = orb6_pair_label(discoverer)
    sbx = norm_catalog(sbx_comp)
    if not sbx:
        return 'unknown_sbx_component'
    if target == sbx:
        return 'agrees'
    known_pairs = {'AB','AC','AD','BC','BD','CD','AE','BE','CE','DE',
                   'AAAB','AABB','BA BB'.replace(' ',''),'BABB', 'A1A2','B1B2'}
    if target in known_pairs and sbx in known_pairs:
        return 'conflict'
    return 'unknown_component_mapping'


def make_position_candidates(first,index,systems,orbits,velocities,sky_index,orb6_index,args):
    pos,source=lookup_target_position(first,orb6_index,not args.no_wds_fallback)
    if pos is None:
        return [],'No usable ORB6 or WDS sky coordinates',[],source
    radius=args.radius_arcsec if source=='ORB6_J2000_precise' else args.wds_radius_arcsec
    neighbors=neighbor_systems(pos,sky_index,radius)
    if not neighbors:
        return [],f'No SBX source within {radius:g} arcsec of {source}',[],source
    period=number(first.get('P_years'))
    ecc=number(first.get('e'))
    if period is None or period<=0 or ecc is None:
        return [],'Missing ORB6 period/eccentricity',[],source
    available_keys=set(target_keys(first))
    matches=[]
    audit=[]
    for seq,sep in neighbors:
        found_keys=[k for k in available_keys if seq in index.get(k,set())]
        comp_check=component_comparison(first.get('discoverer',''),systems[seq].get('component',''))
        for orbit in orbits.get(seq,[]):
            pd=number(orbit['period_d']); e=number(orbit['ecc'])
            if not pd or pd<=0 or e is None:continue
            py=pd/365.242198781
            dp=abs(py-period)/period; de=abs(e-ecc)
            rvs=velocities.get((seq,orbit['orbit']),[])
            passes=(dp<=args.max_period_fraction and de<=args.max_ecc_difference
                    and bool(rvs) and (comp_check!='conflict' or args.allow_component_conflicts))
            audit.append({'file':'','wds':first.get('wds',''),
                'orb6_discoverer':first.get('discoverer',''),
                'sbx_system':seq,'sbx_orbit':orbit['orbit'],
                'separation_arcsec':nfmt(sep,5),'coordinate_source':source,
                'sbx_component':systems[seq].get('component',''),
                'component_check':comp_check,
                'orb6_period_years':nfmt(period,8),'sbx_period_years':nfmt(py,8),
                'fractional_period_difference':nfmt(dp,7),
                'eccentricity_difference':nfmt(de,7),
                'rv_count':len(rvs),'accepted_candidate':'yes' if passes else 'no',
                'sbx_orbit_reference':orbit.get('bibcode',''),
                'shared_identifiers':';'.join(f'{k}:{v}' for k,v in sorted(found_keys))})
            if not passes:continue
            has_hip=any(k[0]=='HIP' for k in found_keys)
            has_hd=any(k[0]=='HD' for k in found_keys)
            quality=3 if (has_hip and has_hd) else 2 if has_hip else 1 if has_hd else 0
            matches.append({'seq':seq,'orbit':orbit,'rvs':rvs,'keys':found_keys,
                'quality':quality,'delta_p':dp,'delta_e':de,'period_yr':py,
                'component':systems[seq].get('component',''),
                'match_method':'RADEC'+('+ALIAS' if found_keys else ''),
                'match_separation_arcsec':sep,'match_coordinate_source':source,
                'match_component_check':comp_check})
    if not matches:
        return [],'Nearby SBX sources found, but orbit/component consistency failed',audit,source
    matches.sort(key=lambda c:(c['match_component_check']!='agrees',
        c['match_separation_arcsec'],c['delta_p'],c['delta_e'],c['seq'],c['orbit']['orbit']))
    return matches,'',audit,source


def match_candidates(first,index,systems,orbits,velocities,sky_index,orb6_index,args):
    """True hybrid union, deduplicated by SBX system+orbit, not fallback.

    Candidate evidence is evaluated on both paths. A unique highest evidence
    tier may be selected, while ties at that tier are always ambiguous.
    Alias-only results inconsistent with precise sky coordinates are excluded;
    their discrepancy is retained in the position audit.
    """
    if args.match_mode == 'coordinates':
        near, why, audit, source = make_position_candidates(
            first,index,systems,orbits,velocities,sky_index,orb6_index,args)
        for c in near:
            c['match_confidence'] = 'medium' if (c['match_component_check']=='agrees' and
                c['match_coordinate_source']=='ORB6_J2000_precise' and
                c['match_separation_arcsec']<=5) else 'review'
            c['match_evidence'] = 'coordinate_candidate'
        return near,why,audit
    aliases, alias_why = make_candidates(first,index,systems,orbits,velocities,args)
    near,sky_why,audit,source = make_position_candidates(
        first,index,systems,orbits,velocities,sky_index,orb6_index,args)
    if args.match_mode == 'aliases':
        near = []
    else:
        for entry in audit:
            if entry.get('accepted_candidate')=='yes':
                entry['accepted_candidate']='passed_position_screen'
    point,source = lookup_target_position(first,orb6_index,not args.no_wds_fallback)
    radius=args.radius_arcsec if source=='ORB6_J2000_precise' else args.wds_radius_arcsec
    by_key={(c['seq'],c['orbit']['orbit']):dict(c) for c in near}
    for c in aliases:
        key=(c['seq'],c['orbit']['orbit'])
        sysrow=systems.get(c['seq'],{})
        comp_check=component_comparison(first.get('discoverer',''),sysrow.get('component',''))
        if comp_check=='conflict' and not args.allow_component_conflicts:
            audit.append({'file':'','wds':first.get('wds',''), 'orb6_discoverer':first.get('discoverer',''),
                'sbx_system':c['seq'],'sbx_orbit':c['orbit']['orbit'],
                'separation_arcsec':'','coordinate_source':source,
                'sbx_component':sysrow.get('component',''), 'component_check':comp_check,
                'orb6_period_years':first.get('P_years',''),
                'sbx_period_years':nfmt(c['period_yr']),
                'fractional_period_difference':nfmt(c['delta_p']),
                'eccentricity_difference':nfmt(c['delta_e']),
                'rv_count':len(c['rvs']),'accepted_candidate':'no',
                'sbx_orbit_reference':c['orbit'].get('bibcode',''),
                'shared_identifiers':';'.join(f'{k}:{v}' for k,v in c['keys']),
                'reject_reason':'alias_component_conflict'})
            continue
        if key in by_key:
            # Merge alias evidence into spatial candidate, including every ID.
            existing=by_key[key]
            existing['keys']=sorted(set(existing['keys'])|set(c['keys']))
            existing['quality']=max(existing.get('quality',0),c.get('quality',0))
            continue
        sep=None
        if point is not None:
            ra=number(sysrow.get('ra_deg'));dec=number(sysrow.get('dec_deg'))
            if ra is not None and dec is not None:
                sep=angular_separation_arcsec(*point,ra,dec)
        if sep is not None and sep>radius:
            audit.append({'file':'','wds':first.get('wds',''), 'orb6_discoverer':first.get('discoverer',''),
                'sbx_system':c['seq'],'sbx_orbit':c['orbit']['orbit'],
                'separation_arcsec':nfmt(sep,4),'coordinate_source':source,
                'sbx_component':sysrow.get('component',''), 'component_check':comp_check,
                'orb6_period_years':first.get('P_years',''),
                'sbx_period_years':nfmt(c['period_yr']),
                'fractional_period_difference':nfmt(c['delta_p']),
                'eccentricity_difference':nfmt(c['delta_e']),
                'rv_count':len(c['rvs']),'accepted_candidate':'no',
                'sbx_orbit_reference':c['orbit'].get('bibcode',''),
                'shared_identifiers':';'.join(f'{k}:{v}' for k,v in c['keys']),
                'reject_reason':'alias_coordinates_disagree'})
            continue
        c=dict(c)
        c.update(match_method='ALIAS_ONLY',match_separation_arcsec=sep,
                 match_coordinate_source=source if sep is not None else 'unavailable',
                 match_component_check=comp_check)
        by_key[key]=c
    matches=list(by_key.values())
    for c in matches:
        keys=c.get('keys',[])
        explicit=any(k[0] in ('HIP','HD') for k in keys)
        sky=c.get('match_method','').startswith('RADEC')
        agree=c.get('match_component_check')=='agrees'
        precise=c.get('match_coordinate_source')=='ORB6_J2000_precise'
        sep=c.get('match_separation_arcsec')
        if sky and explicit and precise and agree:
            tier=4;conf='high';method='RADEC+ALIAS+PAIR'
        elif sky and explicit and precise:
            tier=3;conf='medium';method='RADEC+ALIAS'
        elif sky and precise and agree and sep is not None and sep<=5:
            tier=2;conf='medium';method='RADEC+PAIR'
        elif sky:
            tier=1;conf='review';method='RADEC_ONLY'
        elif explicit:
            tier=1;conf='review';method='ALIAS_ONLY'
        else:
            tier=0;conf='review';method='WDS_ONLY'
        c['match_evidence']=method
        c['match_method']=method
        c['match_confidence']=conf
        c['_tier']=tier
    if not matches:
        return [],f'No accepted candidate: sky={sky_why}; alias={alias_why}',audit
    # Prefer independent evidence, but never break equal-evidence ambiguity
    # using alphabetic references, RV count or numerical proximity alone.
    top=max(c['_tier'] for c in matches)
    winners=[c for c in matches if c['_tier']==top]
    for c in matches:
        audit.append({'file':'','wds':first.get('wds',''), 'orb6_discoverer':first.get('discoverer',''),
            'sbx_system':c['seq'],'sbx_orbit':c['orbit']['orbit'],
            'separation_arcsec':nfmt(c.get('match_separation_arcsec')),
            'coordinate_source':c.get('match_coordinate_source',''),
            'sbx_component':systems.get(c['seq'],{}).get('component',''),
            'component_check':c.get('match_component_check',''),
            'orb6_period_years':first.get('P_years',''),
            'sbx_period_years':nfmt(c['period_yr']),
            'fractional_period_difference':nfmt(c['delta_p']),
            'eccentricity_difference':nfmt(c['delta_e']),
            'rv_count':len(c['rvs']),
            'accepted_candidate':'yes' if c in winners else 'lower_evidence_candidate',
            'sbx_orbit_reference':c['orbit'].get('bibcode',''),
            'shared_identifiers':';'.join(f'{k}:{v}' for k,v in c.get('keys',[])),
            'reject_reason':'' if c in winners else 'lower_evidence'})
    return winners,'',audit


def rv_date(raw: str, time_ref: str, infer_short_jd: bool = True) -> tuple[str,str]:
    """Return sample-file JD-2400000 with explicit conversion provenance.

    SBX labels each orbit JD or MJD. For JD, abbreviated ~5-digit epochs
    follow the SB9/SBX reduced JD convention; the full timescale/heliocentric
    vs barycentric status is NOT established by the catalogue flag alone.
    `infer_short_jd=False` preserves the earlier strict no-short-JD policy.
    """
    epoch=number(raw)
    kind=clean(time_ref).upper()
    if epoch is None or epoch<=0:
        return '', 'invalid_epoch'
    if kind == 'MJD':
        if 10000<=epoch<100000:
            return nfmt(epoch+0.5,8),'MJD_to_JD_minus_2400000_plus_0.5'
        return '', 'MJD_value_outside_supported_range'
    if kind == 'JD':
        if 2_300_000<=epoch<=2_600_000:
            return nfmt(epoch-2_400_000,8),'full_JD_minus_2400000'
        if 10000<=epoch<100000 and infer_short_jd:
            return nfmt(epoch,8),'SBX_abbreviated_JD_assumed_JD_minus_2400000_CHECK_SOURCE'
        return '', 'short_JD_unverified_or_outside_supported_range'
    return '',f'unknown_epoch_reference_{kind or "blank"}'


def map_component(component: str) -> str:
    # SBX documents 'a' and 'b' as RV primary/secondary.  Avoid silently
    # mapping multi-component identifiers (A1, B2, Aa, etc.) to wrong bodies.
    val = clean(component).lower()
    return {'a':'Va','b':'Vb','1':'Va','2':'Vb'}.get(val,'')


def valid_sigma(raw: str) -> float | None:
    val = number(raw)
    return val if val is not None and 0 < val < 1000 else None


def make_inp(meta: dict[str,str], astro: list[dict[str,str]],
             rv: list[dict[str,str]], args) -> tuple[str|None,int,int,str]:
    ast_used = []
    ass_ast = 0
    for rec in astro:
        rho = number(rec.get('rho_arcsec'))
        th = number(rec.get('theta_deg'))
        epoch = number(rec.get('epoch_besselian'))
        if rho is None or rho <= 0 or th is None or epoch is None: continue
        sigma = valid_sigma(rec.get('rho_error_arcsec',''))
        if sigma is None and args.assumed_rho_error is not None:
            sigma = args.assumed_rho_error
            ass_ast += 1
        if sigma is not None:
            ast_used.append((rec,sigma))
    if len(ast_used) < args.min_astrometry_for_inp:
        return None, len(ast_used), 0, 'too few astrometry rows with uncertainties'
    rv_used = []
    assumed_rv = 0
    for rec in rv:
        if rec['rv_component'] not in ('Va','Vb'): continue
        if not rec['rv_epoch_jd_minus_2400000']: continue
        offset = number(rec.get('rv_legacy_offset'))
        # Preserve these in CSV, but cannot fit as if all RV zero points agree.
        if offset is not None and abs(offset)>1e-9 and not args.allow_legacy_offsets:
            continue
        weight = number(rec.get('rv_weight'))
        if weight is not None and weight <= 0: continue
        sigma = valid_sigma(rec.get('rv_error_kms',''))
        if sigma is None and args.assumed_rv_error is not None:
            sigma = args.assumed_rv_error
            assumed_rv += 1
        if sigma is not None: rv_used.append((rec,sigma))
    if not rv_used:
        return None,len(ast_used),0,'no RVs with usable epoch/component/uncertainty'
    # Use the ORB6 geometry, SBX RV parameters. Do not replace ORB6 e/omega/T.
    output = [f"Object: {meta['discoverer']} = WDS {meta['wds']}"]
    ra,dec = number(meta.get('sbx_ra_deg')),number(meta.get('sbx_dec_deg'))
    if ra is not None and 0 <= ra < 360:
        output.append(f'R.A.: {nfmt(ra/15,8)}')
    if dec is not None and -90 <= dec <= 90:
        output.append(f'Dec: {nfmt(dec,8)}')
    px = meta.get('parallax_mas') or (meta.get('sbx_parallax_mas') if (number(meta.get('sbx_parallax_mas')) or 0)>0 else '')
    if px:
        output.append(f'Parallax: {px}    # mas')
    for inpfield, col in [('P','P_years'),('T','T_besselian'),('e','e'),
                          ('a','a_arcsec'),('W','W_deg'),('w','w_deg'),('i','i_deg')]:
        if clean(meta.get(col,'')): output.append(f'{inpfield:<7} {meta[col]}')
    for inpfield, col in [('K1','sbx_K1_kms'),('K2','sbx_K2_kms'),('V0','sbx_V0_kms')]:
        if clean(meta.get(col,'')): output.append(f'{inpfield:<7} {meta[col]}')
    output.append('C Coordinates from SBX (RA decimal hours, Dec decimal degrees)')
    output.append('C Parallax is from '+ ('supplied metadata' if meta.get('parallax_mas') else (meta.get('sbx_parallax_reference') or 'SBX')))
    output.append('C Orbital astrometry from ORB6; RV data/parameters from SBX')
    output.append('C SBX orbit: system '+meta['sbx_system']+', orbit '+meta['sbx_orbit']+'; reference '+meta['sbx_orbit_reference'])
    output.append('C RV epochs represent catalogue JD-2400000, not guaranteed HJD/BJD; verify original references')
    output.append('C RV origin labels and any instrument-dependent zero points are in the merged CSV')
    if ass_ast: output.append(f'C WARNING {ass_ast} astrometric sigma values are explicitly ASSUMED ({args.assumed_rho_error} arcsec)')
    if assumed_rv: output.append(f'C WARNING {assumed_rv} RV sigma values are explicitly ASSUMED ({args.assumed_rv_error} km/s)')
    if args.allow_legacy_offsets:
        output.append('C WARNING: included RVs with non-zero legacy offsets, NOT corrected; check before science use')
    if any('CHECK_SOURCE' in rec.get('rv_epoch_conversion','') for rec,sig in rv_used):
        output.append('C WARNING: abbreviated JD interpreted as JD-2400000; verify original RV publication')
    for comp in ('Va','Vb'):
        comp_rows = sorted(((rec,sigma) for rec,sigma in rv_used if rec['rv_component']==comp),
                           key=lambda pair: float(pair[0]['rv_epoch_jd_minus_2400000']))
        if comp_rows:
            output.append('C Primary RVs' if comp=='Va' else 'C Secondary RVs')
        for rec,sigma in comp_rows:
            ref = re.sub(r'\s+','_',rec.get('rv_source') or 'SBX')[:24]
            output.append(f"{rec['rv_epoch_jd_minus_2400000']} {rec['rv_kms']} {nfmt(sigma,7)} {comp} {ref}")
    output.append('C Visual and interferometric data')
    for rec,sigma in ast_used:
        ref = re.sub(r'\s+','_',rec.get('reference') or 'INT4')[:24]
        output.append(f"{rec['epoch_besselian']} {rec['theta_deg']} {rec['rho_arcsec']} {nfmt(sigma,8)} I1 {ref}")
    return '\n'.join(output)+'\n',len(ast_used),len(rv_used),'ok'


def run(args) -> dict[str,int]:
    if not args.merged_dir.is_dir():
        raise ValueError(f'Missing merged directory: {args.merged_dir}')
    if args.merged_dir.resolve() == args.output_dir.resolve():
        raise ValueError('--output-dir must differ from --merged-dir to protect original measurements')
    if args.min_astrometry_for_inp < 1: raise ValueError('min astro count must be >=1')
    for opt_name in ('max_period_fraction','max_ecc_difference','radius_arcsec','wds_radius_arcsec'):
        val = getattr(args,opt_name)
        if not (math.isfinite(val) and val >= 0): raise ValueError(opt_name+' must be nonnegative finite')
    for opt_name in ('assumed_rho_error','assumed_rv_error'):
        val = getattr(args,opt_name)
        if val is not None and not (math.isfinite(val) and val > 0):
            raise ValueError(opt_name+' must be positive')
    if args.radius_arcsec<=0 or args.wds_radius_arcsec<=0:
        raise ValueError('Sky matching radii must be > 0 arcsec')
    download_files(args.sbx_dir,args.download_method)
    idx, systems, orbits, velocities, stats = index_sbx(args.sbx_dir)
    sky_index=make_sky_index(systems)
    orb6_positions=(load_orb6_positions(args.orb6) if args.match_mode!='aliases' else {})
    print(f'ORB6 precise position groups: {len(orb6_positions):,}',file=sys.stderr)
    print(f'SBX index: {len(systems)} systems; {sum(map(len,orbits.values()))} orbits; '
          f'{stats["velocity_rows"]} velocities', file=sys.stderr)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    paths = sorted(args.merged_dir.glob('*_merged.csv'))
    if not paths:
        raise ValueError('No *_merged.csv in '+str(args.merged_dir)+'. Use the output of merge_int4_orb6_v2.py')
    summary = []
    candidate_audit = []
    position_audit = []
    tally = Counter()
    for path in paths:
        tally['input_files'] += 1
        with path.open(newline='',encoding='utf-8-sig') as f:
            records = list(csv.DictReader(f))
        if not records: continue
        first=records[0]
        # Summaries include all source measurements whether matched or not.
        row = dict.fromkeys(SUMMARY_FIELDS,'')
        row.update({'file':path.name,'wds':first.get('wds',''),
                    'discoverer':first.get('discoverer',''),
                    'astrometry_rows':len(records),
                    'orb6_period_yr':first.get('P_years','')})
        matches,reason,sky_audit = match_candidates(first,idx,systems,orbits,velocities,
            sky_index,orb6_positions,args)
        for entry in sky_audit: entry['file']=path.name
        position_audit.extend(sky_audit)
        if not matches:
            row.update(match_status='unmatched',detail=reason)
            tally['unmatched'] += 1
            summary.append(row)
            if not args.matched_only:
                shutil.copy2(path, args.output_dir/path.name)
                old_inp=args.merged_dir/(path.stem.removesuffix('_merged')+'.inp')
                if old_inp.is_file(): shutil.copy2(old_inp,args.output_dir/old_inp.name)
                tally['preserved_without_sbx'] += 1
            continue
        if len(matches)>1:
            for c in matches:
                candidate_audit.append({'file':path.name,'wds':first.get('wds',''),
                    'seq':c['seq'],'orbit':c['orbit']['orbit'],
                    'period_yr':nfmt(c['period_yr']),
                    'delta_p':nfmt(c['delta_p']), 'delta_e':nfmt(c['delta_e']),
                    'RV_count':len(c['rvs']),
                    'reference':c['orbit']['bibcode'],
                    'component':c['component'],
                    'aliases':';'.join(f'{k}:{v}' for k,v in c['keys']),
                    'match_method':c.get('match_method',''),
                    'separation_arcsec':nfmt(c.get('match_separation_arcsec')),
                    'component_check':c.get('match_component_check',''),
                    'match_confidence':c.get('match_confidence','review')})
            if args.orbit_policy=='skip':
                row.update(match_status='ambiguous',detail=f'{len(matches)} spectroscopic orbits pass thresholds')
                tally['ambiguous_skipped'] += 1
                summary.append(row)
                if not args.matched_only:
                    shutil.copy2(path,args.output_dir/path.name)
                    old_inp=args.merged_dir/(path.stem.removesuffix('_merged')+'.inp')
                    if old_inp.is_file(): shutil.copy2(old_inp,args.output_dir/old_inp.name)
                    tally['preserved_without_sbx'] += 1
                continue
            # Explicit option: prefer most RV points among near-period matches.
            matches.sort(key=lambda c:(-len(c['rvs']),c['delta_p'],c['delta_e']))
        chosen=matches[0]
        seq,orb=chosen['seq'],chosen['orbit']
        sysrow=systems.get(seq,{})
        match_key=';'.join(f'{k}:{v}' for k,v in chosen['keys'])
        if not match_key:match_key='sky coordinates (no shared HIP/HD alias)'
        meta={f:clean(first.get(f,'')) for f in META_FIELDS}
        meta.update({'sbx_system':seq,'sbx_orbit':orb['orbit'],
            'sbx_orbit_reference':orb['bibcode'],'match_identifier':match_key,
            'period_fractional_difference':nfmt(chosen['delta_p']),
            'eccentricity_difference':nfmt(chosen['delta_e']),
            'sbx_parallax_mas':sysrow.get('parallax',''),
            'sbx_ra_deg':sysrow.get('ra_deg',''),
            'sbx_dec_deg':sysrow.get('dec_deg',''),
            'sbx_parallax_error_mas':sysrow.get('parallax_error',''),
            'sbx_parallax_reference':sysrow.get('parallax_ref',''),
            'sbx_K1_kms':orb['K1'],'sbx_K2_kms':orb['K2'],
            'sbx_V0_kms':orb['V0'],
            'match_method':chosen.get('match_method',''),
            'match_separation_arcsec':nfmt(chosen.get('match_separation_arcsec'),5),
            'match_coordinate_source':chosen.get('match_coordinate_source',''),
            'match_component_check':chosen.get('match_component_check',''),
            'match_confidence':chosen.get('match_confidence','review'),
            'match_evidence':chosen.get('match_evidence','')})
        combined=[]
        astro_rows=[]
        for a in records:
            output=dict.fromkeys(OUT_FIELDS,'')
            output.update(meta)
            output['measurement_type']='astrometry'
            for f in ASTRO_FIELDS: output[f]=clean(a.get(f,''))
            combined.append(output)
            astro_rows.append(output)
        rv_rows=[]
        for r in chosen['rvs']:
            output=dict.fromkeys(OUT_FIELDS,'')
            output.update(meta)
            output['measurement_type']='radial_velocity'
            output.update({'rv_epoch_raw':r['epoch'],'rv_epoch_reference':orb['time_ref'],
                           'rv_kms':r['rv'],'rv_error_kms':r['rv_error'],
                           'rv_weight':r['weight'],'rv_component':map_component(r['component']),
                           'rv_component_raw':r['component'],
                           'rv_source':r['vel_ref'],'rv_legacy_offset':r['legacy_offset'],
                           'rv_comment':r['comment']})
            t,tag=rv_date(r['epoch'],orb['time_ref'],not args.strict_epochs)
            output['rv_epoch_jd_minus_2400000']=t
            output['rv_epoch_conversion']=tag
            rv_rows.append(output)
            combined.append(output)
        output_csv=args.output_dir/path.name
        with output_csv.open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=OUT_FIELDS)
            writer.writeheader();writer.writerows(combined)
        inp,nastro,nrv,note=make_inp(meta,astro_rows,rv_rows,args)
        if inp:
            output_inp=args.output_dir/(path.stem.removesuffix('_merged')+'.inp')
            output_inp.write_text(inp,encoding='utf-8')
            tally['inp_written']+=1
        else:
            old_inp = args.merged_dir/(path.stem.removesuffix('_merged')+'.inp')
            if old_inp.is_file():
                shutil.copy2(old_inp,args.output_dir/old_inp.name)
                tally['inp_preserved_without_rvs'] += 1
        tally['matched'] += 1
        tally['rv_measurements'] += len(rv_rows)
        tally['astrometry_measurements'] += len(astro_rows)
        row.update(match_status='matched',detail=note,sbx_system=seq,sbx_orbit=orb['orbit'],
                   match_identifier=match_key,sbx_period_yr=nfmt(chosen['period_yr']),
                   period_fractional_difference=nfmt(chosen['delta_p']),
                   eccentricity_difference=nfmt(chosen['delta_e']),
                   rv_rows=len(rv_rows),
                   rv_a_rows=sum(r['rv_component']=='Va' for r in rv_rows),
                   rv_b_rows=sum(r['rv_component']=='Vb' for r in rv_rows),
                   rv_inp_rows=nrv,astro_inp_rows=nastro,
                   inp_written='yes' if inp else 'no',
                   match_method=chosen.get('match_method',''),
                   match_separation_arcsec=nfmt(chosen.get('match_separation_arcsec'),5),
                   match_coordinate_source=chosen.get('match_coordinate_source',''),
                   match_component_check=chosen.get('match_component_check',''),
                   match_confidence=chosen.get('match_confidence','review'),
                   match_evidence=chosen.get('match_evidence',''))
        summary.append(row)
    with (args.output_dir/'sbx_match_summary.csv').open('w',newline='',encoding='utf-8') as out:
        writer=csv.DictWriter(out,fieldnames=SUMMARY_FIELDS)
        writer.writeheader();writer.writerows(summary)
    if candidate_audit:
        with (args.output_dir/'sbx_candidate_orbits.csv').open('w',newline='',encoding='utf-8') as out:
            writer=csv.DictWriter(out,fieldnames=list(dict.fromkeys(k for rec in candidate_audit for k in rec)))
            writer.writeheader();writer.writerows(candidate_audit)
    if position_audit:
        with (args.output_dir/'sbx_position_candidates.csv').open('w',newline='',encoding='utf-8') as out:
            writer=csv.DictWriter(out,fieldnames=list(dict.fromkeys(k for rec in position_audit for k in rec)))
            writer.writeheader();writer.writerows(position_audit)
    print('Summary:',file=sys.stderr)
    for k,v in tally.items(): print(f'  {k}: {v:,}',file=sys.stderr)
    print('  Input stars unmatched/ambiguous were NOT guessed or merged',file=sys.stderr)
    print('Output:',args.output_dir,file=sys.stderr)
    return dict(tally)


def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--merged-dir',type=Path,required=True,help='Folder made by merge_int4_orb6_v2.py')
    p.add_argument('--sbx-dir',type=Path,default=Path('sbx_catalog'),help='Local SBX .dat directory; missing files auto-downloaded')
    p.add_argument('--download-method', choices=['auto','archive','tap'], default='auto',
                   help='auto: archive then VizieR TAP; tap: query VizieR directly, avoiding blocked FTP archive')
    p.add_argument('--output-dir',type=Path,default=Path('merged_int4_orb6_sbx'))
    p.add_argument('--orb6',type=Path,default=Path('orb6orbits.txt'),
                   help='ORB6 original text with precise J2000 coordinates; downloads if missing')
    p.add_argument('--match-mode',choices=['coordinates','hybrid','aliases'],default='coordinates',
                   help='coordinates (default): RA/Dec spatial candidates + orbit tests; hybrid: alias fallback; aliases: legacy matching')
    p.add_argument('--radius-arcsec',type=float,default=90.0,
                   help='Angular-search radius for precise ORB6 coordinates (default 90 arcsec)')
    p.add_argument('--wds-radius-arcsec',type=float,default=120.0,
                   help='Angular-search radius for coarse WDS designation fallback (default 120 arcsec)')
    p.add_argument('--no-wds-fallback',action='store_true',
                   help='Require precise ORB6 RA/Dec; skip coarse WDS positional fallback')
    p.add_argument('--allow-component-conflicts',action='store_true',
                   help='Unsafe exploratory option: accept contrary binary-component labels')
    p.add_argument('--max-period-fraction',type=float,default=0.15,help='Maximum fractional ORB6-vs-SBX period mismatch (default 0.15)')
    p.add_argument('--max-ecc-difference',type=float,default=0.25,help='Maximum absolute eccentricity difference (default 0.25)')
    p.add_argument('--allow-wds-only',action='store_true',help='Allow weaker WDS-only identifier matches (still require orbit consistency)')
    p.add_argument('--orbit-policy',choices=['skip','most-rvs'],default='skip',help='Skip ambiguous matching SBX orbits (default), or select one with most RV observations')
    p.add_argument('--matched-only',action='store_true',help='Only write stars with secure SBX matches; default preserves unmatched original CSV/INP files')
    p.add_argument('--min-astrometry-for-inp',type=int,default=4,help='Minimum usable astrometric rows for .inp output')
    p.add_argument('--assumed-rho-error',type=float,help='Optional explicitly assumed rho error (arcsec) for .inp; never silently applied')
    p.add_argument('--assumed-rv-error',type=float,help='Optional explicitly assumed RV error (km/s) for .inp; never silently applied')
    p.add_argument('--allow-inferred-short-jd',action='store_true',help='Backward-compatible alias for default SBX abbreviated JD conversion')
    p.add_argument('--strict-epochs',action='store_true', help='Exclude abbreviated JD epochs from .inp; keep raw epochs in CSV')
    p.add_argument('--allow-legacy-offsets',action='store_true',help='Permit .inp RVs with nonzero legacy zero-point offsets without applying correction; NOT recommended')
    args=p.parse_args()
    try: run(args)
    except (ValueError,RuntimeError,OSError) as e: p.exit(2,'ERROR: '+str(e)+'\n')


if __name__=='__main__':main()
