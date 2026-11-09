#!/usr/bin/env python3
"""Audit/integrate SOAR speckle astrometry from VizieR TSV + AAS MRT files.

No measurements are silently added to orbit fitting files. Only confirmed
resolved positions are candidates for future joint orbit fitting. Unresolved
observations are stored as nondetection constraints and NEVER as positions.

Cross-release duplication and INT4 overlap are *screening* decisions, not
claims of independent provenance. The original source values are preserved.
Python >=3.10, standard library only.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

RELEASES = {
    '2015': (1303, '2024', '2016AJ....151..153T'),
    '2016_2017': (2483, '2018', '2018AJ....155..235T'),
    '2018': (3097, '2019', '2019AJ....158...48T'),
    '2019': (2555, '2020', '2020AJ....160....7T'),
    '2020': (1735, '2021', '2021AJ....162...41T'),
    '2021': (2623, '2022', '2022AJ....164...58T'),
    '2022': (2508, '2023', '2023AJ....166..139T'),
    '2023': (1913, '2024', '2024AJ....168...28T'),
    '2024_2025': (5316, '2026', '2026AJ....172....8T'),
}
WDS = re.compile(r'^\d{5}[+-]\d{4}$')
MRT_COL = re.compile(r'^\s*(\d+)\s*-\s*(\d+)\s+\S+\s+\S+\s+(\S+)\s+')


def normalized(s):
    return re.sub(r'[^A-Z0-9]', '', (s or '').upper())


def num(s):
    try:
        x = float(str(s).strip())
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def alias(row, *names):
    # Case/space/underscore normalization; no positional guessing.
    keys = {normalized(k): v for k, v in row.items() if k is not None}
    for n in names:
        val = keys.get(normalized(n))
        if val is not None:
            return str(val).strip()
    return ''


def jyear_to_jd(y):
    return 2451545.0 + (y - 2000.0) * 365.25


def byear_to_jd(y):
    return 2415020.31352 + (y - 1900.0) * 365.242198781


def byear_of_jyear(y):
    return 1900 + (jyear_to_jd(y) - 2415020.31352) / 365.242198781


def angle_diff(a, b):
    return abs((a - b + 180) % 360 - 180)


def tsv(path):
    lines = path.read_text(encoding='utf-8-sig').splitlines()
    if any(line.startswith('#INFO') and 'Error=' in line for line in lines):
        raise ValueError(f'VizieR server returned an error in {path}')
    header = next((i for i, line in enumerate(lines)
                   if 'WDS' in line.split('\t')
                   and any(name in line.split('\t') for name in ('Disc', 'Design', 'Dis', 'DesLet'))
                   and any(name in line.split('\t') for name in ('Epoch', 'Epoch.Y', 'Tobs'))), None)
    if header is None:
        raise ValueError(f'SOAR VizieR table heading not recognized in {path}')
    reader = csv.DictReader(lines[header:], delimiter='\t')
    rows = []
    for row in reader:
        if None in row:
            raise ValueError(f'Inconsistent TSV column count in {path}')
        values = {k: (v or '').strip() for k, v in row.items()}
        # VizieR has a units line and a dashed-separator line after the header.
        if not WDS.fullmatch(alias(values, 'WDS')):
            if alias(values, 'recno').isdigit():
                raise ValueError(f'Record has bad WDS {alias(values, "WDS")!r} in {path}')
            continue
        rows.append(values)
    if not rows:
        raise ValueError(f'No SOAR positions in {path}')
    return rows


def mrt(path):
    lines = path.read_text(encoding='utf-8-sig').splitlines()
    cols = []
    in_desc = False
    for line in lines:
        if 'Byte-by-byte Description of file' in line:
            in_desc = True
            continue
        if not in_desc:
            continue
        match = MRT_COL.match(line)
        if match:
            lo, hi, name = match.groups()
            cols.append((int(lo)-1, int(hi), name))
        elif cols and line.strip().startswith('----'):
            break
    if not cols or not {'WDS', 'Disc'}.issubset(c[2] for c in cols):
        raise ValueError(f'Cannot parse MRT byte positions in {path}')
    width = max(hi for _, hi, _ in cols)
    rows = []
    for line in lines:
        if len(line) < 75 or not WDS.fullmatch(line[:10].strip()):
            continue
        row = {name: line[lo:hi].strip() for lo,hi,name in cols}
        if not alias(row, 'Epoch.Y', 'Epoch'):
            continue
        rows.append(row)
    if not rows:
        raise ValueError(f'No data rows in AAS MRT file {path}')
    return rows


def load_positions(row, release, source_file, ordinal):
    wds = alias(row, 'WDS')
    disc = alias(row, 'Disc', 'Discoverer', 'Design', 'Dis', 'DesLet')
    component = alias(row, 'Comp')
    raw_epoch = alias(row, 'Epoch', 'Epoch.Y', 'Tobs')
    epoch = num(raw_epoch)
    pa = num(alias(row, 'PA', 'theta', 'theta_deg'))
    sep = num(alias(row, 'Sep', 'rho', 'rho_arcsec'))
    if not WDS.fullmatch(wds) or not disc or epoch is None or pa is None or sep is None:
        raise ValueError(f'{source_file} record {ordinal}: required WDS/discoverer/epoch/PA/separation missing; source columns={list(row)}')
    # The 2015 VizieR table publishes epochs as 15.xxxx, i.e. 2015.xxxx.
    # This is a documented representation of the reported observing epoch,
    # not an imputed measurement. Preserve the original token below.
    short_epoch = release == '2015' and 0 <= epoch < 100
    if short_epoch:
        epoch += 2000
    if not 1800 < epoch < 2100 or not 0 <= pa <= 360 or sep <= 0:
        raise ValueError(f'{source_file} record {ordinal}: invalid epoch/PA/separation: '
                         f'epoch={raw_epoch!r} PA={alias(row, "PA")!r} Sep={alias(row, "Sep")!r}')
    # Releases through 2018 store the component separately. Do not merge
    # e.g. CHR 84 Aa,Ab with an unrelated CHR 84 pair.
    if component and not normalized(disc).endswith(normalized(component)):
        pair_designation = disc + ' ' + component
    else:
        pair_designation = disc
    pa_mas = num(alias(row, 'e_PA', 'err'))
    sep_mas = num(alias(row, 'e_Sep'))
    if pa_mas is not None and pa_mas <= 0: pa_mas = None
    if sep_mas is not None and sep_mas <= 0: sep_mas = None
    return dict(
        source='SOAR', release=release, source_file=source_file.name,
        source_epoch_raw=raw_epoch, source_epoch_expanded=('2000+short_epoch' if short_epoch else ''),
        source_discoverer=disc, source_component=component,
        citation=RELEASES[release][2], source_row_number=ordinal,
        source_recno=alias(row, 'recno'), wds=wds, discoverer=pair_designation,
        pair_key=normalized(pair_designation), other_name=alias(row, 'OName'),
        epoch_julian_year=f'{epoch:.7f}', epoch_besselian_year=f'{byear_of_jyear(epoch):.9f}',
        position_angle_deg=f'{pa:.6f}', separation_arcsec=f'{sep:.7f}',
        tangent_sigma_internal_mas='' if pa_mas is None else str(pa_mas),
        separation_sigma_internal_mas='' if sep_mas is None else str(sep_mas),
        separation_sigma_internal_arcsec='' if sep_mas is None else f'{sep_mas/1000:.7f}',
        position_angle_sigma_internal_deg='' if pa_mas is None else f'{math.degrees((pa_mas/1000)/sep):.8f}',
        filter=alias(row, 'Filt', 'Flt'), n_observations=alias(row, 'N', 'Nobs', 'No', 'Num'),
        magnitude_difference=alias(row, 'Delm', 'Dmag', 'Dm'), quadrant_flag=alias(row, 'f_Delm', 'f_Dmag', 'f_Dm'),
        target_tag=alias(row, 'Tag'), orbit_reference=alias(row, 'Ref'),
        source_pa_residual=alias(row, '(O-C)PA', '(O-C)RA', '(O-C)t', 'O-CPA'),
        source_sep_residual=alias(row, '(O-C)Sep', '(O-C)r', 'O-CSep'),
        original_fields_json=json.dumps(row, sort_keys=True, ensure_ascii=False),
        soar_duplicate_status='', soar_duplicate_of='', int4_status='', int4_reference='',
        catalogue_tier='', orbit_pair_found='', fitting_policy='audit_only_not_exported',
    )


def load_unresolved(row, source_file, ordinal):
    wds, disc = alias(row, 'WDS'), alias(row, 'Disc')
    epoch = num(alias(row, 'Epoch.Y', 'Epoch'))
    if not WDS.fullmatch(wds) or epoch is None:
        raise ValueError(f'{source_file} unresolved record {ordinal}: invalid WDS/epoch')
    return dict(source='SOAR', release='2024_2025', source_file=source_file.name,
                citation=RELEASES['2024_2025'][2], source_row_number=ordinal,
                wds=wds, discoverer=disc, other_name=alias(row, 'OName'),
                epoch_julian_year=f'{epoch:.7f}', filter=alias(row, 'Filt'),
                n_observations=alias(row, 'N'),
                min_resolved_separation_arcsec=alias(row, 'MinSep'),
                delta_mag_limit_015_arcsec=alias(row, 'Delm0.15'),
                delta_mag_limit_1_arcsec=alias(row, 'Delm1'),
                original_fields_json=json.dumps(row, sort_keys=True, ensure_ascii=False),
                constraint_type='non_detection_not_position')


def same_position(a,b,days,angle,sep):
    dtime=abs(jyear_to_jd(float(a['epoch_julian_year'])) - jyear_to_jd(float(b['epoch_julian_year'])))
    if dtime > days:
        return ''
    dang=angle_diff(float(a['position_angle_deg']),float(b['position_angle_deg']))
    dsep=abs(float(a['separation_arcsec'])-float(b['separation_arcsec']))
    if dsep <= sep and dang <= angle:
        return 'probable_duplicate'
    if dsep <= sep and abs(dang-180) <= angle:
        return 'quadrant_review'
    return 'same_epoch_discrepancy'


def duplicate_across_soar(records,days,angle,sep):
    grouped=defaultdict(list)
    duplicates=[]
    for r in records:
        key=(r['wds'],r['pair_key'])
        near=[]
        for prior in grouped[key]:
            match=same_position(r,prior,days,angle,sep)
            if match in ('probable_duplicate','quadrant_review'):
                near.append((match,prior))
        if near:
            status, other=near[0]
            r['soar_duplicate_status']=status
            r['soar_duplicate_of']=f"{other['release']}:{other['source_row_number']}"
            duplicates.append(r)
        else:
            r['soar_duplicate_status']='not_identified'
        grouped[key].append(r)
    return duplicates


def int4_index(raw_catalog, relevant_wds):
    if not raw_catalog:
        return {}
    from int4_to_star_csv import read_catalog,make_row
    idx=defaultdict(list)
    for ids, obs in read_catalog(raw_catalog):
        wds=ids.get('wds','')
        if wds not in relevant_wds:
            continue
        pairs={normalized(ids.get('name_1','')),normalized(ids.get('name_2',''))}
        pairs.discard('')
        for measurement in obs:
            row=make_row(ids,measurement,True)
            epoch=num(row.get('epoch_besselian'))
            pa=num(row.get('theta_deg'))
            rho=num(row.get('rho_arcsec'))
            if epoch is None or pa is None or rho is None or rho <= 0:
                continue
            idx[wds].append(dict(pairs=pairs, jd=byear_to_jd(epoch), pa=pa, sep=rho,
                                 reference=row.get('reference','')))
    return idx


def audit_int4(records, index, days, angle, sep):
    probable=[]; review=[]
    for r in records:
        entries=index.get(r['wds'], [])
        pairs=[x for x in entries if r['pair_key'] in x['pairs']]
        if not entries:
            status='wds_not_in_int4'
        elif not pairs:
            status='component_not_confirmed'
        else:
            epoch=jyear_to_jd(float(r['epoch_julian_year']))
            near=[]
            for x in pairs:
                if abs(epoch-x['jd']) > days:continue
                dpa=angle_diff(float(r['position_angle_deg']),x['pa'])
                dsep=abs(float(r['separation_arcsec'])-x['sep'])
                if dpa<=angle and dsep<=sep:cls='probable_duplicate'
                elif abs(dpa-180)<=angle and dsep<=sep:cls='quadrant_review'
                else:cls='same_epoch_discrepancy'
                near.append((cls,abs(epoch-x['jd']),x))
            if near:
                near.sort(key=lambda x:({'probable_duplicate':0,'quadrant_review':1,'same_epoch_discrepancy':2}[x[0]],x[1]))
                status,_,match=near[0]
                r['int4_reference']=match['reference']
            else: status='same_pair_no_epoch_match'
        r['int4_status']=status
        if status == 'probable_duplicate':probable.append(r)
        if status in ('component_not_confirmed','quadrant_review','same_epoch_discrepancy'):review.append(r)
    return probable,review


def build_tier_map(path):
    if path is None:
        return {}
    targets={}
    with path.open(newline='',encoding='utf-8-sig') as stream:
        for r in csv.DictReader(stream):
            key=(r['wds'],normalized(r['discoverer']))
            if key in targets and r['catalogue_tier']!=targets[key]:
                raise ValueError(f'Ambiguous research tier for {key}')
            targets[key]=r['catalogue_tier']
    return targets


def write_csv(path,rows,fields):
    with path.open('w',newline='',encoding='utf-8') as out:
        w=csv.DictWriter(out,fieldnames=fields,extrasaction='ignore')
        w.writeheader();w.writerows(rows)


def execute(source_dir, output_dir, int4_catalog=None, catalogue_index=None, days=2.0, angle=0.5, sep=0.003, require_counts=True):
    output_dir.mkdir(parents=True,exist_ok=True)
    input_status=[]; records=[]
    for release,(expected,_,_) in RELEASES.items():
        filename=(f'soar_{release}.tsv' if release!='2024_2025' else 'soar_2024_2025_positions.txt')
        path=source_dir/filename
        if not path.is_file():
            raise FileNotFoundError(f'Missing SOAR release {release}: {path}')
        original=mrt(path) if release=='2024_2025' else tsv(path)
        if require_counts and len(original)!=expected:
            raise ValueError(f'{release}: parsed {len(original)} measurements; expected {expected} (file {path}). STOP: format/count mismatch')
        for i,row in enumerate(original,1):
            records.append(load_positions(row,release,path,i))
        input_status.append(dict(release=release,source_file=str(path),records=len(original),expected=expected))
    table2=source_dir/'soar_2024_2025_unresolved.txt'
    if not table2.is_file():
        raise FileNotFoundError(f'Missing Table 2 nondetection MRT: {table2}')
    raw_unresolved=mrt(table2)
    if require_counts and len(raw_unresolved)!=1723:
        raise ValueError(f'Table 2: parsed {len(raw_unresolved)} nondetections; expected 1723; STOP')
    unresolved=[load_unresolved(r,table2,i) for i,r in enumerate(raw_unresolved,1)]

    soaring=duplicate_across_soar(records,days,angle,sep)
    index=int4_index(int4_catalog,{r['wds'] for r in records})
    probable_int4, review_int4=audit_int4(records,index,days,angle,sep) if int4_catalog else ([],[])
    tiers=build_tier_map(catalogue_index)
    for r in records:
        tier=tiers.get((r['wds'],r['pair_key']))
        r['catalogue_tier']=tier or 'no_exact_pair_in_selected_index'
        r['orbit_pair_found']='yes' if tier else 'no'
    pairs=defaultdict(list)
    for r in records:
        pairs[(r['wds'],r['pair_key'])].append(r)
    summaries=[]
    for (wds,pair),rows in sorted(pairs.items()):
        summaries.append(dict(wds=wds,pair_key=pair,discoverer=rows[0]['discoverer'],
                              catalogue_tier=rows[0]['catalogue_tier'],soar_measurements=len(rows),
                              earliest_epoch=min(float(r['epoch_julian_year']) for r in rows),
                              latest_epoch=max(float(r['epoch_julian_year']) for r in rows),
                              probable_soar_duplicates=sum(r['soar_duplicate_status']=='probable_duplicate' for r in rows),
                              probable_int4_duplicates=sum(r['int4_status']=='probable_duplicate' for r in rows),
                              unresolved_component_count=sum(r['int4_status']=='component_not_confirmed' for r in rows),
                              eligible_for_fitting='NO - provenance and uncertainty model review required'))
    # A candidate supplemental table only; DO NOT silently merge into .inp or original INT4.
    candidate=[r for r in records if r['int4_status']!='probable_duplicate'
               and r['soar_duplicate_status']=='not_identified'
               and r['orbit_pair_found']=='yes'
               and r['int4_status']!='component_not_confirmed']
    union_review=[r for r in records if r['soar_duplicate_status']!='not_identified'
                  or r['int4_status'] in ('probable_duplicate','component_not_confirmed','quadrant_review','same_epoch_discrepancy')]
    fields=list(records[0]); other_fields=list(unresolved[0]); summary_fields=list(summaries[0])
    for filename,rs,cols in [
        ('all_soar_positions.csv',records,fields),
        ('soar_crossrelease_duplicate_candidates.csv',soaring,fields),
        ('soar_int4_duplicate_candidates.csv',probable_int4,fields),
        ('soar_component_and_position_review.csv',review_int4,fields),
        ('soar_union_review.csv',union_review,fields),
        ('soar_orb6_sbx_supplement_candidates.csv',candidate,fields),
        ('soar_2024_2025_nondetection_constraints.csv',unresolved,other_fields),
        ('soar_pair_summary.csv',summaries,summary_fields),
    ]:write_csv(output_dir/filename,rs,cols)
    stats=dict(
        total_published_resolved_rows=len(records),
        total_published_unresolved_rows=len(unresolved),
        distinct_normalized_pair_keys=len(pairs),
        releases=input_status,
        soar_crossrelease_probable_duplicates=sum(r['soar_duplicate_status']=='probable_duplicate' for r in records),
        soar_crossrelease_quadrant_review=sum(r['soar_duplicate_status']=='quadrant_review' for r in records),
        int4_overlap_status=dict(Counter(r['int4_status'] for r in records)) if int4_catalog else 'not_checked',
        research_tier_counts=dict(Counter(r['catalogue_tier'] for r in records)),
        candidate_supplement_positions=len(candidate),
        notes=[
            'SOAR positions are NOT yet inserted in orbit-fit files.',
            'SOAR undetected stars are separate constraints, NEVER zero-separation detections.',
            'Published-release counts are not necessarily unique observation counts.',
            'All cross-release and INT4 duplicate flags are provisional, not provenance proofs.',
            'Historical TSV error headers may vary; all required WDS/Disc/Epoch/PA/Sep fields are strictly validated.',
            'SOAR Epoch and Epoch.Y are treated as Julian years; converted Besselian epochs stored separately.',
            'e_PA/err is tangential internal uncertainty in mas; e_Sep is radial internal uncertainty in mas.',
            'No external calibration uncertainty is invented or inserted into fitting weights.',
            'Some 2024–2025 release rows have earlier observational epochs; use epoch_julian_year.',
        ])
    (output_dir/'integration_summary.json').write_text(json.dumps(stats,indent=2)+'\n',encoding='utf-8')
    return stats


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--soar-dir',type=Path,default=Path('catalog_run/data/soar/2015_2025'))
    p.add_argument('--output-dir',type=Path,default=Path('catalog_run/results/soar_2015_2025'))
    p.add_argument('--int4-catalog',type=Path,default=Path('catalog_run/data/int4_all_downloaded.dat'))
    p.add_argument('--catalogue-index',type=Path,default=Path('catalog_run/results/research_all/catalogue_index.csv'))
    p.add_argument('--epoch-tolerance-days',type=float,default=2.0)
    p.add_argument('--pa-tolerance-deg',type=float,default=0.5)
    p.add_argument('--sep-tolerance-arcsec',type=float,default=0.003)
    p.add_argument('--skip-count-validation',action='store_true',help='Test fixtures only; not recommended for scientific runs')
    args=p.parse_args()
    for path in (args.int4_catalog,args.catalogue_index):
        if not path.is_file():p.error('Missing required input: '+str(path))
    result=execute(args.soar_dir,args.output_dir,args.int4_catalog,args.catalogue_index,
                   args.epoch_tolerance_days,args.pa_tolerance_deg,args.sep_tolerance_arcsec,
                   not args.skip_count_validation)
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
