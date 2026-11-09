#!/usr/bin/env python3
"""Build 24 traceable PySVOrbit teaching datasets from LOCAL INT4, ORB6, SB9 and SOAR.

No synthetic, assumed, inferred or default uncertainties are inserted. All usable
position measurements go to provenance CSVs; unknown uncertainties remain blank.
Fitting .inp files contain only observations with a usable catalogue-reported
separation error. Rows without a reported error are listed as comments as well
as in the corresponding .csv, because PySVOrbit cannot parse a blank sigma.

Usage:
  python build_24_course_orbits.py --int4 /local/int4_all.txt.gz \
      --orb6 /local/orb6orbits.txt --sb9-dir /local/sb9_catalog \
      --soar-dir /local/SOAR2015_2025 --output course_24
"""
from __future__ import annotations
import argparse
import csv
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "builder_support"))
from int4_to_star_csv import read_catalog, make_row
from merge_int4_orb6_v2 import parse_orb6, designation, radec, num, valid_astrometry
from integrate_soar_2015_2025 import (
    RELEASES, tsv as soar_tsv, mrt as soar_mrt,
    load_positions as normalize_soar, duplicate_across_soar,
    byear_to_jd, jyear_to_jd,
)
from merge_int4_orb6_sbx_radec import (
    index_sbx, make_candidates, angular_separation_arcsec,
    precise_orb6_position, number as sbx_number, field_row, rv_date,
    map_component, valid_sigma as valid_sbx_rv_sigma,
)
from types import SimpleNamespace

PRIORITY = [("21313-0947", "BLA9"), ("19091+3436", "CHR84AAAB")]
DISALLOWED = {"<", ">", "U", "R", "G", "O", "S", "V", "X", "e", ":", "P", "F"}
BAD_SIGMA = {"<", ">", "V", "P", "F", "G", ":"}


def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1<<20), b""): h.update(chunk)
    return h.hexdigest()


def model(orbit, t):
    """Return PA/rho using ORB6's P/T in Besselian years and a in arcsec."""
    P, T, e, a = [float(orbit[k]) for k in ("P_years","T_besselian","eccentricity","a_arcsec")]
    W, w, i = [math.radians(float(orbit[k])) for k in ("node","omega","inclination")]
    M = 2*math.pi*((t-T)/P % 1.0)
    E=M
    for j in range(100):
        denom=1-e*math.cos(E)
        step=(E-e*math.sin(E)-M)/denom
        E -= step
        if abs(step)<1e-13: break
    X = math.cos(E)-e
    Y = math.sqrt((1-e)*(1+e))*math.sin(E)
    aa = a*(math.cos(w)*math.cos(W)-math.sin(w)*math.sin(W)*math.cos(i))
    bb = a*(math.cos(w)*math.sin(W)+math.sin(w)*math.cos(W)*math.cos(i))
    ff = a*(-math.sin(w)*math.cos(W)-math.cos(w)*math.sin(W)*math.cos(i))
    gg = a*(-math.sin(w)*math.sin(W)+math.cos(w)*math.cos(W)*math.cos(i))
    x,y=aa*X+ff*Y,bb*X+gg*Y
    return math.degrees(math.atan2(y,x))%360, math.hypot(x,y)


def get_orbit(index, ids, limit):
    names={designation(ids["name_2"])}
    choices=[o for o in index.get(ids["wds"],[]) if designation(o["discoverer"]) in names
             and o["grade"].isdigit() and 1 <= int(o["grade"]) <= limit]
    choices=[o for o in choices if all(num(o[k]) is not None for k in
             ("P_years","T_besselian","eccentricity","a_arcsec",
              "inclination","node","omega"))]
    choices=[o for o in choices if float(o["P_years"])>0 and float(o["a_arcsec"])>0
             and 0<=float(o["eccentricity"])<1]
    if not choices:return None,"missing_or_unusable_matching_grade"
    choices.sort(key=lambda x:(int(x["grade"]),x["reference"]))
    if len(choices)>1 and choices[0]["grade"]==choices[1]["grade"]:
        return None,"equal_grade_orbit_ambiguity"
    return choices[0],"ok"


def phase_coverage(rows, P, T):
    phase=sorted(set(round((row["epoch"]-T)/P % 1,7) for row in rows))
    if len(phase)<2: return 0,0.0
    gap=max([phase[j+1]-phase[j] for j in range(len(phase)-1)]+[phase[0]+1-phase[-1]])
    bins=len({min(7,int(x*8)) for x in phase})
    return bins,1-gap


def load_soar(master, source_dir, allow_incomplete=False):
    """Read local, published SOAR positions; no downloading or invented values.

    ``master`` is an already normalized all_soar_positions.csv exported by the
    bundled integrator. ``source_dir`` holds its original VizieR TSV / AAS MRT
    releases. Only exact WDS + physical pair-designation matches are allowed.
    """
    if master and source_dir:
        raise ValueError('Use --soar-master OR --soar-dir, not both')
    rows = []
    if master:
        with master.open(newline='', encoding='utf-8-sig') as fh:
            rows = list(csv.DictReader(fh))
        mandatory = {'wds', 'pair_key', 'epoch_besselian_year',
                     'position_angle_deg', 'separation_arcsec',
                     'separation_sigma_internal_arcsec'}
        if not rows or not mandatory.issubset(rows[0]):
            raise ValueError('SOAR master must be the normalized all_soar_positions.csv')
    elif source_dir:
        # Check every expected release before reading any data, so missing
        # 2019/2023 files are reported together rather than as one-at-a-time
        # failures deep inside the parser.
        sources = {
            release: source_dir / (f'soar_{release}.tsv' if release != '2024_2025'
                                    else 'soar_2024_2025_positions.txt')
            for release in RELEASES
        }
        missing = [str(path) for path in sources.values() if not path.is_file()]
        if missing:
            message = ('Missing SOAR releases (2015-2025 coverage is incomplete):\n  '
                       + '\n  '.join(missing))
            if not allow_incomplete:
                raise FileNotFoundError(
                    message + '\nLocate the releases or use --soar-master '
                    'with a verified complete all_soar_positions.csv. '
                    'For an explicitly partial run, add --allow-incomplete-soar.')
            print('WARNING:', message, file=sys.stderr)
            print('WARNING: --allow-incomplete-soar selected: these releases are '
                  'not included, and results must NOT be labelled complete.',
                  file=sys.stderr)
        for release, path in sources.items():
            if not path.is_file():
                continue
            original = soar_mrt(path) if release == '2024_2025' else soar_tsv(path)
            expected = RELEASES[release][0]
            if len(original) != expected:
                raise ValueError(f'SOAR {release}: {len(original)} rows rather than '
                                 f'the expected {expected}; verify the local source file')
            rows.extend(normalize_soar(row, release, path, i)
                        for i, row in enumerate(original, 1))
    if not rows:
        return {}
    # Imported normalized masters may already have these screening flags.
    # Independently recompute cross-release candidates against the local rows.
    for row in rows:
        row['soar_duplicate_status'] = ''
        row['soar_duplicate_of'] = ''
    duplicate_across_soar(rows, days=2.0, angle=0.5, sep=0.003)
    index = defaultdict(list)
    for row in rows:
        key = (row.get('wds', ''), designation(row.get('pair_key','')))
        if key[0] and key[1]:
            index[key].append(row)
    return index


def append_soar(ids, orbit, position_rows, audit, seen_geometry,
                soar_index, allow_soar_internal_sigma):
    """Conservative SOAR supplement; never silently rotate PA or alter sigma."""
    key = (ids['wds'], designation(ids['name_2']))
    records = soar_index.get(key, [])
    for row in records:
        epoch, theta, rho = (num(row.get(n)) for n in
                             ('epoch_besselian_year','position_angle_deg','separation_arcsec'))
        catalog_sigma = num(row.get('separation_sigma_internal_arcsec',''))
        if catalog_sigma is None or catalog_sigma <= 0: catalog_sigma = None
        status = ''
        if row.get('soar_duplicate_status') != 'not_identified':
            status = 'soar_crossrelease_duplicate_or_quadrant_review'
        if (epoch is None or theta is None or rho is None or
                not 0 <= theta <= 360 or rho <= 0):
            status = 'invalid_soar_position'
        if not status:
            # A near-identical epoch and geometry may be one published measure
            # reported both in INT4 and SOAR: do not count it independently.
            for known in position_rows:
                if abs(byear_to_jd(epoch)-byear_to_jd(known['epoch'])) > 2.0:
                    continue
                diff = (theta - known['theta'] + 180) % 360 - 180
                if abs(rho - known['rho']) <= 0.003 and abs(diff) <= 0.5:
                    status = 'soar_probable_int4_duplicate'
                    break
                if abs(rho - known['rho']) <= 0.003 and abs(abs(diff)-180) <= 0.5:
                    status = 'soar_int4_quadrant_review'
                    break
                status = 'soar_same_epoch_discrepancy_review'
                break
        pa_res = rho_res = None
        if not status:
            pa_pred, rho_pred = model(orbit, epoch)
            pa_res = (theta - pa_pred + 180) % 360 - 180
            rho_res = rho - rho_pred
            if abs(abs(pa_res) - 180) < 15:
                status = 'soar_quadrant_review'
            elif abs(pa_res) > 35:
                status = 'soar_angle_mismatch_review'
            elif abs(rho_res) > max(0.03, 0.35 * rho_pred):
                status = 'soar_separation_mismatch_review'
        keygeom = (round(epoch, 5) if epoch is not None else None,
                   round(theta, 4) if theta is not None else None,
                   round(rho, 6) if rho is not None else None)
        if not status and keygeom in seen_geometry:
            status = 'soar_duplicate_geometry_epoch'
        if not status: seen_geometry.add(keygeom)
        if not status:
            if catalog_sigma is None:
                status = 'selected_soar_missing_sigma'
            elif allow_soar_internal_sigma:
                status = 'selected_soar_internal_sigma'
            else:
                status = 'selected_soar_internal_sigma_audit_only'
        provenance = (row.get('citation') or 'SOAR') + ':' + str(row.get('release') or '')
        audit.append(dict(
            wds=ids['wds'],discoverer=ids['name_2'],
            epoch_B=row.get('epoch_besselian_year',''),
            theta_deg=row.get('position_angle_deg',''),rho_arcsec=row.get('separation_arcsec',''),
            catalog_sigma_arcsec='' if catalog_sigma is None else catalog_sigma,
            sigma_flag='',
            sigma_status=('SOAR_REPORTED_INTERNAL_RADIAL' if catalog_sigma is not None
                          else 'MISSING_OR_UNUSABLE'),
            reference=provenance,technique='SOAR',status=status,
            pa_residual_deg=round(pa_res,3) if pa_res is not None else '',
            rho_residual_arcsec=round(rho_res,5) if rho_res is not None else '',
            source='SOAR',release=row.get('release',''),
            source_row_number=row.get('source_row_number','')))
        if status.startswith('selected_'):
            position_rows.append(dict(
                epoch=epoch,theta=theta,rho=rho,
                sigma=(catalog_sigma if allow_soar_internal_sigma else None),
                catalog_sigma=catalog_sigma,
                sigma_status=('SOAR_INTERNAL_RADIAL_AS_SINGLE_SIGMA' if catalog_sigma is not None
                              and allow_soar_internal_sigma else
                              'SOAR_INTERNAL_RADIAL_NOT_USED' if catalog_sigma is not None
                              else 'MISSING_OR_UNUSABLE'),
                ref=provenance,technique='SOAR',source='SOAR',
                release=row.get('release',''),source_row_number=row.get('source_row_number','')))


def screen(ids, obs, orbit, min_points, min_phase_bins, min_coverage,
           min_fit_points, min_fit_phase_bins,
           soar_index=None, allow_soar_internal_sigma=False):
    retained = []
    audit = []
    keys = set()
    for raw in obs:
        row = make_row(ids, raw, True)
        reason = ''
        epoch = num(row['epoch_besselian'])
        theta = num(row['theta_deg'])
        rho = num(row['rho_arcsec'])
        if epoch is None or theta is None or rho is None or not (0 <= theta <= 360) or rho <= 0:
            reason = 'no_resolved_astrometry'
        elif row.get('theta_flag', '').strip() in DISALLOWED or row.get('rho_flag', '').strip() in DISALLOWED:
            reason = 'flagged_nonstandard_or_ambiguous'
        elif not valid_astrometry(row):
            reason = 'invalid_astrometry_flags'
        else:
            key = (round(epoch, 5), round(theta, 4), round(rho, 6))
            if key in keys:
                reason = 'duplicate_geometry_epoch'
            else:
                keys.add(key)
        sigma = num(row['rho_error_arcsec'])
        sigflag = row.get('rho_error_flag', '').strip()
        reported = bool(sigma is not None and sigma > 0 and sigflag not in BAD_SIGMA)
        pa_res, rho_res = None, None
        if not reason:
            pa_pred, rho_pred = model(orbit, epoch)
            pa_res = (theta - pa_pred + 180) % 360 - 180
            rho_res = rho - rho_pred
            if abs(abs(pa_res) - 180) < 15:
                reason = 'quarantined_180_quadrant'
            elif abs(pa_res) > 35:
                reason = 'angle_mismatch_review'
            elif abs(rho_res) > max(0.03, .35 * rho_pred):
                reason = 'separation_mismatch_review'
        status = reason or ('selected_reported_sigma' if reported else 'selected_missing_sigma')
        audit.append(dict(
            wds=ids['wds'], discoverer=ids['name_2'],
            epoch_B=row['epoch_besselian'], theta_deg=row['theta_deg'],
            rho_arcsec=row['rho_arcsec'],
            # Never transform absent / unusable uncertainties into a number.
            catalog_sigma_arcsec=(row['rho_error_arcsec'] if reported else ''),
            sigma_flag=row.get('rho_error_flag', ''),
            sigma_status='CATALOG_REPORTED' if reported else 'MISSING_OR_UNUSABLE',
            reference=row['reference'], technique=row['technique'], status=status,
            pa_residual_deg=round(pa_res, 3) if pa_res is not None else '',
            rho_residual_arcsec=round(rho_res, 5) if rho_res is not None else '',
            source='INT4', release='', source_row_number=''))
        if reason:
            continue
        retained.append(dict(epoch=epoch, theta=theta, rho=rho,
                             sigma=(sigma if reported else None), catalog_sigma=(sigma if reported else None),
                             source='INT4', release='', source_row_number='',
                             sigma_status=('CATALOG_REPORTED' if reported else 'MISSING_OR_UNUSABLE'),
                             ref=row['reference'], technique=row['technique']))
    if soar_index:
        append_soar(ids, orbit, retained, audit, keys, soar_index, allow_soar_internal_sigma)
    if len(retained) < min_points:
        return None, audit, f'only_{len(retained)}_resolved_points'
    P = float(orbit['P_years']); T = float(orbit['T_besselian'])
    bins, coverage = phase_coverage(retained, P, T)
    if bins < min_phase_bins or coverage < min_coverage:
        return None, audit, f'limited_all_position_coverage_bins{bins}_fraction{coverage:.3f}'
    fit_rows = [r for r in retained if r['sigma'] is not None]
    fit_bins, fit_coverage = phase_coverage(fit_rows, P, T)
    if len(fit_rows) < min_fit_points or fit_bins < min_fit_phase_bins:
        return None, audit, (f'insufficient_catalog_uncertainties_{len(fit_rows)}_positions_'
                             f'{fit_bins}_phase_bins')
    deviations = [abs(float(r['pa_residual_deg'])) for r in audit
                  if r['status'].startswith('selected_')]
    if deviations and sorted(deviations)[len(deviations)//2] > 10:
        return None, audit, 'poor_orbit_position_agreement'
    result = dict(ids=ids, orbit=orbit, rows=retained, fit_rows=fit_rows,
                  phase_bins=bins, phase_coverage=coverage,
                  fit_phase_bins=fit_bins, fit_phase_coverage=fit_coverage,
                  n_published=len(fit_rows),
                  median_abs_pa=sorted(deviations)[len(deviations)//2] if deviations else 0)
    return result, audit, 'ok'



def _sbx_pair_suffix(discoverer):
    """Extract component(s) from an ORB6 discoverer designation.

    Examples: CHR 84Aa,Ab -> AAAB, STF 2474 AB -> AB. When the
    designator has no component suffix, the WDS default is AB.
    """
    m = re.match(r'^\s*[A-Za-z]+\s*[0-9]+\s*(.*?)\s*$', discoverer or '')
    if not m:
        return None
    suffix = designation(m.group(1))
    return suffix or 'AB'


def match_sbx_parallax(orbit, alias_index, systems, max_separation_arcsec=30.0):
    """Use SBX system-level Plx only for one uniquely identified object.

    No period/eccentricity/RV match is required, because this is a stellar
    parallax rather than a visual/spectroscopic orbit match. Nevertheless a
    shared HIP/HD designation alone is insufficient in hierarchical multiples:
    check alias uniqueness, position and explicitly conflicting components.
    Return (validated SBX source record or None, source-provenance audit).
    """
    result = dict(status='', sbx_seq='', sbx_component='',
                  orb6_component=_sbx_pair_suffix(orbit.get('discoverer','')) or '',
                  matched_identifiers='', angular_separation_arcsec='',
                  parallax_mas='', parallax_error_mas='',
                  parallax_reference='', source='B/sbx/systems',
                  reason='')
    aliases=[]
    for prefix in ('HIP','HD'):
        raw=str(orbit.get(prefix.lower(),'') or '').strip()
        m=re.search(r'\d+',raw)
        if m:
            aliases.append((prefix,str(int(m.group()))))
    if not aliases:
        result.update(status='unavailable_no_ORB6_HIP_HD',
                      reason='No HIP or HD identifier for reliable SBX crossmatching')
        return None,result

    # Require a single SBX system for the union of all available identifiers;
    # multiple systems remain ambiguous even if only one reports parallax.
    matches={key:alias_index.get(key,set()) for key in aliases}
    seqs=set().union(*matches.values())
    if not seqs:
        result.update(status='unavailable_no_sbx_alias',
                      reason='No SBX HIP/HD alias matches the ORB6 system')
        return None,result
    if len(seqs)!=1:
        result.update(status='review_multiple_sbx_systems',
                      reason='Aliases identify multiple SBX systems: '+','.join(sorted(seqs)))
        return None,result
    seq=next(iter(seqs))
    matched=[f'{k}:{v}' for (k,v),members in matches.items() if seq in members]
    result['sbx_seq']=seq
    result['matched_identifiers']=';'.join(matched)
    src=systems.get(seq)
    if not src:
        result.update(status='unavailable_sbx_system_record_missing',
                      reason='Matched SBX Seq has no systems.dat record')
        return None,result
    result['sbx_component']=src.get('component','')
    result['parallax_mas']=src.get('parallax','')
    result['parallax_error_mas']=src.get('parallax_error','')
    result['parallax_reference']=src.get('parallax_ref','')

    # Explicit incompatible component labels are an unresolved physical-pair
    # mismatch: a valid system-level parallax cannot be silently attributed to
    # the wrong component in a hierarchical multiple.
    orb_comp=result['orb6_component']
    sbx_comp=designation(src.get('component',''))
    if orb_comp and sbx_comp and orb_comp != sbx_comp:
        result.update(status='review_component_mismatch',
                      reason=f'ORB6 pair {orb_comp} differs from SBX component {sbx_comp}')
        return None,result

    pos=precise_orb6_position(orbit.get('coord',''))
    ra=sbx_number(src.get('ra_deg'))
    dec=sbx_number(src.get('dec_deg'))
    if pos is None or ra is None or dec is None or not (0<=ra<360 and -90<=dec<=90):
        result.update(status='review_missing_precise_coordinates',
                      reason='Precise ORB6/SBX coordinates unavailable for matching')
        return None,result
    sep=angular_separation_arcsec(*pos,ra,dec)
    result['angular_separation_arcsec']=f'{sep:.4f}'
    if sep>max_separation_arcsec:
        result.update(status='review_coordinate_disagreement',
                      reason=(f'ORB6/SBX angular separation {sep:.2f} arcsec exceeds '
                              f'{max_separation_arcsec:g} arcsec threshold'))
        return None,result

    px=sbx_number(src.get('parallax'))
    if px is None or px<=0:
        result.update(status='unavailable_nonpositive_or_missing_parallax',
                      reason='SBX Plx is missing, non-finite or nonpositive')
        return None,result
    error_token=str(src.get('parallax_error','') or '').strip()
    error=sbx_number(error_token)
    if error_token and (error is None or error<0):
        result.update(status='review_invalid_parallax_error',
                      reason='SBX e_Plx exists but is not finite and nonnegative')
        return None,result
    result.update(status='accepted_unique_sbx_system',
                  reason='Unique HIP/HD alias, consistent component when specified, '
                         'and precise positional agreement')
    return src,result



# These columns describe *observed* SBX RVs, not reconstructed model velocities.
SBX_RV_FIELDS = [
    'project', 'wds', 'discoverer', 'sbx_seq', 'sbx_orbit',
    'sbx_orbit_reference', 'candidate_status', 'matching_identifiers',
    'orb6_component', 'sbx_component', 'component_match_status',
    'sbx_period_days', 'sbx_eccentricity', 'sbx_K1_kms', 'sbx_K2_kms',
    'sbx_V0_kms', 'source_line_number', 'rv_epoch_original',
    'rv_epoch_reference', 'rv_epoch_jd_minus_2400000_candidate',
    'rv_epoch_conversion_status', 'rv_kms', 'rv_error_kms',
    'rv_weight', 'rv_component_original', 'rv_component_pysvorbit_candidate',
    'rv_source', 'rv_legacy_offset', 'rv_comment',
    'numeric_screening', 'included_in_inp',
]


def sbx_rv_audit_row(candidate, row, source_line_number):
    """Preserve SBX velocities verbatim, never guess missing errors or bodies.

    9"a"/"b" is only a *candidate* PySVOrbit mapping; this method
    does NOT decide if an SBX spectroscopic pair is the ORB6 visual pair.
    """
    orb = candidate['sbx_orbit_record']
    expected_component = candidate['orb6_component']
    observed_component = candidate['sbx_component']
    from_component = designation(observed_component) if observed_component else ''
    pair_status = (
        'explicit_component_agreement_STILL_REVIEW'
        if from_component and from_component == expected_component else
        'explicit_component_disagreement_REVIEW'
        if from_component and expected_component else
        'missing_or_unparsed_component_REVIEW'
    )
    epoch_raw = row.get('epoch', '')
    epoch_ref = orb.get('time_ref', '')
    epoch_converted, date_note = rv_date(epoch_raw, epoch_ref)
    mapped = map_component(row.get('component', ''))
    sigma = valid_sbx_rv_sigma(row.get('rv_error', ''))
    offset = sbx_number(row.get('legacy_offset', ''))
    weight = sbx_number(row.get('weight', ''))
    notes = []
    if sbx_number(epoch_raw) is None or not epoch_converted:
        notes.append('unusable_or_unknown_time_reference')
    if 'CHECK_SOURCE' in date_note:
        notes.append('short_JD_interpretation_requires_publication_review')
    if sbx_number(row.get('rv', '')) is None:
        notes.append('missing_or_invalid_rv')
    if sigma is None:
        notes.append('missing_or_invalid_reported_rv_error')
    if not mapped:
        notes.append('component_not_mapped_to_Va_or_Vb')
    if offset is not None and abs(offset) > 1e-9:
        notes.append('nonzero_uncorrected_legacy_offset')
    if weight is not None and weight <= 0:
        notes.append('nonpositive_catalogue_weight')
    if not notes:
        notes.append('numerically_complete_BUT_identity_and_RV_zero_point_unverified')
    return dict(
        project=candidate['project'], wds=candidate['wds'],
        discoverer=candidate['discoverer'], sbx_seq=candidate['sbx_seq'],
        sbx_orbit=candidate['sbx_orbit'],
        sbx_orbit_reference=orb.get('bibcode', ''),
        candidate_status=candidate['status'],
        matching_identifiers=candidate['matching_identifiers'],
        orb6_component=expected_component, sbx_component=observed_component,
        component_match_status=pair_status,
        sbx_period_days=orb.get('period_d', ''),
        sbx_eccentricity=orb.get('ecc', ''),
        sbx_K1_kms=orb.get('K1', ''),
        sbx_K2_kms=orb.get('K2', ''),
        sbx_V0_kms=orb.get('V0', ''),
        source_line_number=source_line_number,
        rv_epoch_original=epoch_raw,
        rv_epoch_reference=epoch_ref,
        rv_epoch_jd_minus_2400000_candidate=epoch_converted,
        rv_epoch_conversion_status=date_note,
        rv_kms=row.get('rv', ''), rv_error_kms=row.get('rv_error', ''),
        rv_weight=row.get('weight', ''),
        rv_component_original=row.get('component', ''),
        rv_component_pysvorbit_candidate=mapped,
        rv_source=row.get('vel_ref', ''),
        rv_legacy_offset=row.get('legacy_offset', ''),
        rv_comment=row.get('comment', ''),
        numeric_screening=';'.join(notes),
        included_in_inp='no_pending_source_and_component_review',
    )


def export_sbx_rv_measurements(sbx_folder, candidate_records, destination):
    """Export every original SBX velocities.dat row for candidate Seq/O pairs.

    Re-scan the original fixed-width catalogue so even records previously
    skipped by index_sbx (invalid epoch or velocity) remain visible to reviewers.
    No values are populated by defaults or copied between SBX orbits.
    """
    matches = defaultdict(list)
    for candidate in candidate_records:
        matches[(candidate['sbx_seq'], candidate['sbx_orbit'])].append(candidate)
    count = 0
    complete = 0
    with destination.open('w', newline='', encoding='utf-8') as outf:
        writer = csv.DictWriter(outf, fieldnames=SBX_RV_FIELDS)
        writer.writeheader()
        with (sbx_folder/'velocities.dat').open(encoding='latin-1', errors='replace') as inp:
            for source_line_number, line in enumerate(inp, 1):
                row = field_row(line, 'velocities')
                for candidate in matches.get((row['seq'],row['orbit']), ()):
                    record = sbx_rv_audit_row(candidate, row, source_line_number)
                    writer.writerow(record)
                    count += 1
                    if record['numeric_screening'].startswith('numerically_complete'):
                        complete += 1
    return count, complete


def write_input(entry):
    """Valid PySVOrbit file: only published numerical sigma values are active.

    Missing-sigma positions are retained as C comments and in the complete
    observations CSV, never passed to the backend as invented precision.
    """
    ids, orb, rows = entry['ids'], entry['orbit'], entry['rows']
    ra, dec = radec(orb['coord'])
    lines = [
        f"C ORB6 grade {orb['grade']} reference {orb['reference']}; WDS {ids['wds']} {ids['name_2']}",
        'C Astrometry: Fourth Interferometric Catalog (INT4) and optionally SOAR.',
        'C No missing uncertainty is guessed; see per-source CSV provenance.',
        'C SOAR internal radial sigma is only fitted with explicit user opt-in.',
        'C Positions with unknown sigma are comments, not active fitting data.',
        'C Source observations and orbit parameters are unchanged.',
        f"Object: {ids['wds']} {ids['name_2']}",
        f'R.A.: {ra}' if ra else 'C RA unavailable',
        f'Dec: {dec}' if dec else 'C Dec unavailable',
        'Period_Unit: years',
    ]
    # Parallax is read from SBX systems.dat only if the independent identity
    # audit has accepted that exact system. Leave absent otherwise: never 0.
    parallax = entry.get('sbx_parallax_record')
    px_audit = entry.get('sbx_parallax_audit',{})
    if parallax is not None:
        lines.append('C Parallax: SBX B/sbx/systems.dat Seq '
                     + str(px_audit['sbx_seq']))
        lines.append('C Parallax reference: '
                     + (str(parallax.get('parallax_ref') or 'not supplied')))
        lines.append('C Parallax uncertainty (mas): '
                     + (str(parallax.get('parallax_error') or 'not supplied')))
        lines.append('C SBX parallax match: '
                     + str(px_audit.get('matched_identifiers') or ''))
        lines.append('Parallax: ' + str(parallax['parallax']))
    else:
        lines.append('C Parallax unavailable or unmatched in SBX; no value inserted.')
    for name, field in [('P','P_years'),('T','T_besselian'),('e','eccentricity'),
                        ('a','a_arcsec'),('W','node'),('w','omega'),('i','inclination')]:
        lines.append(f'{name} {orb[field]}')
    lines.append('C No RV measurements: K1, K2 and V0 are unspecified, not estimated.')
    for r in rows:
        ref = (r['ref'] or 'INT4').replace(' ', '_')
        if r['sigma'] is None:
            if r.get('catalog_sigma') is None:
                lines.append('C MISSING_SIGMA (not fitted)  epoch theta rho sigma(blank) I1 reference')
                lines.append(f"C {r['epoch']:.5f} {r['theta']:.4f} {r['rho']:.6f} [blank] I1 {ref}")
            else:
                lines.append('C SOAR INTERNAL SIGMA retained in CSV; not used for fitting without explicit opt-in')
                lines.append(f"C {r['epoch']:.5f} {r['theta']:.4f} {r['rho']:.6f} {r['catalog_sigma']:.7f} I1 {ref}")
        else:
            lines.append(f"C {r.get('source','INT4')} reference={ref}; sigma={r['sigma_status']}")
            lines.append(f"{r['epoch']:.5f} {r['theta']:.4f} {r['rho']:.6f} {r['sigma']:.6f} I1 {ref}")
    return '\n'.join(lines) + '\n'


def write_observations_csv(path, entry):
    """Full curated positions; truly missing uncertainties remain blank."""
    fields = ['wds','discoverer','epoch_B','theta_deg','rho_arcsec',
              'rho_sigma_arcsec','sigma_status','source_reference','technique',
              'used_in_inp_fit','source','release','source_row_number']
    records = [dict(wds=entry['ids']['wds'], discoverer=entry['ids']['name_2'],
                    epoch_B=f"{r['epoch']:.5f}",theta_deg=r['theta'],rho_arcsec=r['rho'],
                    rho_sigma_arcsec=r.get('catalog_sigma') if r.get('catalog_sigma') is not None else '',
                    sigma_status=r['sigma_status'],source_reference=r['ref'],
                    technique=r['technique'],
                    used_in_inp_fit='yes' if r['sigma'] is not None else 'no',
                    source=r.get('source','INT4'),release=r.get('release',''),
                    source_row_number=r.get('source_row_number',''))
               for r in entry['rows']]
    write_csv(path,records,fields)


def write_csv(path, rows, fields):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore");w.writeheader();w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--int4',type=Path,required=True,
                    help='Existing local INT4 .txt or .txt.gz; NO downloads')
    ap.add_argument('--orb6',type=Path,required=True,
                    help='Existing local ORB6 orb6orbits.txt; NO downloads')
    ap.add_argument('--soar-master',type=Path,
                    help='Existing normalized all_soar_positions.csv (2015-2025)')
    ap.add_argument('--soar-dir',type=Path,
                    help='Local published SOAR 2015-2025 TSV/MRT release directory')
    ap.add_argument('--allow-incomplete-soar', action='store_true',
                    help='EXPLICIT OPT-IN to omit absent SOAR releases; output is flagged incomplete')
    ap.add_argument('--sb9-dir',type=Path,
                    help='Local SB9/SBX fixed-width folder with alias.dat, systems.dat, orbits.dat, velocities.dat')
    ap.add_argument('--sbx-parallax-max-separation-arcsec',type=float,default=30.0,
                    help='Maximum separation between precise ORB6 and SBX coordinates for parallax crossmatch (default: 30 arcsec)')
    ap.add_argument('--use-soar-internal-radial-as-sigma',action='store_true',
                    help='EXPLICIT OPT-IN: use SOAR e_Sep radial internal error as PySVOrbit single sigma; model approximation')
    ap.add_argument('--output',type=Path,default=Path('course_24'),
                    help='New destination folder; must not exist')
    ap.add_argument('--count',type=int,default=24)
    ap.add_argument('--max-grade',type=int,default=2,choices=[1,2,3])
    ap.add_argument('--min-points',type=int,default=16,
                    help='Minimum independent geometry measurements regardless of sigma')
    ap.add_argument('--min-phase-bins',type=int,default=5)
    ap.add_argument('--min-coverage',type=float,default=0.65)
    ap.add_argument('--min-fit-points',type=int,default=6,
                    help='Minimum observations with actual catalog-reported usable uncertainty')
    ap.add_argument('--min-fit-phase-bins',type=int,default=3)
    opts = ap.parse_args()
    if opts.count < 2: ap.error('--count must be >=2 for the two priority systems')
    if not 0 < opts.sbx_parallax_max_separation_arcsec <= 120:
        ap.error('--sbx-parallax-max-separation-arcsec must be in (0, 120]')
    if opts.min_fit_points < 4: ap.error('--min-fit-points must be at least 4 for seven visual elements')
    if opts.output.exists(): ap.error(f'Destination already exists: {opts.output}')
    for name,path in (('INT4',opts.int4),('ORB6',opts.orb6)):
        if not path.is_file():
            ap.error(f'{name} catalogue not found: {path}. Use your existing local catalogue path.')
    if opts.soar_master and not opts.soar_master.is_file():
        ap.error(f'Missing SOAR master: {opts.soar_master}')
    if opts.soar_dir and not opts.soar_dir.is_dir():
        ap.error(f'Missing SOAR source directory: {opts.soar_dir}')
    if opts.sb9_dir:
        for fname in ('alias.dat','systems.dat','orbits.dat','velocities.dat'):
            if not (opts.sb9_dir/fname).is_file():
                ap.error(f'Missing SB9/SBX source file: {opts.sb9_dir/fname}')
    soar_index = load_soar(opts.soar_master,opts.soar_dir,
                           allow_incomplete=opts.allow_incomplete_soar)
    print('Local SOAR pairs:',len(soar_index))
    print('Reading local ORB6 catalogue...')
    orb6 = parse_orb6(opts.orb6)
    selected_pool = {}
    status = []
    audits = {}
    for ids, observations in read_catalog(opts.int4):
        wds = ids.get('wds','')
        if wds not in orb6: continue
        orbit,reason = get_orbit(orb6,ids,opts.max_grade)
        if not orbit: continue
        result,audit,reason = screen(ids,observations,orbit,opts.min_points,
                                     opts.min_phase_bins,opts.min_coverage,
                                     opts.min_fit_points,opts.min_fit_phase_bins,
                                     soar_index, opts.use_soar_internal_radial_as_sigma)
        key=(wds,designation(ids['name_2']))
        audits[key]=audit
        if result:selected_pool[key]=result
        status.append(dict(wds=wds,designation=ids['name_2'],grade=orbit['grade'],
                           reference=orbit['reference'],raw_measurements=len(observations),
                           retained_all_positions=sum(x['status'].startswith('selected_') for x in audit),
                           retained_reported_sigma=sum(x['status']=='selected_reported_sigma' for x in audit),
                           fit_phase_bins=result['fit_phase_bins'] if result else '',
                           geometry_phase_bins=result['phase_bins'] if result else '',
                           status=reason))
    chosen=[]
    for key in PRIORITY:
        if key in selected_pool:
            chosen.append((key,selected_pool[key]))
        else:
            print(f'WARNING: priority {key} did not qualify. See screening report.')
    remaining=[(k,v) for k,v in selected_pool.items() if k not in set(PRIORITY)]
    remaining.sort(key=lambda kv:(int(kv[1]['orbit']['grade']),
                                   -kv[1]['n_published'],
                                   -kv[1]['fit_phase_bins'],
                                   -len(kv[1]['rows']),
                                   -kv[1]['phase_coverage'],kv[0]))
    chosen += remaining[:max(0,opts.count-len(chosen))]
    if len(chosen) > opts.count: chosen=chosen[:opts.count]
    sbx_data=None
    if opts.sb9_dir:
        print('Indexing SBX once for parallax and advisory RV matching...')
        sbx_data=index_sbx(opts.sb9_dir)
        alias_index, systems, _spec_orbits, _velocities, _counts=sbx_data
        for key, entry in chosen:
            src, px_audit=match_sbx_parallax(entry['orbit'],alias_index,systems,
                max_separation_arcsec=opts.sbx_parallax_max_separation_arcsec)
            entry['sbx_parallax_record']=src
            entry['sbx_parallax_audit']=px_audit
    else:
        for key, entry in chosen:
            entry['sbx_parallax_record']=None
            entry['sbx_parallax_audit']={'status':'not_checked_no_sbx_directory',
                'sbx_seq':'','matched_identifiers':'','angular_separation_arcsec':'',
                'parallax_mas':'','parallax_error_mas':'','parallax_reference':'',
                'reason':'No --sb9-dir supplied'}
    opts.output.mkdir(parents=True)
    if opts.soar_dir:
        coverage = [
            {'release': release,
             'source_file': str(opts.soar_dir / (f'soar_{release}.tsv'
                  if release != '2024_2025' else 'soar_2024_2025_positions.txt')),
             'status': 'PRESENT' if (opts.soar_dir / (f'soar_{release}.tsv'
                  if release != '2024_2025' else 'soar_2024_2025_positions.txt')).is_file()
                  else 'MISSING_NOT_USED'}
            for release in RELEASES
        ]
        write_csv(opts.output/'soar_release_coverage.csv', coverage,
                  ['release','source_file','status'])
    (opts.output/'fit_inputs').mkdir()
    (opts.output/'all_curated_observations').mkdir()
    manifest=[];details=[]
    for number,(key,data) in enumerate(chosen,1):
        orbit=data['orbit']; stem=(f"{number:02d}_{key[0]}_"+
            re.sub(r'[^a-zA-Z0-9]+','_',data['ids']['name_2']).strip('_'))
        inp='fit_inputs/'+stem+'.inp'
        observed='all_curated_observations/'+stem+'.csv'
        (opts.output/inp).write_text(write_input(data),encoding='utf-8')
        write_observations_csv(opts.output/observed,data)
        manifest.append(dict(project=number,wds=key[0],discoverer=data['ids']['name_2'],
                             input_file=inp,all_observations_csv=observed,
                             grade=orbit['grade'],orb6_reference=orbit['reference'],
                             n_curated_positions=len(data['rows']),
                             n_fit_positions_with_catalog_sigma=len(data['fit_rows']),
                             n_missing_sigma=len(data['rows'])-len(data['fit_rows']),
                             all_geometry_phase_bins=data['phase_bins'],
                             fit_phase_bins=data['fit_phase_bins'],
                             phase_coverage_all=round(data['phase_coverage'],3),
                             phase_coverage_fitted=round(data['fit_phase_coverage'],3),
                             initial_median_absolute_pa_residual_deg=data['median_abs_pa'],
                             n_soar_positions=sum(r.get('source')=='SOAR' for r in data['rows']),
                             n_soar_in_fit=sum(r.get('source')=='SOAR' for r in data['fit_rows']),
                             sb9_status='not_checked' if not opts.sb9_dir else 'review_only',
                             sbx_parallax_status=data['sbx_parallax_audit']['status'],
                             sbx_parallax_mas=(data['sbx_parallax_record'] or {}).get('parallax',''),
                             sbx_parallax_error_mas=(data['sbx_parallax_record'] or {}).get('parallax_error',''),
                             sbx_parallax_reference=(data['sbx_parallax_record'] or {}).get('parallax_ref',''),
                             sbx_system_seq=data['sbx_parallax_audit'].get('sbx_seq',''),
                             sbx_parallax_match=data['sbx_parallax_audit'].get('matched_identifiers',''),
                             sbx_parallax_separation_arcsec=data['sbx_parallax_audit'].get('angular_separation_arcsec','')))
        details.extend(audits[key])
    if opts.sb9_dir:
        alias_index, systems, spec_orbits, velocities, _ = sbx_data
        review_rows=[]
        rv_measurement_candidates=[]
        matched_by_key={key: data for key,data in chosen}
        for m in manifest:
            chosen_key=(m['wds'],designation(m['discoverer']))
            entry=matched_by_key[chosen_key]
            v=entry['orbit']
            meta=dict(wds=m['wds'],hip=v.get('hip',''),hd=v.get('hd',''),
                      P_years=v['P_years'],e=v['eccentricity'])
            args=SimpleNamespace(allow_wds_only=False,max_period_fraction=0.15,
                                  max_ecc_difference=0.25)
            candidates,reason = make_candidates(meta,alias_index,systems,
                                                spec_orbits,velocities,args)
            m['sb9_status'] = ('one_preliminary_alias_candidate_REVIEW' if len(candidates)==1
                               else 'ambiguous_alternative_orbits_REVIEW' if candidates
                               else 'no_consistent_hip_hd_alias_match')
            if not candidates:
                review_rows.append(dict(project=m['project'],wds=m['wds'],
                   discoverer=m['discoverer'],sb9_system='',sb9_orbit='',n_rv_rows='',
                   period_yr='',eccentricity='',shared_aliases='',status=m['sb9_status'],
                   note=reason))
            for candidate in candidates:
                rv_measurement_candidates.append(dict(
                    project=m['project'], wds=m['wds'],
                    discoverer=m['discoverer'], sbx_seq=candidate['seq'],
                    sbx_orbit=candidate['orbit']['orbit'],
                    sbx_orbit_record=candidate['orbit'],
                    orb6_component=_sbx_pair_suffix(v.get('discoverer','')) or '',
                    sbx_component=candidate.get('component', ''),
                    status=m['sb9_status'],
                    matching_identifiers=';'.join(f'{k}:{value}' for k,value in candidate['keys']),
                ))
                review_rows.append(dict(project=m['project'],wds=m['wds'],
                   discoverer=m['discoverer'],sb9_system=candidate['seq'],
                   sb9_orbit=candidate['orbit']['orbit'],n_rv_rows=len(candidate['rvs']),
                   period_yr=candidate['period_yr'],
                   eccentricity=candidate['orbit'].get('ecc',''),
                   shared_aliases=';'.join(f'{k}:{v}' for k,v in candidate['keys']),
                   status=m['sb9_status'],
                   note='RV not inserted: component identity, date/velocity convention, offsets and uncertainties need source verification'))
        write_csv(opts.output/'sb9_candidate_orbit_review.csv',review_rows,
                  ['project','wds','discoverer','sb9_system','sb9_orbit','n_rv_rows',
                   'period_yr','eccentricity','shared_aliases','status','note'])
        count_rv, count_numerically_complete = export_sbx_rv_measurements(
            opts.sb9_dir, rv_measurement_candidates,
            opts.output/'sb9_rv_measurement_audit.csv')
        print(f'SBX candidate observed RVs exported: {count_rv}; '
              f'numerically complete (still not verified for fitting): '
              f'{count_numerically_complete}')
    px_rows=[]
    for (key,entry),m in zip(chosen,manifest):
        audit=entry['sbx_parallax_audit']
        px_rows.append(dict(project=m['project'],wds=key[0],
            discoverer=m['discoverer'],status=audit['status'],
            sbx_system_seq=audit.get('sbx_seq',''),
            sbx_component=audit.get('sbx_component',''),
            orb6_component=audit.get('orb6_component',''),
            matched_identifiers=audit.get('matched_identifiers',''),
            separation_arcsec=audit.get('angular_separation_arcsec',''),
            sbx_parallax_mas=audit.get('parallax_mas',''),
            sbx_parallax_error_mas=audit.get('parallax_error_mas',''),
            sbx_parallax_reference=audit.get('parallax_reference',''),
            included_in_inp='yes' if entry['sbx_parallax_record'] else 'no',
            reason=audit.get('reason','')))
    write_csv(opts.output/'sbx_parallax_audit.csv',px_rows,
        ['project','wds','discoverer','status','sbx_system_seq','sbx_component',
         'orb6_component','matched_identifiers','separation_arcsec',
         'sbx_parallax_mas','sbx_parallax_error_mas','sbx_parallax_reference',
         'included_in_inp','reason'])
    if manifest:
        write_csv(opts.output/'selected_systems.csv',manifest,list(manifest[0]))
        write_csv(opts.output/'selected_observation_audit.csv',details,list(details[0]))
    write_csv(opts.output/'candidate_screening.csv',status,
              ['wds','designation','grade','reference','raw_measurements',
               'retained_all_positions','retained_reported_sigma','fit_phase_bins',
               'geometry_phase_bins','status'])
    sonar_missing = (opts.soar_dir and any(not (opts.soar_dir / (f'soar_{release}.tsv'
                    if release != '2024_2025' else 'soar_2024_2025_positions.txt')).is_file()
                    for release in RELEASES))
    (opts.output/'README.txt').write_text(
        f"Requested {opts.count} distinct systems; qualified and exported {len(manifest)}.\n"
        f"This is NOT a claim of 24 if fewer qualified.\n"
        f"SOAR coverage: {'INCOMPLETE (missing releases; see soar_release_coverage.csv)' if sonar_missing else 'no missing local release files reported'}.\n"
        'All positions with valid geometry (including those missing sigma) are in each CSV.\n'
        'Missing/invalid catalog uncertainties are truly blank in CSV; NOT zero.\n'
        'PySVOrbit .inp files contain only measurements with catalog-reported\n'
        'positive uncertainties; other positions are C comments, not fit points.\n'
        'No synthetic astrometry, uncertainties or RVs are generated.\n'
        'SBX catalogue-reported positive parallax is inserted only after unique HIP/HD,\n'
        'compatible explicit component, and positional validation; otherwise absent.\n'
        'sbx_parallax_audit.csv gives the source Seq, uncertainty, reference and any exclusion reason.\n'
        'Published ORB6 starting orbits and original INT4/SOAR theta/rho are untouched.\n'
        'SOAR internal radial sigma is preserved; using it as one PySVOrbit sigma is\n'
        'an explicit modelling approximation only if --use-soar-internal-radial-as-sigma.\n'
        'SB9 is a candidate identity/RV crossmatch audit, never automatically fitted.\n'
        'sb9_rv_measurement_audit.csv preserves actual SBX velocities.dat measurements,\n'
        'including the original time, RV, sigma (or blank), source, component,\n'
        'legacy offset, and original source line. No candidate RV enters the .inp fit\n'
        'until the correct physical orbit and source timing/offsets are verified.\n'
        'SOAR unresolved nondetections are NOT treated as resolved positions.\n'
        'Rejected positions/quadrant ambiguity are preserved in the audit.\n'
        'An ORB6 grade does not guarantee this subset constrains all parameters.\n'
        f'INT4 local file SHA256: {sha256(opts.int4)}\n'
        f'ORB6 local file SHA256: {sha256(opts.orb6)}\n'
        f'Require grade<={opts.max_grade}, all positions>={opts.min_points}, '
        f'fit points>={opts.min_fit_points}, fit phase bins>={opts.min_fit_phase_bins}.\n',
        encoding='utf-8')
    print(f'EXPORTED {len(manifest)}/{opts.count} candidates to {opts.output}')
    if len(manifest)<opts.count:
        print('WARNING: fewer than requested systems qualify under the data-only constraints. '
              'Review candidate_screening.csv; do not invent missing values.')

if __name__=='__main__':
    main()
