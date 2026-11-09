#!/usr/bin/env python3
"""Recreate the REVIEW-ONLY 24-course roster from the standard 48-system pool.
Never fabricates or changes measurements; keeps the source catalogue error model.
"""
import argparse,csv,shutil
from pathlib import Path

REPLACEMENTS={1:25,3:28,7:29,20:38}

def rows(p):
    with p.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def write(p,r):
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(r[0]));w.writeheader();w.writerows(r)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pool',type=Path,required=True,help='48-source builder output')
    ap.add_argument('--output',type=Path,required=True,help='New, nonexistent output directory')
    args=ap.parse_args()
    pool=args.pool.resolve();out=args.output.resolve()
    if out.exists():ap.error('Output exists; refusing to overwrite '+str(out))
    if not (pool/'selected_systems.csv').is_file():ap.error('Missing 48-source pool manifest')
    src=rows(pool/'selected_systems.csv');src={int(z['project']):z for z in src}
    if len(src)<48:ap.error('Pool must contain at least 48 selected systems')
    for name in ['fit_inputs','all_curated_observations']:(out/name).mkdir(parents=True,exist_ok=True)
    selected=[];pmap={};pairs=set()
    for target in range(1,25):
        pool_idx=REPLACEMENTS.get(target,target);record=src[pool_idx].copy()
        pmap[pool_idx]=target;record['project']=str(target);record['source_pool_project']=str(pool_idx)
        for name in ('input_file','all_observations_csv'):
            infile=pool/record[name]
            new_name=f'{target:02d}'+infile.name[2:]
            outfile=out/('fit_inputs' if name=='input_file' else 'all_curated_observations')/new_name
            content=infile.read_bytes()
            if target==12 and name=='input_file':
                doc=content.decode('utf-8')
                if doc.count('\nT ')!=1:raise RuntimeError('Project 12 T marker not found uniquely')
                content=doc.replace('\nT ','\n*T ',1).encode('utf-8')
            outfile.write_bytes(content)
            record[name]=str(outfile.relative_to(out))
        selected.append(record);pairs.add((record['wds'],record['discoverer']))
    write(out/'selected_systems.csv',selected)
    for name in ('sb9_candidate_orbit_review.csv','sb9_rv_measurement_audit.csv','sbx_parallax_audit.csv'):
        filtered=[]
        for z in rows(pool/name):
            idx=int(z['project'])
            if idx in pmap:
                z['project']=str(pmap[idx]);filtered.append(z)
        if filtered:write(out/name,filtered)
    observations=[z for z in rows(pool/'selected_observation_audit.csv') if (z['wds'],z['discoverer']) in pairs]
    write(out/'selected_observation_audit.csv',observations)
    shutil.copyfile(pool/'soar_release_coverage.csv',out/'soar_release_coverage.csv')
    (out/'README_REVIEW.txt').write_text(
        'REVIEW ONLY: 24 catalogue-derived systems; 1<-25, 3<-28, 7<-29, 20<-38. '
        'Project 12 has published T fixed (*T) for conditional nearly circular fit. '
        'Measurements and uncertainties copied verbatim. Do not distribute as publication-verified.\n'
    )
    print('Review-only 24 system roster:',out)

if __name__=='__main__':main()
