#!/usr/bin/env python3
"""Stream SVMU2 CRISPR PAM orthology tables without loading them in memory.

Usage:
    python crispr_cfd_summary.py results/07_summary/conserved_pams.csv \
        results/07_summary/failed_tolerance_pams.csv

Requires: pandas, numpy. Does not modify input files.
"""
import argparse
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd


def summarize(path: Path, chunk_size: int):
    headers = set(pd.read_csv(path, nrows=0).columns)
    required = {'CFD_ref', 'CFD_qry'}
    if not required.issubset(headers):
        raise ValueError(f'{path}: missing {sorted(required - headers)}')
    optional = {'protospacerID', 'syntenic_block_ref', 'syntenic_block_qry',
                'proj_dist', 'pam_ref', 'pam_qry'}
    cols = sorted(required | (optional & headers))
    totals = Counter()

    for chunk in pd.read_csv(path, usecols=cols, chunksize=chunk_size, low_memory=False):
        n = len(chunk)
        totals['rows'] += n
        ref = pd.to_numeric(chunk['CFD_ref'], errors='coerce')
        qry = pd.to_numeric(chunk['CFD_qry'], errors='coerce')
        delta = (qry - ref).to_numpy(dtype=float)
        valid = np.isfinite(delta)
        totals['valid_CFD_pairs'] += int(valid.sum())
        d = delta[valid]
        totals['unchanged_abs_delta_le_1e-9'] += int((np.abs(d) <= 1e-9).sum())
        totals['decreased_delta_lt_minus_1e-9'] += int((d < -1e-9).sum())
        totals['increased_delta_gt_1e-9'] += int((d > 1e-9).sum())
        for t, label in [(0.1, '0.1'), (0.25, '0.25'), (0.5, '0.5'), (0.9, '0.9')]:
            totals[f'delta_le_minus_{label}'] += int((d <= -t).sum())
            totals[f'delta_ge_plus_{label}'] += int((d >= t).sum())
        if {'syntenic_block_ref', 'syntenic_block_qry'} <= headers:
            r = chunk['syntenic_block_ref']
            q = chunk['syntenic_block_qry']
            known = r.notna() & q.notna()
            totals['both_blocks_assigned'] += int(known.sum())
            totals['same_syntenic_block'] += int((known & (r == q)).sum())
            totals['different_syntenic_block'] += int((known & (r != q)).sum())
        if 'proj_dist' in headers:
            dist = pd.to_numeric(chunk['proj_dist'], errors='coerce')
            totals['distance_le_1000_bp'] += int((dist <= 1000).sum())
            totals['distance_gt_1000_bp'] += int((dist > 1000).sum())
        if 'pam_qry' in headers:
            totals['query_noncanonical_NGG'] += int((chunk['pam_qry'].astype(str) != 'GG').sum())

    print(f'\nFILE: {path}')
    for name, value in totals.items():
        denom = totals['valid_CFD_pairs'] if name.startswith(('unchanged_', 'decreased_', 'increased_', 'delta_')) else totals['rows']
        fraction = f' ({value / denom * 100:.3f}%)' if denom and name != 'rows' else ''
        print(f'  {name:35s} {value:>12,d}{fraction}')
    if totals['different_syntenic_block']:
        print('  WARNING: some matches use DIFFERENT syntenic blocks; revisit orthology merge.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('tables', nargs='+', type=Path)
    ap.add_argument('--chunksize', type=int, default=100_000)
    args = ap.parse_args()
    for path in args.tables:
        summarize(path, args.chunksize)
