# -*- coding: utf-8 -*-
"""Batch-annotate a metabolite table with BOND bond codes (generic example).

Usage:
    python scripts/batch_annotate.py input.csv output.tsv --smiles smiles --name name

Input: any CSV/TSV with a SMILES column (encoding auto-detected: utf-8-sig/gbk/utf-8).
Output: TSV with one row per compound: bond_codes, n_bond_types, glycosidic_linkage,
        n_sugar_rings, status.
"""
import argparse
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from bond.judge_core import run_table  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description='BOND batch bond-type annotation')
    ap.add_argument('input', help='input table (CSV/TSV with a SMILES column)')
    ap.add_argument('output', help='output TSV')
    ap.add_argument('--smiles', default='smiles', help='SMILES column name')
    ap.add_argument('--name', default='name', help='compound name column')
    a = ap.parse_args()
    out, cnt, n_ok, n_fail = run_table(a.input, a.output, a.smiles, a.name)
    total = n_ok + n_fail
    hit = (out['n_bond_types'] > 0).sum()
    print(f'rows: {total} | annotated: {n_ok} | parse failed: {n_fail}')
    print(f'coverage (>=1 bond code): {hit} ({hit / total:.1%})')
    print('bond-code counts:', dict(sorted(cnt.items(), key=lambda x: -x[1])))


if __name__ == '__main__':
    main()
